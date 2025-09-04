#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dumb scripts launcher:
- Starts experiment scripts (each script calls the internal runner, which submits many Slurm jobs/arrays)
- Polls runner status for a dashboard
- Optionally, post-launch sweep to auto-release *held* Slurm jobs (JobHeldUser) for each running experiment.
  This uses Slurm's ArrayJobID (%A) so `scontrol release` is called with the correct base ID.

Usage
=====
    python run_all_long_exp.py \
        --scripts-dir scripts/long_shedders_experiments \
        --from-index scripts/long_shedders_experiments/INDEX_exp_scripts.csv \
        --concurrency 4 \
        --poll-seconds 30 \
        [--only 01_generate_data_foo,02_generate_data_bar] \
        [--exp-n 1] \
        [--use-index-order] \
        [--auto-release-held] \
        [--held-threshold-seconds 0] \
        [--debug]
"""

import argparse
import csv
import os
import sys
import time
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import traceback

import simplicity.runners.slurm as slurm
import simplicity.dir_manager as dm

LAUNCHER_FALLBACK_LOGROOT = Path("Data/_launcher_logs")
# broadened checks, all done against lower() of reason
HELD_REASONS_SUBSTR = ("held", "userhold", "jobhelduser", "jobheldadmin")

# ------------------- Slurm helpers (use ArrayJobID for release) -------------------

def _current_user() -> str:
    return os.getenv("USER") or os.getenv("LOGNAME") or "UNKNOWN"

def _squeue_rows_for_user(debug: bool=False) -> List[Tuple[str, str, str, str, str, str]]:
    """
    Returns rows: (array_jobid, jobid_display, name, state, elapsed, reason)
      - %A = ArrayJobID (base numeric id; empty for non-arrays)
      - %i = JOBID (may be like 12345_[1-100] or 12345_7)
      - %.512j = full job name (wide to avoid truncation)
    """
    cmd = ["squeue", "-h", "-u", _current_user(), "-o", "%A|%i|%.512j|%T|%M|%R"]
    try:
        out = subprocess.check_output(cmd, text=True)
        if debug:
            print(f"[launcher][debug] squeue cmd: {' '.join(cmd)}")
            print(f"[launcher][debug] squeue raw lines: {len(out.strip().splitlines())}")
    except Exception as e:
        if debug:
            print(f"[launcher][debug][warn] squeue failed: {e}")
        return []
    rows = []
    for line in out.strip().splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 6:
            if debug:
                print(f"[launcher][debug][warn] skipping malformed squeue line: {line!r}")
            continue
        array_jid, jobid_disp, name, state, elapsed, reason = parts
        rows.append((array_jid, jobid_disp, name, state, elapsed, reason))
    return rows

def _parse_elapsed_to_seconds(elapsed: str) -> int:
    if not elapsed:
        return 0
    days = 0
    if "-" in elapsed:
        d, rest = elapsed.split("-", 1)
        try: days = int(d)
        except: days = 0
        elapsed = rest
    parts = [int(p) for p in elapsed.split(":")]
    if len(parts) == 3: h, m, s = parts
    elif len(parts) == 2: h, m, s = 0, parts[0], parts[1]
    else: h, m, s = 0, 0, parts[0]
    return days*86400 + h*3600 + m*60 + s

def _base_id(array_jid: str, jobid_disp: str) -> str:
    """
    Prefer ArrayJobID (%A). If empty (non-array), derive from %i by stripping _task or [range].
    """
    if array_jid:
        return array_jid
    x = jobid_disp.split("_", 1)[0]
    x = x.split("[", 1)[0]
    return x

def _is_pending(state: str) -> bool:
    """Normalize state so we accept either 'PD' or 'PENDING' (cluster-dependent)."""
    s = (state or "").strip().upper()
    return s == "PD" or s == "PENDING"

def _debug_explain_nonmatches(exp_name: str, rows: List[Tuple[str, str, str, str, str, str]]) -> None:
    """
    When --debug is on, print why pending jobs didn't qualify for release.
    """
    print(f"[launcher][debug] explain non-matches for experiment='{exp_name}':")
    for array_jid, jobid_disp, name, state, elapsed, reason in rows:
        why = []
        if not _is_pending(state):
            why.append(f"state={state}!=PENDING/PD")
        rlow = reason.lower()
        if not any(tag in rlow for tag in HELD_REASONS_SUBSTR):
            why.append(f"reason='{reason}' not held-like")
        if exp_name not in name:
            why.append(f"name does not contain experiment ('{name}')")
        if why:
            print(f"  jobid_disp={jobid_disp} array={array_jid or '-'} name='{name}' -> skip: {', '.join(why)}")

def find_held_jobs_for_experiment(exp_name: str, debug: bool=False) -> List[Tuple[str, str, int]]:
    """
    Return [(base_id, reason, elapsed_sec), ...] for held jobs whose NAME contains the experiment name.
    Filters:
      - state in {PD, PENDING}
      - reason contains HELD_REASONS_SUBSTR
      - name contains the full experiment name (substring)
    """
    rows = _squeue_rows_for_user(debug=debug)
    matches: List[Tuple[str, str, int]] = []
    for array_jid, jobid_disp, name, state, elapsed, reason in rows:
        if not _is_pending(state):
            continue
        rlow = reason.lower()
        if not any(tag in rlow for tag in HELD_REASONS_SUBSTR):
            continue
        if exp_name not in name:
            continue
        base = _base_id(array_jid, jobid_disp)
        matches.append((base, reason, _parse_elapsed_to_seconds(elapsed)))
    if debug and not matches:
        _debug_explain_nonmatches(exp_name, rows)
    return matches

def scontrol_release(job_ids: List[str], debug: bool=False) -> Tuple[bool, str, int]:
    ids = sorted(set(job_ids))
    if not ids:
        return True, "nothing to release", 0
    cmd = ["scontrol", "release"] + ids
    try:
        if debug:
            print(f"[launcher][debug] exec: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        out = (proc.stdout or "").strip()
        err = (proc.stderr or "").strip()
        if debug:
            print(f"[launcher][debug] rc={proc.returncode} stdout={out!r} stderr={err!r}")
        return (proc.returncode == 0), (out + (("\n"+err) if err else "")), proc.returncode
    except Exception as e:
        if debug:
            print(f"[launcher][debug][error] scontrol release exception: {e}")
        return False, str(e), -1

def scancel(job_ids: List[str], debug: bool=False) -> Tuple[bool, str, int]:
    ids = sorted(set(job_ids))
    if not ids:
        return True, "nothing to cancel", 0
    cmd = ["scancel"] + ids
    try:
        if debug:
            print(f"[launcher][debug] exec: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        out = (proc.stdout or "").strip()
        err = (proc.stderr or "").strip()
        if debug:
            print(f"[launcher][debug] rc={proc.returncode} stdout={out!r} stderr={err!r}")
        return (proc.returncode == 0), (out + (("\n"+err) if err else "")), proc.returncode
    except Exception as e:
        if debug:
            print(f"[launcher][debug][error] scancel exception: {e}")
        return False, str(e), -1

# ------------------- Base helpers -------------------

def find_scripts_from_index(index_csv: Path, scripts_dir: Path) -> List[Path]:
    if not index_csv.is_file():
        raise FileNotFoundError(f"Index CSV not found: {index_csv}")
    if not scripts_dir.is_dir():
        raise FileNotFoundError(f"Scripts directory not found: {scripts_dir}")
    scripts: List[Path] = []
    with index_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "filename" not in (reader.fieldnames or []):
            raise ValueError(f"Index {index_csv} missing required 'filename' column; has: {reader.fieldnames}")
        for row in reader:
            fn = (row.get("filename") or "").strip()
            if not fn:
                continue
            p = scripts_dir / fn
            if p.suffix != ".py":
                p = p.with_suffix(".py")
            if not p.is_file():
                raise FileNotFoundError(f"Listed script not found: {p}")
            scripts.append(p)
    return scripts

def filter_only(scripts: List[Path], only_csv: str) -> List[Path]:
    if not only_csv.strip():
        return scripts
    wanted = {s.strip() for s in only_csv.split(",") if s.strip()}
    def match(p: Path) -> bool:
        return (p.stem in wanted) or (p.name in wanted)
    filtered = [p for p in scripts if match(p)]
    missing = [w for w in wanted if not any((p.stem == w or p.name == w) for p in scripts)]
    if missing:
        print(f"[launcher][warn] --only filter not matched: {', '.join(sorted(missing))}")
    return filtered

def experiment_name_from_script(script_path: Path) -> str:
    return script_path.stem

def launcher_log_paths(experiment_name: str) -> Tuple[Path, Path]:
    logs_dir = LAUNCHER_FALLBACK_LOGROOT / experiment_name / "runner_logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return (logs_dir / "launcher_out.log", logs_dir / "launcher_err.log")

def start_script(script_path: Path, exp_n: int, debug: bool = False,
                 extra_env: Optional[Dict[str, str]] = None) -> subprocess.Popen:
    experiment_name = experiment_name_from_script(script_path)
    out_log, err_log = launcher_log_paths(experiment_name)
    out_fp = open(out_log, "ab", buffering=0)
    err_fp = open(err_log, "ab", buffering=0)
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    if extra_env:
        env.update(extra_env)
    cmd = [sys.executable, "-u", str(script_path), "slurm", str(exp_n)]
    if debug:
        print(f"[launcher][debug] exec: {' '.join(cmd)}")
        print(f"[launcher][debug] logs: out={out_log} err={err_log}")
    proc = subprocess.Popen(cmd, stdout=out_fp, stderr=err_fp, env=env)
    return proc

def poll_status(experiment_name: str):
    try:
        return slurm.poll_simulations_status(experiment_name)
    except Exception:
        return None

def compact_status_line(status) -> str:
    return (
        f"tot={status.total} | sub={status.submitted} rel={status.released} "
        f"left={status.left} | pen={status.pending} sta={status.started} "
        f"run={status.running} cmp={status.completed} fail={status.failed}"
    )

# ------------------- Main loop -------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scripts-dir", type=str, default="scripts/long_shedders_experiments")
    parser.add_argument("--from-index", type=str, default="scripts/long_shedders_experiments/INDEX_exp_scripts.csv")
    parser.add_argument("--concurrency", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--only", type=str, default="")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--exp-n", type=int, default=1,
                        help="experiment_number to pass to each script (overridden by --use-index-order)")
    parser.add_argument("--use-index-order", action="store_true",
                        help="Use 'order' column from INDEX as experiment_number per script")
    # Post-launch held handling (no change to runner behavior)
    parser.add_argument("--auto-release-held", action="store_true",
                        help="Periodically scan squeue and scontrol release held jobs for each running experiment.")
    parser.add_argument("--held-threshold-seconds", type=int, default=0,
                        help="Release a held job after it has been seen held for at least this many seconds (default 0 = immediately).")
    args = parser.parse_args()

    scripts_dir = Path(args.scripts_dir)
    index_csv = Path(args.from_index)

    try:
        all_scripts = find_scripts_from_index(index_csv=index_csv, scripts_dir=scripts_dir)
    except Exception:
        print("[launcher][error] Failed reading index or scripts:")
        traceback.print_exc()
        sys.exit(2)

    all_scripts = filter_only(all_scripts, args.only)

    order_map: Dict[str, int] = {}
    try:
        with index_csv.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames and "filename" in reader.fieldnames and "order" in reader.fieldnames:
                for row in reader:
                    fn = (row.get("filename") or "").strip()
                    try:
                        order_map[Path(fn).name] = int(row.get("order"))
                    except Exception:
                        pass
    except Exception:
        order_map = {}

    print(f"[launcher] Found {len(all_scripts)} script(s) from INDEX. Concurrency={args.concurrency}")
    if args.only:
        print(f"[launcher] Filtered with --only: {args.only}")
    if args.use_index_order:
        print("[launcher] Using INDEX 'order' column as experiment_number per script")
    else:
        print(f"[launcher] Using --exp-n={args.exp_n} for all scripts")
    if args.auto_release_held:
        print(f"[launcher] Auto-release held enabled (threshold={args.held_threshold_seconds}s)")
    if args.debug:
        print(f"[launcher][debug] user={_current_user()}")

    queue: List[Path] = list(all_scripts)
    running: Dict[str, Tuple[Path, subprocess.Popen]] = {}
    completed: List[str] = []
    failed: List[str] = []

    held_seen_at: Dict[str, float] = {}  # base_id -> first seen time
    last_print = 0.0

    # Helper closure to do one sweep for one experiment (for reuse)
    def sweep_for_experiment(exp: str):
      if not args.auto_release_held:
          return
  
      rows = _squeue_rows_for_user(debug=args.debug)
      # Gather all matches for this experiment
      matches = []
      for array_jid, jobid_disp, name, state, elapsed, reason in rows:
          if not _is_pending(state):
              continue
          if exp not in name:
              continue
          rlow = reason.lower()
          if not any(tag in rlow for tag in HELD_REASONS_SUBSTR):
              continue
          base = _base_id(array_jid, jobid_disp)
          matches.append((base, reason, _parse_elapsed_to_seconds(elapsed), name))
  
      if args.debug and not matches:
          _debug_explain_nonmatches(exp, rows)
          return
  
      now = time.time()
      to_release: set[str] = set()
      for base_id, reason, elapsed_sec, name in matches:
          first_seen = held_seen_at.get(base_id)
          if first_seen is None:
              held_seen_at[base_id] = now
              first_seen = now
          seen_for = int(now - first_seen)
  
          print(f"    [held] base={base_id} name='{name}' reason='{reason}' seen_for={seen_for}s (elapsed={elapsed_sec}s)")
  
          rlow = reason.lower()
          if "jobheldadmin" in rlow:
              # Don’t cancel automatically; just explain.
              print("    [held] admin-held detected; skipping auto-release (requires admin).")
              continue
  
          # Always allow release for user-held; respect threshold (default 0)
          if seen_for >= args.held_threshold_seconds:
              to_release.add(base_id)
  
      if to_release:
          cmd = ["scontrol", "release"] + sorted(to_release)
          if args.debug:
              print(f"[launcher][debug] exec: {' '.join(cmd)}")
          proc = subprocess.run(cmd, capture_output=True, text=True)
          out = (proc.stdout or "").strip()
          err = (proc.stderr or "").strip()
          rc = proc.returncode
          print(("    [held] released: " if rc == 0 else f"    [held][warn] release failed (rc={rc}): ")
                + (out + (("\n" + err) if err else "") or "(no message)"))
  

    while queue or running:
        # Fill up to concurrency
        while queue and len(running) < args.concurrency:
            script = queue.pop(0)
            exp = experiment_name_from_script(script)
            exp_n = order_map.get(script.name, args.exp_n) if args.use_index_order else args.exp_n
            print(f"[launcher] Starting: {exp}  (experiment_number={exp_n})")
            try:
                proc = start_script(script_path=script, exp_n=exp_n, debug=args.debug)
                running[exp] = (script, proc)
                # Do an immediate sweep right after starting (helps catch fresh holds quickly)
                sweep_for_experiment(exp)
            except Exception as e:
                print(f"[launcher] FAILED to start {exp}: {e}", file=sys.stderr)
                if args.debug:
                    traceback.print_exc()
                failed.append(exp)

        # Check child processes
        to_remove = []
        for exp, (script, proc) in list(running.items()):
            ret = proc.poll()
            if ret is not None:
                if ret == 0:
                    print(f"[launcher] FINISHED: {exp} (exit=0)")
                    completed.append(exp)
                else:
                    print(f"[launcher] FAILED: {exp} (exit={ret})", file=sys.stderr)
                    print(
                        f"[launcher][hint] Check bootstrap logs under {LAUNCHER_FALLBACK_LOGROOT/exp} "
                        f"and Slurm logs under Data/{exp}/slurm/slurm_logs/ if present.",
                        file=sys.stderr,
                    )
                    failed.append(exp)
                to_remove.append(exp)
        for exp in to_remove:
            running.pop(exp, None)

        now = time.time()
        if now - last_print >= max(5, args.poll_seconds):
            last_print = now
            if running:
                print("\n[launcher] ===== Live status (from internal runner) =====")
                for exp in sorted(running.keys()):
                    st = poll_status(exp)
                    if st is None:
                        print(f"  {exp:>32}: (warming up / no seeds yet)")
                    else:
                        print(f"  {exp:>32}: {compact_status_line(st)}")
                    # Post-launch held-job handling (strictly reactive; runner remains unchanged)
                    sweep_for_experiment(exp)
                print("[launcher] =============================================\n")

        time.sleep(1)

    # Summary
    print("\n[launcher] All done.")
    print(f"[launcher] Completed: {len(completed)}")
    if completed:
        for exp in completed:
            print(f"  - {exp}")
    print(f"[launcher] Failed: {len(failed)}")
    if failed:
        for exp in failed:
            print(f"  - {exp}")
    print(
        "[launcher] Tip: bootstrap logs: Data/_launcher_logs/<experiment>/runner_logs/\n"
        "          per-task Slurm logs: Data/<experiment>/slurm/slurm_logs/\n"
        "          Held handling is post-launch only; the internal runner is unchanged."
    )

if __name__ == "__main__":
    main()
