#!/usr/bin/env python3
"""Every figure-preprocessing step after the simulations, from one command.

    python -I figures_gui/prepare.py --stamp S --slurm [--pack] [--yes] [--dry-run] [-- ARGS]
    python -I figures_gui/prepare.py --stamp S --local [--procs 2] [--pack] [--all-builds]
                                     [--only tx,runs] [--skip check] [--dry-run]

--slurm hands the chain to run_precompute.sh: its login-node checks, then the
grid array, tx, runs, explore, warm and, with --pack, the pack job, linked by
afterany. Anything after `--` goes to it unchanged (--limit, --resubmit, ...).

--local runs the same scripts here, one after another, each on --procs cores:

    check    run_precompute.sh --check-only
    tx       export_transmissions.py     -> <store>/tx
    runs     export_clade_explorer.py    -> <store>/runs
    explore  warm_explore.py             daily records, anchors, Figure 5's panel
    grid     compute.py at the operating point (compute.DEFAULTS);
             --all-builds: sweep.py, the whole grid (about a day on 4 cores)
    warm     warm_panels.py              after grid, tx and explore
    pack     pack_cache.sh --here        after the rest, only with --pack;
             --partial unless --all-builds

A step whose prerequisite failed is skipped, not run on half an input. Exit 0
only if every step that ran succeeded.

Figure 1 is not a step: it reads patient data that is not on the cluster and
must never be packed, so the studio builds it locally when its tab is opened.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY = sys.executable
STEPS = ("check", "tx", "runs", "explore", "grid", "warm", "pack")
NEEDS = {"warm": ("grid", "tx", "explore"), "pack": ("tx", "runs", "explore", "grid", "warm")}


def paths(store, stamp, exp_name):
    """sweep.py --print-paths: the one resolution of name, stamp and store."""
    cmd = [PY, "-I", os.path.join(HERE, "sweep.py"), "--print-paths", "--store", store,
           "--stamp", str(stamp)] + (["--exp-name", exp_name] if exp_name else [])
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    return dict(zip(("name", "stamp", "rel", "store", "cache", "root"), out))


def local_steps(a, p):
    """{step: argv} for --local, every script given the experiment explicitly."""
    gui = lambda f: os.path.join("figures_gui", f)
    ident = ["--exp-name", p["name"], "--seeds", str(a.seeds), "--procs", str(a.procs)]
    env_ident = ["--stamp", p["stamp"], "--exp-name", p["name"], "--seeds", str(a.seeds)]
    grid = ([PY, "-I", gui("sweep.py"), "--store", a.store, "--out", p["cache"], "--stamp", p["stamp"]]
            + ident if a.all_builds
            else [PY, "-I", gui("compute.py"), "--exp-num", p["stamp"], "--out", p["cache"]] + ident)
    return {
        "check": ["bash", gui("run_precompute.sh"), "--check-only"] + env_ident,
        "tx": [PY, "-I", gui("export_transmissions.py"), "--exp-num", p["stamp"],
               "--out", os.path.join(p["store"], "tx")] + ident,
        "runs": [PY, "-I", gui("export_clade_explorer.py"), "--stamp", p["stamp"],
                 "--out", os.path.join(p["store"], "runs")] + ident,
        "explore": [PY, "-I", gui("warm_explore.py"), "--stamp", p["stamp"], "--store", p["store"]] + ident,
        "grid": grid,
        "warm": [PY, "-I", gui("warm_panels.py"), "--stamp", p["stamp"], "--exp-name", p["name"],
                 "--base-store", a.store, "--procs", str(a.procs)],
        "pack": ["bash", gui("pack_cache.sh"), "--here"] + ([] if a.all_builds else ["--partial"]),
    }


def run_local(a, p, steps):
    cmds = local_steps(a, p)
    env = dict(os.environ, STORE=a.store, STAMP=p["stamp"], EXP_NAME=p["name"], SEEDS=str(a.seeds))
    done, failed, skipped = [], [], []
    for step in steps:
        miss = [n for n in NEEDS.get(step, ()) if n in steps and n not in done]
        if miss:
            skipped.append(step)
            print(f"[prepare] {step}: SKIPPED, needs {', '.join(miss)}", flush=True)
            continue
        print(f"[prepare] {step}: {' '.join(cmds[step])}", flush=True)
        if a.dry_run:
            done.append(step)
            continue
        t0 = time.time()
        rc = subprocess.run(cmds[step], cwd=ROOT, env=env).returncode
        print(f"[prepare] {step}: exit {rc} after {(time.time() - t0) / 60:.1f} min", flush=True)
        (done if rc == 0 else failed).append(step)
    print(f"[prepare] done: {', '.join(done) or '-'}; failed: {', '.join(failed) or '-'}; "
          f"skipped: {', '.join(skipped) or '-'}")
    return 1 if failed or skipped else 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    extra = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stamp", required=True)
    ap.add_argument("--exp-name", default="")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--slurm", action="store_true")
    mode.add_argument("--local", action="store_true")
    ap.add_argument("--store", default=os.path.join(ROOT, "Data", "figures"))
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=2, help="--local: cores per step (default 2)")
    ap.add_argument("--pack", action="store_true", help="archive into Data_Export/figures/")
    ap.add_argument("--all-builds", action="store_true", help="--local: the whole grid, not the operating point")
    ap.add_argument("--only", default="", help="--local: comma-separated steps")
    ap.add_argument("--skip", default="", help="--local: comma-separated steps")
    ap.add_argument("--yes", action="store_true", help="--slurm: no prompt")
    ap.add_argument("--dry-run", action="store_true", help="print, run nothing")
    a = ap.parse_args(argv)

    if a.slurm:
        if a.only or a.skip or a.all_builds:
            ap.error("--only, --skip and --all-builds are --local options")
        cmd = (["bash", os.path.join("figures_gui", "run_precompute.sh"), "--stamp", a.stamp, "--seeds", str(a.seeds)]
               + (["--exp-name", a.exp_name] if a.exp_name else []) + ([] if a.pack else ["--no-pack"])
               + (["--yes"] if a.yes else []) + (["--dry-run"] if a.dry_run else []) + extra)
        print(f"[prepare] {' '.join(cmd)}", flush=True)
        return subprocess.run(cmd, cwd=ROOT, env=dict(os.environ, STORE=a.store)).returncode

    pick = lambda s: [x for x in s.split(",") if x]
    bad = [x for x in pick(a.only) + pick(a.skip) if x not in STEPS]
    if bad:
        ap.error(f"unknown step(s) {', '.join(bad)}; steps are {', '.join(STEPS)}")
    steps = [s for s in STEPS if (not a.only or s in pick(a.only)) and s not in pick(a.skip)
             and (s != "pack" or a.pack or s in pick(a.only))]
    p = paths(a.store, a.stamp, a.exp_name)
    print(f"[prepare] {p['name']} #{p['stamp']} -> {p['root']}  ({', '.join(steps)}, --procs {a.procs})")
    return run_local(a, p, steps)


if __name__ == "__main__":
    sys.exit(main())
