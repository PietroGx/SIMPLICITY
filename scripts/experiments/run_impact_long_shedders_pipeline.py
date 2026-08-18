#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ============================================================================
# impact_long_shedders -- SINGLE PIPELINE ENTRY POINT
# ----------------------------------------------------------------------------
# Runs the whole pipeline, always in the same order, under one shared
# --exp-num:
#   1. impact_long_shedders_cal_1.py  (isolated long-shedder NSR calibration
#      + its calibration-fit plot, saved into that experiment's 05_Plots)
#   2. impact_long_shedders_cal_2.py  (in-context standard-NSR calibration,
#      one combined SLURM submission for all scenarios -> frozen NSR table)
#   3. impact_long_shedders_exp.py    (production runs, one per scenario)
#   4. submit_sanity_plot.sh          (ONE SLURM job producing a combined
#      4x2 sanity grid -- one row per long-shedder scenario, shared x/y
#      limits per column -- via plot_sot_sanity_regressions.py)
#   5. Once the sanity-plot job has left Slurm's queue, the pipeline log,
#      both calibration-fit plots, and the combined sanity grid are bundled
#      into one zip at Data/pipeline_artifacts/impact_long_shedders_#{exp_num}_
#      artifacts.zip -- a single place to review a run's outputs instead of
#      hunting across per-experiment 05_Plots folders.
#
# Every stage keeps its own independently-runnable CLI; this script is a
# thin sequencer that aborts immediately if any stage fails. All stage
# output is echoed live AND appended to --log-file, except the SLURM
# runner's periodic polling noise (simplicity/runners/slurm.py's
# SimulationsStatus(...) prints and "submitted N seeded simulations" lines),
# which is dropped from the log file only. The old "release up to N seeded
# simulations" ping was removed at the source (slurm.py no longer prints it
# at all -- it fired on every poll cycle once jobs started turning over).
#
# NSR sweep ranges (and per-scenario R_long) are NOT arguments here -- they
# live in impact_long_shedders_config.py (NSR_RANGES / SCENARIOS), reviewed/
# tuned between runs using the calibration-fit and sanity plots this
# pipeline produces.
# ============================================================================

import os
import re
import sys
import time
import zipfile
import argparse
import subprocess
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from impact_long_shedders_config import (
    LONG_NSR_EXP_NAME, STD_NSR_SWEEP_NAME, CAL1_ISOLATED_FIXED_PARAMS,
    USER_FIXED_PARAMS, add_slurm_resource_args,
)

_NOISE_PATTERNS = [
    re.compile(r'^SimulationsStatus\('),
    re.compile(r'^submitted \d+ seeded simulations$'),
]

_SBATCH_JOB_ID_RE = re.compile(r'Submitted batch job (\d+)')
SANITY_PLOT_POLL_INTERVAL_S = 15


def _is_slurm_monitor_noise(line):
    stripped = line.strip()
    return any(p.match(stripped) for p in _NOISE_PATTERNS)


def _log(log_fh, line):
    print(line)
    log_fh.write(line + "\n")
    log_fh.flush()


def run_stage(cmd, log_fh):
    """Run one pipeline stage, streaming output live and to the log file
    (minus SLURM-monitor noise). Raises SystemExit on non-zero exit."""
    header = f"\n$ {' '.join(cmd)}"
    print(header)
    log_fh.write(header + "\n")

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for raw_line in proc.stdout:
        line = raw_line.rstrip("\n")
        print(line)
        # Blank lines (printed liberally as visual spacing, and repeatedly by
        # the SLURM runner's polling wait) carry no diagnostic value in a
        # saved log -- drop them along with the named noise patterns.
        if line.strip() and not _is_slurm_monitor_noise(line):
            log_fh.write(line + "\n")
            log_fh.flush()
    proc.wait()

    if proc.returncode != 0:
        raise SystemExit(f"[FAILED] stage exited {proc.returncode}: {' '.join(cmd)}")


def submit_sanity_plots(exp_num, target_osr_std, target_osr_long, log_fh):
    """Dispatch the ONE combined sanity-plot job (4x2 grid, one row per
    long-shedder scenario, shared axis limits per column -- see
    plot_sot_sanity_regressions.py). Returns a single-entry list so callers
    that iterate over `submitted` (wait_for_sanity_plots, the summary print)
    don't need special-casing versus the old one-job-per-scenario shape."""
    sanity_sh = os.path.join(SCRIPT_DIR, "submit_sanity_plot.sh")
    cmd = ["sbatch", sanity_sh, str(exp_num), str(target_osr_std), str(target_osr_long)]
    header = f"\n$ {' '.join(cmd)}"
    print(header)
    log_fh.write(header + "\n")

    result = subprocess.run(cmd, capture_output=True, text=True)
    for stream in (result.stdout, result.stderr):
        if stream:
            print(stream, end="")
            log_fh.write(stream)
    if result.returncode != 0:
        raise SystemExit("[FAILED] sbatch submission for combined sanity plot")

    return [("combined", result.stdout.strip())]


def wait_for_sanity_plots(submitted, log_fh, poll_interval=SANITY_PLOT_POLL_INTERVAL_S):
    """Block until every submitted sanity-plot job has left Slurm's queue,
    so the plots actually exist on disk before archiving."""
    job_ids = {}
    for scenario, sbatch_out in submitted:
        m = _SBATCH_JOB_ID_RE.search(sbatch_out)
        if m:
            job_ids[m.group(1)] = scenario

    if not job_ids:
        return

    _log(log_fh, f"\nWaiting for {len(job_ids)} sanity plot job(s) to finish: "
                f"{', '.join(job_ids)}")
    while job_ids:
        result = subprocess.run(
            ["squeue", "-h", "-j", ",".join(job_ids), "-o", "%i"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        still_queued = set(result.stdout.split())
        for jid in list(job_ids):
            if jid not in still_queued:
                _log(log_fh, f"  [done] sanity plot job {jid} ({job_ids.pop(jid)})")
        if job_ids:
            time.sleep(poll_interval)


def artifacts_archive_path(exp_num):
    return os.path.join(
        "Data", "pipeline_artifacts", f"impact_long_shedders_#{exp_num}_artifacts.zip")


def write_artifacts_archive(archive_path, exp_num, log_file):
    """Bundle this run's log, both calibration-fit plots, and the combined
    sanity-plot grid into one zip -- a single place to check a run's
    outputs instead of hunting across per-experiment 05_Plots folders.
    A missing file (e.g. the sanity plot job itself failed) is skipped
    with a warning rather than aborting the archive."""
    os.makedirs(os.path.dirname(archive_path), exist_ok=True)

    sanity_exp_name = f"impact_long_shedders_sanity_#{exp_num}"
    candidates = [
        log_file,
        os.path.join(
            "Data", f"{LONG_NSR_EXP_NAME}_#{exp_num}", "05_Plots",
            f"{LONG_NSR_EXP_NAME}_#{exp_num}_long_nsr_calibration_fit.png"),
        os.path.join(
            "Data", f"{STD_NSR_SWEEP_NAME}_#{exp_num}", "05_Plots",
            f"{STD_NSR_SWEEP_NAME}_#{exp_num}_std_nsr_calibration_fit.png"),
        os.path.join(
            "Data", sanity_exp_name, "05_Plots",
            f"{sanity_exp_name}_all_scenarios_global_vs_intrahost.png"),
    ]

    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in candidates:
            if os.path.isfile(path):
                zf.write(path, arcname=os.path.basename(path))
            else:
                print(f"  [skip] artifact not found: {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run the full impact_long_shedders pipeline: cal_1 -> "
                    "cal_2 -> exp -> per-scenario sanity plots.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number, shared by every stage.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--cal-seeds', type=int, default=20,
                        help="Seeds per grid point, cal_1 and cal_2.")
    parser.add_argument('--exp-seeds', type=int, default=20,
                        help="Seeds per scenario, final production stage.")
    parser.add_argument('--skip-cal1', action='store_true',
                        help=f"Reuse an already-computed {LONG_NSR_EXP_NAME}_#{{exp_num}} "
                            "instead of rerunning Stage 1.")
    parser.add_argument('--r-cal1', type=float, default=CAL1_ISOLATED_FIXED_PARAMS['R'],
                        help="R for cal_1's isolated context "
                            f"(default {CAL1_ISOLATED_FIXED_PARAMS['R']}).")
    parser.add_argument('--r-cal2', type=float, default=USER_FIXED_PARAMS['R'],
                        help="R for cal_2 and production "
                            f"(default {USER_FIXED_PARAMS['R']}); unrelated to "
                            "R_long, which is set per scenario in "
                            "impact_long_shedders_config.SCENARIOS.")
    parser.add_argument('--log-file', type=str, default=None,
                        help="Defaults to "
                            "Data/pipeline_logs/impact_long_shedders_pipeline_#{exp_num}.log")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    log_file = args.log_file or os.path.join(
        "Data", "pipeline_logs", f"impact_long_shedders_pipeline_#{args.exp_num}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    py = sys.executable
    cal1_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_cal_1.py")
    cal2_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_cal_2.py")
    exp_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_exp.py")

    slurm_res_args = ["--slurm-mem", args.slurm_mem, "--slurm-time", args.slurm_time]

    with open(log_file, "a") as log_fh:
        _log(log_fh, f"\n===== impact_long_shedders pipeline: exp_num={args.exp_num} "
                    f"@ {datetime.now().isoformat(timespec='seconds')} =====")

        if args.skip_cal1:
            _log(log_fh, f"[skip] cal_1 skipped; reusing {LONG_NSR_EXP_NAME}_#{args.exp_num}")
        else:
            run_stage([py, cal1_path,
                      "--exp-num", str(args.exp_num),
                      "--runner", args.runner,
                      "--seeds", str(args.cal_seeds),
                      "--target-osr-long", str(args.target_osr_long),
                      "--R", str(args.r_cal1),
                      *slurm_res_args], log_fh)

        run_stage([py, cal2_path,
                  "--exp-num", str(args.exp_num),
                  "--runner", args.runner,
                  "--target-osr-std", str(args.target_osr_std),
                  "--target-osr-long", str(args.target_osr_long),
                  "--seeds", str(args.cal_seeds),
                  "--R", str(args.r_cal2),
                  *slurm_res_args], log_fh)

        run_stage([py, exp_path,
                  "--exp-num", str(args.exp_num),
                  "--runner", args.runner,
                  "--seeds", str(args.exp_seeds),
                  *slurm_res_args], log_fh)

        submitted = submit_sanity_plots(
            args.exp_num, args.target_osr_std, args.target_osr_long, log_fh)
        wait_for_sanity_plots(submitted, log_fh)

        table_path = os.path.join(
            "Data", f"impact_long_shedders_calibration_setup_data_#{args.exp_num}",
            "nsr_calibration_table.csv")
        archive_path = artifacts_archive_path(args.exp_num)
        summary_lines = [
            "\n===== Pipeline summary =====",
            f"Frozen calibration table : {table_path}",
            f"Long calibration exp     : {LONG_NSR_EXP_NAME}_#{args.exp_num}",
            f"Standard-NSR sweep exp   : {STD_NSR_SWEEP_NAME}_#{args.exp_num}",
            f"Production experiments   : impact_long_shedders_<scenario>_#{args.exp_num}",
            "Sanity plot jobs submitted:",
        ]
        for scenario, sbatch_out in submitted:
            summary_lines.append(f"  - {scenario}: {sbatch_out}")
        summary_lines.append(f"Full log: {log_file}")
        summary_lines.append(f"Artifacts archive (log + calibration + sanity plots): "
                            f"{archive_path}")

        for line in summary_lines:
            _log(log_fh, line)

        # Archived last, after every summary line above is on disk, so the
        # log copy bundled inside the zip is complete (including this run's
        # own archive path) rather than missing its own tail.
        write_artifacts_archive(archive_path, args.exp_num, log_file)


if __name__ == "__main__":
    main()
