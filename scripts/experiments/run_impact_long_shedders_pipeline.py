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
#   4. submit_sanity_plot.sh          (one SLURM job per long-shedder
#      scenario)
#
# Every stage keeps its own independently-runnable CLI; this script is a
# thin sequencer that aborts immediately if any stage fails. All stage
# output is echoed live AND appended to --log-file, except the SLURM
# runner's periodic polling noise (simplicity/runners/slurm.py's
# SimulationsStatus(...) prints and "submitted/release up to N seeded
# simulations" lines), which is dropped from the log file only.
#
# NSR sweep ranges are NOT arguments here -- they live in the user-editable
# Data/00_Reference_parameters/impact_long_shedders_nsr_ranges.json (see
# impact_long_shedders_config.read_nsr_ranges), reviewed/tuned between runs
# using the calibration-fit and sanity plots this pipeline produces.
# ============================================================================

import os
import re
import sys
import argparse
import subprocess
from datetime import datetime

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from impact_long_shedders_config import (
    LONG_NSR_EXP_NAME, STD_NSR_SWEEP_NAME, add_slurm_resource_args,
)

_NOISE_PATTERNS = [
    re.compile(r'^SimulationsStatus\('),
    re.compile(r'^submitted \d+ seeded simulations$'),
    re.compile(r'^release up to \d+ seeded simulations$'),
]


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
    table_path = os.path.join(
        "Data", f"impact_long_shedders_calibration_setup_data_#{exp_num}",
        "nsr_calibration_table.csv")
    df = pd.read_csv(table_path)

    sanity_sh = os.path.join(SCRIPT_DIR, "submit_sanity_plot.sh")
    submitted = []

    for _, row in df.iterrows():
        if float(row["long_shedders_ratio"]) <= 0.0:
            continue  # control has no long-shedder data to plot
        scenario = row["scenario_name"]
        cmd = ["sbatch", sanity_sh, str(exp_num), scenario,
              str(target_osr_std), str(target_osr_long)]
        header = f"\n$ {' '.join(cmd)}"
        print(header)
        log_fh.write(header + "\n")

        result = subprocess.run(cmd, capture_output=True, text=True)
        for stream in (result.stdout, result.stderr):
            if stream:
                print(stream, end="")
                log_fh.write(stream)
        if result.returncode != 0:
            raise SystemExit(f"[FAILED] sbatch submission for scenario {scenario}")
        submitted.append((scenario, result.stdout.strip()))

    return submitted


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
    parser.add_argument('--cal-seeds', type=int, default=25,
                        help="Seeds per grid point, cal_1 and cal_2.")
    parser.add_argument('--exp-seeds', type=int, default=25,
                        help="Seeds per scenario, final production stage.")
    parser.add_argument('--skip-cal1', action='store_true',
                        help=f"Reuse an already-computed {LONG_NSR_EXP_NAME}_#{{exp_num}} "
                            "instead of rerunning Stage 1.")
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
                      *slurm_res_args], log_fh)

        run_stage([py, cal2_path,
                  "--exp-num", str(args.exp_num),
                  "--runner", args.runner,
                  "--target-osr-std", str(args.target_osr_std),
                  "--target-osr-long", str(args.target_osr_long),
                  "--seeds", str(args.cal_seeds),
                  *slurm_res_args], log_fh)

        run_stage([py, exp_path,
                  "--exp-num", str(args.exp_num),
                  "--runner", args.runner,
                  "--seeds", str(args.exp_seeds),
                  *slurm_res_args], log_fh)

        submitted = submit_sanity_plots(
            args.exp_num, args.target_osr_std, args.target_osr_long, log_fh)

        table_path = os.path.join(
            "Data", f"impact_long_shedders_calibration_setup_data_#{args.exp_num}",
            "nsr_calibration_table.csv")
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

        for line in summary_lines:
            _log(log_fh, line)


if __name__ == "__main__":
    main()
