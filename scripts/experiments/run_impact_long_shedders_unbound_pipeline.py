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
# impact_long_shedders_unbound -- SINGLE PIPELINE ENTRY POINT
# ----------------------------------------------------------------------------
# The second pipeline, run alongside impact_long_shedders (untouched). Each
# cohort is calibrated ALONE and then they are let loose together, so the
# population-level clock is an output rather than a constraint:
#
#   1. cal_1  100% long shedders, sweep NSR_long, intra-host clock
#   2. cal_2  0% long shedders, sweep standard NSR, global clock -> ONE
#             standard rate for every scenario, then the frozen table
#   3. exp    production, both rates applied exactly as calibrated, run ONCE
#             PER CONSENSUS MODE off the SAME frozen table: one production
#             pipeline for the argmax consensus and one for the distributional
#             consensus, differing in nothing but how the distance from the
#             consensus is measured
#   4. sanity combined root-to-tip grid over real output
#   5. artifacts zipped into one file per run
#
# Stages 1 and 2 are independent -- cal_2 never reads cal_1's sweep, only its
# fitted result when assembling the table -- so they could in principle be
# submitted together. They are kept sequential here for the same reason the
# bound pipeline does: one log, one failure point, easy to follow.
#
# Experiment names are all prefixed impact_long_shedders_unbound_*, so this
# pipeline's Data/ output can never collide with the bound one's at the same
# --exp-num.
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
from impact_long_shedders_unbound_config import (
    LONG_NSR_EXP_NAME, STD_NSR_EXP_NAME, PROD_EXP_NAME, SETUP_DIR_TEMPLATE,
    CALIBRATION_SKIPPED_OUTPUTS,
    TABLE_FILENAME, SCENARIOS, CAL1_ISOLATED_FIXED_PARAMS, USER_FIXED_PARAMS,
    add_slurm_resource_args, CONSENSUS_MODES, prod_exp_name,
)

CHECK_SCRIPT = os.path.join(SCRIPT_DIR, os.pardir, "check_completed_simulations.py")

_SBATCH_JOB_ID_RE = re.compile(r'Submitted batch job (\d+)')
SANITY_PLOT_POLL_INTERVAL_S = 15

# cal_2 carries a full population across a 3-year window at up to NSR 1e-3;
# every OOM kill the bound pipeline ever took was in its equivalent stage.
CAL2_SLURM_MEM = "5G"



def _log(log_fh, line):
    print(line)
    log_fh.write(line + "\n")
    log_fh.flush()


def run_stage(cmd, log_fh):
    """Run one stage, streaming output live and to the log. Raises SystemExit
    on non-zero exit."""
    header = f"\n$ {' '.join(cmd)}"
    print(header)
    log_fh.write(header + "\n")

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for raw_line in proc.stdout:
        line = raw_line.rstrip("\n")
        print(line)
        # blank lines are printed liberally as visual spacing and carry no
        # diagnostic value in a saved log; everything else is kept
        if line.strip():
            log_fh.write(line + "\n")
            log_fh.flush()
    proc.wait()

    if proc.returncode != 0:
        raise SystemExit(f"[FAILED] stage exited {proc.returncode}: {' '.join(cmd)}")


def report_simulation_health(exp_names, log_fh, label, skip_outputs=()):
    """Per-grid-point completeness and an extinction/saturation autopsy of
    anything that ended early. Reporting only -- never aborts the run.

    skip_outputs names outputs the stage was told not to write. The checker
    requires a fixed list of files, so without this a calibration stage -- which
    writes four of them on purpose since v2.4.65 -- is reported as every
    simulation incomplete. The run would be fine and the report would say
    otherwise, which is the failure mode this pipeline can least afford.
    """
    _log(log_fh, f"\n===== simulation health: {label} =====")
    for name in exp_names:
        cmd = [sys.executable, CHECK_SCRIPT, name]
        _log(log_fh, f"\n$ {' '.join(cmd)}")
        env = dict(os.environ)
        if skip_outputs:
            env["SIMPLICITY_SKIP_OUTPUTS"] = ",".join(skip_outputs)
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        except Exception as exc:
            _log(log_fh, f"  [skip] health report failed to run: {exc}")
            continue
        for stream in (result.stdout, result.stderr):
            for line in stream.splitlines():
                _log(log_fh, line)
        if result.returncode != 0:
            _log(log_fh, f"  [skip] health report exited {result.returncode} for {name}")


def submit_sanity_plots(exp_num, target_osr_std, target_osr_long, log_fh,
                        consensus_modes):
    """One sanity grid per consensus mode. Each reads its own production
    experiments and writes into its own plots directory, both derived from that
    pipeline's experiment name."""
    sanity_sh = os.path.join(SCRIPT_DIR, "submit_sanity_plot_unbound.sh")
    submitted = []
    for mode in consensus_modes:
        exp_name = prod_exp_name(mode)
        cmd = ["sbatch", sanity_sh, str(exp_num), str(target_osr_std),
               str(target_osr_long), exp_name]
        header = f"\n$ {' '.join(cmd)}"
        print(header)
        log_fh.write(header + "\n")

        result = subprocess.run(cmd, capture_output=True, text=True)
        for stream in (result.stdout, result.stderr):
            if stream:
                print(stream, end="")
                log_fh.write(stream)
        if result.returncode != 0:
            raise SystemExit(
                f"[FAILED] sbatch submission for the {mode} sanity plot")
        submitted.append((mode, result.stdout.strip()))

    return submitted


def wait_for_sanity_plots(submitted, log_fh, poll_interval=SANITY_PLOT_POLL_INTERVAL_S):
    """Block until every sanity-plot job has left Slurm's queue, so the plots
    exist on disk before archiving."""
    job_ids = {}
    for scenario, sbatch_out in submitted:
        m = _SBATCH_JOB_ID_RE.search(sbatch_out)
        if m:
            job_ids[m.group(1)] = scenario

    if not job_ids:
        return

    _log(log_fh, f"\nWaiting for {len(job_ids)} sanity plot job(s): {', '.join(job_ids)}")
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
    return os.path.join("Data", "pipeline_artifacts",
                        f"{PROD_EXP_NAME}_#{exp_num}_artifacts.zip")


def write_artifacts_archive(archive_path, exp_num, log_file):
    os.makedirs(os.path.dirname(archive_path), exist_ok=True)
    candidates = [
        log_file,
        os.path.join("Data", f"{LONG_NSR_EXP_NAME}_#{exp_num}", "05_Plots",
                     f"{LONG_NSR_EXP_NAME}_#{exp_num}_long_nsr_calibration_fit.png"),
        os.path.join("Data", f"{STD_NSR_EXP_NAME}_#{exp_num}", "05_Plots",
                     f"{STD_NSR_EXP_NAME}_#{exp_num}_std_nsr_calibration_fit.png"),
    ]
    # one sanity grid per consensus mode
    for mode in CONSENSUS_MODES:
        sanity_exp = f"{prod_exp_name(mode)}_sanity_#{exp_num}"
        candidates.append(
            os.path.join("Data", sanity_exp, "05_Plots",
                         f"{sanity_exp}_all_scenarios_global_vs_intrahost.png"))
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in candidates:
            if os.path.isfile(path):
                zf.write(path, arcname=os.path.basename(path))
            else:
                print(f"  [skip] artifact not found: {path}")


def compress_results(exp_num, log_fh, keep_slurm_logs=True):
    """Pack everything this run produced into Data_Export/*.tar.xz.

    Hands Data/*_#<exp_num> to package_simplicity_data.sh, which globs every
    experiment directory from this run -- both calibrations, setup data, every
    production scenario of every consensus mode, both sanity directories --
    with no bound/unbound assumption and no run numbers written out anywhere.

    Submitted to Slurm (16 cpus) rather than run here: pixz over several GB
    does not belong on a login node. The pipeline does not wait for it; the
    parts appear in Data_Export/ when the packing job finishes.
    """
    script = os.path.join(SCRIPT_DIR, os.pardir, "package_simplicity_data.sh")
    cmd = [script, "Data", f"_#{exp_num}", "--slurm"]
    env = dict(os.environ,
               SIMPLICITY_KEEP_SLURM_LOGS="1" if keep_slurm_logs else "0")
    header = f"\n$ {' '.join(cmd)}"
    _log(log_fh, header)
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, env=env)
    for line in result.stdout.splitlines():
        _log(log_fh, line)
    if result.returncode != 0:
        # the data is on disk and intact either way, so this is reported and
        # not raised: a packing failure must not read as a failed pipeline
        _log(log_fh, f"[warn] packaging exited {result.returncode}; "
                     f"the run itself is unaffected. Pack it by hand with:\n"
                     f"  {' '.join(cmd)}")
        return None
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(
        description="Run the full impact_long_shedders_unbound pipeline: each "
                    "cohort calibrated alone, population clock left free.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number, shared by every stage.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--cal-seeds', type=int, default=30,
                        help="Seeds per grid point, cal_1 and cal_2.")
    parser.add_argument('--exp-seeds', type=int, default=30,
                        help="Seeds per scenario, production stage.")
    parser.add_argument('--consensus', type=str, nargs='+',
                        choices=list(CONSENSUS_MODES), default=list(CONSENSUS_MODES),
                        help="Consensus modes to run in production, one "
                            "pipeline each, off the same calibration. "
                            "Default: both.")
    parser.add_argument('--skip-cal1', action='store_true',
                        help=f"Reuse an existing {LONG_NSR_EXP_NAME}_#{{exp_num}}.")
    parser.add_argument('--r-cal1', type=float, default=CAL1_ISOLATED_FIXED_PARAMS['R'])
    parser.add_argument('--r-cal2', type=float, default=USER_FIXED_PARAMS['R'],
                        help="R for cal_2 and production; unrelated to R_long.")
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'])
    parser.add_argument('--log-file', type=str, default=None)
    parser.add_argument('--no-compress', action='store_true',
                        help="Skip the final packaging step. By default the "
                            "run's Data/*_#<exp_num> directories are packed "
                            "into Data_Export/*.tar.xz by a Slurm job.")
    parser.add_argument('--drop-slurm-logs', action='store_true',
                        help="Let the packager delete each experiment's slurm/ "
                            "directory (its old behaviour). They are kept by "
                            "default: they are the only record of why a task "
                            "failed, and the deletion is in place on Data/.")
    add_slurm_resource_args(parser)
    parser.add_argument('--slurm-mem-cal2', type=str, default=CAL2_SLURM_MEM,
                        help=f"Per-task SLURM memory for cal_2 only (default "
                            f"{CAL2_SLURM_MEM}).")
    args = parser.parse_args()

    log_file = args.log_file or os.path.join(
        "Data", "pipeline_logs",
        f"{PROD_EXP_NAME}_pipeline_#{args.exp_num}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    py = sys.executable
    cal1_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_unbound_cal_1.py")
    cal2_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_unbound_cal_2.py")
    exp_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_unbound_exp.py")

    slurm_res_args = ["--slurm-mem", args.slurm_mem, "--slurm-time", args.slurm_time]

    with open(log_file, "a") as log_fh:
        _log(log_fh, f"\n===== {PROD_EXP_NAME} pipeline: exp_num={args.exp_num} "
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
                      "--ih-virus-emergence-rate", str(args.ih_virus_emergence_rate),
                      *slurm_res_args], log_fh)

        report_simulation_health([f"{LONG_NSR_EXP_NAME}_#{args.exp_num}"],
                                 log_fh, "stage 1 (long shedders alone)",
                                 skip_outputs=CALIBRATION_SKIPPED_OUTPUTS)

        run_stage([py, cal2_path,
                  "--exp-num", str(args.exp_num),
                  "--runner", args.runner,
                  "--seeds", str(args.cal_seeds),
                  "--target-osr-std", str(args.target_osr_std),
                  "--target-osr-long", str(args.target_osr_long),
                  "--R", str(args.r_cal2),
                  "--ih-virus-emergence-rate", str(args.ih_virus_emergence_rate),
                  "--slurm-mem", args.slurm_mem_cal2,
                  "--slurm-time", args.slurm_time], log_fh)

        report_simulation_health([f"{STD_NSR_EXP_NAME}_#{args.exp_num}"],
                                 log_fh, "stage 2 (standard alone)",
                                 skip_outputs=CALIBRATION_SKIPPED_OUTPUTS)

        # Both pipelines read the SAME frozen table, so the only thing that
        # differs between them is the consensus distance.
        for mode in args.consensus:
            run_stage([py, exp_path,
                      "--exp-num", str(args.exp_num),
                      "--runner", args.runner,
                      "--seeds", str(args.exp_seeds),
                      "--consensus", mode,
                      *slurm_res_args], log_fh)

            report_simulation_health(
                [f"{prod_exp_name(mode)}_{s['name']}_#{args.exp_num}"
                 for s in SCENARIOS],
                log_fh, f"stage 3 (production, consensus={mode})")

        submitted = submit_sanity_plots(
            args.exp_num, args.target_osr_std, args.target_osr_long, log_fh,
            args.consensus)
        wait_for_sanity_plots(submitted, log_fh)

        table_path = os.path.join(
            SETUP_DIR_TEMPLATE.format(exp_num=args.exp_num), TABLE_FILENAME)
        archive_path = artifacts_archive_path(args.exp_num)
        summary_lines = [
            "\n===== Pipeline summary =====",
            f"Frozen calibration table : {table_path}",
            f"Long calibration exp     : {LONG_NSR_EXP_NAME}_#{args.exp_num}",
            f"Standard calibration exp : {STD_NSR_EXP_NAME}_#{args.exp_num}",
            *[f"Production ({m:<12}): {prod_exp_name(m)}_<scenario>_#{args.exp_num}"
              for m in args.consensus],
            "Sanity plot jobs submitted:",
        ]
        for scenario, sbatch_out in submitted:
            summary_lines.append(f"  - {scenario}: {sbatch_out}")
        summary_lines.append(f"Full log: {log_file}")
        summary_lines.append(f"Artifacts archive: {archive_path}")
        if not args.no_compress:
            summary_lines.append("Full results    : Data_Export/*.tar.xz "
                                 "(packing job submitted below)")

        for line in summary_lines:
            _log(log_fh, line)

        write_artifacts_archive(archive_path, args.exp_num, log_file)

        if not args.no_compress:
            _log(log_fh, "\n===== Packing results for download =====")
            compress_results(args.exp_num, log_fh,
                             keep_slurm_logs=not args.drop_slurm_logs)


if __name__ == "__main__":
    main()
