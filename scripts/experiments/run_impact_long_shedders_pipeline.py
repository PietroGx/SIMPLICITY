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
# The second pipeline, run alongside impact_long_shedders (untouched). Each
# cohort is calibrated ALONE and then they are let loose together, so the
# population-level clock is an output rather than a constraint:
#
#   1. cal_1  100% long shedders, sweep NSR_long, intra-host clock
#   2. cal_2  0% long shedders, sweep standard NSR, global clock -> ONE
#             standard rate for every scenario, then the frozen table
#   3. exp    production, both rates applied exactly as calibrated, run ONCE
#             PER CONSENSUS MODE off the SAME frozen table: one for the
#             distributional consensus (the default) and one for argmax,
#             differing in nothing but how the distance from the consensus
#             is measured
#   4. sanity combined root-to-tip grid over real output
#   5. artifacts zipped into one file per run
#
# Stages 1 and 2 are independent -- cal_2 never reads cal_1's sweep, only its
# fitted result when assembling the table -- so they could in principle be
# submitted together. They are kept sequential here for the same reason the
# bound pipeline does: one log, one failure point, easy to follow.
#
# Experiment names are all prefixed impact_long_shedders_*, so this
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
import simplicity.dir_manager as dm

from impact_long_shedders_config import (
    EXP_NAME, SETUP_DIR_TEMPLATE, calib_std_group,
    exp_group, groups_of_phase,
    CALIBRATION_SKIPPED_OUTPUTS,
    TABLE_FILENAME, SCENARIOS, CAL1_ISOLATED_FIXED_PARAMS, USER_FIXED_PARAMS,
    add_slurm_resource_args, CONSENSUS_MODES, DEFAULT_CONSENSUS, exp_name,
    sanity_exp_name, run_export_dirs,
)

CHECK_SCRIPT = os.path.join(SCRIPT_DIR, os.pardir, "check_completed_simulations.py")

_SBATCH_JOB_ID_RE = re.compile(r'Submitted batch job (\d+)')
SANITY_PLOT_POLL_INTERVAL_S = 15

def cal2_memory(slurm_mem, slurm_mem_cal2=None):
    """cal_2's per-task memory: the standard request unless overridden.

    It carried a 5G override until v2.4.106. That was sized when cal_2
    inherited production's 1095-day window and the BOUND pipeline swept it
    inside a mixed population -- where every OOM kill that pipeline ever took
    landed. Neither holds: build_cal2_settings pins CAL2_FINAL_TIME=365 and
    long_shedders_ratio=0, and the stage profiles as the CHEAPEST in the
    pipeline, 301 MB against the default 1G
    (tests/profile_simulation_cost.py). So the one stage asking five times
    the default was the one needing it least.
    """
    return slurm_mem_cal2 or slurm_mem



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


def report_simulation_health(experiment, groups, log_fh, label,
                             skip_outputs=()):
    """Per-grid-point completeness and an extinction/saturation autopsy of
    anything that ended early. Reporting only -- never aborts the run.

    Scoped to ONE phase's groups since 8e: the experiment holds every phase,
    and each group carries its own n_seeds, so an unscoped report would mix
    a calibration sweep with the experiment-phase scenarios and compare both
    against a single seed count.

    skip_outputs names outputs the stage was told not to write. The checker
    requires a fixed list of files, so without this a calibration stage -- which
    writes four of them on purpose since v2.4.65 -- is reported as every
    simulation incomplete. The run would be fine and the report would say
    otherwise, which is the failure mode this pipeline can least afford.
    """
    _log(log_fh, f"\n===== simulation health: {label} =====")
    cmd = [sys.executable, CHECK_SCRIPT, experiment]
    for group in groups:
        cmd += ["--group", group]
    _log(log_fh, f"\n$ {' '.join(cmd)}")
    env = dict(os.environ)
    if skip_outputs:
        env["SIMPLICITY_SKIP_OUTPUTS"] = ",".join(skip_outputs)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    except Exception as exc:
        _log(log_fh, f"  [skip] health report failed to run: {exc}")
        return
    for stream in (result.stdout, result.stderr):
        for line in stream.splitlines():
            _log(log_fh, line)
    if result.returncode != 0:
        _log(log_fh, f"  [skip] health report exited {result.returncode}")


def submit_sanity_plots(exp_num, target_osr_std, target_osr_long, log_fh,
                        consensus_modes):
    """One sanity grid per consensus mode. Each reads its own production
    experiments and writes into its own plots directory, both derived from that
    pipeline's experiment name."""
    sanity_sh = os.path.join(SCRIPT_DIR, "submit_sanity_plot.sh")
    submitted = []
    for mode in consensus_modes:
        name = exp_name(mode)
        cmd = ["sbatch", sanity_sh, str(exp_num), str(target_osr_std),
               str(target_osr_long), name]
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
                        f"{EXP_NAME}_#{exp_num}_artifacts.zip")


def write_artifacts_archive(archive_path, exp_num, log_file, consensus_modes):
    """The grids collected are the ones this run produced.

    This used to walk CONSENSUS_MODES, so a distribution run always reported
    the argmax grid missing and the other way round -- a skip line for a file
    nothing had been asked to write.
    """
    os.makedirs(os.path.dirname(archive_path), exist_ok=True)
    candidates = [
        log_file,
        os.path.join("Data", f"{EXP_NAME}_#{exp_num}", "05_Plots",
                     f"{EXP_NAME}_#{exp_num}_long_nsr_calibration_fit.png"),
        os.path.join("Data", f"{EXP_NAME}_#{exp_num}", "05_Plots",
                     f"{EXP_NAME}_#{exp_num}_std_nsr_calibration_fit.png"),
    ]
    for mode in consensus_modes:
        sanity_exp = sanity_exp_name(mode, exp_num)
        candidates.append(
            os.path.join("Data", sanity_exp, "05_Plots",
                         f"{sanity_exp}_all_scenarios_global_vs_intrahost.png"))
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in candidates:
            if os.path.isfile(path):
                zf.write(path, arcname=os.path.basename(path))
            else:
                print(f"  [skip] artifact not found: {path}")


def compress_results(exp_num, log_fh, consensus, keep_slurm_logs=True):
    """Pack everything this run produced into Data_Export/*.tar.xz.

    Names the directories outright (run_export_dirs) instead of handing
    package_simplicity_data.sh a `Data/*_#<exp_num>` match, which also picked
    up anything else whose name ended in the same number -- a range-finder
    probe, a scratch experiment -- and silently shipped it. Since 8e the list
    is short: the experiment, setup data, the sanity grid.

    The tag puts the experiment name and number in the archive filename, so two
    runs packed on the same day no longer collide on part1.

    Submitted to Slurm (16 cpus) rather than run here: pixz over several GB
    does not belong on a login node. The pipeline does not wait for it; the
    parts appear in Data_Export/ when the packing job finishes.
    """
    script = os.path.join(SCRIPT_DIR, os.pardir, "package_simplicity_data.sh")
    cmd = [script]
    for d in run_export_dirs(exp_num, consensus):
        cmd += ["--dir", d]
    cmd += ["--tag", f"{EXP_NAME}_#{exp_num}", "--slurm"]
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
        description="Run the full impact_long_shedders pipeline: each "
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
                        choices=list(CONSENSUS_MODES),
                        default=[DEFAULT_CONSENSUS],
                        help=f"Consensus mode for the whole run -- both "
                            f"calibration stages and the experiment phase "
                            f"(default: {DEFAULT_CONSENSUS}). Passing two "
                            f"modes is refused: one calibration cannot serve "
                            f"both, so run them under separate --exp-num.")
    parser.add_argument('--skip-cal1', action='store_true',
                        help=f"Reuse the calib_long__* groups already in "
                            f"{EXP_NAME}_#{{exp_num}}.")
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
    parser.add_argument('--slurm-mem-cal2', type=str, default=None,
                        help="Per-task SLURM memory for cal_2 only. Defaults "
                             "to --slurm-mem: cal_2 is the cheapest stage in "
                             "the pipeline and has no reason to ask for more.")
    args = parser.parse_args()

    # One calibration, one consensus mode. Both calibration stages now state
    # the mode they ran under, so serving two experiment modes from one table
    # would mean calibrating under one and running under the other -- the
    # mismatch v2.4.93 removed. Refuse rather than pick, and refuse BEFORE
    # the log file and the Data/ directories exist.
    if len(set(args.consensus)) > 1:
        raise SystemExit(
            "--consensus takes ONE mode here: both calibration stages are run "
            "under it, and a table calibrated under one mode does not transfer "
            f"to another. Got {sorted(set(args.consensus))}. Run each mode "
            "under its own --exp-num.")
    consensus = args.consensus[0]

    log_file = args.log_file or os.path.join(
        "Data", "pipeline_logs",
        f"{EXP_NAME}_pipeline_#{args.exp_num}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    py = sys.executable
    cal1_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_cal_1.py")
    cal2_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_cal_2.py")
    exp_path = os.path.join(SCRIPT_DIR, "impact_long_shedders_exp.py")

    slurm_res_args = ["--slurm-mem", args.slurm_mem, "--slurm-time", args.slurm_time]

    with open(log_file, "a") as log_fh:
        _log(log_fh, f"\n===== {EXP_NAME} pipeline: exp_num={args.exp_num} "
                    f"@ {datetime.now().isoformat(timespec='seconds')} =====")

        if args.skip_cal1:
            _log(log_fh, f"[skip] cal_1 skipped; reusing the "
                         f"calib_long__* groups of {EXP_NAME}_#{args.exp_num}")
        else:
            run_stage([py, cal1_path,
                      "--exp-num", str(args.exp_num),
                      "--runner", args.runner,
                      "--seeds", str(args.cal_seeds),
                      "--target-osr-long", str(args.target_osr_long),
                      "--R", str(args.r_cal1),
                      "--ih-virus-emergence-rate", str(args.ih_virus_emergence_rate),
                      "--consensus", consensus,
                      *slurm_res_args], log_fh)

        numbered = f"{EXP_NAME}_#{args.exp_num}"
        # read the groups off the experiment rather than rebuilding the names
        # here: the phase that wrote them is the authority on what it wrote
        report_simulation_health(
            numbered, groups_of_phase(dm.get_groups(numbered), 'calib_long'),
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
                  "--consensus", consensus,
                  "--slurm-mem", cal2_memory(args.slurm_mem,
                                             args.slurm_mem_cal2),
                  "--slurm-time", args.slurm_time], log_fh)

        report_simulation_health(
            numbered, [calib_std_group()],
            log_fh, "stage 2 (standard alone)",
            skip_outputs=CALIBRATION_SKIPPED_OUTPUTS)

        for mode in args.consensus:
            run_stage([py, exp_path,
                      "--exp-num", str(args.exp_num),
                      "--runner", args.runner,
                      "--seeds", str(args.exp_seeds),
                      "--consensus", mode,
                      *slurm_res_args], log_fh)

            report_simulation_health(
                f"{exp_name(mode)}_#{args.exp_num}",
                [exp_group(s['name']) for s in SCENARIOS],
                log_fh, f"stage 3 (experiment phase, consensus={mode})")

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
            f"Experiment               : {EXP_NAME}_#{args.exp_num}",
            f"  calibration groups     : calib_long__*, calib_std__*",
            *[f"  experiment groups      : exp__* in "
              f"{exp_name(m)}_#{args.exp_num}" for m in args.consensus],
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

        write_artifacts_archive(archive_path, args.exp_num, log_file,
                                args.consensus)

        if not args.no_compress:
            _log(log_fh, "\n===== Packing results for download =====")
            compress_results(args.exp_num, log_fh, consensus,
                             keep_slurm_logs=not args.drop_slurm_logs)


if __name__ == "__main__":
    main()
