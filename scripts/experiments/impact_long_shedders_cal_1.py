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
# impact_long_shedders -- ISOLATED LONG-SHEDDER NSR CALIBRATION (Stage 1)
# ----------------------------------------------------------------------------
# Runs things; builds nothing. All parameter construction lives in
# impact_long_shedders_config.build_cal1_settings.
# ============================================================================

import argparse

from experiment_script_runner import run_experiment_script
from impact_long_shedders_config import (
    LONG_NSR_EXP_NAME, CAL1_R_LONG_WEEKLY_RATE, build_cal1_settings,
    read_nsr_ranges, add_slurm_resource_args, set_slurm_resource_env,
)
from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def run_isolated_calibration(exp_num, runner, seeds, target_osr_long):
    print("\n=========================================================")
    print(f" PHASE 1: ISOLATED LONG-SHEDDER CALIBRATION "
         f"(R_long = {CAL1_R_LONG_WEEKLY_RATE} * tau_3_long/7, 100% LS)")
    print("=========================================================\n")

    exp_name = LONG_NSR_EXP_NAME
    ranges = read_nsr_ranges()['cal1_long_nsr']
    settings_func = build_cal1_settings(seeds, ranges)

    print(f"Dispatching grid experiment: {exp_name}_#{exp_num} to {runner}...")
    run_experiment_script(runner, exp_num, settings_func, exp_name)
    print("\n[Success] Calibration sweep submitted via native runner.")

    numbered = f"{exp_name}_#{exp_num}"
    print(f"\n[Plot] Fitting + plotting long-NSR calibration for {numbered}...")
    plot_and_fit_long_nsr_calibration(numbered, target_osr_long)


def main():
    parser = argparse.ArgumentParser(
        description="Isolated long-shedder NSR calibration (Phase 1 of the "
                    "impact_long_shedders pipeline).")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=25,
                        help="Seeds per grid point.")
    parser.add_argument('--target-osr-long', type=float, default=0.00205,
                        help="Target long OSR, used only for the diagnostic "
                            "calibration-fit plot (not for the grid itself).")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    run_isolated_calibration(args.exp_num, args.runner, args.seeds, args.target_osr_long)


if __name__ == "__main__":
    main()
