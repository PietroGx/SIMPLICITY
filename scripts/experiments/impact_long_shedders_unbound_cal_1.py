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
# impact_long_shedders_unbound -- STAGE 1: LONG SHEDDERS ALONE
# ----------------------------------------------------------------------------
# 100% long-shedder population, sweeps NSR_long, fits each seed's intra-host
# clock and inverts at target_osr_long. Same measurement as the bound
# pipeline's cal_1 -- it calls the same plot_and_fit_long_nsr_calibration, so
# the two can never disagree about what the intra-host clock is.
#
# Independent of stage 2: nothing here reads or waits for the standard side.
#
# Runs things; builds nothing. Parameter construction lives in
# impact_long_shedders_unbound_config.build_cal1_settings.
# ============================================================================

import argparse

from experiment_script_runner import run_experiment_script
from impact_long_shedders_unbound_config import (
    LONG_NSR_EXP_NAME, NSR_RANGES, UNBOUND_CAL1_FINAL_TIME,
    CAL1_ISOLATED_FIXED_PARAMS, USER_FIXED_PARAMS,
    build_cal1_settings, add_slurm_resource_args, set_slurm_resource_env,
    print_fixed_params,
)
from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def run_long_calibration(exp_num, runner, seeds, target_osr_long, R,
                         ih_virus_emergence_rate, model_type='exp',
                         min_seq=30, min_len=100):
    print("\n=========================================================")
    print(f" UNBOUND STAGE 1: LONG SHEDDERS ALONE "
         f"(R={R}, R_long per SCENARIOS, 100% LS, "
         f"k_v={ih_virus_emergence_rate}, final_time={UNBOUND_CAL1_FINAL_TIME})")
    print("=========================================================\n")

    ranges = NSR_RANGES['cal_long']
    settings_func = build_cal1_settings(
        seeds, ranges, R=R, ih_virus_emergence_rate=ih_virus_emergence_rate)

    varying_params, fixed_params, n_seeds = settings_func()
    print_fixed_params(fixed_params)
    print(f"Seeds per grid point: {n_seeds}")
    print("Per (tau_3_long, R_long) group:")
    for g in varying_params['_scenario_groups']:
        nsr_vals = g['nucleotide_substitution_rate_long']
        print(f"  tau_3_long={g['tau_3_long']:7.2f}  R_long={g['R_long']:7.4f}  "
             f"NSR_long grid: {len(nsr_vals)} pts "
             f"[{min(nsr_vals):.6f} .. {max(nsr_vals):.6f}]")

    print(f"\nDispatching {LONG_NSR_EXP_NAME}_#{exp_num} to {runner}...")
    run_experiment_script(runner, exp_num, settings_func, LONG_NSR_EXP_NAME)

    numbered = f"{LONG_NSR_EXP_NAME}_#{exp_num}"
    print(f"\n[Plot] Fitting + plotting long-NSR calibration for {numbered}...")
    plot_and_fit_long_nsr_calibration(numbered, target_osr_long,
                                      model_type=model_type,
                                      min_seq=min_seq, min_len=min_len)


def main():
    parser = argparse.ArgumentParser(
        description="Unbound pipeline stage 1: calibrate NSR_long in a 100% "
                    "long-shedder population.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=30)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--R', type=float, default=CAL1_ISOLATED_FIXED_PARAMS['R'])
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'])
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    run_long_calibration(args.exp_num, args.runner, args.seeds,
                         args.target_osr_long, args.R,
                         args.ih_virus_emergence_rate, args.model,
                         args.min_seq, args.min_len)


if __name__ == "__main__":
    main()
