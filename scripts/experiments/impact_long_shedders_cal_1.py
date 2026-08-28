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
    LONG_NSR_EXP_NAME, CAL1_ISOLATED_FIXED_PARAMS, USER_FIXED_PARAMS, NSR_RANGES,
    build_cal1_settings, add_slurm_resource_args, set_slurm_resource_env,
    print_fixed_params,
)
from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def run_isolated_calibration(exp_num, runner, seeds, target_osr_long, R,
                             ih_virus_emergence_rate, model_type='exp',
                             min_seq=30, min_len=100):
    print("\n=========================================================")
    print(f" PHASE 1: ISOLATED LONG-SHEDDER CALIBRATION "
         f"(R={R}, R_long per SCENARIOS, 100% LS, "
         f"k_v={ih_virus_emergence_rate})")
    print("=========================================================\n")

    exp_name = LONG_NSR_EXP_NAME
    ranges = NSR_RANGES['cal1_long_nsr']
    settings_func = build_cal1_settings(
        seeds, ranges, R=R, ih_virus_emergence_rate=ih_virus_emergence_rate)

    varying_params, fixed_params, n_seeds = settings_func()
    print_fixed_params(fixed_params)
    print(f"Seeds per grid point: {n_seeds}")
    print("Per tau_3_long group:")
    for g in varying_params['_scenario_groups']:
        nsr_vals = g['nucleotide_substitution_rate_long']
        print(f"  tau_3_long={g['tau_3_long']:7.2f}  R_long={g['R_long']:7.4f}  "
             f"NSR_long grid: {len(nsr_vals)} pts [{min(nsr_vals):.6f} .. {max(nsr_vals):.6f}]")

    print(f"\nDispatching grid experiment: {exp_name}_#{exp_num} to {runner}...")
    run_experiment_script(runner, exp_num, settings_func, exp_name)
    print("\n[Success] Calibration sweep submitted via native runner.")

    numbered = f"{exp_name}_#{exp_num}"
    print(f"\n[Plot] Fitting + plotting long-NSR calibration for {numbered}...")
    plot_and_fit_long_nsr_calibration(numbered, target_osr_long,
                                      model_type=model_type,
                                      min_seq=min_seq, min_len=min_len)


def main():
    parser = argparse.ArgumentParser(
        description="Isolated long-shedder NSR calibration (Phase 1 of the "
                    "impact_long_shedders pipeline).")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=20,
                        help="Seeds per grid point.")
    parser.add_argument('--target-osr-long', type=float, default=0.00205,
                        help="Target long OSR, used only for the diagnostic "
                            "calibration-fit plot (not for the grid itself).")
    parser.add_argument('--R', type=float, default=CAL1_ISOLATED_FIXED_PARAMS['R'],
                        help="Reference R for this isolated context's "
                            f"standard individuals (default {CAL1_ISOLATED_FIXED_PARAMS['R']}); "
                            "unrelated to R_long, which is set per scenario in "
                            "impact_long_shedders_config.SCENARIOS.")
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'],
                        help="Intra-host lineage-duplication rate (k_v); "
                            f"default {USER_FIXED_PARAMS['IH_virus_emergence_rate']}, "
                            "matching cal_2/exp so Stage 1 is calibrated under "
                            "the same k_v it will actually run with.")
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'],
                        help="Fit model for this stage's own diagnostic plot. "
                            "Purely diagnostic -- cal_2.py recomputes the "
                            "actual frozen-table long NSR independently with "
                            "its own --model, which must be set to match if "
                            "the two are meant to agree.")
    parser.add_argument('--min-seq', type=int, default=30,
                        help="Min sequences to keep a seed, for this stage's "
                            "own diagnostic plot (same default cal_2.py uses "
                            "for its independent recomputation).")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Min simulation length (days) to keep a seed, "
                            "for this stage's own diagnostic plot (same "
                            "default cal_2.py uses for its independent "
                            "recomputation).")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    run_isolated_calibration(args.exp_num, args.runner, args.seeds, args.target_osr_long, args.R,
                             args.ih_virus_emergence_rate, args.model,
                             args.min_seq, args.min_len)


if __name__ == "__main__":
    main()
