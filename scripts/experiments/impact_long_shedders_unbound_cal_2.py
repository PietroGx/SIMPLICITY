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
# impact_long_shedders_unbound -- STAGE 2: STANDARD INDIVIDUALS ALONE
# ----------------------------------------------------------------------------
# A 0%-long-shedder population at production settings: sweeps the standard
# NSR, fits the global root-to-tip clock, inverts at target_osr_std. ONE
# calibrated standard rate, shared by every scenario.
#
# This is the whole point of the unbound pipeline. The bound one fits this
# rate per scenario INSIDE the mixed population, so the rate absorbs whatever
# divergence long shedders inject -- in run #7 that drove HIV_high to 9.3e-7,
# 150x below control. Here the standard cohort never sees a long shedder while
# being calibrated, so the rate means "the rate standard individuals actually
# mutate at" and the population clock is left free to rise in production.
#
# The sweep is independent of stage 1. Stage 1's result is only read at the
# end, to assemble the frozen table both stages feed to production.
# ============================================================================

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

from experiment_script_runner import run_experiment_script
from impact_long_shedders_unbound_config import (
    SCENARIOS, LONG_NSR_EXP_NAME, STD_NSR_EXP_NAME, NSR_RANGES,
    SETUP_DIR_TEMPLATE, TABLE_FILENAME, TABLE_COLUMNS,
    UNBOUND_CAL2_FINAL_TIME, USER_FIXED_PARAMS,
    derive_scenario_params, lookup_long_nsr, build_cal2_settings,
    add_slurm_resource_args, set_slurm_resource_env, print_fixed_params,
)
from long_nsr_calibration_plot import read_calibrated_long_nsr


def compute_standard_nsr(numbered, target_osr_std, model_type='exp',
                         min_seq=30, min_len=100):
    """Read back the standard-only sweep, fit OSR vs NSR over standard
    individuals, invert at target. Returns the calibrated standard NSR."""
    print(f"\n[Stage 2] Standard-only NSR sweep: {numbered}")

    all_rows = []
    for sod in dm.get_simulation_output_dirs(numbered):
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')

        sod_rows = []
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                # Every individual is standard here; the filter is kept so
                # this measures the same quantity as the bound pipeline's.
                seq_data = seq_data[seq_data['individual_type'] == 'standard']
                if final_time >= min_len and len(seq_data) >= min_seq:
                    sod_rows.append({
                        'nucleotide_substitution_rate': nsr_val,
                        'observed_substitution_rate':
                            er.tempest_regression(seq_data).coef_[0],
                    })
            except Exception:
                continue

        if sod_rows:
            all_rows.append(om.detect_sod_outliers(pd.DataFrame(sod_rows)))

    if not all_rows:
        raise RuntimeError(
            f"[Stage 2] No standard OSR data for {numbered} "
            f"(min_seq={min_seq}, min_len={min_len}).")

    clean_df = pd.concat(all_rows, ignore_index=True)
    clean_df = clean_df[clean_df['is_outlier'] == 0]

    fit_result = er.fit_observed_substitution_rate_regressor(
        numbered, clean_df, model_type,
        parameter_name='nucleotide_substitution_rate',
        experiment_group='standard_only')
    nsr_std = float(er.compute_calibrated_parameter(
        model_type, fit_result, target_osr_std))
    print(f"          calibrated standard NSR = {nsr_std:.8f}")

    plt.figure(figsize=(10, 6))
    plt.scatter(clean_df['nucleotide_substitution_rate'],
                clean_df['observed_substitution_rate'],
                color='#1f77b4', alpha=0.15, s=10)
    x_vals = np.geomspace(clean_df['nucleotide_substitution_rate'].min(),
                          clean_df['nucleotide_substitution_rate'].max(), 100)
    plt.plot(x_vals, fit_result.eval(x=x_vals), color='#1f77b4', linewidth=2,
             label='standard-only fit')
    plt.plot(nsr_std, target_osr_std, marker='*', markersize=14,
             color='#1f77b4', markeredgecolor='black')
    plt.axhline(target_osr_std, color='black', linestyle='--', linewidth=1.5,
                label='Target OSR')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Standard NSR calibration (standard individuals alone)')
    plt.xlabel('Input Nucleotide Substitution Rate (NSR)')
    plt.ylabel('Observed Substitution Rate (OSR, standard individuals)')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f'{numbered}_std_nsr_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Calibration plot saved to: {plot_path}")

    return nsr_std


def write_frozen_table(exp_num, std_nsr, long_nsr_by_group, args,
                       long_calib_exp, std_calib_exp):
    sp = sm.read_standard_parameters_values()
    setup_dir = SETUP_DIR_TEMPLATE.format(exp_num=exp_num)
    os.makedirs(setup_dir, exist_ok=True)

    rows = []
    for scenario in SCENARIOS:
        frozen = derive_scenario_params(scenario, sp)
        long_nsr = (lookup_long_nsr(long_nsr_by_group, frozen["tau_3_long"],
                                    frozen["R_long"])
                    if frozen["is_long"] else None)
        rows.append({
            "scenario_name": frozen["scenario_name"],
            "tau_3_long": frozen["tau_3_long"],
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "susceptibility_long": frozen["susceptibility_long"],
            "R": args.R,
            "IH_virus_emergence_rate": args.ih_virus_emergence_rate,
            "R_long": frozen["R_long"],
            "nucleotide_substitution_rate_long": long_nsr,
            # one standard rate for every scenario -- nothing scenario-specific
            # is fitted once calibration stops compensating for the long side
            "nucleotide_substitution_rate": std_nsr,
            "target_osr_long": args.target_osr_long,
            "model_type": args.model,
            "long_calib_source": long_calib_exp,
            "std_calib_source": std_calib_exp,
        })

    table_path = os.path.join(setup_dir, TABLE_FILENAME)
    pd.DataFrame(rows)[TABLE_COLUMNS].to_csv(table_path, index=False)
    print(f"\n[Success] Frozen calibration table written:\n  {table_path}")
    return table_path


def main():
    parser = argparse.ArgumentParser(
        description="Unbound pipeline stage 2: calibrate the standard NSR in a "
                    "0%-long-shedder population, then freeze the table.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=30)
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--R', type=float, default=USER_FIXED_PARAMS['R'])
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'])
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)

    settings_func = build_cal2_settings(
        args.seeds, NSR_RANGES['cal_std'], args.R, args.ih_virus_emergence_rate)

    print(f"\n[Stage 2] Standard-only sweep (R={args.R}, no long shedders, "
         f"final_time={UNBOUND_CAL2_FINAL_TIME})")
    varying, fixed_params, n_seeds = settings_func()
    print_fixed_params(fixed_params)
    nsr_vals = varying['nucleotide_substitution_rate']
    print(f"Seeds per grid point: {n_seeds}")
    print(f"NSR grid: {len(nsr_vals)} pts "
         f"[{min(nsr_vals):.6f} .. {max(nsr_vals):.6f}]")

    run_experiment_script(args.runner, args.exp_num, settings_func, STD_NSR_EXP_NAME)
    std_calib_exp = f"{STD_NSR_EXP_NAME}_#{args.exp_num}"

    std_nsr = compute_standard_nsr(std_calib_exp, args.target_osr_std,
                                   model_type=args.model,
                                   min_seq=args.min_seq, min_len=args.min_len)

    long_calib_exp = f"{LONG_NSR_EXP_NAME}_#{args.exp_num}"
    long_nsr_by_group = read_calibrated_long_nsr(long_calib_exp, args.target_osr_long)

    write_frozen_table(args.exp_num, std_nsr, long_nsr_by_group, args,
                       long_calib_exp, std_calib_exp)
    print("Next: run impact_long_shedders_unbound_exp.py against this table.")


if __name__ == "__main__":
    main()
