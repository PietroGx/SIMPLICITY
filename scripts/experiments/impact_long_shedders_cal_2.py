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
# impact_long_shedders -- CALIBRATION SCRIPT (File 2 of 2)
# ----------------------------------------------------------------------------
# Produces a FROZEN nsr_calibration_table.csv consumed by the experiment
# runner (impact_long_shedders_exp.py). This script performs NO final
# experiment runs; it only calibrates and writes the table.
#
# Per scenario, two independent absolute rates are derived:
#   1. nucleotide_substitution_rate_long  -- from the EXISTING long-shedder
#      calibration experiment (impact_long_shedders_calibration_lng_nsr_#{exp_num},
#      produced by impact_long_shedders_cal_1.py under the SAME exp_num --
#      this pipeline always runs cal_1 -> cal_2 -> exp sequentially under one
#      shared number), per tau_3_long group, exp-fit + inverted at target
#      long OSR.
#   2. nucleotide_substitution_rate (standard) -- calibrated IN-CONTEXT: a
#      mixed-population NSR sweep with the long side fully frozen, fitting OSR
#      measured in STANDARD individuals only, inverted at target standard OSR.
#
# All scenarios' standard-NSR sweeps are submitted as a SINGLE combined
# experiment (one SLURM array job) via the '_scenario_groups' mechanism in
# simplicity.settings_manager.generate_experiment_settings -- each scenario
# keeps its own NSR sweep range (from the NSR-ranges reference file) and its
# own fixed tau_3_long/ratio/R_long/long-NSR. Per-scenario fitting then reads
# that one merged experiment back and groups by (tau_3_long,
# long_shedders_ratio), mirroring how Stage 1 already group-fits per tau.
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
from impact_long_shedders_config import (
    SCENARIOS, TAU_ROUND, USER_FIXED_PARAMS, DEFAULT_COLORS,
    LONG_NSR_EXP_NAME, STD_NSR_SWEEP_NAME, derive_scenario_params, read_nsr_ranges,
    add_slurm_resource_args, set_slurm_resource_env,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
EXP_NAME = "impact_long_shedders_calibration"


# =============================================================================
# STAGE 1 -- long-shedder NSR per tau_3_long group (from existing calibration)
# =============================================================================
def compute_long_nsr_per_tau(long_calib_exp, target_osr_long,
                             model_type='exp', min_seq=30, min_len=100):
    """
    Replicates the per-group extraction from check_calibration.py:
    extract per-seed OSR grouped by tau_3_long, outlier-filter per grid point,
    fit an 'exp' regressor per group, invert at target_osr_long.

    Returns {round(tau_3_long, TAU_ROUND): calibrated_long_nsr}.
    """
    print(f"\n[Stage 1] Long-shedder NSR calibration from: {long_calib_exp}")
    print(f"          target long OSR = {target_osr_long}")

    simulation_output_dirs = dm.get_simulation_output_dirs(long_calib_exp)

    all_rows = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'tau_3_long')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                if final_time >= min_len and len(seq_data) >= min_seq:
                    osr_val = er.tempest_regression(seq_data).coef_[0]
                    sod_rows.append({
                        'tau_3_long': tau_val,
                        'nucleotide_substitution_rate': nsr_val,
                        'observed_substitution_rate': osr_val,
                    })
            except Exception:
                continue

        if sod_rows:
            sod_df = pd.DataFrame(sod_rows)
            sod_df = om.detect_sod_outliers(sod_df)
            all_rows.append(sod_df)

    if not all_rows:
        raise RuntimeError(
            f"[Stage 1] No valid long-calibration data found in {long_calib_exp}.")

    master_df = pd.concat(all_rows, ignore_index=True)
    clean_df = master_df[master_df['is_outlier'] == 0]

    long_nsr_by_tau = {}
    for tau in sorted(clean_df['tau_3_long'].unique()):
        group_df = clean_df[clean_df['tau_3_long'] == tau]
        fit_result = er.fit_observed_substitution_rate_regressor(
            long_calib_exp,
            group_df,
            model_type=model_type,
            parameter_name='nucleotide_substitution_rate',
            experiment_group=f"tau_{tau}",
        )
        nsr = er.compute_calibrated_parameter(model_type, fit_result, target_osr_long)
        long_nsr_by_tau[round(float(tau), TAU_ROUND)] = float(nsr)
        print(f"          tau_3_long={tau}: long NSR = {nsr:.8f}")

    return long_nsr_by_tau


def lookup_long_nsr(long_nsr_by_tau, tau_3_long):
    """Exact (rounded) match lookup. Raises on miss -- never silently guesses."""
    key = round(float(tau_3_long), TAU_ROUND)
    if key not in long_nsr_by_tau:
        raise KeyError(
            f"No calibrated long NSR for tau_3_long={key}. "
            f"Available keys: {sorted(long_nsr_by_tau.keys())}. "
            f"Check that the scenario tau matches the calibration experiment.")
    return long_nsr_by_tau[key]


# =============================================================================
# STAGE 2 -- standard NSR calibrated in-context (mixed population)
# =============================================================================
def build_scenario_groups(nsr_ranges, long_nsr_by_tau, sp):
    """
    One group dict per scenario: its own nucleotide_substitution_rate sweep
    (list -> the group's local varying axis) plus its fixed tau_3_long/
    long_shedders_ratio/R_long/sequence_long_shedders/long-NSR overrides.
    """
    scenario_groups = []
    scenarios_frozen = []

    for scenario in SCENARIOS:
        frozen = derive_scenario_params(scenario, sp)
        scenarios_frozen.append(frozen)

        name = frozen["scenario_name"]
        r = nsr_ranges[name]
        nsr_list = np.logspace(np.log10(r['min']), np.log10(r['max']), r['steps']).tolist()

        group = {
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "tau_3_long": frozen["tau_3_long"],
            "R_long": frozen["R_long"],
            "sequence_long_shedders": frozen["sequence_long_shedders"],
            "nucleotide_substitution_rate": [float(x) for x in nsr_list],
        }
        if frozen["is_long"]:
            group["nucleotide_substitution_rate_long"] = lookup_long_nsr(
                long_nsr_by_tau, frozen["tau_3_long"])

        scenario_groups.append(group)

    return scenario_groups, scenarios_frozen


def compute_standard_nsr_per_scenario(numbered, scenarios_frozen, target_osr_std,
                                      model_type='exp', min_seq=30, min_len=100):
    """
    Reads back the ONE merged standard-NSR sweep experiment, groups rows by
    matching (tau_3_long, long_shedders_ratio) back to each scenario, and
    fits+inverts per scenario -- mirroring Stage 1's per-tau pattern. OSR is
    computed from STANDARD individuals only (mixed population).

    Returns {scenario_name: calibrated_standard_nsr}.
    """
    print(f"\n[Stage 2] Standard NSR in-context sweep: {numbered}")

    simulation_output_dirs = dm.get_simulation_output_dirs(numbered)

    all_rows = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'tau_3_long')
        ratio_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'long_shedders_ratio')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                seq_data = seq_data[seq_data['individual_type'] == 'standard']
                if final_time >= min_len and len(seq_data) >= min_seq:
                    osr_val = er.tempest_regression(seq_data).coef_[0]
                    sod_rows.append({
                        'tau_3_long': tau_val,
                        'long_shedders_ratio': ratio_val,
                        'nucleotide_substitution_rate': nsr_val,
                        'observed_substitution_rate': osr_val,
                    })
            except Exception:
                continue

        if sod_rows:
            sod_df = pd.DataFrame(sod_rows)
            sod_df = om.detect_sod_outliers(sod_df)
            all_rows.append(sod_df)

    if not all_rows:
        raise RuntimeError(
            f"[Stage 2] No standard-individual OSR data for {numbered}. "
            f"Sweep may have produced too few standard sequences "
            f"(min_seq={min_seq}, min_len={min_len}).")

    master_df = pd.concat(all_rows, ignore_index=True)
    clean_df = master_df[master_df['is_outlier'] == 0]

    plt.figure(figsize=(10, 6))
    color_by_name = {frozen["scenario_name"]: DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                     for i, frozen in enumerate(scenarios_frozen)}

    std_nsr_by_scenario = {}
    for frozen in scenarios_frozen:
        name = frozen["scenario_name"]
        color = color_by_name[name]
        key_tau = round(frozen["tau_3_long"], TAU_ROUND)
        key_ratio = round(frozen["long_shedders_ratio"], TAU_ROUND)

        group_df = clean_df[
            (clean_df['tau_3_long'].round(TAU_ROUND) == key_tau) &
            (clean_df['long_shedders_ratio'].round(TAU_ROUND) == key_ratio)
        ]
        if group_df.empty:
            raise RuntimeError(
                f"[Stage 2] No standard-individual OSR rows matched scenario "
                f"'{name}' (tau_3_long={key_tau}, long_shedders_ratio={key_ratio}) "
                f"in {numbered}.")

        plt.scatter(group_df['nucleotide_substitution_rate'],
                   group_df['observed_substitution_rate'],
                   color=color, alpha=0.15, s=10)

        fit_result = er.fit_observed_substitution_rate_regressor(
            numbered, group_df, model_type,
            parameter_name='nucleotide_substitution_rate',
            experiment_group=name)

        nsr_std = er.compute_calibrated_parameter(model_type, fit_result, target_osr_std)
        std_nsr_by_scenario[name] = float(nsr_std)
        print(f"          {name}: calibrated standard NSR = {nsr_std:.8f}")

        x_vals = np.linspace(group_df['nucleotide_substitution_rate'].min(),
                            group_df['nucleotide_substitution_rate'].max(), 100)
        y_vals = fit_result.eval(x=x_vals)
        plt.plot(x_vals, y_vals, color=color, linewidth=2, label=f"{name} fit")
        plt.plot(nsr_std, target_osr_std, marker='*', markersize=12,
                color=color, markeredgecolor='black')

    plt.axhline(target_osr_std, color='black', linestyle='--', linewidth=1.5,
              label='Target OSR')
    plt.title('Standard NSR Calibration by Scenario (in-context, mixed population)')
    plt.xlabel('Input Nucleotide Substitution Rate (NSR)')
    plt.ylabel('Observed Substitution Rate (OSR, standard individuals)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f'{numbered}_std_nsr_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"\nCalibration plot saved to: {plot_path}")

    return std_nsr_by_scenario


# =============================================================================
# MAIN
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Calibrate impact_long_shedders scenarios -> frozen NSR table.")
    parser.add_argument('--target-osr-std', type=float, default=0.0013,
                        help="Target OSR for standard individuals (default 0.0013).")
    parser.add_argument('--target-osr-long', type=float, required=True,
                        help="Target OSR for long shedders.")
    parser.add_argument('--seeds', type=int, required=True,
                        help="Seeds per NSR point in the standard sweep.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number (shared with cal_1 and exp).")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'], help="Fit model.")
    parser.add_argument('--min-seq', type=int, default=30,
                        help="Min sequences to keep a seed.")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Min simulation length (days) to keep a seed.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)

    sp = sm.read_standard_parameters_values()
    nsr_ranges = read_nsr_ranges()['cal2_standard_nsr']

    # This pipeline always runs cal_1 -> cal_2 sequentially under the same
    # exp_num, so the long calibration to read from is fixed, not passed in.
    long_calib_exp = f"{LONG_NSR_EXP_NAME}_#{args.exp_num}"

    setup_dir = f"Data/{EXP_NAME}_setup_data_#{args.exp_num}"
    os.makedirs(setup_dir, exist_ok=True)

    # --- Stage 1: long NSR per tau group (once, reused across scenarios) ---
    long_nsr_by_tau = compute_long_nsr_per_tau(
        long_calib_exp, args.target_osr_long,
        model_type=args.model, min_seq=args.min_seq, min_len=args.min_len)

    # --- Build one group per scenario and submit ONE combined experiment ---
    scenario_groups, scenarios_frozen = build_scenario_groups(nsr_ranges, long_nsr_by_tau, sp)

    def settings_func():
        varying_params = {'_scenario_groups': scenario_groups}
        return (varying_params, USER_FIXED_PARAMS.copy(), args.seeds)

    run_experiment_script(args.runner, args.exp_num, settings_func, STD_NSR_SWEEP_NAME)
    numbered = f"{STD_NSR_SWEEP_NAME}_#{args.exp_num}"

    # --- Per-scenario fitting (local, no further submissions) ---
    std_nsr_by_scenario = compute_standard_nsr_per_scenario(
        numbered, scenarios_frozen, args.target_osr_std,
        model_type=args.model, min_seq=args.min_seq, min_len=args.min_len)

    # --- Write frozen table ---
    rows = []
    for frozen in scenarios_frozen:
        name = frozen["scenario_name"]
        long_nsr = (lookup_long_nsr(long_nsr_by_tau, frozen["tau_3_long"])
                   if frozen["is_long"] else None)
        rows.append({
            "scenario_name": name,
            "tau_3_long": frozen["tau_3_long"],
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "R_long": frozen["R_long"],
            "nucleotide_substitution_rate_long": long_nsr,
            "nucleotide_substitution_rate": std_nsr_by_scenario[name],
            "target_osr_std": args.target_osr_std,
            "target_osr_long": args.target_osr_long,
            "model_type": args.model,
            "long_calib_source": long_calib_exp,
        })

    table_path = os.path.join(setup_dir, "nsr_calibration_table.csv")
    cols = ["scenario_name", "tau_3_long", "long_shedders_ratio", "R_long",
            "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
            "target_osr_std", "target_osr_long", "model_type", "long_calib_source"]
    pd.DataFrame(rows)[cols].to_csv(table_path, index=False)

    print(f"\n[Success] Frozen calibration table written:\n  {table_path}")
    print("Next: run impact_long_shedders_exp.py against this table.")


if __name__ == "__main__":
    main()
