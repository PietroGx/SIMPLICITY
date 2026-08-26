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
    TAU_ROUND, DEFAULT_COLORS, LONG_NSR_EXP_NAME, STD_NSR_SWEEP_NAME,
    USER_FIXED_PARAMS, NSR_RANGES, lookup_long_nsr, build_cal2_scenario_groups,
    build_cal2_settings, add_slurm_resource_args,
    set_slurm_resource_env, print_fixed_params,
)
from long_nsr_calibration_plot import read_calibrated_long_nsr

# =============================================================================
# CONFIGURATION
# =============================================================================
EXP_NAME = "impact_long_shedders_calibration"


# =============================================================================
# STAGE 1 -- long-shedder NSR per tau_3_long group
# ----------------------------------------------------------------------------
# Not computed here: cal_1.py already fit it and persisted the result (see
# long_nsr_calibration_plot.py's write_calibrated_long_nsr /
# read_calibrated_long_nsr). This used to be a second, independent
# re-derivation from raw simulation output -- see BACKLOG for how that let
# it silently disagree with cal_1's own diagnostic fit.
# =============================================================================


# =============================================================================
# STAGE 2 -- standard NSR calibrated in-context (mixed population)
# ----------------------------------------------------------------------------
# Scenario-group construction (build_cal2_scenario_groups/build_cal2_settings)
# lives in impact_long_shedders_config.py; this script only reads back what
# it submitted.
# =============================================================================
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
    plt.xscale('log')
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
    parser.add_argument('--seeds', type=int, default=20,
                        help="Seeds per NSR point in the standard sweep.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number (shared with cal_1 and exp).")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--R', type=float, default=USER_FIXED_PARAMS['R'],
                        help="Standard-population R for this stage and exp "
                            f"(default {USER_FIXED_PARAMS['R']}); unrelated to "
                            "R_long, which is set per scenario in "
                            "impact_long_shedders_config.SCENARIOS.")
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'],
                        help="Intra-host lineage-duplication rate (k_v) for "
                            "this stage and exp (default "
                            f"{USER_FIXED_PARAMS['IH_virus_emergence_rate']}). "
                            "Recorded in the frozen table so exp.py's "
                            "production run automatically matches whatever "
                            "value was calibrated under.")
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'],
                        help="Fit model for this stage's OWN standard-NSR "
                            "fit. Does not affect the long-NSR value (read "
                            "back from cal_1's already-computed calibration, "
                            "see read_calibrated_long_nsr) -- to change that "
                            "fit's model, rerun cal_1.py with a different "
                            "--model.")
    parser.add_argument('--min-seq', type=int, default=30,
                        help="Min sequences to keep a seed, for this stage's "
                            "own standard-NSR fit.")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Min simulation length (days) to keep a seed, "
                            "for this stage's own standard-NSR fit.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)

    nsr_ranges = NSR_RANGES['cal2_standard_nsr']

    # This pipeline always runs cal_1 -> cal_2 sequentially under the same
    # exp_num, so the long calibration to read from is fixed, not passed in.
    long_calib_exp = f"{LONG_NSR_EXP_NAME}_#{args.exp_num}"

    setup_dir = f"Data/{EXP_NAME}_setup_data_#{args.exp_num}"
    os.makedirs(setup_dir, exist_ok=True)

    # --- Stage 1: long NSR per tau group, already computed by cal_1.py ---
    long_nsr_by_tau = read_calibrated_long_nsr(long_calib_exp, args.target_osr_long)

    # --- Build one group per scenario and submit ONE combined experiment ---
    scenario_groups, scenarios_frozen = build_cal2_scenario_groups(
        nsr_ranges, long_nsr_by_tau)
    settings_func = build_cal2_settings(scenario_groups, args.seeds, args.R,
                                        args.ih_virus_emergence_rate)

    print(f"\n[Stage 2] Submitting standard-NSR sweep (R={args.R}, "
         f"R_long per SCENARIOS)")
    _, fixed_params, n_seeds = settings_func()
    print_fixed_params(fixed_params)
    print(f"Seeds per grid point: {n_seeds}")
    print("Per scenario:")
    for g, frozen in zip(scenario_groups, scenarios_frozen):
        nsr_vals = g['nucleotide_substitution_rate']
        print(f"  {frozen['scenario_name']:10s} ratio={frozen['long_shedders_ratio']:.3f}  "
             f"tau_3_long={frozen['tau_3_long']:7.2f}  R_long={frozen['R_long']:7.4f}  "
             f"NSR grid: {len(nsr_vals)} pts [{min(nsr_vals):.6f} .. {max(nsr_vals):.6f}]")

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
        long_nsr = (lookup_long_nsr(long_nsr_by_tau, frozen["tau_3_long"], frozen["R_long"])
                   if frozen["is_long"] else None)
        rows.append({
            "scenario_name": name,
            "tau_3_long": frozen["tau_3_long"],
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "R": args.R,
            "IH_virus_emergence_rate": args.ih_virus_emergence_rate,
            "R_long": frozen["R_long"],
            "nucleotide_substitution_rate_long": long_nsr,
            "nucleotide_substitution_rate": std_nsr_by_scenario[name],
            "target_osr_std": args.target_osr_std,
            "target_osr_long": args.target_osr_long,
            "model_type": args.model,
            "long_calib_source": long_calib_exp,
        })

    table_path = os.path.join(setup_dir, "nsr_calibration_table.csv")
    cols = ["scenario_name", "tau_3_long", "long_shedders_ratio", "R",
            "IH_virus_emergence_rate", "R_long",
            "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
            "target_osr_std", "target_osr_long", "model_type", "long_calib_source"]
    pd.DataFrame(rows)[cols].to_csv(table_path, index=False)

    print(f"\n[Success] Frozen calibration table written:\n  {table_path}")
    print("Next: run impact_long_shedders_exp.py against this table.")


if __name__ == "__main__":
    main()
