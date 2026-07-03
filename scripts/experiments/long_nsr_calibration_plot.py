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
# Shared long-shedder NSR calibration fit + plot.
# ----------------------------------------------------------------------------
# Per-tau_3_long-group outlier-filtered regression of observed substitution
# rate (OSR) vs. input nucleotide substitution rate (NSR), inverted at the
# target long OSR, plus a multi-curve calibration figure.
#
# Used by impact_long_shedders_cal_1.py right after its long-NSR calibration
# grid finishes (saved into that experiment's own 05_Plots folder, so it can
# be reviewed at the end of a pipeline run and used to retune
# cal1_long_nsr in the NSR-ranges reference file), and by
# scripts/check_calibration.py as a standalone manual-rerun tool.
# ============================================================================

import os
import sys
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

# This module is imported both as a plain sibling (from
# impact_long_shedders_cal_1.py, run directly from this directory) and as a
# dotted submodule (experiments.long_nsr_calibration_plot, from
# scripts/check_calibration.py, run directly from the parent directory) --
# those two entry points put different directories on sys.path[0]. Ensure
# this file's own directory is always importable so the sibling import below
# resolves in both cases.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from impact_long_shedders_config import SCENARIOS, TAU_ROUND, derive_scenario_params

DEFAULT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]


def _tau_labels_and_colors(sp):
    """Map each corrected tau_3_long value back to its scenario name/color."""
    labels, colors = {}, {}
    color_cycle = iter(DEFAULT_COLORS)
    for scenario in SCENARIOS:
        if scenario["long_shedders_ratio"] <= 0.0:
            continue
        frozen = derive_scenario_params(scenario, sp)
        tau = round(frozen["tau_3_long"], TAU_ROUND)
        if tau not in labels:
            labels[tau] = scenario["name"]
            colors[tau] = next(color_cycle, "black")
    return labels, colors


def plot_and_fit_long_nsr_calibration(experiment_name, target_osr_long,
                                      model_type='exp', min_seq=30, min_len=100):
    """
    Extracts per-seed OSR for every (tau_3_long, NSR) grid point of
    `experiment_name`, outlier-filters per grid point, fits a regressor per
    tau_3_long group, inverts at `target_osr_long`, and saves a multi-curve
    calibration figure into the experiment's 05_Plots folder.

    Returns {tau_3_long: calibrated_nsr}.
    """
    sp = sm.read_standard_parameters_values()
    labels, colors = _tau_labels_and_colors(sp)

    print(f"--- Analyzing long-shedder NSR calibration: {experiment_name} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)

    all_sod_results = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'nucleotide_substitution_rate')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')

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
            all_sod_results.append(sod_df)

    if not all_sod_results:
        raise RuntimeError(
            f"No valid long-calibration data found in {experiment_name}. "
            f"Check if simulations finished and lengths/seqs meet minimums.")

    master_df = pd.concat(all_sod_results, ignore_index=True)
    clean_df = master_df[master_df['is_outlier'] == 0]

    plt.figure(figsize=(10, 6))
    tau_groups = sorted(clean_df['tau_3_long'].unique())
    results = {}

    print(f"\n--- Calibration results (target long OSR = {target_osr_long}) ---")

    for tau in tau_groups:
        group_df = clean_df[clean_df['tau_3_long'] == tau]
        key = round(float(tau), TAU_ROUND)
        color = colors.get(key, 'black')
        label = labels.get(key, f'tau={tau}')

        plt.scatter(group_df['nucleotide_substitution_rate'],
                   group_df['observed_substitution_rate'],
                   color=color, alpha=0.15, s=10)

        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name,
            group_df,
            model_type=model_type,
            parameter_name='nucleotide_substitution_rate',
            experiment_group=label,
        )

        try:
            calibrated_nsr = er.compute_calibrated_parameter(model_type, fit_result, target_osr_long)
            results[tau] = float(calibrated_nsr)
            print(f"{label:<17}: NSR = {calibrated_nsr:.6f}")
        except Exception as e:
            print(f"{label:<17}: Could not compute calibration ({e})")
            continue

        x_vals = np.linspace(group_df['nucleotide_substitution_rate'].min(),
                            group_df['nucleotide_substitution_rate'].max(), 100)
        y_vals = fit_result.eval(x=x_vals)
        plt.plot(x_vals, y_vals, color=color, linewidth=2, label=f"{label} fit")
        plt.plot(calibrated_nsr, target_osr_long, marker='*', markersize=12,
                color=color, markeredgecolor='black')

    plt.axhline(target_osr_long, color='black', linestyle='--', linewidth=1.5,
              label='Target OSR')
    plt.title('Long-Shedder OSR Calibration by Scenario')
    plt.xlabel('Input Nucleotide Substitution Rate (NSR)')
    plt.ylabel('Observed Substitution Rate (OSR)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(experiment_name),
                             f'{experiment_name}_long_nsr_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"\nCalibration plot saved to: {plot_path}")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Standalone rerun of the long-shedder NSR calibration fit/plot.")
    parser.add_argument('--experiment-name', type=str, required=True,
                        help="Numbered long-calibration experiment, e.g. calibrate_long_nsr_#3.")
    parser.add_argument('--target-osr', type=float, default=0.00205)
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    args = parser.parse_args()

    try:
        plot_and_fit_long_nsr_calibration(
            args.experiment_name, args.target_osr, args.model, args.min_seq, args.min_len)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
