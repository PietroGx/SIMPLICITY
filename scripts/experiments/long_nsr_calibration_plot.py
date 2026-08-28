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
# grid finishes (plot saved into that experiment's own 05_Plots folder, so it
# can be reviewed at the end of a pipeline run and used to retune
# cal1_long_nsr in impact_long_shedders_config.py), and by
# check_calibration.py as a standalone manual-rerun tool.
#
# The computed {(tau_3_long, R_long): calibrated_nsr} result is ALSO
# persisted to a small CSV (write_calibrated_long_nsr, in the experiment's
# Fit_results dir) -- the single source of truth impact_long_shedders_cal_2.py
# reads back (read_calibrated_long_nsr) instead of independently re-deriving
# the same fit from raw simulation output a second time. Re-running this
# module against a past experiment (e.g. via check_calibration.py with a new
# --target-osr) overwrites that file too, so cal_2's next run picks up the
# new calibration -- not a side-effect-free plot refresh anymore.
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

# Imported as a plain sibling by both impact_long_shedders_cal_1.py and
# check_calibration.py (both live in this same directory). Ensure this
# file's own directory is on sys.path regardless of the caller's cwd, so the
# sibling import below always resolves.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from impact_long_shedders_config import (
    SCENARIOS, TAU_ROUND, DEFAULT_COLORS, derive_tau_3_long,
)


def _tau_r_long_labels_and_colors(sp):
    """Map each (tau_3_long, R_long) pair back to a label and color. When
    multiple scenarios share a (tau, R_long) pair (e.g. HIV_low/HIV_high --
    same duration and R_long, only ratio differs, which Stage 1's isolated
    context doesn't use), the label combines every scenario name sharing it,
    joined with '+' -- not '/', which write_fit_results_csv's file-saving
    reads as a directory separator."""
    names_by_key = {}
    for scenario in SCENARIOS:
        if scenario["long_shedders_ratio"] <= 0.0:
            continue
        key = (round(derive_tau_3_long(scenario, sp), TAU_ROUND),
               round(float(scenario["R_long"]), TAU_ROUND))
        names_by_key.setdefault(key, []).append(scenario["name"])

    labels, colors = {}, {}
    color_cycle = iter(DEFAULT_COLORS)
    for key, names in names_by_key.items():
        labels[key] = "+".join(names)
        colors[key] = next(color_cycle, "black")
    return labels, colors


def _calibrated_long_nsr_path(experiment_name):
    return os.path.join(dm.get_experiment_fit_result_dir(experiment_name),
                        f'{experiment_name}_calibrated_long_nsr.csv')


def write_calibrated_long_nsr(experiment_name, results, target_osr_long, model_type):
    """Persist Stage 1's calibrated long NSR per (tau_3_long, R_long) group
    to a small CSV -- the single computed result, so cal_2.py can read it
    straight back instead of independently re-extracting/re-fitting from
    raw simulation output a second time (two implementations of the same
    fit is exactly what let cal_1's and cal_2's methodology drift apart
    before; see BACKLOG)."""
    rows = [
        {'tau_3_long': tau, 'R_long': r_long,
         'nucleotide_substitution_rate_long': nsr,
         'target_osr_long': target_osr_long, 'model_type': model_type}
        for (tau, r_long), nsr in results.items()
    ]
    path = _calibrated_long_nsr_path(experiment_name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f"Saved calibrated long NSR per group to: {path}")


def read_calibrated_long_nsr(experiment_name, target_osr_long):
    """Read back Stage 1's calibrated long NSR per (tau_3_long, R_long)
    group, written by write_calibrated_long_nsr. Raises rather than
    silently guessing if the file is missing (cal_1 never ran under this
    exp_num) or if target_osr_long doesn't match what the file was
    actually calibrated against -- a caller wanting a different target
    must re-run cal_1.py with it, not just pass a different flag here and
    reuse a stale calibration.

    Returns {(tau_3_long, R_long): calibrated_nsr}.
    """
    path = _calibrated_long_nsr_path(experiment_name)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"No calibrated long NSR file at {path}. Run "
            f"impact_long_shedders_cal_1.py for this exp_num first.")

    df = pd.read_csv(path)
    mismatched = df[df['target_osr_long'] != target_osr_long]
    if not mismatched.empty:
        raise ValueError(
            f"{path} was calibrated against target_osr_long="
            f"{mismatched['target_osr_long'].iloc[0]}, but this run requested "
            f"target_osr_long={target_osr_long}. Re-run cal_1.py with "
            f"--target-osr-long {target_osr_long} first, or pass the matching "
            f"target here.")

    return {
        (round(float(row['tau_3_long']), TAU_ROUND),
         round(float(row['R_long']), TAU_ROUND)): float(row['nucleotide_substitution_rate_long'])
        for _, row in df.iterrows()
    }


def plot_and_fit_long_nsr_calibration(experiment_name, target_osr_long,
                                      model_type='exp', min_seq=30, min_len=100):
    """
    Extracts per-seed OSR for every (tau_3_long, R_long, NSR) grid point of
    `experiment_name`, outlier-filters per grid point, fits a regressor per
    (tau_3_long, R_long) group -- two scenarios can share a duration but
    have different R_long, so grouping by tau alone would conflate them --
    inverts at `target_osr_long`, and saves a multi-curve calibration
    figure into the experiment's 05_Plots folder.

    Returns {(tau_3_long, R_long): calibrated_nsr}.
    """
    sp = sm.read_standard_parameters_values()
    labels, colors = _tau_r_long_labels_and_colors(sp)

    print(f"--- Analyzing long-shedder NSR calibration: {experiment_name} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)

    all_sod_results = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate_long')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
        r_long_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'R_long')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                # Same final_time/sequence-count gate as everywhere else that
                # fits an OSR regressor -- decides whether this seed ran long
                # enough / sequenced enough to trust at all. What actually
                # gets fit below is NOT seq_data (population-pooled, global-
                # reference, absolute-time -- see evolutionary_rate.py's
                # extract_ih_regression_data docstring): cal_1 is meant to
                # calibrate against each long shedder's OWN intra-host clock,
                # not the population root-to-tip rate.
                if final_time >= min_len and len(seq_data) >= min_seq:
                    ih_points = er.extract_ih_regression_data(ssod)
                    if ih_points.empty:
                        continue
                    osr_val = er.tempest_regression(ih_points).coef_[0]
                    sod_rows.append({
                        'tau_3_long': tau_val,
                        'R_long': r_long_val,
                        'nucleotide_substitution_rate_long': nsr_val,
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
    clean_df['_group_key'] = list(zip(clean_df['tau_3_long'], clean_df['R_long']))
    tau_r_long_groups = sorted(clean_df['_group_key'].unique())
    results = {}

    print(f"\n--- Calibration results (target long OSR = {target_osr_long}) ---")

    for tau, r_long in tau_r_long_groups:
        group_df = clean_df[clean_df['_group_key'] == (tau, r_long)]
        key = (round(float(tau), TAU_ROUND), round(float(r_long), TAU_ROUND))
        color = colors.get(key, 'black')
        label = labels.get(key, f'tau={tau},R_long={r_long}')

        plt.scatter(group_df['nucleotide_substitution_rate_long'],
                   group_df['observed_substitution_rate'],
                   color=color, alpha=0.15, s=10)

        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name,
            group_df,
            model_type=model_type,
            parameter_name='nucleotide_substitution_rate_long',
            experiment_group=label,
        )

        try:
            calibrated_nsr = er.compute_calibrated_parameter(model_type, fit_result, target_osr_long)
            results[key] = float(calibrated_nsr)
            print(f"{label:<17}: NSR = {calibrated_nsr:.6f}")
        except Exception as e:
            print(f"{label:<17}: Could not compute calibration ({e})")
            continue

        x_vals = np.linspace(group_df['nucleotide_substitution_rate_long'].min(),
                            group_df['nucleotide_substitution_rate_long'].max(), 100)
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

    write_calibrated_long_nsr(experiment_name, results, target_osr_long, model_type)

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Standalone rerun of the long-shedder NSR calibration fit/plot.")
    parser.add_argument('--experiment-name', type=str, required=True,
                        help="Numbered long-calibration experiment, e.g. "
                            "impact_long_shedders_calibration_lng_nsr_#3.")
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
