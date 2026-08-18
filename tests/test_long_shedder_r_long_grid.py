# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ============================================================================
# HPC test: isolated long-shedder NSR calibration over a grid of
# R_long x infected_individuals_at_start values, for a single tau_3_long
# (borrowed from one of impact_long_shedders' 5 production scenarios), at
# final_time = 3 years. Not part of the impact_long_shedders pipeline itself
# -- a standalone exploratory tool to visually compare calibration curves
# across R_long/starting-cohort-size before picking defaults for
# scripts/experiments/impact_long_shedders_config.py.
#
# Reuses NSR_RANGES['cal1_long_nsr'] from impact_long_shedders_config as the
# NSR sweep range (same range the real pipeline uses), and
# CAL1_ISOLATED_FIXED_PARAMS['R'] as the standard-population reference R.
# population_size is fixed (not swept).
# ============================================================================

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

from simplicity.runme import run_experiment
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm
import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'scripts', 'experiments'))
from impact_long_shedders_config import (
    SCENARIOS, derive_tau_3_long, CAL1_ISOLATED_FIXED_PARAMS, NSR_RANGES,
)

EXP_NAME = "long_shedder_r_long_grid_test"

R_LONG_VALUES = [1.0, 1.5, 2.0, 4.0]
INFECTED_START_VALUES = [10, 50, 100]
POPULATION_SIZE = 1000        # fixed -- not swept
FINAL_TIME = 365 * 3          # 3 years
TARGET_OSR_LONG = 0.00205     # reference line only (same default as cal_1)

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}


def build_settings(tau_3_long, n_seeds):
    """One scenario group per (R_long, infected_individuals_at_start) cell,
    each sweeping the SAME nucleotide_substitution_rate grid (reusing
    impact_long_shedders_config.NSR_RANGES['cal1_long_nsr']). tau_3_long is
    shared across every cell -- only R_long and the starting cohort size
    vary between groups."""
    ranges = NSR_RANGES['cal1_long_nsr']
    nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()

    scenario_groups = [
        {
            'nucleotide_substitution_rate': nsr_values,
            'R_long': r_long,
            'infected_individuals_at_start': infected_start,
        }
        for r_long in R_LONG_VALUES
        for infected_start in INFECTED_START_VALUES
    ]

    fixed_params = {
        'long_shedders_ratio': 1.0,
        'R': CAL1_ISOLATED_FIXED_PARAMS['R'],
        'population_size': POPULATION_SIZE,
        'tau_3_long': tau_3_long,
        'final_time': FINAL_TIME,
        'sequence_long_shedders': True,
    }

    def make_settings():
        return ({'_scenario_groups': scenario_groups}, fixed_params, n_seeds)

    return make_settings


def analyze_and_plot(numbered, tau_3_long, min_seq=30, min_len=100, model_type='exp'):
    """Extracts per-seed OSR for every (R_long, infected_start, NSR) point,
    outlier-filters per grid point, fits a regressor per (R_long,
    infected_start) cell, and saves ONE combined figure: R_long as rows
    (ascending), starting-infected count as columns (ascending), shared x/y
    axes so every cell is directly comparable at a glance.
    """
    print(f"--- Analyzing R_long/starting-infected grid: {numbered} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(numbered)

    all_rows = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')
        r_long_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'R_long')
        infected_start_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'infected_individuals_at_start')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        final_times = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                final_times.append(final_time)
                seq_data = om.read_sequencing_data_regression(ssod)
                if final_time >= min_len and len(seq_data) >= min_seq:
                    osr_val = er.tempest_regression(seq_data).coef_[0]
                    sod_rows.append({
                        'R_long': r_long_val,
                        'infected_individuals_at_start': infected_start_val,
                        'nucleotide_substitution_rate': nsr_val,
                        'observed_substitution_rate': osr_val,
                    })
            except Exception:
                continue

        if sod_rows:
            sod_df = pd.DataFrame(sod_rows)
            sod_df = om.detect_sod_outliers(sod_df)
            sod_df['mean_final_time'] = np.mean(final_times) if final_times else np.nan
            sod_df['n_seeds_total'] = len(seeded_dirs)
            all_rows.append(sod_df)

    if not all_rows:
        raise RuntimeError(
            f"No valid data found in {numbered}. "
            f"Check if simulations finished and lengths/seqs meet minimums.")

    master_df = pd.concat(all_rows, ignore_index=True)
    clean_df = master_df[master_df['is_outlier'] == 0]

    n_rows, n_cols = len(R_LONG_VALUES), len(INFECTED_START_VALUES)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.5 * n_rows),
                             sharex=True, sharey=True)

    print(f"\n--- Fit results (tau_3_long={tau_3_long:.2f}) ---")
    for i, r_long in enumerate(R_LONG_VALUES):
        for j, infected_start in enumerate(INFECTED_START_VALUES):
            ax = axes[i, j]
            cell_df = clean_df[(clean_df['R_long'] == r_long) &
                               (clean_df['infected_individuals_at_start'] == infected_start)]

            label = f"R_long={r_long}, start={infected_start}"
            if cell_df.empty:
                ax.text(0.5, 0.5, "no valid data", ha='center', va='center',
                       transform=ax.transAxes)
                print(f"{label:30s}: no valid data")
            else:
                ax.scatter(cell_df['nucleotide_substitution_rate'],
                          cell_df['observed_substitution_rate'],
                          alpha=0.2, s=10, color='tab:blue')
                try:
                    fit_result = er.fit_observed_substitution_rate_regressor(
                        numbered, cell_df, model_type=model_type,
                        parameter_name='nucleotide_substitution_rate',
                        experiment_group=label)
                    x_vals = np.linspace(cell_df['nucleotide_substitution_rate'].min(),
                                        cell_df['nucleotide_substitution_rate'].max(), 100)
                    y_vals = fit_result.eval(x=x_vals)
                    ax.plot(x_vals, y_vals, color='tab:red', linewidth=2)
                    r2 = fit_result.rsquared
                    n_pts = len(cell_df)
                    mean_ft = cell_df['mean_final_time'].iloc[0]
                    ax.text(0.03, 0.97, f"R2={r2:.3f}\nn={n_pts}\nmean t={mean_ft:.0f}d",
                           transform=ax.transAxes, va='top', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                    print(f"{label:30s}: R2={r2:.4f}  n_points={n_pts}  "
                         f"mean_final_time={mean_ft:.1f}")
                except Exception as e:
                    ax.text(0.5, 0.5, "fit failed", ha='center', va='center',
                           transform=ax.transAxes)
                    print(f"{label:30s}: fit failed ({e})")

            ax.axhline(TARGET_OSR_LONG, color='black', linestyle='--',
                      linewidth=1, alpha=0.6)
            if i == 0:
                ax.set_title(f"start={infected_start}")
            if j == 0:
                ax.set_ylabel(f"R_long={r_long}\nOSR")
            if i == n_rows - 1:
                ax.set_xlabel("NSR")

    fig.suptitle(f"Long-shedder isolated calibration grid "
                f"({numbered}, tau_3_long={tau_3_long:.2f})")
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f"{numbered}_r_long_grid.png")
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"\nGrid plot saved to: {plot_path}")
    return plot_path


def main():
    parser = argparse.ArgumentParser(
        description="HPC test: sweep the isolated long-shedder NSR calibration "
                    "over a grid of R_long x starting-infected-count values, "
                    "for a single tau_3_long (borrowed from one of the 5 "
                    "production scenarios), at final_time=3 years. Produces "
                    "one combined figure to visually compare calibration "
                    "curves across the grid.")
    parser.add_argument('test_n', type=int,
                        help="Test number (numbers this test's experiment).")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--n-seeds', type=int, default=10,
                        help="Seeds per (R_long, start, NSR) grid point (default 10).")
    parser.add_argument('--scenario', type=str, default='edge_case',
                        choices=[s['name'] for s in SCENARIOS],
                        help="Which of the 5 production scenarios' tau_3_long "
                            "to use (default 'edge_case' -- the longest-"
                            "duration, hardest-to-calibrate scenario).")
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()
    scenario = next(s for s in SCENARIOS if s['name'] == args.scenario)
    tau_3_long = derive_tau_3_long(scenario, sp)

    experiment_name = f"{EXP_NAME}_#{args.test_n}"
    print("\n=========================================================")
    print(f" R_long / starting-infected grid test #{args.test_n}")
    print(f" scenario={args.scenario}  tau_3_long={tau_3_long:.2f}")
    print(f" R_long values: {R_LONG_VALUES}")
    print(f" starting-infected values: {INFECTED_START_VALUES}")
    print(f" population_size={POPULATION_SIZE} (fixed)  final_time={FINAL_TIME} (3 years)")
    print("=========================================================\n")

    settings_func = build_settings(tau_3_long, args.n_seeds)
    run_experiment(experiment_name, settings_func,
                   simplicity_runner=RUNNERS[args.runner], archive_experiment=False)

    analyze_and_plot(experiment_name, tau_3_long)


if __name__ == "__main__":
    main()
