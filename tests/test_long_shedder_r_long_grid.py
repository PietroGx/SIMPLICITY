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
# HPC test: isolated long-shedder NSR calibration across EVERY long-shedder
# production scenario (by unique tau_3_long -- HIV_low/HIV_high share a
# duration in this isolated 100%-long-shedder context, so are tested once
# and reported together, not duplicated) x R_long in {1.0, 1.3}, sampled at
# only 3 NSR anchor points (min/mid/max, overridable via --nsr-min/--nsr-max)
# -- this is a bracketing check (does the sampled range reach target OSR?),
# not a full curve-shape scan. final_time = 3 years.
# infected_individuals_at_start is fixed (not swept). Not part of the
# impact_long_shedders pipeline itself -- a standalone exploratory tool, and
# deliberately decoupled from NSR_RANGES['cal1_long_nsr'] in
# impact_long_shedders_config.py (exploring here never changes that shared
# production value).
#
# Mirrors cal_1's methodology exactly, and must keep doing so -- this tool is
# what gets reached for when cal_1's grid needs retuning, so any divergence
# produces confident numbers about a different quantity (see BACKLOG):
#   - sweeps nucleotide_substitution_rate_long, the rate mutations.py applies
#     to long_shedder_i -- the standard rate governs nobody at ratio 1.0;
#   - takes fixed_params from CAL1_ISOLATED_FIXED_PARAMS (ratio=1.0, so every
#     individual carries the long-shedder trait from creation, founders
#     included) rather than rebuilding them;
#   - fits each seed's genuine intra-host clock via
#     er.extract_ih_regression_data, not the population-pooled root-to-tip
#     regression.
#
# Writes a plain-text recap (Data/tests/test_cal_1_out_#{test_n}.txt) with
# per-cell fit stats AND raw per-NSR-node mean OSR (so bracketing can be
# checked directly, without trusting a fit that may only have 3 points to
# work with), plus the combined grid figure, both zipped into a single file
# at Data/{experiment_name}_artifacts.zip.
# ============================================================================

import os
import sys
import zipfile
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
    SCENARIOS, TAU_ROUND, derive_tau_3_long, CAL1_ISOLATED_FIXED_PARAMS,
    USER_FIXED_PARAMS,
)

EXP_NAME = "long_shedder_scenario_r_long_test"

R_LONG_VALUES = [1.0, 1.3]
INFECTED_START = 10           # fixed -- not swept
POPULATION_SIZE = 1000        # fixed -- not swept
FINAL_TIME = 365 * 3          # 3 years
TARGET_OSR_LONG = 0.00205     # reference line only (same default as cal_1)

# NSR_long bracket for this test. The old 0.001-0.1 default was chosen under
# the pre-2.4.24 measurement (standard rate, global clock) and does not carry
# over: on the intra-host clock, run #5 calibrated NSR_long to 5.5e-4 (SOT)
# and 1.0e-3 (HIV), so 0.001 starts at the answer and 0.1 is far into
# raw-Hamming saturation. Brackets those values instead; override via
# --nsr-min/--nsr-max.
DEFAULT_NSR_MIN = 1e-05
DEFAULT_NSR_MAX = 1e-02

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}


def unique_tau_labels(sp):
    """{rounded tau_3_long: 'SOT' / 'HIV_low+HIV_high'} across every
    long-shedder scenario, deduplicated by tau (HIV_low and HIV_high
    share a duration, hence a label, in this isolated context). '+' not '/'
    -- this label is also used as fit_observed_substitution_rate_regressor's
    experiment_group, which becomes part of a saved-results file path, and
    '/' would be read as a directory separator there."""
    labels = {}
    for s in SCENARIOS:
        if s["long_shedders_ratio"] <= 0.0:
            continue
        tau = round(derive_tau_3_long(s, sp), TAU_ROUND)
        labels.setdefault(tau, []).append(s["name"])
    return {tau: "+".join(names) for tau, names in labels.items()}


def build_settings(nsr_values, tau_labels, n_seeds):
    """One scenario group per (tau_3_long, R_long) cell, each sweeping the
    same 3-point NSR grid."""
    scenario_groups = [
        {
            'nucleotide_substitution_rate_long': list(nsr_values),
            'tau_3_long': tau,
            'R_long': r_long,
        }
        for tau in tau_labels
        for r_long in R_LONG_VALUES
    ]

    # cal_1's own isolated context, with only this test's deliberate
    # overrides on top. k_v is not in CAL1_ISOLATED_FIXED_PARAMS --
    # build_cal1_settings injects it, so it has to be set here too or it
    # silently defaults to standard_values.json's 0.
    fixed_params = CAL1_ISOLATED_FIXED_PARAMS.copy()
    fixed_params.update({
        'population_size': POPULATION_SIZE,
        'infected_individuals_at_start': INFECTED_START,
        'final_time': FINAL_TIME,
        'IH_virus_emergence_rate': USER_FIXED_PARAMS['IH_virus_emergence_rate'],
    })

    def make_settings():
        return ({'_scenario_groups': scenario_groups}, fixed_params, n_seeds)

    return make_settings


def analyze_and_report(numbered, nsr_values, tau_labels, test_n,
                       min_seq=30, min_len=100, model_type='exp'):
    """Extracts per-seed OSR for every (tau_3_long, R_long, NSR) point,
    outlier-filters per grid point, fits a regressor per (tau, R_long) cell,
    saves ONE combined grid figure (scenario rows x R_long columns, shared
    axes), writes a plain-text recap with fit stats AND raw per-NSR-node
    mean OSR, and zips the recap + figure into one file in Data/.
    """
    print(f"--- Analyzing scenario/R_long grid: {numbered} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(numbered)

    all_rows = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate_long')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
        r_long_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'R_long')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        final_times = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                final_times.append(final_time)
                seq_data = om.read_sequencing_data_regression(ssod)
                # seq_data is the eligibility gate only; the fit is on the
                # intra-host clock, same as long_nsr_calibration_plot.py.
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
            sod_df['mean_final_time'] = np.mean(final_times) if final_times else np.nan
            sod_df['n_seeds_total'] = len(seeded_dirs)
            all_rows.append(sod_df)

    master_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    clean_df = (master_df[master_df['is_outlier'] == 0]
               if not master_df.empty else master_df)

    taus = sorted(tau_labels.keys())
    n_rows, n_cols = len(taus), len(R_LONG_VALUES)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3.5 * n_rows),
                             sharex=True, sharey=True)
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    recap = []
    recap.append(f"=== Long-shedder scenario/R_long calibration test: {numbered} ===")
    recap.append("Swept: nucleotide_substitution_rate_long; "
                 "measured: intra-host clock (extract_ih_regression_data)")
    recap.append(f"NSR_long anchor points (3): {nsr_values}")
    recap.append(f"R_long values: {R_LONG_VALUES}")
    recap.append("Scenarios (by unique tau_3_long): " +
                 ", ".join(f"{lbl} (tau_3_long={tau:.2f})" for tau, lbl in tau_labels.items()))
    recap.append(f"population_size={POPULATION_SIZE} (fixed), "
                 f"infected_individuals_at_start={INFECTED_START} (fixed), "
                 f"final_time={FINAL_TIME} (3 years)")
    recap.append(f"Target OSR = {TARGET_OSR_LONG}")
    recap.append("")

    for i, tau in enumerate(taus):
        label = tau_labels[tau]
        for j, r_long in enumerate(R_LONG_VALUES):
            ax = axes[i, j]
            cell_label = f"{label} (tau={tau:.2f}), R_long={r_long}"
            recap.append(f"--- {cell_label} ---")

            cell_df = (clean_df[(clean_df['tau_3_long'].round(TAU_ROUND) == tau) &
                                (clean_df['R_long'] == r_long)]
                      if not clean_df.empty else clean_df)

            if cell_df.empty:
                ax.text(0.5, 0.5, "no valid data", ha='center', va='center',
                       transform=ax.transAxes)
                recap.append("  no valid data")
                print(f"{cell_label:40s}: no valid data")
            else:
                ax.scatter(cell_df['nucleotide_substitution_rate_long'],
                          cell_df['observed_substitution_rate'],
                          alpha=0.3, s=15, color='tab:blue')

                for nsr in nsr_values:
                    sub = cell_df[np.isclose(
                        cell_df['nucleotide_substitution_rate_long'], nsr)]
                    if not sub.empty:
                        recap.append(f"  NSR={nsr:.6f} -> mean OSR="
                                    f"{sub['observed_substitution_rate'].mean():.6f} "
                                    f"(n={len(sub)})")
                    else:
                        recap.append(f"  NSR={nsr:.6f} -> no valid points")

                try:
                    fit_result = er.fit_observed_substitution_rate_regressor(
                        numbered, cell_df, model_type=model_type,
                        parameter_name='nucleotide_substitution_rate_long',
                        experiment_group=cell_label)
                    x_vals = np.geomspace(
                        cell_df['nucleotide_substitution_rate_long'].min(),
                        cell_df['nucleotide_substitution_rate_long'].max(), 100)
                    y_vals = fit_result.eval(x=x_vals)
                    ax.plot(x_vals, y_vals, color='tab:red', linewidth=2)
                    r2 = fit_result.rsquared
                    n_pts = len(cell_df)
                    mean_ft = cell_df['mean_final_time'].iloc[0]
                    fit_params = {name: p.value for name, p in fit_result.params.items()}
                    ax.text(0.03, 0.97, f"R2={r2:.3f}\nn={n_pts}\nmean t={mean_ft:.0f}d",
                           transform=ax.transAxes, va='top', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                    recap.append(f"  fit: {fit_params}  R2={r2:.4f}  n={n_pts}  "
                                f"mean_final_time={mean_ft:.1f}")
                    print(f"{cell_label:40s}: R2={r2:.4f}  n={n_pts}  "
                         f"mean_final_time={mean_ft:.1f}")
                except Exception as e:
                    ax.text(0.5, 0.5, "fit failed", ha='center', va='center',
                           transform=ax.transAxes)
                    recap.append(f"  fit failed: {e}")
                    print(f"{cell_label:40s}: fit failed ({e})")

            ax.axhline(TARGET_OSR_LONG, color='black', linestyle='--',
                      linewidth=1, alpha=0.6)
            # Anchors are geometric over ~3 decades; linear x crushes them.
            ax.set_xscale('log')
            if i == 0:
                ax.set_title(f"R_long={r_long}")
            if j == 0:
                ax.set_ylabel(f"{label}\ntau={tau:.1f}\nintra-host OSR")
            if i == n_rows - 1:
                ax.set_xlabel("NSR_long")
            recap.append("")

    fig.suptitle(f"Long-shedder scenario/R_long calibration grid ({numbered})")
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f"{numbered}_scenario_r_long_grid.png")
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"\nGrid plot saved to: {plot_path}")

    recap_dir = os.path.join(dm.get_data_dir(), "tests")
    os.makedirs(recap_dir, exist_ok=True)
    recap_path = os.path.join(recap_dir, f"test_cal_1_out_#{test_n}.txt")
    with open(recap_path, "w") as f:
        f.write("\n".join(recap))
    print(f"Recap written to: {recap_path}")

    zip_path = os.path.join(dm.get_data_dir(), f"{numbered}_artifacts.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(recap_path, arcname=os.path.basename(recap_path))
        zf.write(plot_path, arcname=os.path.basename(plot_path))
    print(f"Artifacts zipped to: {zip_path}")

    return recap_path, plot_path, zip_path


def main():
    parser = argparse.ArgumentParser(
        description="HPC test: sweep the isolated long-shedder NSR calibration "
                    "across every long-shedder production scenario (by unique "
                    "tau_3_long) x R_long in {1.0, 1.3}, at 3 NSR anchor points "
                    "(min/mid/max), final_time=3 years, infected_individuals_"
                    "at_start fixed. Writes a recap file + combined figure, "
                    "zipped together into Data/.")
    parser.add_argument('test_n', type=int,
                        help="Test number (numbers this test's experiment).")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--n-seeds', type=int, default=10,
                        help="Seeds per (tau, R_long, NSR) grid point (default 10).")
    parser.add_argument('--nsr-min', type=float, default=DEFAULT_NSR_MIN,
                        help=f"Min NSR anchor point (default {DEFAULT_NSR_MIN}). "
                            "Deliberately decoupled from NSR_RANGES['cal1_long_nsr'] "
                            "in impact_long_shedders_config.py -- exploring here "
                            "never changes that shared production value.")
    parser.add_argument('--nsr-max', type=float, default=DEFAULT_NSR_MAX,
                        help=f"Max NSR anchor point (default {DEFAULT_NSR_MAX}).")
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()
    tau_labels = unique_tau_labels(sp)

    nsr_min, nsr_max = args.nsr_min, args.nsr_max
    nsr_mid = float(np.sqrt(nsr_min * nsr_max))
    nsr_values = [nsr_min, nsr_mid, nsr_max]

    experiment_name = f"{EXP_NAME}_#{args.test_n}"
    print("\n=========================================================")
    print(f" Scenario / R_long calibration test #{args.test_n}")
    print(f" scenarios (by tau): {tau_labels}")
    print(f" R_long values: {R_LONG_VALUES}")
    print(f" NSR_long anchor points: {nsr_values}")
    print(f" population_size={POPULATION_SIZE} (fixed)  "
         f"infected_individuals_at_start={INFECTED_START} (fixed)  "
         f"final_time={FINAL_TIME} (3 years)")
    print("=========================================================\n")

    settings_func = build_settings(nsr_values, tau_labels, args.n_seeds)
    run_experiment(experiment_name, settings_func,
                   simplicity_runner=RUNNERS[args.runner], archive_experiment=False)

    analyze_and_report(experiment_name, nsr_values, tau_labels, args.test_n)


if __name__ == "__main__":
    main()
