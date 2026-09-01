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
# HPC test: in-context (mixed-population) STANDARD-individual NSR
# calibration -- the Stage 2 (cal_2) counterpart to
# test_long_shedder_r_long_grid.py's Stage 1 (cal_1) check.
#
# Real pipeline run #1 (exp-num 1, 2026-08) got past Stage 1 cleanly (its
# widened NSR_RANGES['cal1_long_nsr']=0.01-0.5 worked exactly as validated),
# then Stage 2 aborted: SOT's standard-NSR fit came back essentially flat
# (R2=0.0003) within the OLD NSR_RANGES['cal2_standard_nsr']=1e-5-3e-4,
# triggering the NSR_SANITY_MAX guard. `control` (no long shedders) fit
# fine in that same run (0.00019, within range) -- so the problem is
# specific to scenarios WITH long shedders: their now much higher
# calibrated long-NSR (0.06-0.12, from Stage 1) injects far more diversity
# into the mixed population than before, apparently swamping the standard
# individuals' own signal within the old narrow range. This test checks
# whether a wider NSR_RANGES['cal2_standard_nsr'] fixes that, the same way
# widening cal1_long_nsr fixed Stage 1.
#
# Uses the REAL calibrated long-NSR values Stage 1 produced in run #1
# (LONG_NSR_BY_SCENARIO below) as fixed inputs -- these would need
# re-deriving if cal_1's R_long or NSR_RANGES['cal1_long_nsr'] ever change
# again. R_long is NOT swept here (already decided, fixed at whatever
# SCENARIOS currently has) -- only standard NSR varies, per scenario
# (including `control`, tested for comparison even though it already
# works, as a sanity reference). Not part of the impact_long_shedders
# pipeline itself -- a standalone exploratory tool, decoupled from
# NSR_RANGES['cal2_standard_nsr'] (exploring here never changes that
# shared production value).
#
# Reuses USER_FIXED_PARAMS from impact_long_shedders_config rather than
# duplicating it.
#
# Writes a plain-text recap (Data/tests/test_cal_2_out_#{test_n}.txt) with
# per-scenario fit stats AND raw per-NSR-node mean OSR, plus the combined
# grid figure, both zipped into Data/{experiment_name}_artifacts.zip.
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
    SCENARIOS, TAU_ROUND, derive_tau_3_long, USER_FIXED_PARAMS,
)

EXP_NAME = "long_shedder_std_nsr_test"

TARGET_OSR_STD = 0.0013       # reference line only (same default as cal_2)
DEFAULT_NSR_MIN = 0.00001
DEFAULT_NSR_MAX = 0.05

# Stage 1's real fitted long-NSR per scenario, from run #1 (R_long=1.1,
# NSR_RANGES['cal1_long_nsr']=0.01-0.5). HIV_low and HIV_high share
# tau_3_long=94.23 so share a fit. `control` has no long side.
LONG_NSR_BY_SCENARIO = {
    "SOT": 0.086046,
    "HIV_low": 0.062089,
    "HIV_high": 0.062089,
}

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}


def build_settings(nsr_values, sp, n_seeds):
    """One scenario group per SCENARIOS entry, each sweeping the same NSR
    anchor points over `nucleotide_substitution_rate` (standard). R_long/
    tau_3_long/long-NSR/ratio are fixed per group, not swept -- already
    decided by Stage 1."""
    scenario_groups = []
    for scenario in SCENARIOS:
        name = scenario["name"]
        is_control = scenario["long_shedders_ratio"] <= 0.0
        group = {
            'nucleotide_substitution_rate': list(nsr_values),
            'long_shedders_ratio': float(scenario["long_shedders_ratio"]),
            'tau_3_long': derive_tau_3_long(scenario, sp),
            'sequence_long_shedders': not is_control,
        }
        if not is_control:
            group['R_long'] = float(scenario["R_long"])
            group['nucleotide_substitution_rate_long'] = LONG_NSR_BY_SCENARIO[name]
        scenario_groups.append(group)

    fixed_params = USER_FIXED_PARAMS.copy()

    def make_settings():
        return ({'_scenario_groups': scenario_groups}, fixed_params, n_seeds)

    return make_settings


def analyze_and_report(numbered, nsr_values, sp, test_n, min_seq=30, min_len=100, model_type='exp'):
    """Extracts per-seed OSR (STANDARD individuals only) for every
    (scenario, NSR) point, outlier-filters per grid point, fits a
    regressor per scenario, saves ONE combined figure (one row per
    scenario, shared axes), writes a plain-text recap with fit stats AND
    raw per-NSR-node mean OSR, and zips the recap + figure into Data/.
    """
    print(f"--- Analyzing standard-NSR grid: {numbered} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(numbered)

    all_rows = []
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
        ratio_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'long_shedders_ratio')

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_rows = []
        final_times = []
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                final_times.append(final_time)
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
            sod_df['mean_final_time'] = np.mean(final_times) if final_times else np.nan
            sod_df['n_seeds_total'] = len(seeded_dirs)
            all_rows.append(sod_df)

    master_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    clean_df = (master_df[master_df['is_outlier'] == 0]
               if not master_df.empty else master_df)

    n_rows = len(SCENARIOS)
    fig, axes = plt.subplots(n_rows, 1, figsize=(6, 3 * n_rows), sharex=True, sharey=True)
    if n_rows == 1:
        axes = [axes]

    recap = []
    recap.append(f"=== Long-shedder standard-NSR calibration test: {numbered} ===")
    recap.append(f"NSR anchor points: {nsr_values}")
    recap.append("Scenarios: " + ", ".join(s["name"] for s in SCENARIOS))
    recap.append(f"Fixed long-NSR per scenario (from Stage 1 run #1): {LONG_NSR_BY_SCENARIO}")
    recap.append(f"USER_FIXED_PARAMS: {USER_FIXED_PARAMS}")
    recap.append(f"Target standard OSR = {TARGET_OSR_STD}")
    recap.append("")

    for i, scenario in enumerate(SCENARIOS):
        name = scenario["name"]
        ax = axes[i]
        recap.append(f"--- {name} ---")

        # long_shedders_ratio alone is NOT unique across scenarios (SOT,
        # HIV_low both share ratio=0.01) -- match on the
        # (tau_3_long, ratio) pair together, same as elsewhere in this repo.
        key_tau = round(derive_tau_3_long(scenario, sp), TAU_ROUND)
        key_ratio = round(float(scenario["long_shedders_ratio"]), TAU_ROUND)
        cell_df = (clean_df[(clean_df['tau_3_long'].round(TAU_ROUND) == key_tau) &
                            (clean_df['long_shedders_ratio'].round(TAU_ROUND) == key_ratio)]
                  if not clean_df.empty else clean_df)

        if cell_df.empty:
            ax.text(0.5, 0.5, "no valid data", ha='center', va='center', transform=ax.transAxes)
            recap.append("  no valid data")
            print(f"{name:12s}: no valid data")
        else:
            ax.scatter(cell_df['nucleotide_substitution_rate'],
                      cell_df['observed_substitution_rate'],
                      alpha=0.3, s=15, color='tab:blue')
            for nsr in nsr_values:
                sub = cell_df[np.isclose(cell_df['nucleotide_substitution_rate'], nsr)]
                if not sub.empty:
                    recap.append(f"  NSR={nsr:.6f} -> mean OSR="
                                f"{sub['observed_substitution_rate'].mean():.6f} (n={len(sub)})")
                else:
                    recap.append(f"  NSR={nsr:.6f} -> no valid points")
            try:
                fit_result = er.fit_observed_substitution_rate_regressor(
                    numbered, cell_df, model_type=model_type,
                    parameter_name='nucleotide_substitution_rate',
                    experiment_group=name)
                x_vals = np.linspace(cell_df['nucleotide_substitution_rate'].min(),
                                    cell_df['nucleotide_substitution_rate'].max(), 100)
                y_vals = fit_result.eval(x=x_vals)
                ax.plot(x_vals, y_vals, color='tab:red', linewidth=2)
                r2 = fit_result.rsquared
                n_pts = len(cell_df)
                mean_ft = cell_df['mean_final_time'].iloc[0]
                fit_params = {p_name: p.value for p_name, p in fit_result.params.items()}
                ax.text(0.03, 0.97, f"R2={r2:.3f}\nn={n_pts}\nmean t={mean_ft:.0f}d",
                       transform=ax.transAxes, va='top', fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                recap.append(f"  fit: {fit_params}  R2={r2:.4f}  n={n_pts}  "
                            f"mean_final_time={mean_ft:.1f}")
                print(f"{name:12s}: R2={r2:.4f}  n={n_pts}  mean_final_time={mean_ft:.1f}")
            except Exception as e:
                ax.text(0.5, 0.5, "fit failed", ha='center', va='center', transform=ax.transAxes)
                recap.append(f"  fit failed: {e}")
                print(f"{name:12s}: fit failed ({e})")

        ax.axhline(TARGET_OSR_STD, color='black', linestyle='--', linewidth=1, alpha=0.6)
        ax.set_ylabel(f"{name}\nOSR")
        recap.append("")

    axes[-1].set_xlabel("NSR (standard)")
    fig.suptitle(f"Standard-NSR calibration by scenario ({numbered})")
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f"{numbered}_std_nsr_grid.png")
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"\nGrid plot saved to: {plot_path}")

    recap_dir = os.path.join(dm.get_data_dir(), "tests")
    os.makedirs(recap_dir, exist_ok=True)
    recap_path = os.path.join(recap_dir, f"test_cal_2_out_#{test_n}.txt")
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
        description="HPC test: sweep the in-context (mixed-population) "
                    "standard-individual NSR calibration across every "
                    "production scenario (control included), using Stage "
                    "1's real fitted long-NSR values as fixed inputs, at "
                    "3 NSR anchor points (min/mid/max). Writes a recap "
                    "file + combined figure, zipped together into Data/.")
    parser.add_argument('test_n', type=int,
                        help="Test number (numbers this test's experiment).")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--n-seeds', type=int, default=10,
                        help="Seeds per (scenario, NSR) grid point (default 10).")
    parser.add_argument('--nsr-min', type=float, default=DEFAULT_NSR_MIN,
                        help=f"Min NSR anchor point (default {DEFAULT_NSR_MIN}). "
                            "Deliberately decoupled from "
                            "NSR_RANGES['cal2_standard_nsr'] in "
                            "impact_long_shedders_config.py -- exploring "
                            "here never changes that shared production value.")
    parser.add_argument('--nsr-max', type=float, default=DEFAULT_NSR_MAX,
                        help=f"Max NSR anchor point (default {DEFAULT_NSR_MAX}).")
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()

    nsr_min, nsr_max = args.nsr_min, args.nsr_max
    nsr_mid = float(np.sqrt(nsr_min * nsr_max))
    nsr_values = [nsr_min, nsr_mid, nsr_max]

    experiment_name = f"{EXP_NAME}_#{args.test_n}"
    print("\n=========================================================")
    print(f" Standard-NSR calibration test #{args.test_n}")
    print(f" scenarios: {[s['name'] for s in SCENARIOS]}")
    print(f" NSR anchor points: {nsr_values}")
    print(f" fixed long-NSR: {LONG_NSR_BY_SCENARIO}")
    print(f" USER_FIXED_PARAMS: {USER_FIXED_PARAMS}")
    print("=========================================================\n")

    settings_func = build_settings(nsr_values, sp, args.n_seeds)
    run_experiment(experiment_name, settings_func,
                   simplicity_runner=RUNNERS[args.runner], archive_experiment=False)

    analyze_and_report(experiment_name, nsr_values, sp, args.test_n)


if __name__ == "__main__":
    main()
