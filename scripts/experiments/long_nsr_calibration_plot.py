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
# Constants only: both pipelines share this module, and the scenario list it
# used to import from here is what put a bound-pipeline label on an unbound
# group (see _cal1_groups).
from impact_long_shedders_config import TAU_ROUND, DEFAULT_COLORS


def _cal1_groups(experiment_name):
    """Every group of a cal_1 experiment, in declared order, with the
    (tau_3_long, R_long) cohort it covers and a plot colour.

    The cohort comes from the group's OWN simulation records, and a group
    holding more than one pair raises rather than silently taking the first:
    a cal_1 group IS one long-shedder cohort swept over NSR_long, and anything
    else means the config and this fit disagree about what a group is.

    This used to be a {(tau, R_long) -> name} map, with the partition itself
    done by re-reading those two floats off every simulation output directory.
    Before that the name was re-derived by matching the pair against a scenario
    list passed in by the caller, defaulting to the BOUND pipeline's SCENARIOS
    -- the unbound pipeline's edge_case is not in that list, so its group fell
    through to a raw 'tau=350.23,R_long=1.1' in unbound run #1, and that string
    landed in the saved fit-results filename.

    Names are joined with '+' by the config, never '/', which
    write_fit_results_csv's file-saving reads as a directory separator.
    """
    pairs = {}
    for record in sm.read_simulations(experiment_name):
        params = record['parameters']
        pairs.setdefault(record['group'], set()).add(
            (round(float(params['tau_3_long']), TAU_ROUND),
             round(float(params['R_long']), TAU_ROUND)))

    color_cycle = iter(DEFAULT_COLORS)
    groups = []
    for name in dm.get_groups(experiment_name):
        cohorts = pairs.get(name)
        if not cohorts:
            continue
        if len(cohorts) > 1:
            raise RuntimeError(
                f"Group '{name}' of {experiment_name} covers more than one "
                f"(tau_3_long, R_long) cohort: {sorted(cohorts)}. A cal_1 "
                f"group is one cohort swept over NSR_long.")
        tau, r_long = cohorts.pop()
        groups.append({'name': name, 'tau_3_long': tau, 'R_long': r_long,
                       'color': next(color_cycle, 'black')})

    # Two groups on one cohort is a config error -- unique_long_tau_r_long_pairs
    # returns a set, so neither pipeline can produce it -- but it has to be
    # loud rather than silent: the frozen CSV this stage hands cal_2 is keyed
    # on the pair, so the second group would overwrite the first and cal_2
    # would calibrate one cohort against the other's rate.
    by_key = {}
    for group in groups:
        by_key.setdefault((group['tau_3_long'], group['R_long']), []).append(
            group['name'])
    clashing = {key: names for key, names in by_key.items() if len(names) > 1}
    if clashing:
        raise RuntimeError(
            f"{experiment_name} has groups sharing one (tau_3_long, R_long) "
            f"cohort at {TAU_ROUND} decimals: {clashing}. cal_2 looks the "
            f"calibrated long NSR up by that pair and could not tell them "
            f"apart.")
    return groups


def _calibrated_long_nsr_path(experiment_name):
    return os.path.join(dm.get_experiment_fit_result_dir(experiment_name),
                        f'{experiment_name}_calibrated_long_nsr.csv')


def write_calibrated_long_nsr(experiment_name, results, target_osr_long,
                              model_type, group_by_key=None):
    """Persist Stage 1's calibrated long NSR per (tau_3_long, R_long) group
    to a small CSV -- the single computed result, so cal_2.py can read it
    straight back instead of independently re-extracting/re-fitting from
    raw simulation output a second time (two implementations of the same
    fit is exactly what let cal_1's and cal_2's methodology drift apart
    before; see BACKLOG).

    Keyed on the pair, not on the group name, and deliberately so: cal_1 and
    cal_2 partition differently. cal_1 fits one curve per long-shedder cohort,
    merging HIV_low and HIV_high because they share one; cal_2 then looks each
    of its own scenarios up against that cohort. The group name rides along as
    a column so the handoff can be audited, never matched on.
    """
    group_by_key = group_by_key or {}
    rows = [
        {'group': group_by_key.get((tau, r_long), ''),
         'tau_3_long': tau, 'R_long': r_long,
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
                                      model_type='exp', min_seq=30, min_len=100,
                                      ):
    """
    Extracts per-seed OSR for every (tau_3_long, R_long, NSR) grid point of
    `experiment_name`, outlier-filters per grid point, fits a regressor per
    DECLARED GROUP, inverts at `target_osr_long`, and saves a multi-curve
    calibration figure into the experiment's 05_Plots folder.

    One group is one long-shedder cohort: two scenarios can share a duration
    but differ in R_long, so the cohort is the (tau_3_long, R_long) pair and
    never tau alone. Which simulations belong to it is read from the record,
    not recomputed from their parameter values -- see _cal1_groups.

    Group names for the legend and the saved fit-results filename come from the
    experiment's own settings, where the config that defined the group wrote
    them. There used to be a `scenarios` argument here for the caller to say
    which scenario list to match pairs against -- needed only because this
    module re-derived the names instead of reading them.

    Returns {(tau_3_long, R_long): calibrated_nsr}.
    """
    groups = _cal1_groups(experiment_name)

    print(f"--- Analyzing long-shedder NSR calibration: {experiment_name} ---")

    # Partitioned by the group each simulation was ASSIGNED to, not by reading
    # (tau_3_long, R_long) back off every output directory and bucketing on the
    # float pair. The two agree while every group carries a distinct rounded
    # pair, which is what unique_long_tau_r_long_pairs builds; they stop
    # agreeing the moment one does not, and nothing would say so.
    per_group = {}
    for group in groups:
        sod_frames = []
        for sod in dm.get_simulation_output_dirs(experiment_name,
                                                 group=group['name']):
            nsr_val = sm.get_parameter_value_from_simulation_output_dir(
                sod, 'nucleotide_substitution_rate_long')

            sod_rows = []
            for ssod in dm.get_seeded_simulation_output_dirs(sod):
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
                            'group': group['name'],
                            'nucleotide_substitution_rate_long': nsr_val,
                            'observed_substitution_rate': osr_val,
                        })
                except Exception:
                    continue

            if sod_rows:
                sod_df = om.detect_sod_outliers(pd.DataFrame(sod_rows))
                sod_frames.append(sod_df[sod_df['is_outlier'] == 0])

        if sod_frames:
            per_group[group['name']] = pd.concat(sod_frames, ignore_index=True)

    if not per_group:
        raise RuntimeError(
            f"No valid long-calibration data found in {experiment_name}. "
            f"Check if simulations finished and lengths/seqs meet minimums.")

    plt.figure(figsize=(10, 6))
    results = {}
    group_by_key = {}

    print(f"\n--- Calibration results (target long OSR = {target_osr_long}) ---")

    for group in groups:
        label = group['name']
        group_df = per_group.get(label)
        if group_df is None or group_df.empty:
            print(f"{label:<17}: no usable data, skipped")
            continue
        key = (group['tau_3_long'], group['R_long'])
        group_by_key[key] = label
        color = group['color']

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

        # geomspace, not linspace: the sweep is geometric and the axis is log,
        # so evenly-spaced-in-linear points would put almost none of the curve
        # in the lowest decade -- which is where the target is crossed
        x_vals = np.geomspace(group_df['nucleotide_substitution_rate_long'].min(),
                              group_df['nucleotide_substitution_rate_long'].max(), 200)
        y_vals = fit_result.eval(x=x_vals)
        plt.plot(x_vals, y_vals, color=color, linewidth=2, label=f"{label} fit")
        plt.plot(calibrated_nsr, target_osr_long, marker='*', markersize=12,
                color=color, markeredgecolor='black')

    plt.axhline(target_osr_long, color='black', linestyle='--', linewidth=1.5,
              label='Target OSR')
    # Log x, matching cal_2's plot, which has had it all along. The sweep is
    # np.geomspace(1e-4, 1e-2, 10), so on a linear axis five of the ten points
    # fall below 7.8e-4 and pile into the leftmost ~8% of the panel -- which is
    # exactly where the target line is crossed and the calibrated value is read
    # off. The figure was compressing the one region it exists to show.
    plt.xscale('log')
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

    write_calibrated_long_nsr(experiment_name, results, target_osr_long,
                              model_type, group_by_key)

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
