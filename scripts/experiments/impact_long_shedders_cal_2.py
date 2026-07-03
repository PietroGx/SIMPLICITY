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
# impact_long_shedders -- CALIBRATION SCRIPT (File 1 of 2)
# ----------------------------------------------------------------------------
# Produces a FROZEN nsr_calibration_table.csv consumed by the experiment
# runner (impact_long_shedders_exp.py). This script performs NO final
# experiment runs; it only calibrates and writes the table.
#
# Per scenario, two independent absolute rates are derived:
#   1. nucleotide_substitution_rate_long  -- from the EXISTING long-shedder
#      calibration experiment (default: calibrate_long_nsr_#2), per tau_3_long
#      group, exp-fit + inverted at target long OSR.
#   2. nucleotide_substitution_rate (standard) -- calibrated IN-CONTEXT: a
#      mixed-population NSR sweep with the long side fully frozen, fitting OSR
#      measured in STANDARD individuals only, inverted at target standard OSR.
#
# NOTE: This is a TEST pipeline (low seeds). Rerun at high seeds once the
# mechanics are confirmed.
# ============================================================================

import os
import argparse
import numpy as np
import pandas as pd

import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

from experiments.experiment_script_runner import run_experiment_script

# =============================================================================
# CONFIGURATION
# =============================================================================
EXP_NAME = "impact_long_shedders_calibration"

# Fixed experiment overrides (shared by all scenarios). Mirrors the population
# context of the impact experiment.
USER_FIXED_PARAMS = {
    "population_size": 1000,
    "infected_individuals_at_start": 100,
    "R": 1.05,
    "final_time": 1095,
    "IH_virus_emergence_rate": 0.1,
}

# Scenario definitions: (name, D, long_shedders_ratio)
#   - 'freq' multiplier removed: R_long is now tau_3_long / 7.
#   - control is a placeholder (ratio=0): no long side, uses standard defaults.
#
# TODO(temp-hack): The long-shedder calibration (calibrate_long_nsr_#2) stored
# tau_3_long DIRECTLY as 63/109/365 (no phase subtraction). To make the
# stage-1 join match those keys, we treat D AS tau_3_long here (no subtraction
# of tau_1+tau_2+tau_4). This is a known inconsistency with the intended
# derivation tau_3_long = D - (tau_1+tau_2+tau_4); fix the calibration and this
# derivation together in a future pass, then rerun everything.
SCENARIOS = [
    ("control",   1,     0.00),
    ("SOT",       63.0,  0.01),
    ("HIV_low",   109.0, 0.01),
    ("HIV_high",  109.0, 0.12),
    ("edge_case", 365.0, 0.01),
]

TAU_ROUND = 3  # decimals for tau_3_long dict-key matching


# =============================================================================
# STAGE 0 -- scenario parameter derivation
# =============================================================================
def derive_scenario_params(name, D, ratio, sp):
    """
    Derive the frozen long-side parameters for one scenario.

    Returns a dict with: scenario_name, long_shedders_ratio, tau_3_long,
    R_long, sequence_long_shedders, is_long.
    """
    if ratio == 0.0:
        # Control: no long shedders. Use standard defaults; long NSR unused.
        return {
            "scenario_name": name,
            "long_shedders_ratio": 0.0,
            "tau_3_long": float(sp["tau_3_long"]),
            "R_long": float(sp["R_long"]),
            "sequence_long_shedders": False,
            "is_long": False,
        }

    # TODO(temp-hack): no subtraction -- tau_3_long = D (see SCENARIOS note).
    tau_3_long = float(D)
    R_long = tau_3_long / 7.0
    return {
        "scenario_name": name,
        "long_shedders_ratio": float(ratio),
        "tau_3_long": tau_3_long,
        "R_long": R_long,
        "sequence_long_shedders": True,
        "is_long": True,
    }


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
def build_standard_sweep_settings(frozen, long_nsr, nsr_min, nsr_max, steps, n_seeds):
    """
    Returns a make_settings closure for a standard-NSR sweep with the long
    side fully frozen. Long NSR is set as the new absolute parameter.
    """
    nsr_values = np.logspace(np.log10(nsr_min), np.log10(nsr_max), num=steps)
    nsr_list = [float(x) for x in nsr_values]

    fixed = USER_FIXED_PARAMS.copy()
    fixed.update({
        "long_shedders_ratio": frozen["long_shedders_ratio"],
        "tau_3_long": frozen["tau_3_long"],
        "R_long": frozen["R_long"],
        "sequence_long_shedders": frozen["sequence_long_shedders"],
    })
    # Only set long NSR when there is a long side to calibrate against.
    if frozen["is_long"]:
        fixed["nucleotide_substitution_rate_long"] = float(long_nsr)

    def make_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = fixed.copy()
        return (varying_params, fixed_params, n_seeds)

    return make_settings


def calibrate_standard_nsr_in_context(sweep_exp_name, exp_num, runner,
                                      settings_func, target_osr_std,
                                      model_type='exp', min_seq=30, min_len=100):
    """
    Runs the in-context standard sweep, fits OSR measured in STANDARD
    individuals only, and inverts at target_osr_std.

    Relies on the write->read CSV path: outliers are auto-dropped by
    read_OSR_vs_parameter_csv (include_outliers=False default), so no manual
    detect_sod_outliers call here.
    """
    print(f"\n[Stage 2] Standard NSR in-context sweep: {sweep_exp_name}")
    run_experiment_script(runner, exp_num, settings_func, sweep_exp_name)

    numbered = f"{sweep_exp_name}_#{exp_num}"

    # Type-filtered extraction (Edit 7 makes 'standard' filtering effective).
    om.write_OSR_vs_parameter_csv(
        numbered, 'nucleotide_substitution_rate',
        min_seq, min_len, individual_type='standard')

    osr_df = om.read_OSR_vs_parameter_csv(
        numbered, 'nucleotide_substitution_rate',
        min_seq, min_len, individual_type='standard')

    if osr_df.empty:
        raise RuntimeError(
            f"[Stage 2] No standard-individual OSR data for {numbered}. "
            f"Sweep may have produced too few standard sequences "
            f"(min_seq={min_seq}, min_len={min_len}).")

    fit_result = er.fit_observed_substitution_rate_regressor(
        numbered, osr_df, model_type,
        parameter_name='nucleotide_substitution_rate')

    nsr_std = er.compute_calibrated_parameter(model_type, fit_result, target_osr_std)
    print(f"          calibrated standard NSR = {nsr_std:.8f}")
    return float(nsr_std)


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
                        help="Seeds per NSR point in the standard sweeps.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--long-calib-exp', type=str, default='calibrate_long_nsr_#2',
                        help="Existing long-shedder calibration experiment (numbered).")
    parser.add_argument('--min-nsr', type=float, default=3e-5,
                        help="Min standard NSR for sweep.")
    parser.add_argument('--max-nsr', type=float, default=5e-4,
                        help="Max standard NSR for sweep.")
    parser.add_argument('--steps', type=int, default=8,
                        help="NSR points per standard sweep.")
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'], help="Fit model.")
    parser.add_argument('--min-seq', type=int, default=30,
                        help="Min sequences to keep a seed.")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Min simulation length (days) to keep a seed.")
    args = parser.parse_args()

    # Apply pop overrides that come from defaults (kept for parity / future args).
    sp = sm.read_standard_parameters_values()

    setup_dir = f"Data/{EXP_NAME}_setup_data_#{args.exp_num}"
    os.makedirs(setup_dir, exist_ok=True)

    # --- Stage 1: long NSR per tau group (once, reused across scenarios) ---
    long_nsr_by_tau = compute_long_nsr_per_tau(
        args.long_calib_exp, args.target_osr_long,
        model_type=args.model, min_seq=args.min_seq, min_len=args.min_len)

    # --- Per scenario: derive -> lookup long NSR -> calibrate standard NSR ---
    rows = []
    for name, D, ratio in SCENARIOS:
        print(f"\n{'='*60}\n[Scenario] {name}\n{'='*60}")
        frozen = derive_scenario_params(name, D, ratio, sp)

        if frozen["is_long"]:
            long_nsr = lookup_long_nsr(long_nsr_by_tau, frozen["tau_3_long"])
        else:
            long_nsr = None  # control: no long side

        settings_func = build_standard_sweep_settings(
            frozen, long_nsr, args.min_nsr, args.max_nsr, args.steps, args.seeds)

        sweep_name = f"imp_ls_cal_{name}"
        std_nsr = calibrate_standard_nsr_in_context(
            sweep_name, args.exp_num, args.runner, settings_func,
            args.target_osr_std, model_type=args.model,
            min_seq=args.min_seq, min_len=args.min_len)

        rows.append({
            "scenario_name": name,
            "tau_3_long": frozen["tau_3_long"],
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "R_long": frozen["R_long"],
            "nucleotide_substitution_rate_long": long_nsr,
            "nucleotide_substitution_rate": std_nsr,
            "target_osr_std": args.target_osr_std,
            "target_osr_long": args.target_osr_long,
            "model_type": args.model,
            "long_calib_source": args.long_calib_exp,
        })

    # --- Write frozen table ---
    table_path = os.path.join(setup_dir, "nsr_calibration_table.csv")
    cols = ["scenario_name", "tau_3_long", "long_shedders_ratio", "R_long",
            "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
            "target_osr_std", "target_osr_long", "model_type", "long_calib_source"]
    pd.DataFrame(rows)[cols].to_csv(table_path, index=False)

    print(f"\n[Success] Frozen calibration table written:\n  {table_path}")
    print("Next: run impact_long_shedders_exp.py against this table.")


if __name__ == "__main__":
    main()
