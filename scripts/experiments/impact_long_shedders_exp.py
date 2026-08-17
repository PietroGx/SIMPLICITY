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
# impact_long_shedders -- EXPERIMENT RUNNER (final stage)
# ----------------------------------------------------------------------------
# Reads the FROZEN nsr_calibration_table.csv produced by
# impact_long_shedders_cal_2.py and runs one simulation scenario per row.
#
# Both rates are set as INDEPENDENT ABSOLUTE parameters:
#   - nucleotide_substitution_rate        (standard, calibrated in-context)
#   - nucleotide_substitution_rate_long   (long shedders, per tau group)
#
# This script performs NO calibration, NO sweeps, NO fitting. It only reads
# frozen values and dispatches. If the table is missing a scenario or a value,
# it fails loudly rather than guessing.
#
# The calibration table path is derived internally from --exp-num.
# ============================================================================

import os
import argparse
import pandas as pd

from experiment_script_runner import run_experiment_script
from impact_long_shedders_config import (
    build_exp_scenario_settings, add_slurm_resource_args, set_slurm_resource_env,
)

EXP_NAME = "impact_long_shedders"

# Location of the frozen calibration table (written by cal_2). Derived from
# --exp-num; not passed as an argument.
SETUP_DIR_TEMPLATE = "Data/impact_long_shedders_calibration_setup_data_#{exp_num}"
TABLE_FILENAME = "nsr_calibration_table.csv"

# Columns the frozen table must contain.
REQUIRED_COLUMNS = [
    "scenario_name", "tau_3_long", "long_shedders_ratio", "R", "R_long",
    "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
]


def default_table_path(exp_num):
    """Derive the frozen calibration table path from the experiment number."""
    return os.path.join(SETUP_DIR_TEMPLATE.format(exp_num=exp_num), TABLE_FILENAME)


def load_calibration_table(path):
    """Load and validate the frozen NSR table."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Calibration table not found: {path}\n"
            f"Run impact_long_shedders_cal_2.py first.")

    df = pd.read_csv(path)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Calibration table {path} is missing required columns: {missing}")

    if df.empty:
        raise ValueError(f"Calibration table {path} contains no scenarios.")

    return df


def dispatch_scenario(row, exp_num, runner, n_seeds):
    """Dispatch a single scenario run. Blocks until the runner completes."""
    name = row["scenario_name"]
    settings_func = build_exp_scenario_settings(row, n_seeds)

    print(f"\n{'='*60}")
    print(f"[Dispatch] {EXP_NAME}_{name}")
    print(f"   R            : {float(row['R'])}")
    print(f"   standard NSR : {float(row['nucleotide_substitution_rate']):.8f}")
    if float(row["long_shedders_ratio"]) > 0.0:
        print(f"   long NSR     : {float(row['nucleotide_substitution_rate_long']):.8f}")
        print(f"   tau_3_long   : {float(row['tau_3_long'])}")
        print(f"   R_long       : {float(row['R_long']):.4f}")
        print(f"   LS ratio     : {float(row['long_shedders_ratio'])}")
    else:
        print(f"   (control: no long shedders)")
    print(f"{'='*60}")

    run_experiment_script(runner, exp_num, settings_func, f"{EXP_NAME}_{name}")


def main():
    parser = argparse.ArgumentParser(
        description="Run impact_long_shedders scenarios from a frozen NSR table.")
    parser.add_argument('--seeds', type=int, required=True,
                        help="Number of seeds per scenario.")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number (also selects the calibration table).")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--only', type=str, default=None,
                        help="Optional: run only this scenario_name.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)

    table_path = default_table_path(args.exp_num)
    df = load_calibration_table(table_path)

    if args.only is not None:
        df = df[df["scenario_name"] == args.only]
        if df.empty:
            raise ValueError(
                f"Scenario '{args.only}' not found in the calibration table.")

    print(f"\n[Runner] Dispatching {len(df)} scenario(s) from {table_path}")

    for _, row in df.iterrows():
        dispatch_scenario(row, args.exp_num, args.runner, args.seeds)

    print(f"\n[Success] All scenarios dispatched.")


if __name__ == "__main__":
    main()
