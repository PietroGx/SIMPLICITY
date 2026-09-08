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
# impact_long_shedders_unbound -- STAGE 3: PRODUCTION
# ----------------------------------------------------------------------------
# Reads the frozen table and runs one scenario per row. Both rates are applied
# exactly as calibrated -- NSR_long from the 100%-long-shedder population,
# the standard NSR from the 0%-long-shedder one, the same standard rate for
# every scenario. Nothing compensates for their interaction, so the
# population-level clock is an OUTPUT of these runs, not a target they were
# tuned to hit.
#
# No calibration, no sweeps, no fitting. Fails loudly on a missing value.
# ============================================================================

import os
import argparse
import pandas as pd

from experiment_script_runner import run_experiment_script
from impact_long_shedders_unbound_config import (
    PROD_EXP_NAME, SETUP_DIR_TEMPLATE, TABLE_FILENAME, USER_FIXED_PARAMS,
    build_exp_scenario_settings, add_slurm_resource_args,
    set_slurm_resource_env, print_fixed_params,
)

REQUIRED_COLUMNS = [
    "scenario_name", "tau_3_long", "long_shedders_ratio", "susceptibility_long",
    "R", "IH_virus_emergence_rate", "R_long",
    "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
]


def default_table_path(exp_num):
    return os.path.join(SETUP_DIR_TEMPLATE.format(exp_num=exp_num), TABLE_FILENAME)


def load_calibration_table(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Calibration table not found: {path}\n"
            f"Run impact_long_shedders_unbound_cal_2.py first.")
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Calibration table {path} is missing columns: {missing}")
    if df.empty:
        raise ValueError(f"Calibration table {path} contains no scenarios.")
    return df


def dispatch_scenario(row, exp_num, runner, n_seeds):
    name = row["scenario_name"]
    settings_func = build_exp_scenario_settings(row, n_seeds)

    print(f"\n{'='*60}")
    print(f"[Dispatch] {PROD_EXP_NAME}_{name}")
    print(f"   R            : {float(row['R'])}")
    print(f"   k_v          : {float(row['IH_virus_emergence_rate'])}")
    print(f"   standard NSR : {float(row['nucleotide_substitution_rate']):.8f}")
    if float(row["long_shedders_ratio"]) > 0.0:
        print(f"   long NSR     : {float(row['nucleotide_substitution_rate_long']):.8f}")
        print(f"   tau_3_long   : {float(row['tau_3_long'])}")
        print(f"   R_long       : {float(row['R_long']):.4f}")
        print(f"   LS ratio     : {float(row['long_shedders_ratio'])}")
    else:
        print(f"   (control: no long shedders)")
    print(f"{'='*60}")

    run_experiment_script(runner, exp_num, settings_func, f"{PROD_EXP_NAME}_{name}")


def main():
    parser = argparse.ArgumentParser(
        description="Unbound pipeline stage 3: production runs from the frozen "
                    "table, both rates applied as calibrated.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=50)
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
            raise ValueError(f"Scenario '{args.only}' not found in the table.")

    print(f"\n[Runner] Dispatching {len(df)} scenario(s) from {table_path}")
    shared = {k: v for k, v in USER_FIXED_PARAMS.items() if k != "R"}
    print_fixed_params(shared, label="Shared parameters (R shown per scenario below)")

    for _, row in df.iterrows():
        dispatch_scenario(row, args.exp_num, args.runner, args.seeds)

    print(f"\n[Success] All scenarios dispatched.")


if __name__ == "__main__":
    main()
