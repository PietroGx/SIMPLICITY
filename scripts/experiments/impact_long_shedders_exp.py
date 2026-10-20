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
# impact_long_shedders -- STAGE 3: PRODUCTION
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

from experiment_script_runner import run_experiment_phase_script
from impact_long_shedders_config import (
    EXP_NAME, table_path as default_table_path, USER_FIXED_PARAMS,
    build_exp_phase_settings, exp_group, add_slurm_resource_args,
    set_slurm_resource_env, print_fixed_params,
    set_experiment_output_env,
    CONSENSUS_MODES, DEFAULT_CONSENSUS, exp_name,
)

REQUIRED_COLUMNS = [
    "scenario_name", "tau_3_long", "long_shedders_ratio", "susceptibility_long",
    "R", "IH_virus_emergence_rate", "R_long",
    "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
]


def load_calibration_table(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Calibration table not found: {path}\n"
            f"Run impact_long_shedders_cal_2.py first.")
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Calibration table {path} is missing columns: {missing}")
    if df.empty:
        raise ValueError(f"Calibration table {path} contains no scenarios.")
    return df


def describe(row):
    """What this scenario will run with, before anything is submitted."""
    name = row["scenario_name"]
    print(f"\n{'='*60}")
    print(f"[Group] {exp_group(name)}")
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


def dispatch_all(rows, exp_num, runner, n_seeds, consensus):
    """The whole experiment phase, as ONE phase of ONE experiment.

    Each scenario is a named group appended to impact_long_shedders_#N, and
    the runner is given all of them at once -- one Slurm array for the phase.

    Until 8e each scenario was its own experiment, so this submitted five of
    them from five threads and divided SIMPLICITY_MAX_PARALLEL_SEEDED_
    SIMULATIONS_SLURM between them to keep the queue bounded. One submission
    needs neither: the cap applies to the single array, whole.
    """
    for row in rows:
        describe(row)

    settings_func = build_exp_phase_settings(rows, n_seeds, consensus)
    return run_experiment_phase_script(runner, exp_num, settings_func,
                                       exp_name(consensus), append=True)


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3: experiment-phase runs from the frozen "
                    "table, both rates applied as calibrated.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=30)
    parser.add_argument('--consensus', type=str, choices=list(CONSENSUS_MODES),
                        default=DEFAULT_CONSENSUS,
                        help=f"How the consensus distance is measured "
                            f"(default: {DEFAULT_CONSENSUS}). 'argmax' writes "
                            f"to the _amx experiment names, so both modes can "
                            f"share one --exp-num.")
    parser.add_argument('--only', type=str, default=None,
                        help="Optional: run only this scenario_name.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    set_experiment_output_env()

    table_path = default_table_path(args.exp_num, args.consensus)
    df = load_calibration_table(table_path)

    if args.only is not None:
        df = df[df["scenario_name"] == args.only]
        if df.empty:
            raise ValueError(f"Scenario '{args.only}' not found in the table.")

    print(f"\n[Runner] Dispatching {len(df)} scenario(s) from {table_path} "
          f"with consensus={args.consensus} -> "
          f"{exp_name(args.consensus)}_#{args.exp_num}, one group each")
    shared = {k: v for k, v in USER_FIXED_PARAMS.items() if k != "R"}
    print_fixed_params(shared, label="Shared parameters (R shown per scenario below)")

    groups = dispatch_all([row for _, row in df.iterrows()], args.exp_num,
                          args.runner, args.seeds, args.consensus)

    print(f"\n[Success] All scenarios dispatched: {groups}")


if __name__ == "__main__":
    main()
