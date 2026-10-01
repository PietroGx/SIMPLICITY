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
    CONSENSUS_MODES, prod_exp_name,
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


def dispatch_scenario(row, exp_num, runner, n_seeds, consensus="argmax"):
    name = row["scenario_name"]
    settings_func = build_exp_scenario_settings(row, n_seeds, consensus)
    prefix = prod_exp_name(consensus)

    print(f"\n{'='*60}")
    print(f"[Dispatch] {prefix}_{name}   (consensus: {consensus})")
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

    run_experiment_script(runner, exp_num, settings_func, f"{prefix}_{name}")


# Scenarios are independent: each submits its own Slurm array and each array's
# tasks share one cap. Run them one after another and only one scenario's seeds
# are ever in flight -- 30 at a time, five times over, against a cap of 200.
# Submitting them together puts every scenario's seeds in the queue at once and
# finishes the stage in roughly the time of its slowest scenario rather than the
# sum of all five.
#
# Threads, not processes: run_seeded_simulations is entirely subprocess calls
# and sleeps, so the GIL is never the constraint, and threads share the Data
# directory and signal files without any extra coordination.
def dispatch_all(rows, exp_num, runner, n_seeds, consensus, parallel=True):
    if runner != "slurm" or not parallel or len(rows) <= 1:
        for row in rows:
            dispatch_scenario(row, exp_num, runner, n_seeds, consensus)
        return

    import threading

    # Split the global cap between the scenarios rather than letting each take
    # it in full: five scenarios at the 200 default would otherwise put 1000
    # tasks in the queue. run_seeded_simulations reads this per call, so it has
    # to be set before any thread starts.
    key = "SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"
    cap = int(os.environ.get(key, 200))
    budget = max(1, cap // len(rows))
    print(f"\n[Runner] Submitting {len(rows)} scenarios together: "
          f"{budget} concurrent seeds each, {budget * len(rows)} of a {cap} cap.")

    errors = []
    lock = threading.Lock()

    def worker(row):
        try:
            dispatch_scenario(row, exp_num, runner, n_seeds, consensus)
        except BaseException as exc:            # noqa: BLE001 -- re-raised below
            with lock:
                errors.append((row["scenario_name"], exc))

    previous = os.environ.get(key)
    os.environ[key] = str(budget)
    try:
        threads = [threading.Thread(target=worker, args=(row,),
                                    name=f"scenario-{row['scenario_name']}")
                   for row in rows]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    finally:
        # restore it: a second dispatch_all in the same process would otherwise
        # divide the already-divided budget again
        if previous is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = previous

    if errors:
        # every scenario is reported, then the first failure is re-raised so
        # the pipeline still stops here rather than calibrating against data
        # that was never produced
        for name, exc in errors:
            print(f"[FAILED] scenario {name}: {type(exc).__name__}: {exc}")
        raise errors[0][1]


def main():
    parser = argparse.ArgumentParser(
        description="Unbound pipeline stage 3: production runs from the frozen "
                    "table, both rates applied as calibrated.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=30)
    parser.add_argument('--consensus', type=str, choices=list(CONSENSUS_MODES),
                        default='argmax',
                        help="How the consensus distance is measured. "
                            "'distribution' writes to the _dist experiment "
                            "names, so both pipelines can share one --exp-num.")
    parser.add_argument('--only', type=str, default=None,
                        help="Optional: run only this scenario_name.")
    parser.add_argument('--sequential', action='store_true',
                        help="Submit one scenario at a time and wait for each, "
                            "as before v2.4.49. The default submits all "
                            "scenarios together, splitting the concurrency cap "
                            "between them.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)

    table_path = default_table_path(args.exp_num)
    df = load_calibration_table(table_path)

    if args.only is not None:
        df = df[df["scenario_name"] == args.only]
        if df.empty:
            raise ValueError(f"Scenario '{args.only}' not found in the table.")

    print(f"\n[Runner] Dispatching {len(df)} scenario(s) from {table_path} "
          f"with consensus={args.consensus} -> {prod_exp_name(args.consensus)}_<scenario>")
    shared = {k: v for k, v in USER_FIXED_PARAMS.items() if k != "R"}
    print_fixed_params(shared, label="Shared parameters (R shown per scenario below)")

    dispatch_all([row for _, row in df.iterrows()], args.exp_num, args.runner,
                 args.seeds, args.consensus, parallel=not args.sequential)

    print(f"\n[Success] All scenarios dispatched.")


if __name__ == "__main__":
    main()
