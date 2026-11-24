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
# impact_long_shedders -- STAGE 2: STANDARD INDIVIDUALS ALONE
# ----------------------------------------------------------------------------
# A 0%-long-shedder population at production settings: sweeps the standard
# NSR, fits the global root-to-tip clock, inverts at target_osr_std. ONE
# calibrated standard rate, shared by every scenario.
#
# This is the whole point of the design. The deleted bound pipeline fit this
# rate per scenario INSIDE the mixed population, so the rate absorbs whatever
# divergence long shedders inject -- in run #7 that drove HIV_high to 9.3e-7,
# 150x below control. Here the standard cohort never sees a long shedder while
# being calibrated, so the rate means "the rate standard individuals actually
# mutate at" and the population clock is left free to rise in production.
#
# The sweep is independent of stage 1. Stage 1's result is only read at the
# end, to assemble the frozen table both stages feed to production.
# ============================================================================

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.phenotype.distance as dis

from experiment_script_runner import run_experiment_phase_script
from impact_long_shedders_config import (
    SETTINGS, check_stamp,
    SCENARIOS, EXP_NAME, NSR_RANGES, calib_std_group,
    setup_dir, TABLE_FILENAME, TABLE_COLUMNS,
    CAL2_FINAL_TIME, USER_FIXED_PARAMS,
    derive_scenario_params, lookup_long_nsr, build_cal2_settings,
    CONSENSUS_MODES, DEFAULT_CONSENSUS,
    add_slurm_resource_args, set_slurm_resource_env, print_fixed_params,
    set_calibration_output_env,
)
from long_nsr_calibration_plot import (read_calibrated_long_nsr,
                                       write_osr_per_seed,
                                       EXPECTED_SEED_SKIPS, record_seed_skip,
                                       report_seed_usage, merge_seed_skips,
                                       format_seed_skips,
                                       MAX_CALIBRATION_WORKERS,
                                       calibration_workers, seed_task_mapper,
                                       read_seed_outcomes)

_REF_LEN = len(dis.reference)


def _std_seed_osr_row(task):
    """One seeded simulation output dir's standard-cohort OSR row, run here or
    in a pool worker. Tagged exactly as cal_1's _seed_osr_row: ('row', dict),
    ('empty', None) for a seed that fails the gate, ('skip', exception name)
    for an expected read failure. Anything structural is re-raised naming its
    seeded dir, the traceback not surviving the pickle."""
    ssod, nsr_val, sim_name, min_seq, min_len = task
    try:
        final_time = om.read_final_time(ssod)
        seq_data = om.read_sequencing_data_regression(ssod)
        # Every individual is standard here; the filter is kept so this
        # measures the same quantity as the bound pipeline's.
        seq_data = seq_data[seq_data['individual_type'] == 'standard']
        if not (final_time >= min_len and len(seq_data) >= min_seq):
            return 'empty', None
        return 'row', {
            'nucleotide_substitution_rate': nsr_val,
            'observed_substitution_rate':
                er.tempest_regression(seq_data).coef_[0],
            'simulation': sim_name,
            'seed': os.path.basename(ssod),
            'n_points': len(seq_data),
            'muts_median': float(
                seq_data['Distance_from_root'].median()) * _REF_LEN,
            'muts_max': float(
                seq_data['Distance_from_root'].max()) * _REF_LEN,
        }
    except EXPECTED_SEED_SKIPS as exc:
        return 'skip', type(exc).__name__
    except Exception as exc:
        raise RuntimeError(f"{type(exc).__name__} in {ssod}: {exc}") from exc


def compute_standard_nsr(numbered, target_osr_std, model_type='exp',
                         min_seq=30, min_len=100, groups=None, workers=None):
    """Read back the standard-only sweep, fit OSR vs NSR over standard
    individuals, invert at target. Returns the calibrated standard NSR.

    `groups` scopes the read. It matters since 8e: the experiment also holds
    the long-calibration cohorts and the experiment-phase scenarios, and
    every one of them has a nucleotide_substitution_rate, so an unscoped
    listing would pool a 100%-long-shedder sweep and single-point scenario
    runs into this fit and still produce a healthy-looking R^2. None means
    the calib_std group the config names.

    `workers` is the process count for the per-seed extraction; None asks
    calibration_workers. It changes only how long this takes -- outcomes are
    read back in task order, so the fit is fed the identical rows.
    """
    wanted = list(groups) if groups is not None else [calib_std_group()]
    print(f"\n[Stage 2] Standard-only NSR sweep: {numbered}  groups={wanted}")

    all_rows = []
    skipped_all = {}
    sods = [sod for group in wanted
            for sod in dm.get_simulation_output_dirs(numbered, group=group)]
    if not sods:
        raise RuntimeError(
            f"{numbered} has no simulation output in group(s) {wanted}. "
            f"Groups on disk: {dm.get_groups(numbered)}.")
    # Flat task list over every (simulation output dir, seed), so one pool
    # covers the stage; `plan` holds the structure the outcomes are read back
    # with -- they arrive in task order, so each directory's slice is
    # consumed, reported and outlier-filtered in its turn, as before.
    plan = []
    tasks = []
    for sod in sods:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'nucleotide_substitution_rate')
        ssods = dm.get_seeded_simulation_output_dirs(sod)
        plan.append((sod, len(ssods)))
        tasks.extend((ssod, nsr_val, os.path.basename(sod), min_seq, min_len)
                     for ssod in ssods)

    with seed_task_mapper(calibration_workers(workers), len(tasks)) as mapper:
        outcomes = mapper(_std_seed_osr_row, tasks)
        for sod, n_seeds in plan:
            sod_rows, skips = read_seed_outcomes(outcomes, n_seeds)
            report_seed_usage('  [Stage 2]', sod, len(sod_rows), n_seeds,
                              skips)
            merge_seed_skips(skipped_all, skips)

            if sod_rows:
                all_rows.append(om.detect_sod_outliers(pd.DataFrame(sod_rows)))

    if not all_rows:
        raise RuntimeError(
            f"[Stage 2] No standard OSR data for {numbered} "
            f"(min_seq={min_seq}, min_len={min_len}): "
            f"{format_seed_skips(skipped_all)}.")

    clean_df = pd.concat(all_rows, ignore_index=True)
    write_osr_per_seed(numbered, clean_df, 'std')
    clean_df = clean_df[clean_df['is_outlier'] == 0]

    # The group name comes from the config, via `wanted` -- not a string typed
    # here, and NOT whatever is on disk. The fit-results filename carries it.
    #
    # This read dm.get_groups(numbered) until v2.4.103 and required exactly
    # one. That was true while cal_2 had an experiment to itself; since 8e the
    # experiment also holds cal_1's cohorts, so the stage died here -- AFTER
    # its 300 simulations had run -- on an experiment that was perfectly fine.
    if len(wanted) != 1:
        raise RuntimeError(
            f"[Stage 2] one fit needs one group, got {wanted}. "
            f"Groups on disk: {dm.get_groups(numbered)}.")
    fit_result = er.fit_observed_substitution_rate_regressor(
        clean_df, model_type,
        parameter_name='nucleotide_substitution_rate')
    om.write_fit_results_csv(numbered, model_type, fit_result,
                             experiment_group=wanted[0])
    nsr_std = float(er.compute_calibrated_parameter(
        model_type, fit_result, target_osr_std))
    print(f"          calibrated standard NSR = {nsr_std:.8f}")

    plt.figure(figsize=(10, 6))
    plt.scatter(clean_df['nucleotide_substitution_rate'],
                clean_df['observed_substitution_rate'],
                color='#1f77b4', alpha=0.15, s=10)
    x_vals = np.geomspace(clean_df['nucleotide_substitution_rate'].min(),
                          clean_df['nucleotide_substitution_rate'].max(), 100)
    plt.plot(x_vals, fit_result.eval(x=x_vals), color='#1f77b4', linewidth=2,
             label='standard-only fit')
    plt.plot(nsr_std, target_osr_std, marker='*', markersize=14,
             color='#1f77b4', markeredgecolor='black')
    plt.axhline(target_osr_std, color='black', linestyle='--', linewidth=1.5,
                label='Target OSR')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Standard NSR calibration (standard individuals alone)')
    plt.xlabel('Input Nucleotide Substitution Rate (NSR)')
    plt.ylabel('Observed Substitution Rate (OSR, standard individuals)')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f'{numbered}_std_nsr_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Calibration plot saved to: {plot_path}")

    return nsr_std


def write_frozen_table(stamp, std_nsr, long_nsr_by_group, args,
                       long_calib_exp, std_calib_exp):
    sp = sm.read_standard_parameters_values()
    target_dir = setup_dir(stamp, args.consensus)
    os.makedirs(target_dir, exist_ok=True)

    rows = []
    for scenario in SCENARIOS:
        frozen = derive_scenario_params(scenario, sp)
        long_nsr = (lookup_long_nsr(long_nsr_by_group, frozen["tau_3_long"],
                                    frozen["R_long"])
                    if frozen["is_long"] else None)
        rows.append({
            "scenario_name": frozen["scenario_name"],
            "tau_3_long": frozen["tau_3_long"],
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "susceptibility_long": frozen["susceptibility_long"],
            "R": args.R,
            "IH_virus_emergence_rate": args.ih_virus_emergence_rate,
            "R_long": frozen["R_long"],
            "nucleotide_substitution_rate_long": long_nsr,
            # one standard rate for every scenario -- nothing scenario-specific
            # is fitted once calibration stops compensating for the long side
            "nucleotide_substitution_rate": std_nsr,
            "target_osr_long": args.target_osr_long,
            "model_type": args.model,
            "long_calib_source": long_calib_exp,
            "std_calib_source": std_calib_exp,
        })

    table_path = os.path.join(target_dir, TABLE_FILENAME)
    pd.DataFrame(rows)[TABLE_COLUMNS].to_csv(table_path, index=False)
    print(f"\n[Success] Frozen calibration table written:\n  {table_path}")
    return table_path


def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: calibrate the standard NSR in a "
                    "0%-long-shedder population, then freeze the table.")
    parser.add_argument('--stamp', type=check_stamp, required=True)
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=SETTINGS['seeds']['cal'])
    parser.add_argument('--target-osr-std', type=float,
                        default=SETTINGS['targets']['osr_std'])
    parser.add_argument('--target-osr-long', type=float,
                        default=SETTINGS['targets']['osr_long'])
    parser.add_argument('--R', type=float, default=USER_FIXED_PARAMS['R'])
    parser.add_argument('--consensus', type=str, choices=list(CONSENSUS_MODES),
                        default=SETTINGS['consensus'],
                        help=f"Consensus mode to CALIBRATE under (default: "
                            f"{SETTINGS['consensus']}). Must match the mode the "
                            f"experiment phase will run.")
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'])
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    parser.add_argument('--workers', type=int, default=None,
                        help="Processes for the per-seed extraction "
                             "(default: SLURM_CPUS_PER_TASK, else the core "
                             f"count, capped at {MAX_CALIBRATION_WORKERS}). "
                             "1 runs it inline.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    # calibration reads four files; do not write the other three
    set_calibration_output_env()

    settings_func = build_cal2_settings(
        args.seeds, NSR_RANGES['cal_std'], args.R,
        args.ih_virus_emergence_rate, consensus=args.consensus)

    print(f"\n[Stage 2] Standard-only sweep (R={args.R}, no long shedders, "
         f"final_time={CAL2_FINAL_TIME}, consensus={args.consensus})")
    varying, fixed_params, n_seeds = settings_func()
    print_fixed_params(fixed_params)
    print(f"Seeds per grid point: {n_seeds}")
    # Over the groups, like bound cal_2's banner: 7a moved the sweep axis
    # inside the named group and this still read varying['nucleotide_
    # substitution_rate'], so the stage died here, before dispatching.
    for g in varying['_scenario_groups']:
        nsr_vals = g['nucleotide_substitution_rate']
        print(f"  {g['name']:14s} NSR grid: {len(nsr_vals)} pts "
             f"[{min(nsr_vals):.6f} .. {max(nsr_vals):.6f}]")

    numbered = f"{EXP_NAME}_#{args.stamp}"
    dm.mark_experiment_started(numbered)
    # append: cal_1 wrote this experiment, and its groups are finished. Only
    # the groups below are run.
    groups = run_experiment_phase_script(args.runner, args.stamp,
                                         settings_func, EXP_NAME, append=True)

    std_nsr = compute_standard_nsr(numbered, args.target_osr_std,
                                   model_type=args.model,
                                   min_seq=args.min_seq, min_len=args.min_len,
                                   groups=groups, workers=args.workers)

    long_nsr_by_group = read_calibrated_long_nsr(numbered, args.target_osr_long)

    write_frozen_table(args.stamp, std_nsr, long_nsr_by_group, args,
                       numbered, numbered)
    print("Next: run impact_long_shedders_exp.py against this table.")


if __name__ == "__main__":
    main()
