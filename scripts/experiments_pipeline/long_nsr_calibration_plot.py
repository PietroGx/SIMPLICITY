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
# the same fit from raw simulation output a second time.
#
# Nothing may write into a finished experiment, by two mechanisms that do not
# overlap. LAYERING scopes a refit: the fit computes and returns, the stage
# persists, and passing refit_out_dir puts plot_and_fit_long_nsr_calibration in
# refit mode -- figure to that directory, not one byte into the experiment.
# The LIFECYCLE MARKER seals a run: dir_manager.assert_experiment_writable,
# which the two write_* below call, refuses a completed experiment.
# ============================================================================

import contextlib
import multiprocessing as mp
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
import simplicity.phenotype.distance as dis

_REF_LEN = len(dis.reference)

# A per-seed read fails for two unrelated reasons and only one is a reason to
# skip: the repeat produced no output, or its output cannot be read/parsed.
# Anything structural -- a missing column, a type error -- is a defect of this
# code, and absorbing it reported "no data" for an experiment that had plenty.
EXPECTED_SEED_SKIPS = (OSError,                      # absent / unreadable file
                       pd.errors.EmptyDataError,     # file exists, no content
                       pd.errors.ParserError)        # truncated mid-write


def count_seed_skip(skips, name):
    skips[name] = skips.get(name, 0) + 1


def record_seed_skip(skips, exc):
    count_seed_skip(skips, type(exc).__name__)


def report_seed_usage(prefix, sod, n_used, n_seeds, skips):
    """One line per simulation output directory, printed whether or not
    anything was skipped: an empty result has to show its counts instead of
    reading as missing data."""
    why = ('; ' + ', '.join(f'{name} x{n}' for name, n in sorted(skips.items()))
           if skips else '')
    print(f"{prefix} {os.path.basename(sod)}: {n_used}/{n_seeds} seeds used, "
          f"{sum(skips.values())} skipped{why}")


def merge_seed_skips(total, skips):
    for name, n in skips.items():
        total[name] = total.get(name, 0) + n
    return total


def format_seed_skips(total):
    if not total:
        return 'no seed was skipped'
    return (f"{sum(total.values())} seed(s) skipped: "
            + ', '.join(f'{name} x{n}' for name, n in sorted(total.items())))


# ---------------------------------------------------------------------------
# Parallel per-seed extraction, shared with cal_2's standard-NSR loop.
#
# er.extract_ih_regression_data is ~99.8% of the per-seed cost and one seeded
# output dir never touches another, so the loop runs in a process pool --
# processes, not threads, the cost being CPU-bound Python. Three things the
# serial loop had for free must be rebuilt across the boundary: outcomes are
# collected BY POSITION, so the fit sees the same seeds in the same order at
# any worker count; a skip is RETURNED and counted in the parent, a worker's
# counter dying with it; and a structural error is re-raised naming its
# seeded dir, the traceback not surviving the pickle.
# ---------------------------------------------------------------------------

# Measured on run #3's calib_long__edge_case: one seed of its top NSR grid
# point carries ~180k intra-host lineages and peaks at 1.6 GB, and a worker
# holds one seed at a time. 4 x 1.6 GB + the parent is ~6.5 GB, inside the 8G
# an analysis job asks for (submit_sanity_plot.sh); 8 would be ~13 GB. A node
# with more memory takes an explicit --workers, honoured above this cap.
MAX_CALIBRATION_WORKERS = 4


def calibration_workers(override=None):
    """How many processes the per-seed loops may use. 1 means no pool at all.

    `override`, else SIMPLICITY_CALIBRATION_WORKERS, is taken as given: the
    caller knows its node, so the cap only warns. Auto-detection prefers
    SLURM_CPUS_PER_TASK and is capped, because os.cpu_count() on a login node
    reports the whole machine and each worker holds one seed's intra-host
    points at once -- 1.6 GB of them on edge_case data.
    """
    explicit = (override if override is not None
                else os.environ.get('SIMPLICITY_CALIBRATION_WORKERS'))
    if explicit not in (None, ''):
        n = max(1, int(explicit))
        if n > MAX_CALIBRATION_WORKERS:
            print(f"  [workers] {n} asked for explicitly, above the "
                  f"{MAX_CALIBRATION_WORKERS} cap -- honoured. Budget "
                  f"~1.6 GB per worker: edge_case seeds peak there.")
        return n
    detected = (int(os.environ.get('SLURM_CPUS_PER_TASK') or 0)
                or os.cpu_count() or 1)
    return max(1, min(detected, MAX_CALIBRATION_WORKERS))


@contextlib.contextmanager
def seed_task_mapper(workers, n_tasks):
    """A `map(fn, tasks)` over `workers` processes that yields IN TASK ORDER.

    Ordered, never as-completed: outlier detection is per simulation output
    dir and the fit is fed the concatenation, so the same seeds have to arrive
    in the same places whatever the worker count.

    A context manager so ONE pool covers every group of a run, rather than a
    fork per simulation output dir.
    """
    n = min(workers, n_tasks) if n_tasks else 1
    if n <= 1:
        yield lambda fn, tasks: iter(map(fn, tasks))
        return
    # fork, so a worker resolves om./er. through the module it inherited and a
    # caller that monkeypatched either still measures what it installed.
    with mp.get_context('fork').Pool(processes=n) as pool:
        print(f"  [workers] {n} process(es) over {n_tasks} seed task(s)")
        yield lambda fn, tasks: pool.imap(fn, tasks)


def read_seed_outcomes(outcomes, n_seeds):
    """The next `n_seeds` outcomes of one simulation output dir: its rows, and
    its skips counted by exception name."""
    rows, skips = [], {}
    for _ in range(n_seeds):
        kind, payload = next(outcomes)
        if kind == 'row':
            rows.append(payload)
        elif kind == 'skip':
            count_seed_skip(skips, payload)
    return rows, skips


def _seed_osr_row(task):
    """One seeded simulation output dir's intra-host OSR row.

    ('row', dict) when the seed is used; ('empty', None) when it fails the
    final_time/sequence gate or carries no intra-host lineage -- neither used
    nor skipped, as the serial loop counted it; ('skip', exception name) for
    an EXPECTED_SEED_SKIPS read failure.
    """
    ssod, group_name, nsr_val, sim_name, min_seq, min_len = task
    try:
        final_time = om.read_final_time(ssod)
        seq_data = om.read_sequencing_data_regression(ssod)
        # Same final_time/sequence-count gate as everywhere else that fits an
        # OSR regressor. What gets fit is NOT seq_data (population-pooled,
        # global-reference, absolute-time -- see evolutionary_rate.py's
        # extract_ih_regression_data docstring): cal_1 calibrates against each
        # long shedder's OWN intra-host clock, not the population
        # root-to-tip rate.
        if not (final_time >= min_len and len(seq_data) >= min_seq):
            return 'empty', None
        ih_points = er.extract_ih_regression_data(ssod)
        if ih_points.empty:
            return 'empty', None
        return 'row', {
            'group': group_name,
            'nucleotide_substitution_rate_long': nsr_val,
            'observed_substitution_rate':
                er.tempest_regression(ih_points).coef_[0],
            'simulation': sim_name,
            'seed': os.path.basename(ssod),
            'n_points': len(ih_points),
            'muts_median': float(
                ih_points['Distance_from_root'].median()) * _REF_LEN,
            'muts_max': float(
                ih_points['Distance_from_root'].max()) * _REF_LEN,
        }
    except EXPECTED_SEED_SKIPS as exc:
        return 'skip', type(exc).__name__
    except Exception as exc:
        raise RuntimeError(f"{type(exc).__name__} in {ssod}: {exc}") from exc


# Imported as a plain sibling by impact_long_shedders_cal_1.py and, from
# scripts/checks/, by check_calibration.py. Ensure this file's own directory is
# on sys.path regardless of the caller's cwd, so the sibling import below
# always resolves.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Constants only: importing the scenario LIST from here is what once put one
# pipeline's label on another's group (see _cal1_groups).
from impact_long_shedders_config import (TAU_ROUND, DEFAULT_COLORS,
                                         groups_of_phase)


def cal1_group_names(experiment_name):
    """The groups of `experiment_name` that belong to the long calibration.

    Stage 8e put every phase into ONE experiment, so dm.get_groups now also
    returns the standard-NSR sweep and the experiment-phase scenarios. Those
    carry a tau_3_long and an R_long like everything else, so handing them to
    the fit below would have it treat each of them as another long-shedder
    cohort -- a fit over rows that never swept NSR_long, reported under a
    cohort name. Scope by the phase prefix the config stamps on the name.

    An experiment with no prefixed group at all is one that predates the
    convention or was built outside it; all of its groups are returned, which
    is what this function did for every caller before 8e.
    """
    named = groups_of_phase(dm.get_groups(experiment_name), 'calib_long')
    return named or dm.get_groups(experiment_name)


def _cal1_groups(experiment_name, groups=None):
    """Every long-calibration group, in declared order, with the
    (tau_3_long, R_long) cohort it covers and a plot colour.

    `groups` names them explicitly; None asks cal1_group_names.

    The cohort comes from the group's OWN simulation records, and a group
    holding more than one pair raises rather than silently taking the first:
    a cal_1 group IS one long-shedder cohort swept over NSR_long, and anything
    else means the config and this fit disagree about what a group is.

    This used to be a {(tau, R_long) -> name} map, with the partition itself
    done by re-reading those two floats off every simulation output directory.
    Before that the name was re-derived by matching the pair against a scenario
    list passed in by the caller, defaulting to the deleted pipeline's
    SCENARIOS -- edge_case was not in that list, so its group fell through to a
    raw 'tau=350.23,R_long=1.1' and that string landed in the saved
    fit-results filename.

    Names are joined with '+' by the config, never '/', which
    write_fit_results_csv's file-saving reads as a directory separator.
    """
    wanted = (list(groups) if groups is not None
              else cal1_group_names(experiment_name))
    wanted_set = set(wanted)

    # Scoped to `wanted`: one 8e experiment also holds cal_2's standard-only
    # sweep, whose records have no tau_3_long to subscript.
    present = set()
    pairs = {}
    for record in sm.read_simulations(experiment_name):
        present.add(record['group'])
        if record['group'] not in wanted_set:
            continue
        params = record['parameters']
        pairs.setdefault(record['group'], set()).add(
            (round(float(params['tau_3_long']), TAU_ROUND),
             round(float(params['R_long']), TAU_ROUND)))

    missing = [name for name in wanted if name not in pairs]
    if missing:
        raise RuntimeError(
            f"{experiment_name} has no simulations in group(s) {missing}. "
            f"Present: {sorted(present)}.")

    color_cycle = iter(DEFAULT_COLORS)
    groups = []
    for name in wanted:
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


CALIBRATED_LONG_NSR_COLUMNS = ['group', 'tau_3_long', 'R_long',
                               'nucleotide_substitution_rate_long',
                               'target_osr_long', 'model_type']


def _cohort_key(tau_3_long, r_long):
    """The (tau_3_long, R_long) cohort key, rounded exactly as
    read_calibrated_long_nsr and config.lookup_long_nsr round it."""
    return (round(float(tau_3_long), TAU_ROUND),
            round(float(r_long), TAU_ROUND))


def _preserved_cohorts(path, rows, target_osr_long, model_type):
    """The cohorts already in `path` that this write does not cover, keyed as
    `rows` is. A preserved row calibrated against another target or fitted
    with another model is refused: one table cannot hold two calibrations,
    and read_calibrated_long_nsr raises on any target mismatch anyway. A
    write covering EVERY cohort preserves nothing, so re-running the whole
    stage at a new target stays possible."""
    if not os.path.isfile(path):
        return {}
    kept = {}
    # round_trip: the default parser drops the last digits of a float, so a
    # preserved row would come back rewritten.
    for _, row in pd.read_csv(path, float_precision='round_trip').iterrows():
        key = _cohort_key(row['tau_3_long'], row['R_long'])
        if key in rows:
            continue
        for field, writing in (('target_osr_long', target_osr_long),
                               ('model_type', model_type)):
            if row[field] != writing:
                raise RuntimeError(
                    f"{path} holds cohort {row['group']!r} "
                    f"{key} with {field}={row[field]!r}, and this write "
                    f"carries {field}={writing!r}. Merging them would leave "
                    f"one table holding two different calibrations. Re-run "
                    f"the whole stage at one {field}, or write elsewhere.")
        kept[key] = {name: row[name] for name in CALIBRATED_LONG_NSR_COLUMNS}
    return kept


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

    MERGES on that key: cohorts in `results` are replaced, cohorts already in
    the file and not in `results` are kept. A stage scoped to one group used
    to rebuild the file from that group alone, leaving one row where cal_2
    reads three -- and the seal does not cover it, since that happens while
    the experiment is still legitimately writable.
    """
    dm.assert_experiment_writable(experiment_name)
    path = _calibrated_long_nsr_path(experiment_name)
    group_by_key = group_by_key or {}
    rows = {
        _cohort_key(tau, r_long):
            {'group': group_by_key.get((tau, r_long), ''),
             'tau_3_long': tau, 'R_long': r_long,
             'nucleotide_substitution_rate_long': nsr,
             'target_osr_long': target_osr_long, 'model_type': model_type}
        for (tau, r_long), nsr in results.items()
    }
    rows.update(_preserved_cohorts(path, rows, target_osr_long, model_type))

    ordered = [rows[key] for key in sorted(rows)]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(ordered)[CALIBRATED_LONG_NSR_COLUMNS].to_csv(path, index=False)
    print(f"Saved calibrated long NSR for {len(ordered)} cohort(s) to: {path}")


def _osr_per_seed_path(experiment_name, stage):
    return os.path.join(dm.get_experiment_fit_result_dir(experiment_name),
                        f'{experiment_name}_osr_per_seed_{stage}.csv')


def _merge_osr_per_seed(old, new, path):
    """`new`'s groups replace `old`'s; `old`'s other groups are kept, in the
    order the file already had them."""
    if list(old.columns) != list(new.columns):
        raise RuntimeError(
            f"{path} has columns {list(old.columns)} and this write carries "
            f"{list(new.columns)}: the file was written by another version of "
            f"this stage. Merging two layouts into one table would put two "
            f"different measurements in one column.")
    if 'group' not in old.columns:
        return new
    replaced = set(new['group'])
    order = list(dict.fromkeys(list(old['group']) + list(new['group'])))
    return pd.concat(
        [new[new['group'] == g] if g in replaced else old[old['group'] == g]
         for g in order], ignore_index=True)


def write_osr_per_seed(experiment_name, df, stage):
    """Persist the per-seed OSR frame the calibration fit is derived from --
    one row per seed, outliers included and flagged.

    Only the fitted A/B and the inverted NSR were kept before, so no run on
    disk let its own calibration be re-examined. The mutation counts are here
    because a median of 0 means that grid point is reporting the smallest
    distance the reference can express, not a rate (see BACKLOG).

    MERGES on `group`, like the calibrated table: a scoped run replaces its
    own groups and leaves the rest. The 'std' frame has no group column and
    only ever covers cal_2's single standard group, so it is replaced whole.
    """
    dm.assert_experiment_writable(experiment_name)
    path = _osr_per_seed_path(experiment_name, stage)
    if 'group' in df.columns and os.path.isfile(path):
        df = _merge_osr_per_seed(
            pd.read_csv(path, float_precision='round_trip'), df, path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    print(f"Saved per-seed OSR frame ({len(df)} rows) to: {path}")


def read_calibrated_long_nsr(experiment_name, target_osr_long):
    """Read back Stage 1's calibrated long NSR per (tau_3_long, R_long)
    group, written by write_calibrated_long_nsr. Raises rather than
    silently guessing if the file is missing (cal_1 never ran under this
    stamp) or if target_osr_long doesn't match what the file was
    actually calibrated against -- a caller wanting a different target
    must re-run cal_1.py with it, not just pass a different flag here and
    reuse a stale calibration.

    Returns {(tau_3_long, R_long): calibrated_nsr}.
    """
    path = _calibrated_long_nsr_path(experiment_name)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"No calibrated long NSR file at {path}. Run "
            f"impact_long_shedders_cal_1.py for this stamp first.")

    # round_trip: the default parser drops the last digits, so the rate cal_2
    # hands production would not be the rate cal_1 calibrated.
    df = pd.read_csv(path, float_precision='round_trip')
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
                                      groups=None, refit_out_dir=None,
                                      workers=None):
    """
    Extracts per-seed OSR for every (tau_3_long, R_long, NSR) grid point of
    `experiment_name`, outlier-filters per grid point, fits a regressor per
    DECLARED GROUP, inverts at `target_osr_long`, and saves a multi-curve
    calibration figure into the experiment's 05_Plots folder.

    refit_out_dir switches the whole function from calibration stage to refit
    tool: the figure goes there and NOTHING is written into the experiment --
    no fit-results CSV, no calibrated-NSR table, no per-seed frame. One
    argument rather than a flag per write, so a half-persisting refit of a
    finished experiment cannot be expressed.

    One group is one long-shedder cohort: two scenarios can share a duration
    but differ in R_long, so the cohort is the (tau_3_long, R_long) pair and
    never tau alone. Which simulations belong to it is read from the record,
    not recomputed from their parameter values -- see _cal1_groups.

    Group names for the legend and the saved fit-results filename come from the
    experiment's own settings, where the config that defined the group wrote
    them. There used to be a `scenarios` argument here for the caller to say
    which scenario list to match pairs against -- needed only because this
    module re-derived the names instead of reading them.

    `workers` is the process count for the per-seed extraction; None asks
    calibration_workers. It changes only how long this takes -- outcomes are
    read back in task order, so the fit is fed the identical rows.

    Returns {(tau_3_long, R_long): calibrated_nsr}.
    """
    groups = _cal1_groups(experiment_name, groups)
    persist = refit_out_dir is None

    print(f"--- Analyzing long-shedder NSR calibration: {experiment_name} ---")
    print(f"    groups: {[g['name'] for g in groups]}")

    # Partitioned by the group each simulation was ASSIGNED to, not by reading
    # (tau_3_long, R_long) back off every output directory and bucketing on the
    # float pair. The two agree while every group carries a distinct rounded
    # pair, which is what unique_long_tau_r_long_pairs builds; they stop
    # agreeing the moment one does not, and nothing would say so.
    # One flat task list over every (simulation output dir, seed) of every
    # group, so one pool covers the whole run and no directory waits on the
    # previous one's slowest seed. `plan` holds the structure the outcomes are
    # read back with -- they arrive in task order, so each directory's slice is
    # consumed, reported and outlier-filtered in its turn, exactly as before.
    plan = []
    tasks = []
    for group in groups:
        for sod in dm.get_simulation_output_dirs(experiment_name,
                                                 group=group['name']):
            nsr_val = sm.get_parameter_value_from_simulation_output_dir(
                sod, 'nucleotide_substitution_rate_long')
            ssods = dm.get_seeded_simulation_output_dirs(sod)
            plan.append((group['name'], sod, len(ssods)))
            tasks.extend((ssod, group['name'], nsr_val, os.path.basename(sod),
                          min_seq, min_len) for ssod in ssods)

    group_frames = {}
    audit_frames = []
    skipped_all = {}
    with seed_task_mapper(calibration_workers(workers), len(tasks)) as mapper:
        outcomes = mapper(_seed_osr_row, tasks)
        for group_name, sod, n_seeds in plan:
            sod_rows, skips = read_seed_outcomes(outcomes, n_seeds)
            report_seed_usage(f"  [{group_name}]", sod, len(sod_rows),
                              n_seeds, skips)
            merge_seed_skips(skipped_all, skips)

            if sod_rows:
                sod_df = om.detect_sod_outliers(pd.DataFrame(sod_rows))
                audit_frames.append(sod_df)
                group_frames.setdefault(group_name, []).append(
                    sod_df[sod_df['is_outlier'] == 0])

    per_group = {name: pd.concat(frames, ignore_index=True)
                 for name, frames in group_frames.items()}

    if not per_group:
        raise RuntimeError(
            f"No valid long-calibration data found in {experiment_name}: "
            f"{format_seed_skips(skipped_all)}. Check if simulations finished "
            f"and lengths/seqs meet minimums.")

    if audit_frames and persist:
        write_osr_per_seed(experiment_name,
                           pd.concat(audit_frames, ignore_index=True), 'long')

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
            group_df,
            model_type=model_type,
            parameter_name='nucleotide_substitution_rate_long',
        )
        if persist:
            om.write_fit_results_csv(experiment_name, model_type, fit_result,
                                     experiment_group=label)

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

    plots_dir = (dm.get_experiment_plots_dir(experiment_name) if persist
                 else refit_out_dir)
    os.makedirs(plots_dir, exist_ok=True)
    plot_path = os.path.join(
        plots_dir, f'{experiment_name}_long_nsr_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"\nCalibration plot saved to: {plot_path}")

    if persist:
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
    # Same refit tool as check_calibration.py, so the same rule: a rerun
    # against an experiment on disk writes its figure outside it.
    parser.add_argument('--out', type=str, required=True,
                        help="Directory for the figure. Required: a rerun "
                             "writes nothing into the experiment.")
    parser.add_argument('--workers', type=int, default=None,
                        help="Processes for the per-seed extraction "
                             "(default: SLURM_CPUS_PER_TASK, else the core "
                             f"count, capped at {MAX_CALIBRATION_WORKERS}). "
                             "1 runs it inline.")
    args = parser.parse_args()

    try:
        plot_and_fit_long_nsr_calibration(
            args.experiment_name, args.target_osr, args.model, args.min_seq,
            args.min_len, refit_out_dir=args.out, workers=args.workers)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
