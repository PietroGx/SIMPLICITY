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
# NSR SENSITIVITY ANALYSIS -- how does the calibrated NSR_long depend on
# infection duration and on k_v?
# ----------------------------------------------------------------------------
# Motivation. Unbound run #1 calibrated NSR_long per duration group and got a
# monotone rise: 5.56e-4 (63 d), 1.07e-3 (109 d), 1.70e-3 (365 d). Fit quality
# rose with it (R2 0.906 -> 0.967 -> 0.996; B 0.574 -> 0.730 -> 0.904). Two
# explanations fit that pattern and they imply opposite things for the paper:
#
#   QUANTISATION (artifact). Hamming distance is an integer count. At the
#   target rate a 63-day infection accumulates <2 substitutions in total, so
#   most intra-host points sit at 0 or 1. A through-origin slope fitted to
#   mostly-zero integers is biased low, so the inversion demands a HIGHER NSR
#   to reach the target. If this is the cause, the three NSR_long values are a
#   measurement artifact and the underlying rate is duration-independent.
#
#   REAL SUBLINEARITY. Something in the model genuinely slows per-day
#   accumulation in long infections -- the IH_lineages_max cap, or lineage
#   turnover via k_v. If so, the per-duration calibration is correct and long
#   infections really do mutate more slowly per day.
#
# k_v is in the grid because it plausibly drives BOTH: sub_events draws
# Poisson(NSR * dt * L * IH_lineages), so mutation opportunity scales with the
# number of intra-host lineages, and k_v (rate k_v * infected) is what adds
# them. Duration and k_v may be two routes to one variable -- lineage count.
#
# ----------------------------------------------------------------------------
# TWO STAGES, so figures are replottable without re-simulating
# ----------------------------------------------------------------------------
# TWO COMMANDS, not three.
#
#   --run    ONE command for the whole HPC job: simulate, extract the RAW
#            per-point intra-host data to CSV, fit, draw every figure, and
#            compress the full output via scripts/package_simplicity_data.sh.
#            It finishes by printing exactly where everything landed.
#
#   --plot   replot, later, from the datasets alone -- no simulation output
#            needed, so it runs in seconds on a laptop against the small
#            artifacts zip downloaded from the HPC (unzip it and point
#            --dataset-dir at it). Changing the fit model or a figure never
#            costs a rerun.
#
# The run writes TWO things worth knowing apart:
#   * the ARTIFACTS ZIP  -- small: log, figures, recap, and the tidy datasets.
#                           This is the one to download, share, and replot from.
#   * the FULL ARCHIVE   -- large: every simulation output, pixz-compressed,
#                           for the record.
#
# Methodology is NOT re-implemented: the intra-host clock comes from
# simplicity.tuning.evolutionary_rate.extract_ih_regression_data and the
# isolated 100%-long-shedder context from
# impact_long_shedders_config.CAL1_ISOLATED_FIXED_PARAMS, exactly as cal_1
# uses them. A diagnostic that drifts from the stage it diagnoses gives
# confident numbers about a different quantity (see BACKLOG).
# ============================================================================

import os
import sys
import glob
import json
import zipfile
import argparse
import subprocess
from datetime import datetime

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
                                'experiments'))
from impact_long_shedders_config import (
    CAL1_ISOLATED_FIXED_PARAMS, USER_FIXED_PARAMS, phase_offset,
    add_slurm_resource_args, set_slurm_resource_env,
)

EXP_NAME = "nsr_sensitivity"

# =============================================================================
# GRID
# =============================================================================
# Clinical infection durations (days). The three production values plus two
# fillers, so duration is sampled as a CURVE rather than three points.
DURATIONS = [63.0, 109.0, 180.0, 270.0, 365.0]

# Intra-host lineage emergence rate. Log-spaced around production's 0.1.
K_V_VALUES = [0.001, 0.01, 0.1]

# NSR_long sweep, shared by every cell.
NSR_MIN, NSR_MAX, NSR_STEPS = 1e-05, 1e-02, 6

TARGET_OSR_LONG = 0.00205
POPULATION_SIZE = 1000
INFECTED_START = 10

# Simulated window: at least a year, and at least 3 infection cycles, rounded
# up to whole years. 63 d and 109 d -> 365; 180 d -> 730; 270 d and 365 d ->
# 1095. Long infections need the room to complete cycles at all; short ones
# still get a full year so no cell is measured over a stub.
WINDOW_BASE_DAYS = 365
MIN_CYCLES = 3

# Local development preset: enough to exercise every path, not to answer
# anything. 2 durations x 2 k_v x 3 NSR x 3 seeds = 36 simulations.
QUICK_DURATIONS = [63.0, 180.0]
QUICK_K_V = [0.001, 0.1]
QUICK_NSR_STEPS = 3
QUICK_POPULATION = 300
QUICK_INFECTED_START = 10

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}

# --- figure palette -----------------------------------------------------------
# Both series axes are ORDERED magnitudes, not categories, so each gets a
# single-hue ordinal ramp (light -> dark) rather than categorical hues. Both
# ramps are validated: monotone lightness, adjacent dL >= 0.06, light end
# clears the surface, single hue.
DURATION_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
K_V_RAMP = ["#86b6ef", "#2a78d6", "#104281"]
INK, INK_MUTED, GRID = "#1a1a19", "#5c5c58", "#dcdcd8"
ACCENT = "#eb6834"   # target lines only -- never a series


def set_paper_rcparams():
    plt.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 300,
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": INK_MUTED, "axes.linewidth": 0.6,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
        "figure.facecolor": "white", "savefig.bbox": "tight",
    })


def ramp_for(values, ramp):
    """Map an ordered list of values onto an ordinal ramp, low -> light."""
    if len(values) == 1:
        return {values[0]: ramp[len(ramp) // 2]}
    idx = np.linspace(0, len(ramp) - 1, len(values)).round().astype(int)
    return {v: ramp[i] for v, i in zip(values, idx)}


# =============================================================================
# GRID HELPERS
# =============================================================================
def window_for(inf_duration):
    """Simulated window in days: >= 1 year and >= MIN_CYCLES infection cycles,
    rounded up to whole years."""
    needed = MIN_CYCLES * inf_duration
    years = max(1, int(np.ceil(needed / WINDOW_BASE_DAYS)))
    return years * WINDOW_BASE_DAYS


def tau_3_long_for(inf_duration, sp):
    """Model parameter from the clinical duration, same derivation the
    pipeline uses."""
    return float(inf_duration) - phase_offset(sp)


def grid_cells(durations, k_vs, sp):
    """(inf_duration, tau_3_long, k_v, window) for every cell."""
    return [(d, round(tau_3_long_for(d, sp), 3), k, window_for(d))
            for d in durations for k in k_vs]


def nsr_nodes(steps):
    return np.geomspace(NSR_MIN, NSR_MAX, steps).tolist()


# =============================================================================
# STAGE: RUN
# =============================================================================
def build_settings(durations, k_vs, nsr_values, n_seeds, population, infected_start, sp):
    """One scenario group per (duration, k_v) cell, each sweeping the same
    NSR_long grid, each with its own tau_3_long, k_v and window.

    Base context is CAL1_ISOLATED_FIXED_PARAMS -- the same 100%-long-shedder
    isolated population cal_1 calibrates in -- so this measures the quantity
    the pipeline actually calibrates.
    """
    fixed_params = CAL1_ISOLATED_FIXED_PARAMS.copy()
    fixed_params.update({
        'population_size': population,
        'infected_individuals_at_start': infected_start,
    })
    fixed_params.pop('final_time', None)   # per-group below

    scenario_groups = [
        {
            'nucleotide_substitution_rate_long': list(nsr_values),
            'tau_3_long': tau,
            'R_long': 1.1,
            'IH_virus_emergence_rate': k_v,
            'final_time': window,
        }
        for _, tau, k_v, window in grid_cells(durations, k_vs, sp)
    ]

    def make_settings():
        return ({'_scenario_groups': scenario_groups}, fixed_params, n_seeds)

    return make_settings


def extract_datasets(numbered, sp, min_seq=30, min_len=100):
    """Walk the simulation output once and write the RAW data two ways.

    points: one row per intra-host lineage observation -- the (time since own
            infection, distance from own founder) pairs the clock is fitted
            to, via extract_ih_regression_data.
    hosts:  one row per long-shedder host -- lineage count and realised
            infection length, for the lineage-count hypothesis.

    Neither is summarised or fitted here; that is --plot's job, so the fit
    model can change without re-simulating.
    """
    point_rows, host_rows = [], []
    kept = dropped = 0

    for sod in dm.get_simulation_output_dirs(numbered):
        def p(name):
            return sm.get_parameter_value_from_simulation_output_dir(sod, name)

        tau = float(p('tau_3_long'))
        cell = {
            'inf_duration': round(tau + phase_offset(sp), 3),
            'tau_3_long': round(tau, 3),
            'k_v': float(p('IH_virus_emergence_rate')),
            'nsr_long': float(p('nucleotide_substitution_rate_long')),
            'final_time_setting': float(p('final_time')),
        }

        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            seed = dm.get_seed_from_SSOD(ssod)
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
            except Exception:
                dropped += 1
                continue
            if final_time < min_len or len(seq_data) < min_seq:
                dropped += 1
                continue

            try:
                ih = er.extract_ih_regression_data(ssod)
            except Exception:
                dropped += 1
                continue
            if ih.empty:
                dropped += 1
                continue
            kept += 1

            for x, y in zip(ih['Sequencing_time'].to_numpy(),
                            ih['Distance_from_root'].to_numpy()):
                point_rows.append({**cell, 'seed': seed,
                                   'final_time': final_time,
                                   't_since_infection': float(x),
                                   'distance': float(y)})

            try:
                ind = om.read_individuals_data(ssod)
                ind = ind[ind['type'] == 'long_shedder']
                for _, row in ind.iterrows():
                    t_inf = row.get('t_infection')
                    t_end = row.get('t_not_infected')
                    host_rows.append({
                        **cell, 'seed': seed,
                        'n_lineages': int(row.get('IH_lineages_number') or 0),
                        'n_unique_lineages': int(row.get('IH_unique_lineages_number') or 0),
                        'realised_duration':
                            float(t_end) - float(t_inf)
                            if pd.notna(t_inf) and pd.notna(t_end) else np.nan,
                    })
            except Exception:
                pass

    print(f"  seeds kept: {kept}   dropped (gate or no data): {dropped}")
    return pd.DataFrame(point_rows), pd.DataFrame(host_rows)


def stage_run(args, sp):
    durations = QUICK_DURATIONS if args.quick else DURATIONS
    k_vs = QUICK_K_V if args.quick else K_V_VALUES
    steps = QUICK_NSR_STEPS if args.quick else NSR_STEPS
    population = QUICK_POPULATION if args.quick else POPULATION_SIZE
    infected_start = QUICK_INFECTED_START if args.quick else INFECTED_START
    nsr_values = nsr_nodes(steps)

    cells = grid_cells(durations, k_vs, sp)
    n_sims = len(cells) * steps * args.n_seeds

    print("\n=========================================================")
    print(f" NSR SENSITIVITY -- run #{args.run_n}{'  [QUICK]' if args.quick else ''}")
    print(f" durations (d): {durations}")
    print(f" k_v          : {k_vs}")
    print(f" NSR_long     : {steps} nodes "
          f"[{min(nsr_values):.2e} .. {max(nsr_values):.2e}]")
    print(f" population={population}  infected_at_start={infected_start}  "
          f"seeds={args.n_seeds}")
    print(f" cells={len(cells)}  simulations={n_sims}")
    print(" windows per duration:")
    for d in durations:
        print(f"   {d:6.0f} d infection -> final_time {window_for(d)}")
    print("=========================================================\n")

    settings_func = build_settings(durations, k_vs, nsr_values, args.n_seeds,
                                   population, infected_start, sp)
    numbered = f"{EXP_NAME}_#{args.run_n}"
    run_experiment(numbered, settings_func,
                   simplicity_runner=RUNNERS[args.runner], archive_experiment=False)

    print(f"\n[extract] reading {numbered} into tidy datasets...")
    points, hosts = extract_datasets(numbered, sp,
                                     min_seq=args.min_seq, min_len=args.min_len)
    out_dir = dataset_dir(args.run_n)
    os.makedirs(out_dir, exist_ok=True)

    points_path = os.path.join(out_dir, "intra_host_points.csv.gz")
    hosts_path = os.path.join(out_dir, "long_shedder_hosts.csv.gz")
    points.to_csv(points_path, index=False, compression="gzip")
    hosts.to_csv(hosts_path, index=False, compression="gzip")

    meta = {
        "run_n": args.run_n, "quick": bool(args.quick),
        "experiment": numbered,
        "durations": durations, "k_v": k_vs,
        "nsr_nodes": nsr_values, "n_seeds": args.n_seeds,
        "population_size": population,
        "infected_individuals_at_start": infected_start,
        "windows": {str(d): window_for(d) for d in durations},
        "target_osr_long": TARGET_OSR_LONG,
        "min_seq": args.min_seq, "min_len": args.min_len,
        "written": datetime.now().isoformat(timespec="seconds"),
    }
    with open(os.path.join(out_dir, "metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"  points: {len(points):>8d} rows -> {points_path}")
    print(f"  hosts : {len(hosts):>8d} rows -> {hosts_path}")

    if points.empty:
        raise SystemExit(
            "[FAILED] no intra-host observations were extracted -- nothing to "
            "plot. Check the simulation-health output above.")

    # figures + artifacts zip, from the datasets just written
    print("\n[plot] building figures from the datasets...")
    zip_path = stage_plot(args, dataset_override=out_dir)

    archive_paths = []
    if not args.no_package:
        print("\n[package] compressing the full simulation output...")
        archive_paths = package_raw_output(args.run_n, args.slurm)

    report_locations(args.run_n, zip_path, out_dir, archive_paths)
    return points_path


def report_locations(run_n, zip_path, dataset_dir_path, archive_paths):
    """Final block: say plainly where everything is, so it can be found and
    downloaded without hunting."""
    bar = "=" * 74
    print(f"\n{bar}\n WHERE THE OUTPUT IS\n{bar}")
    print("\n  DOWNLOAD THIS -- artifacts (log, figures, recap, tidy datasets):")
    print(f"    {os.path.abspath(zip_path)}")
    if os.path.isfile(zip_path):
        print(f"    ({os.path.getsize(zip_path) / 1e6:.1f} MB)")
    print("\n  Replot from it later, anywhere:")
    print(f"    unzip -d <dir> {os.path.basename(zip_path)}")
    print(f"    python scripts/nsr_sensitivity_analysis.py --plot "
          f"--run-n {run_n} --dataset-dir <dir>")
    print("\n  Datasets on this machine (what --plot reads by default):")
    print(f"    {os.path.abspath(dataset_dir_path)}")
    if archive_paths:
        total = sum(os.path.getsize(p) for p in archive_paths if os.path.isfile(p))
        print(f"\n  Full compressed simulation output ({total / 1e6:.0f} MB, "
              f"{len(archive_paths)} part(s)) -- for the record, not needed to replot:")
        for p in archive_paths:
            print(f"    {os.path.abspath(p)}")
    print(f"\n{bar}")


# =============================================================================
# STAGE: PLOT -- reads only the datasets
# =============================================================================
def dataset_dir(run_n):
    return os.path.join(dm.get_data_dir(), "nsr_sensitivity", f"run_{run_n}")


def load_datasets(run_n, dataset_override=None):
    d = dataset_override or dataset_dir(run_n)
    points_path = os.path.join(d, "intra_host_points.csv.gz")
    if not os.path.isfile(points_path):
        raise SystemExit(
            f"No dataset at {points_path}.\n"
            f"Either run:  --run --run-n {run_n}\n"
            f"or unzip the artifacts from the HPC and point at them:\n"
            f"  --plot --run-n {run_n} --dataset-dir <unzipped dir>")
    points = pd.read_csv(points_path)
    hosts_path = os.path.join(d, "long_shedder_hosts.csv.gz")
    hosts = pd.read_csv(hosts_path) if os.path.isfile(hosts_path) else pd.DataFrame()
    meta_path = os.path.join(d, "metadata.json")
    meta = json.load(open(meta_path)) if os.path.isfile(meta_path) else {}
    return points, hosts, meta


CELL = ['inf_duration', 'k_v', 'nsr_long']


def per_seed_osr(points):
    """Per (cell, seed) intra-host clock: through-origin slope of distance on
    time since infection -- the same tempest_regression the pipeline fits."""
    rows = []
    for key, g in points.groupby(CELL + ['seed']):
        if len(g) < 2:
            continue
        df = g.rename(columns={'t_since_infection': 'Sequencing_time',
                               'distance': 'Distance_from_root'})
        rows.append({
            'inf_duration': key[0], 'k_v': key[1], 'nsr_long': key[2],
            'seed': key[3], 'n_points': len(g),
            'osr': float(er.tempest_regression(df).coef_[0]),
        })
    return pd.DataFrame(rows)


def fit_cells(osr_df, model_type='exp', target=TARGET_OSR_LONG):
    """Per (duration, k_v): fit mean OSR vs NSR_long and invert at target --
    the same fit/invert the calibration performs."""
    rows = []
    for (dur, k_v), g in osr_df.groupby(['inf_duration', 'k_v']):
        fit_df = g.rename(columns={'nsr_long': 'nucleotide_substitution_rate',
                                   'osr': 'observed_substitution_rate'})
        rec = {'inf_duration': dur, 'k_v': k_v, 'n_points': len(g)}
        try:
            # Same model and same inversion the pipeline uses
            # (factory_model_lmfit + compute_calibrated_parameter), but NOT
            # via fit_observed_substitution_rate_regressor: that writes a CSV
            # into the experiment directory, which would tie --plot to the
            # simulation output it is meant to be independent of.
            model, init = er.factory_model_lmfit(model_type)
            fit = model.fit(fit_df['observed_substitution_rate'], init,
                            x=fit_df['nucleotide_substitution_rate'])
            params = fit.params.valuesdict()
            rec.update({
                'A': params.get('A'), 'B': params.get('B'),
                'r_squared': getattr(fit, 'rsquared', np.nan),
                'calibrated_nsr_long': float(
                    er.compute_calibrated_parameter(model_type, fit, target)),
            })
            lo, hi = g['nsr_long'].min(), g['nsr_long'].max()
            rec['extrapolated'] = not (lo <= rec['calibrated_nsr_long'] <= hi)
        except Exception as exc:
            rec.update({'A': np.nan, 'B': np.nan, 'r_squared': np.nan,
                        'calibrated_nsr_long': np.nan, 'extrapolated': True,
                        'error': str(exc)})
        rows.append(rec)
    return pd.DataFrame(rows)


def _legend(ax, title=None):
    ax.legend(title=title, frameon=False, fontsize=7, title_fontsize=7,
              labelcolor=INK)


def fig_accumulation(points, outdir, meta):
    """THE DECISIVE FIGURE. Mean divergence vs time since infection, one line
    per duration, at a fixed NSR_long. Lines lying on a common straight line
    => accumulation is linear and the duration dependence of the calibration
    is a measurement artifact. Long durations bending below => the
    sublinearity is real."""
    durations = sorted(points['inf_duration'].unique())
    colors = ramp_for(durations, DURATION_RAMP)
    k_vs = sorted(points['k_v'].unique())
    nsrs = sorted(points['nsr_long'].unique())
    nsr_mid = nsrs[len(nsrs) // 2]

    fig, axes = plt.subplots(1, len(k_vs), figsize=(3.4 * len(k_vs), 3.4),
                             sharey=True, squeeze=False, layout="constrained")
    for j, k_v in enumerate(k_vs):
        ax = axes[0, j]
        sub = points[(points['k_v'] == k_v) & (points['nsr_long'] == nsr_mid)]
        for dur in durations:
            d = sub[sub['inf_duration'] == dur]
            if d.empty:
                continue
            # bin on time so the line is a mean accumulation curve, not a cloud
            bins = np.linspace(0, d['t_since_infection'].max(), 12)
            d = d.assign(_b=pd.cut(d['t_since_infection'], bins))
            m = d.groupby('_b', observed=True).agg(
                x=('t_since_infection', 'mean'), y=('distance', 'mean'),
                n=('distance', 'size')).dropna()
            m = m[m['n'] >= 3]
            ax.plot(m['x'], m['y'], color=colors[dur], linewidth=2,
                    marker='o', markersize=3.5, label=f"{dur:g}")
        ax.set_title(f"$k_v$ = {k_v:g}", color=INK)
        ax.set_xlabel("Time since infection (years)")
        if j == 0:
            ax.set_ylabel("Mean distance from founder (subs/site)")
            _legend(ax, "Infection (d)")
    fig.suptitle(f"Intra-host accumulation at NSR$_{{long}}$ = {nsr_mid:.2e}"
                 "   (straight and overlapping = linear accumulation)",
                 fontsize=9, color=INK)
    path = os.path.join(outdir, "fig1_accumulation_curves.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_zero_inflation(points, outdir):
    """Share of intra-host observations sitting at distance zero, by duration.
    High zero fractions are the quantisation signature: a through-origin fit
    over mostly-zero integers is biased low."""
    durations = sorted(points['inf_duration'].unique())
    k_vs = sorted(points['k_v'].unique())
    colors = ramp_for(k_vs, K_V_RAMP)

    z = (points.assign(is_zero=points['distance'] <= 0)
                .groupby(['inf_duration', 'k_v'])['is_zero'].mean()
                .reset_index())

    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for k_v in k_vs:
        d = z[z['k_v'] == k_v].sort_values('inf_duration')
        ax.plot(d['inf_duration'], d['is_zero'], color=colors[k_v],
                linewidth=2, marker='o', markersize=5, label=f"{k_v:g}")
    ax.set_xlabel("Infection duration (days)")
    ax.set_ylabel("Fraction of observations at distance 0")
    ax.set_ylim(0, 1)
    ax.set_title("Zero-inflation of the intra-host clock", color=INK)
    _legend(ax, "$k_v$")
    path = os.path.join(outdir, "fig2_zero_inflation.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_calibration_surface(fits, outdir):
    """The headline: calibrated NSR_long against duration, one line per k_v.
    Flat lines would mean the calibration is duration-independent."""
    k_vs = sorted(fits['k_v'].unique())
    colors = ramp_for(k_vs, K_V_RAMP)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for k_v in k_vs:
        d = (fits[fits['k_v'] == k_v]
             .dropna(subset=['calibrated_nsr_long'])
             .sort_values('inf_duration'))
        if d.empty:
            continue
        ax.plot(d['inf_duration'], d['calibrated_nsr_long'], color=colors[k_v],
                linewidth=2, marker='o', markersize=5, label=f"{k_v:g}")
        ext = d[(d['extrapolated'] == True) & d['calibrated_nsr_long'].notna()]  # noqa: E712
        if not ext.empty:
            ax.scatter(ext['inf_duration'], ext['calibrated_nsr_long'],
                       facecolors='none', edgecolors=ACCENT, s=90, zorder=4,
                       linewidths=1.2)
    usable = fits['calibrated_nsr_long'].dropna()
    if (usable > 0).any():
        ax.set_yscale('log')
    else:
        ax.text(0.5, 0.5, "no cell produced a usable calibration",
                ha='center', va='center', transform=ax.transAxes,
                fontsize=7, color=INK_MUTED)
    ax.set_xlabel("Infection duration (days)")
    ax.set_ylabel("Calibrated NSR$_{long}$")
    ax.set_title("Calibrated NSR$_{long}$ vs infection duration\n"
                 "(ringed = inverted outside the swept range)", color=INK)
    _legend(ax, "$k_v$")
    path = os.path.join(outdir, "fig3_calibration_surface.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_fit_quality(fits, outdir):
    """Fit exponent B and R2 against duration. Both rising with duration is
    what better resolution looks like: OSR becomes closer to proportional to
    NSR (B -> 1) and the fit tightens."""
    k_vs = sorted(fits['k_v'].unique())
    colors = ramp_for(k_vs, K_V_RAMP)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), layout="constrained")
    for ax, col, lab in ((axes[0], 'B', "Fit exponent $B$  (OSR $\\propto$ NSR$^B$)"),
                         (axes[1], 'r_squared', "Fit $R^2$")):
        for k_v in k_vs:
            d = (fits[fits['k_v'] == k_v].dropna(subset=[col])
                 .sort_values('inf_duration'))
            if d.empty:
                continue
            ax.plot(d['inf_duration'], d[col], color=colors[k_v],
                    linewidth=2, marker='o', markersize=5, label=f"{k_v:g}")
        ax.set_xlabel("Infection duration (days)")
        ax.set_ylabel(lab)
    axes[0].axhline(1.0, color=ACCENT, linestyle='--', linewidth=1)
    axes[0].annotate("proportional", xy=(0.02, 0.94), xycoords='axes fraction',
                     fontsize=6.5, color=ACCENT)
    _legend(axes[1], "$k_v$")
    fig.suptitle("Conditioning of the calibration fit", fontsize=9, color=INK)
    path = os.path.join(outdir, "fig4_fit_quality.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_lineages(hosts, outdir):
    """Mean intra-host lineages per host against duration and k_v -- tests
    whether duration and k_v act through one shared variable, lineage count
    (sub_events scales mutation opportunity with it)."""
    if hosts.empty or 'n_lineages' not in hosts.columns:
        return None
    k_vs = sorted(hosts['k_v'].unique())
    colors = ramp_for(k_vs, K_V_RAMP)

    m = hosts.groupby(['inf_duration', 'k_v'])['n_lineages'].mean().reset_index()
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for k_v in k_vs:
        d = m[m['k_v'] == k_v].sort_values('inf_duration')
        ax.plot(d['inf_duration'], d['n_lineages'], color=colors[k_v],
                linewidth=2, marker='o', markersize=5, label=f"{k_v:g}")
    ax.set_xlabel("Infection duration (days)")
    ax.set_ylabel("Mean intra-host lineages per host")
    ax.set_title("Lineage load", color=INK)
    _legend(ax, "$k_v$")
    path = os.path.join(outdir, "fig5_lineage_load.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_sampling_time(points, outdir):
    """Mean time-since-infection of the observations each cell's clock is
    fitted to, against NSR_long, one line per duration.

    This is the confound the other figures hint at. The clock is a
    through-origin slope over whatever intra-host observations exist, and
    WHICH observations exist depends on the parameters: at low NSR_long and
    low k_v few lineages are born, and those that are appear early, so the
    whole sample sits at small x. A slope measured over a sample concentrated
    near the origin is not comparable to one measured across a full infection,
    and the difference tracks duration and k_v exactly as the calibrated
    NSR_long does. Flat, overlapping lines here would mean the clocks are
    being compared over the same time window and the duration effect is real.
    """
    durations = sorted(points['inf_duration'].unique())
    colors = ramp_for(durations, DURATION_RAMP)
    k_vs = sorted(points['k_v'].unique())

    m = (points.groupby(['inf_duration', 'k_v', 'nsr_long'])['t_since_infection']
               .mean().reset_index())

    fig, axes = plt.subplots(1, len(k_vs), figsize=(3.4 * len(k_vs), 3.2),
                             sharey=True, squeeze=False, layout="constrained")
    for j, k_v in enumerate(k_vs):
        ax = axes[0, j]
        for dur in durations:
            d = m[(m['k_v'] == k_v) & (m['inf_duration'] == dur)].sort_values('nsr_long')
            if d.empty:
                continue
            ax.plot(d['nsr_long'], d['t_since_infection'], color=colors[dur],
                    linewidth=2, marker='o', markersize=4, label=f"{dur:g}")
        ax.set_xscale('log')
        ax.set_title(f"$k_v$ = {k_v:g}", color=INK)
        ax.set_xlabel("NSR$_{long}$")
        if j == 0:
            ax.set_ylabel("Mean time since infection\nof fitted observations (years)")
            _legend(ax, "Infection (d)")
    fig.suptitle("What time window is each clock actually measured over?",
                 fontsize=9, color=INK)
    path = os.path.join(outdir, "fig6_sampling_time.png")
    fig.savefig(path)
    plt.close(fig)
    return path


def write_recap(outdir, meta, osr_df, fits, points, hosts):
    lines = [
        "=== NSR sensitivity analysis ===",
        f"run: {meta.get('experiment', '?')}   written: {meta.get('written', '?')}",
        f"quick preset: {meta.get('quick')}",
        f"durations (d): {meta.get('durations')}",
        f"k_v: {meta.get('k_v')}",
        f"seeds/cell: {meta.get('n_seeds')}   population: {meta.get('population_size')}",
        f"windows: {meta.get('windows')}",
        f"target intra-host OSR: {meta.get('target_osr_long', TARGET_OSR_LONG)}",
        "",
        f"intra-host observations: {len(points)}",
        f"long-shedder host records: {len(hosts)}",
        f"per-seed clocks fitted: {len(osr_df)}",
        "",
        "--- per-cell calibration ---",
        f"{'duration':>9} {'k_v':>7} {'cal NSR_long':>14} {'B':>7} {'R2':>7} {'n':>5}",
    ]
    for _, r in fits.sort_values(['k_v', 'inf_duration']).iterrows():
        flag = "  EXTRAPOLATED" if r.get('extrapolated') else ""
        nsr = "-" if pd.isna(r['calibrated_nsr_long']) else f"{r['calibrated_nsr_long']:.4e}"
        b = "-" if pd.isna(r['B']) else f"{r['B']:.3f}"
        r2 = "-" if pd.isna(r['r_squared']) else f"{r['r_squared']:.3f}"
        lines.append(f"{r['inf_duration']:9.0f} {r['k_v']:7.3g} {nsr:>14} "
                     f"{b:>7} {r2:>7} {int(r['n_points']):5d}{flag}")

    z = (points.assign(is_zero=points['distance'] <= 0)
                .groupby('inf_duration')['is_zero'].mean())
    lines += ["", "--- zero-inflation by duration (quantisation signature) ---"]
    for dur, frac in z.items():
        lines.append(f"{dur:9.0f} d: {frac*100:5.1f}% of observations at distance 0")

    lines += ["", "--- time window each clock is measured over ---",
              "(a slope fitted near the origin is not comparable to one fitted",
              " across a whole infection; this is what fig6 shows)",
              f"{'duration':>9} {'k_v':>7} {'NSR_long':>11} {'mean t (y)':>11} {'max t (y)':>10}"]
    w = points.groupby(['inf_duration', 'k_v', 'nsr_long'])['t_since_infection'].agg(['mean', 'max'])
    for (dur, k_v, nsr), r in w.iterrows():
        lines.append(f"{dur:9.0f} {k_v:7.3g} {nsr:11.2e} {r['mean']:11.4f} {r['max']:10.4f}")

    path = os.path.join(outdir, "recap.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return path


def stage_plot(args, dataset_override=None):
    """Figures from the datasets alone. Returns the artifacts zip path."""
    set_paper_rcparams()
    base = dataset_override or getattr(args, 'dataset_dir', None) or dataset_dir(args.run_n)
    points, hosts, meta = load_datasets(args.run_n, base)
    outdir = os.path.join(base, "figures")
    os.makedirs(outdir, exist_ok=True)

    print(f"[plot] {len(points)} intra-host observations, "
          f"{len(hosts)} host records")
    osr_df = per_seed_osr(points)
    fits = fit_cells(osr_df, model_type=args.model)

    osr_path = os.path.join(base, "per_seed_osr.csv")
    fits_path = os.path.join(base, "per_cell_fits.csv")
    osr_df.to_csv(osr_path, index=False)
    fits.to_csv(fits_path, index=False)

    made = [
        fig_accumulation(points, outdir, meta),
        fig_zero_inflation(points, outdir),
        fig_calibration_surface(fits, outdir),
        fig_fit_quality(fits, outdir),
        fig_lineages(hosts, outdir),
        fig_sampling_time(points, outdir),
    ]
    made = [m for m in made if m]
    recap = write_recap(outdir, meta, osr_df, fits, points, hosts)

    # Alongside the datasets being plotted: replotting an archive unzipped
    # anywhere must not write back into this checkout's Data/.
    explicit = dataset_override or getattr(args, 'dataset_dir', None)
    zip_dir = os.path.dirname(os.path.abspath(base)) if explicit else dm.get_data_dir()
    zip_path = os.path.join(zip_dir, f"{EXP_NAME}_#{args.run_n}_artifacts.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in made + [recap, osr_path, fits_path,
                         os.path.join(base, "metadata.json"),
                         os.path.join(base, "intra_host_points.csv.gz"),
                         os.path.join(base, "long_shedder_hosts.csv.gz")]:
            if p and os.path.isfile(p):
                zf.write(p, arcname=os.path.basename(p))

    print(f"\n[figures] {len(made)} written to {outdir}")
    print(f"[artifacts] {os.path.abspath(zip_path)}")
    return zip_path


# =============================================================================
# STAGE: PACKAGE
# =============================================================================
def package_raw_output(run_n, use_slurm=False):
    """Hand the raw simulation output to the repo's existing pixz packager.
    Returns the archive paths it produced. Never fatal: the artifacts zip is
    already written by this point, and a packaging failure must not cost the
    run."""
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "package_simplicity_data.sh")
    suffix = f"{EXP_NAME}_#{run_n}"
    cmd = [script, dm.get_data_dir(), suffix]
    if use_slurm:
        cmd.append("--slurm")

    before = set(glob.glob(os.path.join("Data_Export", "*.xz")))
    print(f"$ {' '.join(cmd)}")
    rc = subprocess.call(cmd)
    if rc != 0:
        print(f"  [warn] packaging exited {rc}; the artifacts zip is unaffected.")
        return []
    return sorted(set(glob.glob(os.path.join("Data_Export", "*.xz"))) - before)


# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="NSR sensitivity to infection duration and k_v. "
                    "--run is the whole HPC job (simulate, extract, plot, "
                    "compress) and prints where the output landed. --plot "
                    "replots later from the datasets alone, so changing a "
                    "figure or the fit model never costs a rerun.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true',
                      help="The whole job: simulate, extract, plot and compress. "
                          "Prints where everything landed.")
    mode.add_argument('--plot', action='store_true',
                      help="Replot from the datasets alone -- no simulation "
                          "output needed.")

    parser.add_argument('--run-n', type=int, required=True, help="Run number.")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--n-seeds', type=int, default=50,
                        help="Seeds per grid point (default 50, for the HPC).")
    parser.add_argument('--quick', action='store_true',
                        help="Local development preset: a small corner of the "
                            "grid, enough to exercise every path.")
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'],
                        help="Fit model, --plot only -- changing it needs no rerun.")
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    parser.add_argument('--dataset-dir', type=str, default=None,
                        help="--plot only: read the datasets from here instead "
                            "of Data/ (e.g. an unzipped artifacts archive).")
    parser.add_argument('--no-package', action='store_true',
                        help="--run only: skip compressing the full output.")
    parser.add_argument('--slurm', action='store_true',
                        help="--run only: submit the packaging step to SLURM.")
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()

    if args.run:
        set_slurm_resource_env(args.slurm_mem, args.slurm_time)
        stage_run(args, sp)
    else:
        stage_plot(args)


if __name__ == "__main__":
    main()
