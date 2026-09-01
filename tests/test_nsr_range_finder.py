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
# HPC test: find the NSR sweep ranges for cal_1 and cal_2, i.e. what
# NSR_RANGES in impact_long_shedders_config.py should be set to before a
# production run.
#
# Runs the two stages SEQUENTIALLY under one test number, the way the
# pipeline does: stage 1 probes NSR_long, fits and inverts its own result,
# and hands that frozen long-side rate to stage 2 -- so stage 2 is probed in
# the population context stage 1 actually implies, not against a hardcoded
# or borrowed long NSR. --stage 1 or --stage 2 runs just one (stage 2 alone
# then needs --long-calib-exp or --nsr-long).
#
# Coarse, wide and cheap: a few seeds over a few nodes, to locate the target
# rather than to fit it precisely. The two existing grid tests answer "does
# this range work?"; this one answers "which range should I use?", and
# prints a ready-to-paste NSR_RANGES entry per stage.
#
# Both stages go through impact_long_shedders_config's own builders
# (build_cal1_settings / build_cal2_scenario_groups + build_cal2_settings),
# so the population context and swept parameter are cal_1's and cal_2's by
# construction and cannot drift from them. What differs per stage is only:
#
#   stage 1  sweeps nucleotide_substitution_rate_long in the isolated 100%
#            long-shedder context, measures each seed's INTRA-HOST clock
#            (er.extract_ih_regression_data), targets --target-osr-long,
#            and groups by (tau_3_long, R_long);
#   stage 2  sweeps nucleotide_substitution_rate in the mixed production
#            population with the long side frozen, measures the GLOBAL clock
#            over STANDARD individuals only, targets --target-osr-std, and
#            groups by scenario.
#
# Three things are reported per node, and all three constrain the range:
#
#   OSR    -- mean observed rate; the target must be bracketed, and the
#             recommendation is the bracketing pair.
#   yield  -- share of seeds surviving the min_len/min_seq gate. Run #5's
#             Stage 2 kept 3.5-15% of its seeds and fit HIV_high on 7
#             points; a node that brackets the target but yields almost
#             nothing is not usable, so the recommended max is clipped to
#             the highest node still meeting --min-yield.
#   d      -- substitutions/site accumulated over the window, from the
#             MEASURED rate (mean OSR x window / 365.25). Raw Hamming
#             distance is uncorrected for multiple hits, so nodes past ~0.1
#             measure a saturated quantity; these are flagged and kept out
#             of the recommendation. Deliberately not NSR x window_days:
#             although sub_events draws Poisson(NSR*dt*L*lineages), so NSR
#             is nominally per site per day, realized divergence comes out
#             roughly two orders of magnitude below that -- the input rate
#             would put the saturation boundary in the wrong place.
#
# Not part of the impact_long_shedders pipeline; writes nothing it reads.
#
# Writes Data/tests/test_nsr_range_out_#{test_n}.txt covering every stage
# run, plus one figure per stage, zipped into
# Data/{EXP_NAME}_#{test_n}_artifacts.zip.
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
    SCENARIOS, TAU_ROUND, USER_FIXED_PARAMS, CAL1_ISOLATED_FIXED_PARAMS,
    CAL2_FINAL_TIME, phase_offset, derive_scenario_params,
    unique_long_tau_r_long_pairs, build_cal1_settings,
    build_cal2_scenario_groups, build_cal2_settings,
)
from long_nsr_calibration_plot import read_calibrated_long_nsr

EXP_NAME = "nsr_range_finder_test"

DEFAULT_SCENARIO = "SOT"
DEFAULT_NSR_MIN = 1e-04
DEFAULT_NSR_MAX = 1e-03
DEFAULT_STEPS = 3
DEFAULT_SEEDS = 10
SATURATION_D = 0.1        # subs/site past which raw Hamming saturates

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}


# =============================================================================
# SETTINGS -- both stages via the config's own builders
# =============================================================================
def scenario_pair(scenario_name, sp):
    """(tau_3_long, R_long) of one scenario -- cal_1's grouping key."""
    scenario = next(s for s in SCENARIOS if s["name"] == scenario_name)
    frozen = derive_scenario_params(scenario, sp)
    return (round(frozen["tau_3_long"], TAU_ROUND),
            round(frozen["R_long"], TAU_ROUND))


def build_stage1(ranges, n_seeds, final_time, k_v, scenario_name, sp):
    """cal_1's own builder, subset to one scenario's (tau_3_long, R_long)
    group -- a range probe does not need the whole grid."""
    varying, fixed, seeds = build_cal1_settings(
        n_seeds, ranges, R=CAL1_ISOLATED_FIXED_PARAMS['R'],
        ih_virus_emergence_rate=k_v)()

    want = scenario_pair(scenario_name, sp)
    kept = [g for g in varying['_scenario_groups']
            if (round(g['tau_3_long'], TAU_ROUND),
                round(g['R_long'], TAU_ROUND)) == want]
    if not kept:
        raise SystemExit(f"no cal_1 group matches scenario '{scenario_name}' {want}")

    fixed = dict(fixed, final_time=final_time)
    return lambda: ({'_scenario_groups': kept}, fixed, seeds)


def build_stage2(ranges, n_seeds, final_time, k_v, long_nsr_by_group,
                 scenario_name):
    """cal_2's own builder, subset to one scenario's group."""
    nsr_ranges = {s["name"]: ranges for s in SCENARIOS}
    scenario_groups, scenarios_frozen = build_cal2_scenario_groups(
        nsr_ranges, long_nsr_by_group)

    kept = [g for g, f in zip(scenario_groups, scenarios_frozen)
            if f["scenario_name"] == scenario_name]
    if not kept:
        raise SystemExit(f"no cal_2 group matches scenario '{scenario_name}'")

    varying, fixed, seeds = build_cal2_settings(
        kept, n_seeds, USER_FIXED_PARAMS['R'], k_v)()
    fixed = dict(fixed, final_time=final_time)
    return lambda: (varying, fixed, seeds)


# =============================================================================
# MEASUREMENT -- each stage measures its own clock
# =============================================================================
def measure_stage1(ssod, min_seq, min_len):
    """Intra-host clock, as cal_1 fits it."""
    final_time = om.read_final_time(ssod)
    seq_data = om.read_sequencing_data_regression(ssod)
    if final_time < min_len or len(seq_data) < min_seq:
        return None
    ih_points = er.extract_ih_regression_data(ssod)
    if ih_points.empty:
        return None
    return float(er.tempest_regression(ih_points).coef_[0])


def measure_stage2(ssod, min_seq, min_len):
    """Global clock over standard individuals only, as cal_2 fits it."""
    final_time = om.read_final_time(ssod)
    seq_data = om.read_sequencing_data_regression(ssod)
    seq_data = seq_data[seq_data['individual_type'] == 'standard']
    if final_time < min_len or len(seq_data) < min_seq:
        return None
    return float(er.tempest_regression(seq_data).coef_[0])


def collect_points(numbered, stage, nsr_param, group_keys, min_seq, min_len):
    """Per-surviving-seed OSR rows, plus seeds attempted per node (the
    denominator of `yield`, which is only knowable here)."""
    measure = measure_stage1 if stage == 1 else measure_stage2
    points, attempts = [], []

    for sod in dm.get_simulation_output_dirs(numbered):
        nsr_val = float(sm.get_parameter_value_from_simulation_output_dir(
            sod, nsr_param))
        key = tuple(round(float(sm.get_parameter_value_from_simulation_output_dir(
            sod, k)), TAU_ROUND) for k in group_keys)

        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        for ssod in seeded_dirs:
            try:
                osr = measure(ssod, min_seq, min_len)
            except Exception:
                continue
            if osr is not None and np.isfinite(osr) and osr > 0:
                points.append({'group_key': key, 'nsr': nsr_val,
                               'observed_substitution_rate': osr})
        attempts.append({'group_key': key, 'nsr': nsr_val,
                         'n_seeds': len(seeded_dirs)})

    return pd.DataFrame(points), pd.DataFrame(attempts)


def node_table(points_df, attempts_df, group_key):
    """Mean OSR / surviving seeds / seeds attempted, per NSR node."""
    att = attempts_df[attempts_df['group_key'].apply(lambda k: k == group_key)]
    if att.empty:
        return pd.DataFrame()

    pts = (points_df[points_df['group_key'].apply(lambda k: k == group_key)]
           if not points_df.empty else points_df)

    rows = []
    for _, a in att.iterrows():
        sub = (pts[np.isclose(pts['nsr'], a['nsr'])]
               if not pts.empty else pts)
        rows.append({
            'nsr': a['nsr'],
            'mean_osr': float(sub['observed_substitution_rate'].mean())
                        if len(sub) else np.nan,
            'n_ok': len(sub),
            'n_seeds': int(a['n_seeds']),
        })
    return pd.DataFrame(rows)


# =============================================================================
# RECOMMENDATION
# =============================================================================
def recommend_range(nodes, target, window_days, min_yield):
    """Bracketing NSR pair for `target`, clipped to nodes that are neither
    starved of surviving seeds nor saturated."""
    df = nodes.sort_values('nsr').reset_index(drop=True)
    df['yield'] = df['n_ok'] / df['n_seeds'].replace(0, np.nan)
    # Divergence over the window, from the MEASURED rate (subs/site/year), not
    # from NSR x days: the two disagree by ~100x in practice, so deriving it
    # from the input rate would flag saturation in the wrong place entirely.
    df['d'] = df['mean_osr'] * window_days / 365.25

    usable = df[(df['yield'] >= min_yield) & df['mean_osr'].notna()
                & (df['d'] <= SATURATION_D)]
    if usable.empty:
        return {'ok': False, 'table': df,
                'reason': 'no node had both usable yield and an unsaturated, '
                          'measurable OSR'}

    below = usable[usable['mean_osr'] <= target]
    above = usable[usable['mean_osr'] >= target]
    if below.empty:
        return {'ok': False, 'table': df,
                'reason': f"every usable node is ABOVE target {target:g} "
                          f"(lowest mean OSR {usable['mean_osr'].min():.6f} at "
                          f"NSR={usable['nsr'].iloc[0]:g}) -- extend the sweep down"}
    if above.empty:
        return {'ok': False, 'table': df,
                'reason': f"every usable node is BELOW target {target:g} "
                          f"(highest mean OSR {usable['mean_osr'].max():.6f} at "
                          f"NSR={usable['nsr'].iloc[-1]:g}) -- extend the sweep up"}

    return {
        'ok': True,
        'lo': float(below['nsr'].max()),
        'hi': float(above['nsr'].min()),
        'ceiling': float(usable['nsr'].max()),
        'osr_lo': float(below.loc[below['nsr'].idxmax(), 'mean_osr']),
        'osr_hi': float(above.loc[above['nsr'].idxmin(), 'mean_osr']),
        'table': df,
    }


def calibrated_nsr_from_probe(numbered, group_points, nsr_param, target, rec,
                              model_type, label):
    """Stage 1's own calibrated NSR_long, to hand to stage 2. Fit+invert
    exactly as cal_1 does, but a probe has few nodes, so the result is
    sanity-checked and falls back to the bracketing pair rather than
    handing stage 2 an extrapolated number. Returns (value, how)."""
    lo_bound, hi_bound = group_points['nsr'].min(), group_points['nsr'].max()
    # fit_observed_substitution_rate_regressor looks the x column up by
    # parameter name, and falls back to column 0 with only a warning if it
    # is missing -- which silently fits against group_key.
    fit_df = group_points.rename(columns={'nsr': nsr_param})

    try:
        fit = er.fit_observed_substitution_rate_regressor(
            numbered, fit_df, model_type=model_type,
            parameter_name=nsr_param, experiment_group=label)
        nsr = float(er.compute_calibrated_parameter(model_type, fit, target))
        if np.isfinite(nsr) and nsr > 0:
            if lo_bound <= nsr <= hi_bound:
                return nsr, f"fit+invert ({model_type})"
            return nsr, (f"fit+invert ({model_type}), EXTRAPOLATED outside the "
                         f"probed range [{lo_bound:g}, {hi_bound:g}]")
    except Exception as e:
        print(f"  [probe fit failed: {e}]")

    if rec['ok']:
        nsr = float(np.sqrt(rec['lo'] * rec['hi']))
        return nsr, "geometric mean of the bracketing nodes (fit unusable)"

    closest = group_points.groupby('nsr')['observed_substitution_rate'].mean()
    nsr = float((closest - target).abs().idxmin())
    return nsr, "node with mean OSR closest to target (no bracket, fit unusable)"


def format_node_table(df):
    lines = [f"    {'NSR':>12s} {'mean OSR':>12s} {'yield':>7s} {'n':>7s} "
             f"{'subs/site':>10s}"]
    for _, r in df.iterrows():
        osr = "-" if pd.isna(r['mean_osr']) else f"{r['mean_osr']:.6f}"
        d = "-" if pd.isna(r['d']) else f"{r['d']:.3g}"
        flag = "  SATURATED" if pd.notna(r['d']) and r['d'] > SATURATION_D else ""
        lines.append(f"    {r['nsr']:12.6g} {osr:>12s} "
                     f"{r['yield']:7.2f} {int(r['n_ok']):3d}/{int(r['n_seeds']):<3d} "
                     f"{d:>10s}{flag}")
    return lines


# =============================================================================
# GROUPS -- match key and the window each stage's saturation uses
# =============================================================================
def stage_groups(stage, sp, scenario_name, final_time):
    """({label: match key}, {label: window in days}) for one stage."""
    if stage == 1:
        key = scenario_pair(scenario_name, sp)
        # Intra-host divergence accumulates over one infection, not the run.
        return {scenario_name: key}, {scenario_name: key[0] + phase_offset(sp)}

    scenario = next(s for s in SCENARIOS if s["name"] == scenario_name)
    frozen = derive_scenario_params(scenario, sp)
    key = (round(frozen["tau_3_long"], TAU_ROUND),
           round(frozen["long_shedders_ratio"], TAU_ROUND))
    return {scenario_name: key}, {scenario_name: float(final_time)}


# =============================================================================
# ONE STAGE, END TO END
# =============================================================================
def run_stage(stage, args, sp, long_nsr_by_group=None):
    """Submit, measure, recommend and plot one stage. Returns a dict with
    the recap lines, the plot path, and (stage 1) the calibrated NSR_long to
    carry into stage 2."""
    nsr_param = ('nucleotide_substitution_rate_long' if stage == 1
                 else 'nucleotide_substitution_rate')
    group_keys = (('tau_3_long', 'R_long') if stage == 1
                  else ('tau_3_long', 'long_shedders_ratio'))
    target = args.target_osr_long if stage == 1 else args.target_osr_std
    final_time = (args.final_time if args.final_time is not None else
                  (CAL1_ISOLATED_FIXED_PARAMS['final_time'] if stage == 1
                   else CAL2_FINAL_TIME))

    ranges = {'min': args.nsr_min, 'max': args.nsr_max, 'steps': args.steps}
    if stage == 1:
        settings_func = build_stage1(ranges, args.n_seeds, final_time,
                                     args.ih_virus_emergence_rate,
                                     args.scenario, sp)
    else:
        settings_func = build_stage2(ranges, args.n_seeds, final_time,
                                     args.ih_virus_emergence_rate,
                                     long_nsr_by_group, args.scenario)

    groups, windows = stage_groups(stage, sp, args.scenario, final_time)
    numbered = f"{EXP_NAME}_s{stage}_#{args.test_n}"
    _, fixed_params, _ = settings_func()

    print("\n=========================================================")
    print(f" STAGE {stage} -- {'cal_1, intra-host clock' if stage == 1 else 'cal_2, global clock (standard only)'}")
    print(f" experiment: {numbered}   scenario: {args.scenario}")
    print(f" sweep {nsr_param}: {args.nsr_min:g} .. {args.nsr_max:g} "
          f"({args.steps} nodes), {args.n_seeds} seeds/node")
    for k, v in fixed_params.items():
        print(f"   {k}: {v}")
    print("=========================================================\n")

    run_experiment(numbered, settings_func,
                   simplicity_runner=RUNNERS[args.runner],
                   archive_experiment=False)

    points_df, attempts_df = collect_points(
        numbered, stage, nsr_param, group_keys, args.min_seq, args.min_len)

    recap = [
        f"=== STAGE {stage}: {numbered} ===",
        f"swept: {nsr_param}",
        f"measured: {'intra-host clock' if stage == 1 else 'global clock, standard individuals only'}",
        f"target OSR = {target}",
        f"sweep: {args.nsr_min:g} .. {args.nsr_max:g} ({args.steps} nodes), "
        f"{args.n_seeds} seeds/node, final_time={final_time}",
        f"gates: min_seq={args.min_seq}, min_len={args.min_len}, "
        f"min_yield={args.min_yield}, saturation at d>{SATURATION_D}",
    ]
    if stage == 2:
        recap.append(f"frozen long NSR: {long_nsr_by_group}")
    recap.append("")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    result = {'recap': recap, 'plot_path': None, 'nsr_long': None,
              'recommendations': {}}

    for label, key in groups.items():
        window = windows[label]
        recap.append(f"--- {label} (window {window:.1f} d) ---")

        nodes = node_table(points_df, attempts_df, key)
        if nodes.empty:
            recap.extend(["  no simulation output for this group", ""])
            print(f"{label:20s}: no simulation output")
            continue

        rec = recommend_range(nodes, target, window, args.min_yield)
        recap.extend(format_node_table(rec['table']))

        if rec['ok']:
            recap.append(f"  brackets target: OSR {rec['osr_lo']:.6f} @ NSR "
                         f"{rec['lo']:g}  ->  {rec['osr_hi']:.6f} @ NSR {rec['hi']:g}")
            recap.append(f"  highest usable node (yield+saturation): {rec['ceiling']:g}")
            recap.append(f'  RECOMMENDED: {{"min": {rec["lo"]:g}, '
                         f'"max": {rec["hi"]:g}, "steps": _NSR_STEPS}}')
            result['recommendations'][label] = (rec['lo'], rec['hi'])
            print(f"{label:20s}: brackets target in [{rec['lo']:g}, {rec['hi']:g}]"
                  f"  (usable ceiling {rec['ceiling']:g})")
        else:
            recap.append(f"  NO RECOMMENDATION: {rec['reason']}")
            print(f"{label:20s}: {rec['reason']}")

        group_points = points_df[points_df['group_key'].apply(lambda k: k == key)] \
            if not points_df.empty else points_df
        if stage == 1 and not group_points.empty:
            nsr_long, how = calibrated_nsr_from_probe(
                numbered, group_points, nsr_param, target, rec,
                args.model, label)
            result['nsr_long'] = nsr_long
            recap.append(f"  calibrated NSR_long for stage 2: {nsr_long:.8f}  [{how}]")
            print(f"{label:20s}: NSR_long for stage 2 = {nsr_long:.8f} ({how})")
        recap.append("")

        plot_df = rec['table'].dropna(subset=['mean_osr'])
        line, = ax.plot(plot_df['nsr'], plot_df['mean_osr'], marker='o',
                        linewidth=1.5, label=label)
        starved = plot_df[plot_df['yield'] < args.min_yield]
        ax.scatter(starved['nsr'], starved['mean_osr'], marker='x', s=90,
                   color=line.get_color(), zorder=3)
        saturated = plot_df[plot_df['d'] > SATURATION_D]
        ax.scatter(saturated['nsr'], saturated['mean_osr'], facecolors='none',
                   edgecolors=line.get_color(), s=140, zorder=3)

    ax.axhline(target, color='black', linestyle='--', linewidth=1.5,
               label=f'target {target:g}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(nsr_param)
    ax.set_ylabel("mean observed substitution rate")
    ax.set_title(f"NSR range finder -- stage {stage} ({numbered})\n"
                 "x = below min yield, o = saturated")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f"{numbered}_nsr_range_finder.png")
    plt.savefig(plot_path, dpi=200)
    plt.close()
    result['plot_path'] = plot_path
    print(f"Plot saved to: {plot_path}")

    return result


def write_artifacts(test_n, stage_results):
    """One recap covering every stage run, zipped with each stage's plot."""
    recap = []
    for r in stage_results:
        recap.extend(r['recap'])

    paste = [(label, rng) for r in stage_results
             for label, rng in r['recommendations'].items()]
    if paste:
        recap.append("--- paste into impact_long_shedders_config.NSR_RANGES ---")
        for label, (lo, hi) in paste:
            recap.append(f'  {label}: {{"min": {lo:g}, "max": {hi:g}, '
                         f'"steps": _NSR_STEPS}}')
        recap.append("")

    recap_dir = os.path.join(dm.get_data_dir(), "tests")
    os.makedirs(recap_dir, exist_ok=True)
    recap_path = os.path.join(recap_dir, f"test_nsr_range_out_#{test_n}.txt")
    with open(recap_path, "w") as f:
        f.write("\n".join(recap))
    print(f"\nRecap written to: {recap_path}")

    zip_path = os.path.join(dm.get_data_dir(),
                            f"{EXP_NAME}_#{test_n}_artifacts.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(recap_path, arcname=os.path.basename(recap_path))
        for r in stage_results:
            if r['plot_path'] and os.path.isfile(r['plot_path']):
                zf.write(r['plot_path'], arcname=os.path.basename(r['plot_path']))
    print(f"Artifacts zipped to: {zip_path}")

    return recap_path, zip_path


def stage2_long_nsr(args, sp, stage1_result):
    """Frozen long-side NSR for stage 2: stage 1's own probe when the two ran
    together, else an explicit override or a real cal_1 experiment."""
    if args.nsr_long is not None:
        return {pair: args.nsr_long for pair in unique_long_tau_r_long_pairs(sp)}
    if stage1_result is not None and stage1_result['nsr_long'] is not None:
        return {pair: stage1_result['nsr_long']
                for pair in unique_long_tau_r_long_pairs(sp)}
    if args.long_calib_exp:
        return read_calibrated_long_nsr(args.long_calib_exp, args.target_osr_long)
    raise SystemExit(
        "stage 2 needs a frozen long-side NSR. Run both stages together "
        "(the default), or pass --long-calib-exp "
        "<impact_long_shedders_calibration_lng_nsr_#N> or --nsr-long <value>.")


def main():
    parser = argparse.ArgumentParser(
        description="HPC test: find the NSR sweep ranges for cal_1 and cal_2. "
                    "Runs both stages sequentially by default, stage 1's "
                    "calibrated NSR_long feeding stage 2.")
    parser.add_argument('test_n', type=int)
    parser.add_argument('--stage', type=str, default='both',
                        choices=['both', '1', '2'],
                        help="Default 'both': cal_1 then cal_2, sequentially.")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--scenario', type=str, default=DEFAULT_SCENARIO,
                        choices=[s["name"] for s in SCENARIOS],
                        help=f"Scenario to probe (default {DEFAULT_SCENARIO}). "
                            "One group per run keeps the probe cheap; stage 1 "
                            "resolves it to its (tau_3_long, R_long) pair.")
    parser.add_argument('--n-seeds', type=int, default=DEFAULT_SEEDS,
                        help="Seeds per node; this is reconnaissance, not a fit.")
    parser.add_argument('--nsr-min', type=float, default=DEFAULT_NSR_MIN)
    parser.add_argument('--nsr-max', type=float, default=DEFAULT_NSR_MAX)
    parser.add_argument('--steps', type=int, default=DEFAULT_STEPS)
    parser.add_argument('--final-time', type=int, default=None,
                        help="Override both stages' window. Default: each "
                            f"stage's own (cal_1 "
                            f"{CAL1_ISOLATED_FIXED_PARAMS['final_time']}, "
                            f"cal_2 {CAL2_FINAL_TIME}).")
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'],
                        help="Fit model used to invert stage 1's probe into "
                            "the NSR_long handed to stage 2.")
    parser.add_argument('--ih-virus-emergence-rate', type=float,
                        default=USER_FIXED_PARAMS['IH_virus_emergence_rate'])
    parser.add_argument('--long-calib-exp', type=str, default=None,
                        help="--stage 2 only: completed cal_1 experiment to "
                            "read the frozen long-side NSR from.")
    parser.add_argument('--nsr-long', type=float, default=None,
                        help="Override the frozen long-side NSR stage 2 runs "
                            "with, instead of stage 1's own result.")
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    parser.add_argument('--min-yield', type=float, default=0.25,
                        help="Min share of seeds a node must keep to be "
                            "eligible for the recommended range.")
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()
    stages = [1, 2] if args.stage == 'both' else [int(args.stage)]

    print("\n#########################################################")
    print(f" NSR range finder, test #{args.test_n} -- stages {stages}")
    print(f" scenario {args.scenario}, {args.n_seeds} seeds/node, "
          f"{args.steps} nodes over {args.nsr_min:g}..{args.nsr_max:g}")
    print("#########################################################")

    stage_results, stage1_result = [], None
    for stage in stages:
        if stage == 1:
            result = run_stage(1, args, sp)
            stage1_result = result
        else:
            long_nsr = stage2_long_nsr(args, sp, stage1_result)
            result = run_stage(2, args, sp, long_nsr_by_group=long_nsr)
        stage_results.append(result)

    write_artifacts(args.test_n, stage_results)


if __name__ == "__main__":
    main()
