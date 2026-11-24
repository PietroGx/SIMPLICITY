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
# Is the long-shedder NSR calibration's inversion biased by the width of the
# swept range? Reads an existing calibration's own output; simulates nothing.
#
# cal_1 fits OSR = A * NSR_long**B over the WHOLE swept range (100x in run #3)
# and inverts it at target_osr_long. A power law is only locally true here, so
# the single global fit is inverted at a level set by the far ends of the
# sweep. This prints, per calib_long group, the four numbers needed to judge
# that:
#
#   1. the global fit the run itself made, READ from its own output -- A and B
#      from Fit_results/<exp>_<model>_fit_results_<group>.csv and the inverted
#      NSR_long from <exp>_calibrated_long_nsr.csv. Not refitted here: the
#      baseline is the number the production run used, not a re-derivation of
#      it that could silently differ;
#   2. the local root: interpolation of the measured per-grid-point means
#      across the two points that bracket the target, with a bootstrap
#      interval over the seeds;
#   3. narrow-window power-law inversions (+/-1.2x .. +/-3x of the local
#      root), which show the failure mode: the fitted exponent barely moves
#      while the inverted level drifts by tens of percent. R^2 does not see
#      it, so a good fit is no evidence here;
#   4. the conditioning number 1/B_local -- how a 1% error in measured OSR is
#      amplified into NSR_long -- and the regression's resolution floor,
#      (1/len(reference)) / (tau_3_long/365.25): the smallest non-zero slope
#      one lineage can produce over one shedding window.
#
# It then prints the proposed stage-2 grid: NSR_long values geometric within
# +/-1.5x of the local root, as an NSR_RANGES-shaped dict to paste.
#
# MEASUREMENT PROVENANCE. The per-seed OSR comes from the pipeline's own
# functions and nothing else: er.extract_ih_regression_data for the intra-host
# points, er.tempest_regression for the slope, om.detect_sod_outliers per grid
# point, under the same min_seq/min_len gates -- the extraction loop here is a
# mirror of plot_and_fit_long_nsr_calibration's. If the run already stored that
# frame (Fit_results/*_osr_per_seed_long.csv) it is read instead of
# re-extracted; extraction is otherwise cached under Data/tests/.
#
# THIS TOOL WRITES NOTHING INTO THE EXPERIMENT. It does not call
# plot_and_fit_long_nsr_calibration, whose write_calibrated_long_nsr rebuilds
# <exp>_calibrated_long_nsr.csv -- the frozen table cal_2 reads back as
# Stage 1's result -- from only the groups just fitted, so a one-group run
# through that path would drop the other two cohorts' rows. The baseline is
# read, and the narrow-window and local fits persist nothing: those fit
# subsets the experiment never fitted, so they are new quantities rather than
# re-derivations, and they go through er's own fitter to stay the pipeline's
# arithmetic.
#
# Usage, from the repo root:
#   python tests/analyze_long_nsr_local_fit.py
#   python tests/analyze_long_nsr_local_fit.py --experiment impact_long_shedders_#3
#   python tests/analyze_long_nsr_local_fit.py --group calib_long__edge_case
#   python tests/analyze_long_nsr_local_fit.py --all-groups   # ~41 min
#   python tests/analyze_long_nsr_local_fit.py --target-osr 0.00205 --boot 4000
#
# --target-osr must be the target the run was calibrated against: the stored
# baseline is an inversion at that target, and read_calibrated_long_nsr
# refuses to hand back a value calibrated against another one.
#
# Default is calib_long__SOT alone: the cheapest group and the binding case
# (its target sits closest to the resolution floor). --all-groups re-extracts
# every group that has no stored frame; calib_long__edge_case alone is ~35 min
# at 6.9 s per seed.
# ============================================================================

import io
import os
import sys
import glob
import argparse
import contextlib

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'scripts', 'experiments_pipeline'))

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.phenotype.distance as dis
from impact_long_shedders_config import calib_long_group
import long_nsr_calibration_plot as c1

REF_LEN = len(dis.reference)
NSR_COL = 'nucleotide_substitution_rate_long'
OSR_COL = 'observed_substitution_rate'
DEFAULT_GROUP = calib_long_group('SOT')
DEFAULT_TARGET = 0.00205
WINDOW_FACTORS = (1.2, 1.5, 2.0, 3.0)
STAGE2_SPAN = 1.5
STAGE2_STEPS = 6
UNRESOLVED_AT = 3.0       # x the resolution floor, below which a point is noise


# =============================================================================
# GROUPS AND COHORTS -- through cal_1's own discovery, never re-derived
# =============================================================================
def long_groups(experiment):
    """The calib_long groups of `experiment`, through cal_1's own discovery
    (long_nsr_calibration_plot._cal1_groups), so this tool and the stage it
    diagnoses cannot disagree about which groups exist."""
    return [g['name'] for g in c1._cal1_groups(experiment, None)]


def cohort_of(experiment, group):
    """(tau_3_long, R_long) of one calib_long group, as cal_1 reads it.
    _cal1_groups raises itself if a group covers more than one cohort."""
    g = c1._cal1_groups(experiment, [group])[0]
    return g['tau_3_long'], g['R_long']


# =============================================================================
# THE PER-SEED OSR FRAME
# =============================================================================
def _stored_frame(experiment, group):
    """The run's own per-seed frame, if cal_1 wrote one covering `group`."""
    pattern = os.path.join(dm.get_experiment_fit_result_dir(experiment),
                           '*_osr_per_seed_long.csv')
    for path in sorted(glob.glob(pattern)):
        df = pd.read_csv(path)
        if 'group' not in df.columns or NSR_COL not in df.columns:
            continue
        sub = df[df['group'] == group]
        if not sub.empty:
            return sub.reset_index(drop=True), path
    return None, None


def _cache_path(experiment, group):
    return os.path.join(dm.get_data_dir(), 'tests',
                        f'{experiment}_{group}_osr_per_seed.csv')


def extract_per_seed(experiment, group, min_seq=30, min_len=100):
    """Per-seed intra-host OSR for one calib_long group.

    A mirror of plot_and_fit_long_nsr_calibration's extraction loop: the same
    min_len/min_seq gate on the sequencing record, the same
    er.extract_ih_regression_data intra-host points (NOT the population
    root-to-tip record), the same er.tempest_regression slope, the same
    om.detect_sod_outliers applied per grid point.
    """
    frames = []
    for sod in dm.get_simulation_output_dirs(experiment, group=group):
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(sod, NSR_COL)
        rows = []
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                if final_time >= min_len and len(seq_data) >= min_seq:
                    ih_points = er.extract_ih_regression_data(ssod)
                    if ih_points.empty:
                        continue
                    rows.append({
                        'group': group,
                        NSR_COL: nsr_val,
                        OSR_COL: er.tempest_regression(ih_points).coef_[0],
                        'simulation': os.path.basename(sod),
                        'seed': os.path.basename(ssod),
                        'n_points': len(ih_points),
                        'muts_median': float(
                            ih_points['Distance_from_root'].median()) * REF_LEN,
                        'muts_max': float(
                            ih_points['Distance_from_root'].max()) * REF_LEN,
                    })
            except Exception:
                continue
        if rows:
            frames.append(om.detect_sod_outliers(pd.DataFrame(rows)))

    if not frames:
        raise RuntimeError(
            f"no usable intra-host data in {experiment} group '{group}' "
            f"(min_seq={min_seq}, min_len={min_len}).")
    return pd.concat(frames, ignore_index=True)


def stored_or_cached_frame(experiment, group):
    """(frame, source), or (None, None) -- never extracts. For callers that
    must stay cheap, e.g. a check running on a login node."""
    df, path = _stored_frame(experiment, group)
    if df is not None:
        return df, path
    cache = _cache_path(experiment, group)
    if os.path.isfile(cache):
        return pd.read_csv(cache), cache
    return None, None


def per_seed_frame(experiment, group, min_seq=30, min_len=100, use_cache=True):
    """(frame, source). Prefers the run's stored frame, then the local cache,
    then extraction -- which is cached so a second run is fast."""
    if use_cache:
        df, path = stored_or_cached_frame(experiment, group)
        if df is not None:
            return df, path
    else:
        df, path = _stored_frame(experiment, group)
        if df is not None:
            return df, path

    cache = _cache_path(experiment, group)
    df = extract_per_seed(experiment, group, min_seq, min_len)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    df.to_csv(cache, index=False)
    return df, f'extracted -> cached at {cache}'


def grid_summary(df):
    """Per-grid-point mean / SD / n / SEM / CV over the non-outlier seeds."""
    ok = df[df['is_outlier'] == 0]
    out = (ok.groupby(NSR_COL)[OSR_COL]
             .agg(mean_osr='mean', sd_osr='std', n='count')
             .reset_index().sort_values(NSR_COL).reset_index(drop=True))
    out['sem'] = out['sd_osr'] / np.sqrt(out['n'])
    out['cv_pct'] = 100.0 * out['sd_osr'] / out['mean_osr']
    return out


# =============================================================================
# THE GLOBAL FIT -- read back from the run, never re-derived
# =============================================================================
def stored_global_fit(experiment, group, model_type, target):
    """(calibrated NSR_long, {A, B}) as the calibration stage left them on
    disk: A/B from Fit_results/<exp>_<model>_fit_results_<group>.csv and the
    inversion from <exp>_calibrated_long_nsr.csv, keyed on the group's
    (tau_3_long, R_long) cohort the way cal_2 reads it.

    Refitting here would compare the local root against a re-derivation
    rather than against the number the run actually used."""
    ab_path = om.get_fit_results_filepath(experiment, model_type, None, group)
    ab = om.read_fit_results_csv(experiment, model_type,
                                 experiment_group=group)
    if ab is None:
        raise SystemExit(
            f"no stored {model_type} fit for group '{group}': {ab_path} is "
            f"missing. This tool reads the run's own global fit and does not "
            f"refit it.")

    cal_path = c1._calibrated_long_nsr_path(experiment)
    try:
        by_cohort = c1.read_calibrated_long_nsr(experiment, target)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc))
    key = cohort_of(experiment, group)
    if key not in by_cohort:
        raise SystemExit(
            f"{cal_path} has no row for cohort (tau_3_long, R_long)={key}, "
            f"which is group '{group}'. Rows present: {sorted(by_cohort)}.")

    return float(by_cohort[key]), {k: float(v) for k, v in ab.items()}


# =============================================================================
# THE LOCAL ROOT
# =============================================================================
def bracket(x, y, target):
    """Index i with y[i] <= target <= y[i+1], over x-sorted arrays."""
    for i in range(len(y) - 1):
        if y[i] <= target <= y[i + 1]:
            return i
    return None


def local_interpolation(x, y, target):
    """The target crossing read off the measured means themselves.

    Two readings of the same two bracketing points: log-log (a local power
    law, the model cal_1 fits, and the one the geometric sweep is spaced for)
    and straight linear. Their spread is the interpolation's own ambiguity,
    and it is far smaller than the gap to the global fit.
    """
    i = bracket(x, y, target)
    if i is None:
        raise RuntimeError(
            f"target OSR {target} is not bracketed by the measured means "
            f"[{y.min():.6f} .. {y.max():.6f}]; the sweep does not cross it.")
    b = np.log(y[i + 1] / y[i]) / np.log(x[i + 1] / x[i])
    return {
        'lo': x[i], 'hi': x[i + 1], 'osr_lo': y[i], 'osr_hi': y[i + 1],
        'exponent': b,
        'nsr_loglog': x[i] * (target / y[i]) ** (1.0 / b),
        'nsr_linear': x[i] + (target - y[i]) * (x[i + 1] - x[i])
                             / (y[i + 1] - y[i]),
    }


def bootstrap_local_root(df, target, n_boot=2000, seed=0):
    """95% interval on the log-log local root, resampling the SEEDS at each
    grid point. The bracket is re-chosen each draw, so the interval covers
    the grid point the target crosses between as well as the mean's noise."""
    ok = df[df['is_outlier'] == 0]
    groups = [g[OSR_COL].values for _, g in ok.groupby(NSR_COL)]
    x = np.array(sorted(ok[NSR_COL].unique()), dtype=float)
    rng = np.random.default_rng(seed)
    roots = []
    for _ in range(n_boot):
        y = np.array([vals[rng.integers(0, len(vals), len(vals))].mean()
                      for vals in groups])
        try:
            roots.append(local_interpolation(x, y, target)['nsr_loglog'])
        except RuntimeError:
            continue
    if not roots:
        return None
    roots = np.array(roots)
    return {'lo': float(np.percentile(roots, 2.5)),
            'hi': float(np.percentile(roots, 97.5)),
            'median': float(np.median(roots)), 'n': len(roots)}


# =============================================================================
# NARROW-WINDOW FITS -- the invisibility of the bias
# =============================================================================
def window_table(df, centre, target, factors=WINDOW_FACTORS, model_type='exp'):
    """Fit A*x**B on the grid points within each factor of `centre` and invert
    at `target`. Goes through the pipeline's own fitter and inverse dispatcher,
    so the arithmetic is cal_1's; only the row subset differs. These are
    subsets the experiment never fitted, so they are new quantities, not
    re-derivations -- and the fitter persists nothing, so nothing is written.
    Its lmfit report is swallowed: four per group would bury the table."""
    ok = df[df['is_outlier'] == 0]
    rows = []
    for f in factors:
        w = ok[(ok[NSR_COL] >= centre / f) & (ok[NSR_COL] <= centre * f)]
        n_pts = w[NSR_COL].nunique()
        row = {'factor': f, 'lo': centre / f, 'hi': centre * f,
               'n_grid_points': n_pts, 'n_seeds': len(w),
               'A': np.nan, 'B': np.nan, 'rsq': np.nan, 'nsr': np.nan,
               'pct_vs_centre': np.nan, 'note': ''}
        if n_pts < 2:
            row['note'] = 'under-determined: 2 free parameters, '\
                          f'{n_pts} distinct NSR_long in window'
            rows.append(row)
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            fit = er.fit_observed_substitution_rate_regressor(
                w, model_type, parameter_name=NSR_COL)
        nsr = float(er.compute_calibrated_parameter(model_type, fit, target))
        row.update(A=fit.params['A'].value, B=fit.params['B'].value,
                   rsq=1.0 - fit.residual.var() / np.var(w[OSR_COL].values),
                   nsr=nsr, pct_vs_centre=100.0 * (nsr / centre - 1.0))
        if n_pts == 2:
            row['note'] = 'exactly determined: 2 points, 2 parameters'
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# RESOLUTION FLOOR AND THE PROPOSED STAGE-2 GRID
# =============================================================================
def resolution_floor(tau_3_long):
    """Smallest non-zero slope the intra-host regression can report: one
    substitution over the reference (1/4967 subs/site) accumulated across one
    shedding window (tau_3_long days, expressed in years, as
    extract_ih_regression_data's x axis is)."""
    return (1.0 / REF_LEN) / (tau_3_long / 365.25)


def stage2_grid(root, span=STAGE2_SPAN, steps=STAGE2_STEPS):
    return np.geomspace(root / span, root * span, steps)


# =============================================================================
# CRITERION [2] -- DIRECTIONAL
# =============================================================================
# The predicted separation a null must beat before it refutes the prediction,
# in units of the SEM of the mean OSR. A bootstrap 95% interval is +/-1.96
# SEM, so a true displacement of d SEM is missed with probability ~Phi(1.96-d):
# 20% at d=2.8, 0.1% at d=5, nil at d=9.
MIN_SEPARATION_SEM = 5.0

SIDE = {1: 'above', -1: 'below', 0: 'on'}


def directional_verdict(target, mean, lo, hi, sem, nsr_disp_pct, pred_pct,
                        sep_pred=float('nan')):
    """Criterion [2]: is the measured mean OSR displaced from `target` in the
    direction the local exponent predicts, by a resolvable amount?

    `nsr_disp_pct` is the cohort's own stored-NSR-vs-local-root displacement
    and is the only source of the predicted direction. `pred_pct` is the
    predicted OSR gap (~B x that displacement) and `sep_pred` that gap in
    units of the source run's SEM. Returns the state and the lines to print.
    No cohort is named or branched on here.
    """
    sign = 0 if nsr_disp_pct == 0 else (1 if nsr_disp_pct > 0 else -1)
    want = SIDE[sign]
    got = SIDE[0 if mean == target else (1 if mean > target else -1)]
    up, dn = lo > target, hi < target
    as_pred = (sign > 0 and up) or (sign < 0 and dn)
    opposite = (sign > 0 and dn) or (sign < 0 and up)
    null = not up and not dn
    meas_pct = 100 * (mean / target - 1)
    sep_run = (abs(pred_pct) / (100 * sem / target)
               if sem and np.isfinite(sem) and np.isfinite(pred_pct)
               else float('nan'))
    weak = [(s, n) for s, n in ((sep_pred, 'predicted, from the source run'),
                                (sep_run, 'predicted, at this run\'s own SEM'))
            if np.isfinite(s) and s < MIN_SEPARATION_SEM]
    outcome = ('no-prediction' if sign == 0 else 'as-predicted' if as_pred
               else 'opposite' if opposite else 'null')

    L = [f'[2] {"PASS" if as_pred else "FAIL"}  mean OSR at the old '
         f'global-fit value is displaced from the target in the direction '
         f'the local exponent predicts',
         f'      predicted  {pred_pct:+.2f}% of target ({want} it) -- the '
         f'sign and size of this cohort\'s own {nsr_disp_pct:+.2f}% NSR '
         f'displacement through its local exponent',
         f'      measured   {meas_pct:+.2f}% of target '
         f'({got + " it" if sign or mean != target else "exactly on it"}) -- '
         f'{mean:.6f}, boot 95% CI {lo:.6f} .. {hi:.6f}'
         + (f', {abs(mean - target) / sem:.1f} SEM from the target'
            if sem and np.isfinite(sem) else '')]
    if np.isfinite(pred_pct) and abs(pred_pct) > 1e-9:
        ratio = meas_pct / pred_pct
        L.append(f'      agreement  measured gap / predicted gap = '
                 f'{0.0 if abs(ratio) < 5e-3 else ratio:+.2f}x  '
                 f'(measurement first, then the model of the error it tests)')

    if outcome == 'as-predicted':
        L += ['      OUTCOME  displaced in the predicted direction and the CI '
              'clears the target.',
              '               The old stored NSR_long is genuinely wrong for '
              'this cohort.']
    elif outcome == 'opposite':
        L += [f'      OUTCOME  displaced the OPPOSITE way to the prediction: '
              f'predicted {pred_pct:+.2f}%, measured {meas_pct:+.2f}%, and '
              f'the CI clears',
              '               the target on the wrong side.',
              '               THIS IS THE THIRD OUTCOME AND IT MATTERS MORE '
              'THAN EITHER OF THE OTHER TWO.',
              '               The model of the error -- stored NSR_long off '
              'its local root, carried into',
              '               OSR by the local exponent -- does not describe '
              'what the simulation does, so',
              '               neither rate is trustworthy and the stage-2 '
              'argument needs explaining first.']
    elif outcome == 'null':
        L += ['      OUTCOME  lands ON the target: the boot 95% CI contains '
              'it, so no displacement is',
              '               resolved in either direction.']
        if weak:
            L += [f'               LOW POWER -- predicted separation '
                  f'{min(s for s, _ in weak):.1f} SEM '
                  f'({min(weak, key=lambda w: w[0])[1]}), under the',
                  f'               {MIN_SEPARATION_SEM:.1f} SEM this test '
                  f'needs. A null here is INCONCLUSIVE, not a refutation: a '
                  f'displacement',
                  '               of the predicted size would often be missed '
                  'at this noise level.']
        else:
            L += ['               THIS IS THE INFORMATIVE NULL: this cohort '
                  'had no error worth fixing, and',
                  '               that is the result that kills the premise '
                  'for it. Do not recalibrate it.']
    else:
        L += ['      OUTCOME  the stored NSR_long sits exactly on its local '
              'root, so the local exponent',
              '               predicts no displacement and [2] has no '
              'direction to test. [1] is the result.']

    return {'sign': sign, 'outcome': outcome, 'pass2': as_pred,
            'null': null, 'opposite': opposite, 'low_power': bool(weak),
            'want': want, 'got': got, 'meas_pct': meas_pct,
            'sep_run': sep_run, 'lines': L}


# =============================================================================
# REPORT
# =============================================================================
def report_group(experiment, group, target, args):
    print()
    print('=' * 78)
    print(f' {experiment}   group {group}')
    print('=' * 78)

    tau_3_long, r_long = cohort_of(experiment, group)
    df, source = per_seed_frame(experiment, group, args.min_seq, args.min_len,
                                use_cache=not args.no_cache)
    print(f'cohort            tau_3_long={tau_3_long} d, R_long={r_long}')
    print(f'per-seed OSR from {source}')
    print(f'rows              {len(df)}  '
          f'({int(df["is_outlier"].sum())} flagged outliers, excluded below)')
    print(f'target OSR        {target}')

    grid = grid_summary(df)
    floor = resolution_floor(tau_3_long)
    grid['x_floor'] = grid['mean_osr'] / floor

    print()
    print('--- measured OSR per grid point -------------------------------------')
    print(f'{"NSR_long":>12} {"mean OSR":>11} {"SD":>10} {"CV%":>6} '
          f'{"SEM%":>6} {"n":>4} {"x floor":>8}  flag')
    for _, r in grid.iterrows():
        flag = 'UNRESOLVED' if r['x_floor'] < UNRESOLVED_AT else ''
        print(f'{r[NSR_COL]:12.8f} {r["mean_osr"]:11.6f} {r["sd_osr"]:10.6f} '
              f'{r["cv_pct"]:6.2f} {100*r["sem"]/r["mean_osr"]:6.2f} '
              f'{int(r["n"]):4d} {r["x_floor"]:8.2f}  {flag}')
    print(f'resolution floor  {floor:.8f} subs/site/y  = (1/{REF_LEN}) over '
          f'{tau_3_long} d')
    print(f'target / floor    {target/floor:.2f}x   '
          f'(a point under {UNRESOLVED_AT:g}x is flagged: its mean is set by '
          f'how often a lineage shows its first substitution, not by a rate)')

    x = grid[NSR_COL].values.astype(float)
    y = grid['mean_osr'].values.astype(float)

    print()
    print('--- inversions ------------------------------------------------------')
    g_nsr, g_params = stored_global_fit(experiment, group, args.model, target)
    print(f'stored global fit A={g_params["A"]!r}  B={g_params["B"]!r}'
          f'   over the whole swept range')
    print(f'  -> NSR_long     {g_nsr!r}      '
          f'[stored inversion: what the run calibrated to]')

    loc = local_interpolation(x, y, target)
    print(f'local bracket     {loc["lo"]:.8f} (OSR {loc["osr_lo"]:.6f})  ..  '
          f'{loc["hi"]:.8f} (OSR {loc["osr_hi"]:.6f})')
    print(f'  -> NSR_long     {loc["nsr_loglog"]:.9f}      '
          f'[log-log interpolation of the measured means]')
    print(f'  -> NSR_long     {loc["nsr_linear"]:.9f}      '
          f'[linear interpolation, same two points]')
    print(f'stored vs local   {100*(g_nsr/loc["nsr_loglog"]-1):+.1f}%  '
          f'(stored global over log-log local)')

    boot = None
    if args.boot > 0:
        boot = bootstrap_local_root(df, target, args.boot, args.boot_seed)
        if boot:
            print(f'bootstrap 95% CI  {boot["lo"]:.9f} .. {boot["hi"]:.9f}   '
                  f'({boot["n"]}/{args.boot} draws bracketed the target)')
            inside = boot['lo'] <= g_nsr <= boot['hi']
            print(f'  stored value is  {"INSIDE" if inside else "OUTSIDE"} '
                  f'that interval')

    print()
    print('--- conditioning ----------------------------------------------------')
    b = loc['exponent']
    print(f'local log-log exponent B   {b:.4f}   (stored fit B={g_params["B"]:.4f})')
    print(f'conditioning 1/B           {1.0/b:.2f}x   -- a 1% error in measured '
          f'OSR becomes {100*((1.01)**(1.0/b)-1):.1f}% in NSR_long')
    worst = grid.loc[(grid[NSR_COL] >= loc['lo']) & (grid[NSR_COL] <= loc['hi'])]
    if not worst.empty:
        sem_pct = (100 * worst['sem'] / worst['mean_osr']).max()
        print(f'worst SEM near target      {sem_pct:.2f}% of the mean over '
              f'{int(worst["n"].max())} seeds '
              f'-> {sem_pct/b:.2f}% in NSR_long')

    print()
    print('--- narrow-window power-law inversions ------------------------------')
    print('centre = log-log local root. B and R^2 are what a reader would')
    print('check; watch them hold -- R^2 can even improve -- while the')
    print('inverted level walks away. That is why a good fit is not evidence.')
    wt = window_table(df, loc['nsr_loglog'], target, model_type=args.model)
    print(f'{"window":>8} {"pts":>4} {"seeds":>6} {"B":>8} {"R^2":>8} '
          f'{"NSR_long":>12} {"vs centre":>10}  note')
    for _, r in wt.iterrows():
        print(f'{r["factor"]:7.2f}x {int(r["n_grid_points"]):4d} '
              f'{int(r["n_seeds"]):6d} '
              f'{r["B"]:8.4f} {r["rsq"]:8.5f} {r["nsr"]:12.8f} '
              f'{r["pct_vs_centre"]:+9.1f}%  {r["note"]}'
              if np.isfinite(r['B']) else
              f'{r["factor"]:7.2f}x {int(r["n_grid_points"]):4d} '
              f'{int(r["n_seeds"]):6d} {"":8} {"":8} {"":12} {"":10}  '
              f'{r["note"]}')

    print()
    print('--- proposed stage-2 grid -------------------------------------------')
    root = loc['nsr_loglog']
    pts = stage2_grid(root, args.span, args.steps)
    print(f'provisional root  {root:.9f}   '
          f'(log-log local interpolation, stage 1)')
    print(f'window            +/-{args.span:g}x  -> '
          f'{pts[0]:.9f} .. {pts[-1]:.9f}')
    print(f'{args.steps} NSR_long values (geometric):')
    print('  [' + ', '.join(f'{v:.9f}' for v in pts) + ']')
    print('NSR_RANGES-shaped, to paste:')
    print(f'    {{"min": {pts[0]:.3e}, "max": {pts[-1]:.3e}, '
          f'"steps": {args.steps}}}')

    return {'group': group, 'tau_3_long': tau_3_long, 'global': g_nsr,
            'local': root, 'boot': boot, 'exponent': b, 'floor': floor,
            'grid': pts}


def main():
    p = argparse.ArgumentParser(
        description="Global-fit vs local inversion of an existing "
                    "long-shedder NSR calibration. Simulates nothing.")
    p.add_argument('--experiment', default='impact_long_shedders_#3')
    p.add_argument('--target-osr', type=float, default=DEFAULT_TARGET)
    p.add_argument('--group', action='append', default=None,
                   help=f"calib_long group; repeatable. Default {DEFAULT_GROUP}.")
    p.add_argument('--all-groups', action='store_true',
                   help="every calib_long group (~41 min if nothing is cached)")
    p.add_argument('--model', default='exp',
                   choices=['lin', 'log', 'exp', 'tan'])
    p.add_argument('--min-seq', type=int, default=30)
    p.add_argument('--min-len', type=int, default=100)
    p.add_argument('--boot', type=int, default=2000,
                   help="bootstrap draws for the local root (0 to skip)")
    p.add_argument('--boot-seed', type=int, default=0)
    p.add_argument('--span', type=float, default=STAGE2_SPAN,
                   help="half-width of the proposed stage-2 grid, as a factor")
    p.add_argument('--steps', type=int, default=STAGE2_STEPS)
    p.add_argument('--no-cache', action='store_true')
    args = p.parse_args()

    present = long_groups(args.experiment)
    if not present:
        raise SystemExit(f"{args.experiment} has no calib_long group. "
                         f"Groups on disk: {dm.get_groups(args.experiment)}")
    if args.all_groups:
        wanted = present
    elif args.group:
        wanted = args.group
    else:
        wanted = [DEFAULT_GROUP]
    missing = [g for g in wanted if g not in present]
    if missing:
        raise SystemExit(f"{args.experiment} has no group(s) {missing}. "
                         f"calib_long groups present: {present}")

    results = [report_group(args.experiment, g, args.target_osr, args)
               for g in wanted]

    print()
    print('=' * 78)
    print(' summary')
    print('=' * 78)
    print(f'{"group":<30} {"stored":>12} {"local root":>12} {"stored-local":>13}')
    for r in results:
        print(f'{r["group"]:<30} {r["global"]:12.8f} {r["local"]:12.8f} '
              f'{100*(r["global"]/r["local"]-1):+12.1f}%')


if __name__ == '__main__':
    main()
