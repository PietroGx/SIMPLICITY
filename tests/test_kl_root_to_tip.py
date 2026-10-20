#!/usr/bin/env python3
"""KL divergence between the patient and model root-to-tip distributions.

The question panel D cannot answer by eye: given the times real samples were
taken, does the model produce the observed spread of root-to-tip divergence --
and is any gap bigger than the model's own seed-to-seed noise?

THREE CORRECTIONS, none optional:

  TIME      divergence grows with time, so the raw divergence marginal is a
            mixture over each side's sampling schedule. The model is resampled
            to the PATIENT's time distribution (stratified by the patient's own
            time deciles) before anything is compared. Without this the KL
            measures sampling design.

  LENGTH    the patient overlap is ~3821 nt, the model reference 4967. A model
            sequence has 1.3x more sites to mutate. Both sides are therefore
            compared as divergence PER SITE, never as raw substitution counts.

  COVERAGE  a half-covered patient sequence can only show half its mutations.
            Sequences below MIN_COVERAGE are dropped (it costs ~20% of them).

SUPPORT: both sides are discrete -- the patient carries integer substitution
counts over a fixed overlap, the model integer counts over a fixed reference.
So there is no bandwidth to choose: the histogram is binned at the PATIENT's
own resolution, 1/overlap_length. That is the one analysis choice here, and
--bin-scale sweeps it.

NULL: the nats value alone is uninterpretable. For every seed, that seed is
treated as "the observed cohort" (drawn down to the patient's n, time-matched
the same way) and scored against the OTHER seeds. That gives the divergence
the model produces against itself at the patient's sample size -- the only
yardstick the patient value can be read against.

TWO READINGS, reported side by side:
  raw        includes the rate offset; the standard clock runs ~15% low, and
             that is a real discrepancy the figure should not hide.
  detrended  each side's own zero-intercept fit removed first, so only the
             SHAPE of the scatter is compared.

Run from the repo root:
    python tests/test_kl_root_to_tip.py --exp-num 3
"""
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, 'scripts', 'nature_plots'))
sys.path.insert(0, HERE)

import simplicity.dir_manager as dm
import simplicity.sequencing as sq
import simplicity.phenotype.distance as dis
import fig1_plots as plots
import fig1_preprocess_data as p1
from _scenarios import EXP_NAME
from test_global_clock_window import seed_points

MIN_COVERAGE = 0.95
N_TIME_BINS = 10
SMOOTH = 0.5            # Krichevsky-Trofimov; keeps KL finite on empty bins


# ---------------------------------------------------------------- loading ---
def patient_frame():
    """Patient root-to-tip: time in days from the first sample, divergence per
    site, restricted to well-covered sequences."""
    rs, _ = p1.get_patient_clock_data()
    r = rs.dropna(subset=['days_since_ref', 'hamming_per_overlap']).copy()
    before = len(r)
    if 'overlap_coverage' in r.columns:
        r = r[r['overlap_coverage'] >= MIN_COVERAGE]
    t = r['days_since_ref'].astype(float).values
    out = pd.DataFrame({'t': t - t.min(),
                        'd': r['hamming_per_overlap'].astype(float).values})
    res = 1.0 / float(np.median(r['overlap_length'])) if 'overlap_length' \
        in r.columns else 1.0 / 3821.0
    print(f"patient : {len(out)} of {before} sequences at coverage "
          f">= {MIN_COVERAGE}, t 0-{out['t'].max():.0f} d, "
          f"resolution 1/{1 / res:.0f}")
    return out, res


def model_seeds(exp_num):
    """[(seed, rows)] for the control scenario's standard individuals, with the
    host index attached so a host-level draw is possible."""
    experiment = f"{EXP_NAME}_#{exp_num}"
    out = []
    for sod in dm.get_simulation_output_dirs(experiment, group='exp__control'):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            reg = seed_points(ssod)
            if reg is None:
                continue
            reg = reg[reg['individual_type'] == 'standard']
            if reg.empty:
                continue
            out.append((os.path.basename(ssod),
                        [{'individual_index': int(i), 't': float(tt) * 365.25,
                          'd': float(dd)}
                         for i, tt, dd in zip(reg['individual_index'],
                                              reg['Sequencing_time'],
                                              reg['Distance_from_root'])]))
    print(f"model   : {len(out)} seeds, "
          f"{sum(len(r) for _, r in out)} sequences at 100%")
    return out


# ------------------------------------------------------------- the machine --
def time_match(rows, edges, want, rng):
    """Resample `rows` so its time histogram matches `want` over `edges`.

    With replacement: a bin the model under-fills would otherwise drop the
    comparison to whatever that bin happens to hold.
    """
    t = np.array([r['t'] for r in rows])
    d = np.array([r['d'] for r in rows])
    idx = np.digitize(t, edges[1:-1])
    picked = []
    for b, n_want in enumerate(want):
        pool = np.flatnonzero(idx == b)
        if len(pool) == 0 or n_want == 0:
            continue
        picked.append(rng.choice(pool, size=int(n_want), replace=True))
    if not picked:
        return np.empty(0), np.empty(0)
    sel = np.concatenate(picked)
    return t[sel], d[sel]


def pmf(values, res, bin_scale, lo, hi):
    """Histogram on the patient's own resolution, smoothed so KL stays finite."""
    width = res * bin_scale
    edges = np.arange(lo, hi + 2 * width, width)
    counts, _ = np.histogram(values, bins=edges)
    p = counts + SMOOTH
    return p / p.sum(), edges


def kl(p, q):
    """D(p||q) in nats. Both already smoothed and normalised."""
    return float(np.sum(p * np.log(p / q)))


def detrend(t, d):
    """Residuals about this sample's own zero-intercept fit."""
    if len(t) < 2 or not np.any(t > 0):
        return d
    slope = float((t @ d) / (t @ t))
    return d - slope * t


def compare(pat, res, seeds, ratio, bin_scale, reps, rng_seed):
    """KL(patient || model) at one sampling ratio, with the model's own null."""
    rng = np.random.default_rng(rng_seed)
    edges_t = np.quantile(pat['t'].values, np.linspace(0, 1, N_TIME_BINS + 1))
    edges_t[0], edges_t[-1] = -np.inf, np.inf
    want = np.histogram(pat['t'].values,
                        bins=np.quantile(pat['t'].values,
                                         np.linspace(0, 1, N_TIME_BINS + 1)))[0]
    n_pat = len(pat)

    drawn = {}
    for s, rows in seeds:
        r = rows if ratio >= 1.0 else sq.subsample(rows, ratio, seed=rng_seed)
        if len(r) >= 2:
            drawn[s] = r
    if len(drawn) < 3:
        return None

    def sample_of(keys, n):
        rows = [r for k in keys for r in drawn[k]]
        return time_match(rows, edges_t,
                          np.round(want * n / n_pat).astype(int), rng)

    keys = list(drawn)
    out = {'ratio': ratio, 'seeds': len(keys),
           'per_seed_n': int(np.mean([len(drawn[k]) for k in keys]))}

    for mode in ('raw', 'detrended'):
        prep = (lambda t, d: d) if mode == 'raw' else detrend
        # reference: every seed pooled, time-matched, at full size
        tm, dm_ = sample_of(keys, sum(len(drawn[k]) for k in keys))
        ref_vals = prep(tm, dm_)
        pv = prep(pat['t'].values, pat['d'].values)
        lo = min(ref_vals.min(), pv.min())
        hi = max(ref_vals.max(), pv.max())
        q, _ = pmf(ref_vals, res, bin_scale, lo, hi)
        p, _ = pmf(pv, res, bin_scale, lo, hi)
        observed = kl(p, q)

        null = []
        for k in keys:                       # leave-one-seed-out
            others = [x for x in keys if x != k]
            for rep in range(reps):
                th, dh = sample_of([k], n_pat)
                if len(th) < 10:
                    continue
                to, do = sample_of(others,
                                   sum(len(drawn[x]) for x in others))
                ph, _ = pmf(prep(th, dh), res, bin_scale, lo, hi)
                qo, _ = pmf(prep(to, do), res, bin_scale, lo, hi)
                null.append(kl(ph, qo))
        null = np.array(null)
        pct = float((null < observed).mean() * 100) if len(null) else np.nan
        out[mode] = {'kl': observed, 'null_median': float(np.median(null)),
                     'null_p95': float(np.percentile(null, 95)),
                     'percentile': pct, 'n_null': len(null)}
        out.setdefault('_arrays', {})[mode] = {
            'p': p, 'q': q, 'lo': lo, 'hi': hi, 'null': null,
            'width': res * bin_scale, 'ref': ref_vals, 'pat': pv}
    return out


def make_plot(results, path):
    """Left: the two distributions at full sampling. Right: the null at each
    sampling ratio, with the patient value on it."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    byr = {r['ratio']: r for r in results if '_arrays' in r}
    if not byr:
        return
    full = byr[max(byr)]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))

    # Each side on ITS OWN resolution. Binning both at the patient's width
    # put roughly every third bin over two model lattice points and gave the
    # model a bimodality it does not have -- which an earlier version of this
    # panel, and the page built from it, reported as a finding.
    a = full['_arrays']['raw']
    hi = float(max(a['ref'].max(), a['pat'].max()))
    grid = np.linspace(0.0, hi, 600)
    lat = plots._lattice(a['ref'])
    for vals, step_w, colour, style, name, lt in (
            (a['ref'], (lat or a['width']) * 2, '#333333', 'solid',
             f"model, time-matched (n/seed {full['per_seed_n']})", lat),
            (a['pat'], a['width'] * 2, '#1b9e77', (0, (4, 2)),
             'patient cohort', None)):
        c, y = plots._density(vals, step_w, 0.0, hi)
        if c is not None:
            axes[0].step(c * 1e3, y / 1e3, where='mid', color=colour, lw=0.6,
                         alpha=0.30, zorder=3)
        sm = plots._smoothed(vals, grid, lattice=lt)
        if sm is None:
            continue
        if lt is not None:
            axes[0].fill_between(grid * 1e3, 0, sm / 1e3, color=colour,
                                 alpha=0.12, lw=0, zorder=2)
        axes[0].plot(grid * 1e3, sm / 1e3, color=colour, lw=1.7,
                     linestyle=style, label=name, zorder=6)
    axes[0].set_xlim(0, hi * 1e3)
    axes[0].set_ylim(bottom=0)
    axes[0].set_xlabel('root-to-tip divergence (substitutions per 1000 sites)')
    axes[0].set_ylabel('density (per 1000 sites)')
    axes[0].set_title(f"The two distributions "
                      f"(KL = {full['raw']['kl']:.3f} nats)", fontsize=10)
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].spines[['top', 'right']].set_visible(False)

    colours = {0.10: '#D55E00', 0.50: '#E69F00', 1.0: '#0072B2'}
    for ratio in sorted(byr):
        r = byr[ratio]
        null = r['_arrays']['raw']['null']
        c = colours.get(ratio, '#888888')
        axes[1].hist(null, bins=24, density=True, color=c, alpha=0.35,
                     label=f"null, {ratio*100:.0f}% sampling "
                           f"(n/seed {r['per_seed_n']})")
        axes[1].axvline(r['raw']['kl'], color=c, lw=1.8, ls='--')
    axes[1].set_xlabel('KL (nats)')
    axes[1].set_ylabel('density')
    axes[1].set_title('Model against itself (filled) vs the patient value '
                      '(dashed)', fontsize=10)
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"[saved] {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--ratios', type=float, nargs='+', default=[0.10, 0.50, 1.0])
    p.add_argument('--bin-scales', type=float, nargs='+', default=[1.0])
    p.add_argument('--reps', type=int, default=2)
    p.add_argument('--rng-seed', type=int, default=4242)
    p.add_argument('--json-out', default='')
    p.add_argument('--plot-out', default='')
    args = p.parse_args()

    pat, res = patient_frame()
    seeds = model_seeds(args.exp_num)
    if not seeds:
        raise SystemExit("no model seeds")

    results = []
    for bs in args.bin_scales:
        print(f"\n=== bin width = {bs:g} x patient resolution ===")
        print(f"  {'ratio':>6s} {'seqs/seed':>10s} {'mode':>10s} {'KL':>9s} "
              f"{'null med':>9s} {'null p95':>9s} {'pctile':>7s}")
        for ratio in args.ratios:
            got = compare(pat, res, seeds, ratio, bs, args.reps, args.rng_seed)
            if got is None:
                print(f"  {ratio:6.2f}  too few usable seeds")
                continue
            got['bin_scale'] = bs
            results.append(got)
            for mode in ('raw', 'detrended'):
                m = got[mode]
                print(f"  {ratio:6.2f} {got['per_seed_n']:10d} {mode:>10s} "
                      f"{m['kl']:9.4f} {m['null_median']:9.4f} "
                      f"{m['null_p95']:9.4f} {m['percentile']:6.1f}%")
    if args.plot_out:
        make_plot([r for r in results if r.get('bin_scale') == 1.0]
                  or results, args.plot_out)
    if args.json_out:
        slim = [{k: v for k, v in r.items() if k != '_arrays'}
                for r in results]
        with open(args.json_out, 'w') as fh:
            json.dump({'patient_n': len(pat), 'resolution': res,
                       'results': slim}, fh, indent=2)
        print(f"\n[saved] {args.json_out}")


if __name__ == '__main__':
    main()
