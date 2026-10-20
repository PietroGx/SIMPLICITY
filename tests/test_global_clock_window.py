#!/usr/bin/env python3
"""Does the standard global clock depend on the window, and on how much of the
sequencing record goes into the fit?

cal_2 fits the standard NSR over CAL2_FINAL_TIME (730 d) and inverts at
target_osr_std; the experiment phase then runs to 1095 d. Run #3's control --
long_shedders_ratio 0.0, the same composition cal_2 calibrates in -- realizes
0.00107 against the 0.0013 it was calibrated to. Truncating a 1095-day run at
day 730 is the same record a 730-day run would have written (stopping does not
change the past), so the same seeds answer both halves:

    window  does the gap close when production is cut to cal_2's window?
    ratio   does the slope move when only 50/25/10% of the HOSTS are
            sequenced? production writes the COMPLETE diagnosed record and the
            ratio is applied afterwards (simplicity/sequencing.py).

Both fits go through er.tempest_regression -- the zero-intercept regression of
Distance_from_root on Sequencing_time that cal_2 and the sanity plot both use.
Nothing is re-derived here.

Subsampling is per HOST, via sq.subsample: a host is sequenced or it is not,
and then all of its distinct genomes are read. sequencing_data_regression.csv
carries no individual index, so it is joined positionally to
sequencing_data.csv, which does -- the two are written from the same loop in
output_manager.save_sequencing_dataset, in the same order. The join is
asserted, not assumed.

Run from the repo root:
    python tests/test_global_clock_window.py --exp-num 3
"""
import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
import simplicity.sequencing as sq
import simplicity.tuning.evolutionary_rate as er
import impact_long_shedders_config as cfg

DAYS_PER_YEAR = 365.25
SCENARIOS = ['control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case']


def seed_points(ssod):
    """One seed's regression points with the host index attached.

    Returns None when the two files disagree, rather than fitting a join that
    silently slipped by a row.
    """
    reg_path = os.path.join(ssod, 'sequencing_data_regression.csv')
    if not os.path.isfile(reg_path):
        return None
    reg = pd.read_csv(reg_path)

    # Runs from v2.4.114 on carry the host index in the regression file
    # itself, and the experiment phase no longer writes sequencing_data.csv
    # (EXPERIMENT_SKIPPED_OUTPUTS), so the join below is dead for them.
    if 'individual_index' in reg.columns:
        return reg

    # Older runs: join positionally to sequencing_data.csv, which does carry
    # it. Both files are written from the same loop in
    # output_manager.save_sequencing_dataset, in the same order -- asserted,
    # not assumed.
    seq_path = os.path.join(ssod, 'sequencing_data.csv')
    if not os.path.isfile(seq_path):
        return None
    seq = pd.read_csv(seq_path)
    if len(seq) != len(reg):
        return None
    if not (seq['individual_type'].values == reg['individual_type'].values).all():
        return None
    reg = reg.copy()
    reg['individual_index'] = seq['individual_index'].values
    return reg


def load_scenario(exp_name, exp_num, scenario, individual_type):
    """Every seed of one experiment-phase scenario, as lists of row dicts --
    the shape sq.subsample consumes."""
    experiment = f"{cfg.exp_name(cfg.DEFAULT_CONSENSUS)}_#{exp_num}"
    group = cfg.exp_group(scenario)
    seeds = []
    for sod in dm.get_simulation_output_dirs(experiment, group=group):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            reg = seed_points(ssod)
            if reg is None:
                continue
            if individual_type is not None:
                reg = reg[reg['individual_type'] == individual_type]
            if reg.empty:
                continue
            seeds.append([
                {'individual_index': int(i), 't': float(t), 'd': float(d)}
                for i, t, d in zip(reg['individual_index'],
                                   reg['Sequencing_time'],
                                   reg['Distance_from_root'])])
    return seeds


def fit(rows):
    """er.tempest_regression on pooled rows -> (slope, r2, n)."""
    df = pd.DataFrame({'Sequencing_time': [r['t'] for r in rows],
                       'Distance_from_root': [r['d'] for r in rows]})
    fitted = er.tempest_regression(df)
    x = df['Sequencing_time'].values.reshape(-1, 1)
    return (float(fitted.coef_[0]),
            float(fitted.score(x, df['Distance_from_root'].values)),
            len(df))


def one_draw(seeds, window_days, ratio, rep, min_seq):
    """Subsample hosts, cut at the window, drop seeds under the sequence gate,
    pool, fit. The gate is applied AFTER subsampling on purpose: a thinner
    record really does cost a seed its place in the fit, and hiding that would
    report only half of what a lower ratio does."""
    cut = window_days / DAYS_PER_YEAR
    pooled, kept = [], 0
    for s, rows in enumerate(seeds):
        drawn = rows if ratio >= 1.0 else sq.subsample(
            rows, ratio, seed=hash((s, rep, int(ratio * 1000))) % (2 ** 32))
        inside = [r for r in drawn if r['t'] <= cut]
        if len(inside) < min_seq:
            continue
        pooled.extend(inside)
        kept += 1
    if not pooled:
        return None
    slope, r2, n = fit(pooled)
    return {'slope': slope, 'r2': r2, 'points': n, 'seeds': kept}


def sweep(seeds, windows, ratios, reps, min_seq):
    out = {}
    for w in windows:
        for ratio in ratios:
            draws = [one_draw(seeds, w, ratio, rep, min_seq)
                     for rep in (range(1) if ratio >= 1.0 else range(reps))]
            draws = [d for d in draws if d]
            if not draws:
                out[(w, ratio)] = None
                continue
            out[(w, ratio)] = {
                'slope': float(np.mean([d['slope'] for d in draws])),
                'sd': float(np.std([d['slope'] for d in draws])),
                'r2': float(np.mean([d['r2'] for d in draws])),
                'points': float(np.mean([d['points'] for d in draws])),
                'seeds': float(np.mean([d['seeds'] for d in draws])),
            }
    return out


def scatter_and_fits(ax, seeds, window_days, ratios, target, min_seq):
    """One scenario, one window: the pooled cloud and a fitted line per ratio."""
    cut = window_days / DAYS_PER_YEAR
    full = [r for rows in seeds for r in rows if r['t'] <= cut]
    ax.scatter([r['t'] for r in full], [r['d'] for r in full],
               s=3, alpha=0.07, color='#3b528b', linewidths=0, zorder=1)
    xs = np.linspace(0, cut, 50)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(ratios)))
    for ratio, color in zip(ratios, colors):
        got = one_draw(seeds, window_days, ratio, 0, min_seq)
        if not got:
            continue
        ax.plot(xs, got['slope'] * xs, color=color, lw=1.6, zorder=3,
                label=f"{ratio * 100:.0f}%  {got['slope']:.5f}")
    ax.plot(xs, target * xs, color='k', ls='--', lw=1.2, zorder=4,
            label=f"target {target:g}")
    ax.set_xlabel('Sequencing time (years)')
    ax.set_ylabel('Distance from root (subs/site)')
    ax.set_title(f"window {window_days:g} d", fontsize=9)
    ax.legend(fontsize=6, loc='upper left', framealpha=0.9)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--scenarios', nargs='+', default=SCENARIOS)
    p.add_argument('--windows', type=float, nargs='+',
                   default=[cfg.CAL2_FINAL_TIME, 1095.0])
    p.add_argument('--ratios', type=float, nargs='+',
                   default=[1.0, 0.5, 0.25, 0.10])
    p.add_argument('--reps', type=int, default=10,
                   help="Random draws per ratio below 1.0.")
    p.add_argument('--min-seq', type=int, default=30,
                   help="Per-seed sequence gate, as the calibration fits use.")
    p.add_argument('--individual-type', default='standard',
                   help="'standard' is what cal_2 measures; 'all' pools "
                        "every diagnosed host.")
    p.add_argument('--target', type=float, default=0.0013)
    p.add_argument('--out', default=os.path.join('Data', 'figures',
                                                 'global_clock_window.png'))
    args = p.parse_args()

    itype = None if args.individual_type == 'all' else args.individual_type
    print(f"experiment  : {cfg.exp_name(cfg.DEFAULT_CONSENSUS)}_#{args.exp_num}")
    print(f"individuals : {args.individual_type}")
    print(f"windows     : {args.windows} d    ratios: {args.ratios}")
    print(f"gate        : >= {args.min_seq} sequences per seed, after "
          f"subsampling\n")

    results, loaded = {}, {}
    for scenario in args.scenarios:
        seeds = load_scenario(cfg.EXP_NAME, args.exp_num, scenario, itype)
        if not seeds:
            print(f"{scenario}: no usable seeds")
            continue
        loaded[scenario] = seeds
        results[scenario] = sweep(seeds, args.windows, args.ratios,
                                  args.reps, args.min_seq)

        print(f"--- {scenario}  ({len(seeds)} seeds on disk) ---")
        print(f"  {'window':>7s} {'ratio':>6s} {'OSR':>10s} {'sd':>9s} "
              f"{'vs target':>10s} {'R2':>6s} {'points':>8s} {'seeds':>6s}")
        for w in args.windows:
            for ratio in args.ratios:
                r = results[scenario][(w, ratio)]
                if not r:
                    print(f"  {w:7.0f} {ratio:6.2f}   (no data)")
                    continue
                print(f"  {w:7.0f} {ratio:6.2f} {r['slope']:10.6f} "
                      f"{r['sd']:9.6f} "
                      f"{100 * (r['slope'] / args.target - 1):+9.1f}% "
                      f"{r['r2']:6.2f} {r['points']:8.0f} {r['seeds']:6.1f}")
        print()

    if not results:
        raise SystemExit("nothing to plot")

    n_win = len(args.windows)
    fig, axes = plt.subplots(2, n_win, figsize=(5.5 * n_win, 9))
    axes = np.atleast_2d(axes)
    if n_win == 1:
        axes = axes.reshape(2, 1)

    head = args.scenarios[0] if args.scenarios[0] in loaded else \
        next(iter(loaded))
    for col, w in enumerate(args.windows):
        scatter_and_fits(axes[0, col], loaded[head], w, args.ratios,
                         args.target, args.min_seq)
        axes[0, col].set_title(f"{head}: window {w:g} d", fontsize=9)

        ax = axes[1, col]
        for scenario in loaded:
            ys = [results[scenario][(w, r)]['slope'] if results[scenario][(w, r)]
                  else np.nan for r in args.ratios]
            es = [results[scenario][(w, r)]['sd'] if results[scenario][(w, r)]
                  else np.nan for r in args.ratios]
            ax.errorbar([r * 100 for r in args.ratios], ys, yerr=es,
                        marker='o', ms=4, lw=1.3, capsize=2, label=scenario)
        ax.axhline(args.target, color='k', ls='--', lw=1.2,
                   label=f"target {args.target:g}")
        ax.set_xscale('log')
        ax.set_xticks([r * 100 for r in args.ratios])
        ax.set_xticklabels([f"{r * 100:.0f}%" for r in args.ratios])
        ax.set_xlabel('fraction of hosts sequenced')
        ax.set_ylabel('standard global clock (subs/site/year)')
        ax.set_title(f"window {w:g} d", fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    fig.suptitle(f"Standard global clock: analysis window and sequencing "
                 f"effort  ({cfg.exp_name(cfg.DEFAULT_CONSENSUS)}"
                 f"_#{args.exp_num}, {args.individual_type} individuals)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches='tight')
    print(f"[saved] {args.out}")


if __name__ == '__main__':
    main()
