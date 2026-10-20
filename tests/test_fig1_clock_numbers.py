#!/usr/bin/env python3
"""Print figure 1 panel D/E's fitted numbers to full precision.

The panels render their slopes as five decimals inside a text box, which is
enough to read and not enough to quote. This prints what the panel computed:
the mean of the per-seed slopes, the SD ACROSS seeds that the wedge draws, the
honest SE of that mean, and the gap against the patient cohort.

It mirrors fig1 call for call -- same preprocessing functions, same
_through_origin, same per-seed grouping as fig1_plots.plot_clock_overlay --
so a number here cannot drift from the number in the figure. Nothing is
re-derived; in particular the intra-host clock still comes from
evolutionary_rate.extract_ih_regression_data via the preprocessor.

    python tests/test_fig1_clock_numbers.py --exp-num 1
"""
import argparse
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import fig1_preprocess_data as preproc          # noqa: E402
import fig1_plots as plots                      # noqa: E402
from _scenarios import EXP_NAME                 # noqa: E402

YEAR = 365.25


def per_seed_slopes(df, xcol='Sequencing_time', ycol='Distance_from_root'):
    """The per-seed slopes plot_clock_overlay averages, in subs/site/day."""
    out = []
    for _, g in df.groupby('seed'):
        s, _ = plots._through_origin(g[xcol].astype(float).values * YEAR,
                                     g[ycol].astype(float).values)
        if np.isfinite(s):
            out.append(s)
    return np.array(out)


def report(label, df, patient_rate):
    if df is None or df.empty:
        print(f"  {label:12s} -- no data")
        return
    slopes = per_seed_slopes(df)
    if not len(slopes):
        print(f"  {label:12s} -- no finite per-seed slope")
        return
    rate = slopes.mean() * YEAR                 # subs/site/year, as plotted
    sd = slopes.std(ddof=1) * YEAR if len(slopes) > 1 else 0.0
    se = sd / np.sqrt(len(slopes))
    pooled, r2 = plots._through_origin(
        df['Sequencing_time'].astype(float).values * YEAR,
        df['Distance_from_root'].astype(float).values)
    gap = 100 * (rate - patient_rate) / patient_rate if patient_rate else float('nan')
    inside = abs(rate - patient_rate) <= sd
    print(f"  {label:12s} n_seeds={len(slopes):3d}  points={len(df):7d}")
    print(f"    mean of per-seed slopes  {rate:.6f} s/s/y   (pooled {pooled * YEAR:.6f})")
    print(f"    SD across seeds (wedge)  {sd:.6f}          SE of mean {se:.6f}")
    print(f"    vs patient {patient_rate:.6f}         gap {gap:+.2f}%  "
          f"patient {'INSIDE' if inside else 'outside'} +/-1 SD")


def patient_rate(df, xcol, ycol):
    """The patient slope exactly as the panel draws it: x re-zeroed on the
    cohort's own first observation, then one fit through the origin."""
    r = df.dropna(subset=[xcol, ycol])
    x = r[xcol].astype(float).values
    y = r[ycol].astype(float).values
    s, r2 = plots._through_origin(x - x.min(), y)
    return s * YEAR, r2, len(r)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--exp-name', default=EXP_NAME)
    args = p.parse_args()

    df_real_std, df_real_long = preproc.get_patient_clock_data()

    print("\n=== D : global clock (standard) ===")
    pr, pr2, pn = patient_rate(df_real_std, 'days_since_ref',
                               'hamming_per_overlap')
    print(f"  patient      n={pn}  {pr:.6f} s/s/y (R2={pr2:.3f})")
    g = preproc.get_model_global_clock(exp_num=args.exp_num,
                                       exp_name=args.exp_name)
    for cohort in sorted(set(g.get('cohort', []))):
        report(cohort, g[g['cohort'] == cohort], pr)

    print("\n=== F : intra-host clock (long shedders) ===")
    pr, pr2, pn = patient_rate(df_real_long, 'time_point',
                               'hamming_distance_per_overlap')
    print(f"  patient      n={pn}  {pr:.6f} s/s/y (R2={pr2:.3f})")
    ih = preproc.get_model_intrahost_clock(exp_num=args.exp_num,
                                    exp_name=args.exp_name)
    for cohort in sorted(set(ih.get('cohort', []))):
        report(cohort, ih[ih['cohort'] == cohort], pr)

    print(f"\nratios: global {preproc.GLOBAL_CLOCK_RATIO} "
          f"intra-host {preproc.INTRAHOST_CLOCK_RATIO}")


if __name__ == '__main__':
    main()
