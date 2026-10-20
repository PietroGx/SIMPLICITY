#!/usr/bin/env python3
"""Panel E claims the model is unimodal and the cohort is not. Check it.

Panel E and tests/test_kl_root_to_tip.py read the same two frames and must
stay the same comparison: the KL number in the caption is computed by the
test, the picture by the figure. This exercises the FIGURE's own functions --
fig1_plots._lattice / _time_match / _density -- rather than reimplementing
them, so the two cannot drift apart silently.

    python tests/test_fig1_panel_e.py --exp-num 1
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, 'scripts', 'nature_plots'))

import fig1_plots as plots              # noqa: E402
import fig1_preprocess_data as preproc  # noqa: E402


def modes(centres, y, min_prom=0.15):
    """Local maxima of a density, keeping only those with real prominence."""
    peak = y.max()
    out = []
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] >= min_prom * peak:
            lo = y[:i].min() if i else y[i]
            hi = y[i + 1:].min() if i + 1 < len(y) else y[i]
            if y[i] - max(lo, hi) >= min_prom * peak:
                out.append((centres[i], y[i]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    args = p.parse_args()

    real, _ = preproc.get_patient_clock_data()
    model = preproc.get_model_global_clock(exp_num=args.exp_num)
    real = real.dropna(subset=['days_since_ref', 'hamming_per_overlap'])

    tp = real['days_since_ref'].astype(float).values
    tp = tp - tp.min()
    dp = real['hamming_per_overlap'].astype(float).values
    res = 1.0 / float(np.median(real['overlap_length']))

    rng = np.random.default_rng(4242)
    tm = model['Sequencing_time'].astype(float).values * 365.25
    dmv = model['Distance_from_root'].astype(float).values
    tm, dmv = plots._time_match(tm, dmv, tp, 10, rng)

    step = plots._lattice(dmv)
    print(f"model lattice : {step:.8f}  (1/{1/step:.0f}), "
          f"{len(np.unique(np.round(dmv, 12)))} distinct values, n={len(dmv)}")
    print(f"patient res   : {res:.8f}  (1/{1/res:.0f}), "
          f"{len(np.unique(np.round(dp, 12)))} distinct values, n={len(dp)}")
    assert plots._lattice(dp) is None, "patient data should have no lattice"

    # time-matching must leave the model on the cohort's schedule
    qs = [10, 25, 50, 75, 90]
    print("\nsampling time deciles after matching (days):")
    print("   patient " + " ".join(f"{np.percentile(tp, q):7.0f}" for q in qs))
    print("   model   " + " ".join(f"{np.percentile(tm, q):7.0f}" for q in qs))
    worst = max(abs(np.percentile(tp, q) - np.percentile(tm, q)) for q in qs)
    assert worst < 20, f"time matching off by {worst:.0f} d"

    hi = max(dp.max(), dmv.max())
    print("\nmodes of each density (subs/1000 sites):")
    found = {}
    for name, vals, width in (('patient', dp, res * 2), ('model', dmv, step * 2)):
        c, y = plots._density(vals, width, 0.0, hi)
        g = np.linspace(0, hi, 600)
        sm = plots._smoothed(vals, g, lattice=(step if name == 'model' else None))
        found[name] = modes(g, sm)
        peaks = ", ".join(f"{x*1e3:.2f}" for x, _ in found[name])
        print(f"   {name:>7}: {len(found[name])} mode(s) at [{peaks}]")

    assert len(found['model']) == 1, \
        f"model should be unimodal, found {len(found['model'])}"
    assert len(found['patient']) == 2, \
        f"cohort should be bimodal, found {len(found['patient'])}"

    lo_pk, hi_pk = sorted(x for x, _ in found['patient'])
    print(f"\ncohort modes {lo_pk*1e3:.2f} and {hi_pk*1e3:.2f} /1000; "
          f"model single mode {found['model'][0][0]*1e3:.2f} /1000")
    print("OK")


if __name__ == '__main__':
    main()
