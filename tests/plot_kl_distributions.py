#!/usr/bin/env python3
"""Figure-1-panel-D's two root-to-tip distributions, drawn honestly.

The left panel of the KL figure binned both sides at the PATIENT's median
resolution, 1/3822. The model steps by 1/4967, so roughly every third bin
swallowed two model lattice points and the model acquired a zigzag that looked
like bimodality. It is an artefact of the binning, not a property of the model.

This draws the same two distributions without that artefact:

  model     binned on ITS OWN lattice (width 1/4967, or a multiple), so every
            bin holds exactly the same number of lattice points
  patient   binned freely -- it has no lattice. Each patient sequence has its
            own overlap length, so `hamming_per_overlap` takes 231 distinct
            values over 1470 sequences; 1/3822 is only the median denominator.

Both are drawn as DENSITY, per unit divergence, which is what makes two
different bin widths comparable at all. The KDE overlays are a visual guide to
shape; the KL is computed on the histograms, never on them.

Everything upstream of the plot -- loading, coverage filter, time-matching --
is imported from test_kl_root_to_tip so the panel cannot drift from the number.

    python tests/plot_kl_distributions.py --exp-num 1 --out kl_panel.png
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))

import test_kl_root_to_tip as K            # noqa: E402

MODEL_RES = 1.0 / 4967.0                   # the reference genome length
PATIENT_COLOR = '#1b9e77'
MODEL_COLOR = '#333333'


def density(values, width, lo, hi):
    """Histogram as a density, so two different bin widths are comparable."""
    edges = np.arange(lo, hi + width, width)
    counts, _ = np.histogram(values, bins=edges)
    centres = edges[:-1] + width / 2
    return centres, counts / (counts.sum() * width)


def kde(values, grid, lattice=None):
    """Gaussian KDE, for the eye only. Returns None if scipy is absent.

    The bandwidth is floored at 1.2x the lattice spacing where the variable is
    discrete. Scott's rule on 145k model points gives a bandwidth SMALLER than
    the 1/4967 step, so the curve ripples once per lattice point -- drawing the
    discreteness as if it were structure, which is the same mistake the
    patient-width histogram made.
    """
    try:
        from scipy.stats import gaussian_kde
    except ImportError:
        return None
    sigma = float(np.std(values))
    bw = len(values) ** (-1.0 / 5.0)                     # Scott
    if lattice and sigma > 0:
        bw = max(bw, 1.2 * lattice / sigma)
    k = gaussian_kde(values, bw_method=bw)
    return k(grid)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--xmax', type=float, default=4.5,
                   help="x limit in substitutions per 1000 sites")
    p.add_argument('--model-bins', type=int, default=2,
                   help="model bin width, in multiples of 1/4967")
    p.add_argument('--patient-bins', type=int, default=2,
                   help="patient bin width, in multiples of its median 1/overlap")
    p.add_argument('--rng-seed', type=int, default=4242)
    p.add_argument('--no-kde', action='store_true')
    args = p.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    pat, res = K.patient_frame()
    seeds = K.model_seeds(args.exp_num)
    if not seeds:
        raise SystemExit("no model seeds")

    # time-match the pooled model to the patient schedule, as compare() does
    rng = np.random.default_rng(args.rng_seed)
    q = np.linspace(0, 1, K.N_TIME_BINS + 1)
    edges_t = np.quantile(pat['t'].values, q)
    want = np.histogram(pat['t'].values, bins=edges_t)[0]
    edges_t[0], edges_t[-1] = -np.inf, np.inf
    rows = [r for _, rr in seeds for r in rr]
    _, dm = K.time_match(rows, edges_t, want * (len(rows) // len(pat)), rng)
    dp = pat['d'].values

    print(f"model   : {len(dm)} time-matched, {len(np.unique(np.round(dm, 12)))} "
          f"distinct values on a 1/4967 lattice")
    print(f"patient : {len(dp)}, {len(np.unique(np.round(dp, 12)))} distinct "
          f"values -- no lattice")

    lo, hi = 0.0, max(dm.max(), dp.max()) + MODEL_RES
    wm = MODEL_RES * args.model_bins
    wp = res * args.patient_bins
    cm, ym = density(dm, wm, lo, hi)
    cp, yp = density(dp, wp, lo, hi)

    fig, ax = plt.subplots(figsize=(7.4, 4.6))

    ax.fill_between(cm * 1e3, 0, ym / 1e3, step='mid',
                    color=MODEL_COLOR, alpha=0.13, lw=0)
    ax.step(cm * 1e3, ym / 1e3, where='mid', color=MODEL_COLOR, lw=1.4,
            label=f"model, time-matched (n={len(dm):,}, bin {args.model_bins}/4967)")
    ax.step(cp * 1e3, yp / 1e3, where='mid', color=PATIENT_COLOR, lw=1.6,
            ls='--',
            label=f"patient cohort (n={len(dp):,}, bin {args.patient_bins}/3822)")

    if not args.no_kde:
        grid = np.linspace(lo, hi, 600)
        for vals, colour, name, lat in ((dm, MODEL_COLOR, 'model', MODEL_RES),
                                        (dp, PATIENT_COLOR, 'patient', res)):
            y = kde(vals, grid, lattice=lat)
            if y is None:
                print("[warn] scipy missing; KDE overlay skipped")
                break
            ax.plot(grid * 1e3, y / 1e3, color=colour, lw=2.2, alpha=0.55,
                    label=f"{name}, smoothed (guide only)")

    ax.set_xlim(0, args.xmax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('root-to-tip divergence (substitutions per 1000 sites)')
    ax.set_ylabel('density (per 1000 sites)')
    ax.set_title('Root-to-tip divergence, patients vs model, 100% sampling',
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches='tight')
    print(f"[saved] {args.out}")


if __name__ == '__main__':
    main()
