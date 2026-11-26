#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does the emergence rise time measure the trajectory or the smoothing window?

The sweep showed the median emergence rise time moving 1 -> 25 -> 34 -> 37 days
as the smoothing window goes 1 -> 7 -> 15 -> 21 days, at every k. Two different
things could produce that, and the sweep cannot tell them apart:

  smoothing   a wider window flattens the approach to the peak, so the 10% and
              90% crossings of one episode move further apart
  merging     a wider window fuses neighbouring episodes (2.03 -> 1.37 episodes
              per clade over the same range), so the surviving episodes are
              bigger objects than the ones a narrow window measured

This separates them by FREEZING the episode set:

  free     detect episodes at each window and measure the emergence rise --
           what the panels do, and what the sweep reported
  frozen   detect episodes ONCE at --ref-win, keep that emergence episode's
           index window, and re-measure the rise inside it on trajectories
           smoothed at each other window

If frozen is flat, the window only matters because it merges episodes, and the
rise time of a given episode is a property of the trajectory. If frozen rises
like free, the window itself sets the timescale and the number is partly an
artefact of the choice.

`raw` is the same frozen window measured with no smoothing at all, which is the
shortest defensible answer for that episode.

Reads the daily lineage_frequency record through compute.run_context, the same
path the clade build uses, so these are the panels' own trajectories. Clades
are those rule B keeps at --mt. Read-only. Run from the repo root:

    python -I tests/rise_window_effect.py --seeds 10
"""
import os
import sys
import argparse

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in ("", "tests", "figures_gui",
           os.path.join("figures_gui", "prep")):
    sys.path.insert(0, os.path.join(ROOT, _p))

import simplicity.clustering as cl
import compute

SCEN = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
WINS = [1.0, 7.0, 15.0, 21.0]


def emergence(ys, t, prom, lo, hi):
    """(episode window, rise days) for the emergence episode, or (None, None)."""
    h, eps = compute.episodes_of(ys, prom)
    rs = [(ab, r) for ab, r in
          ((ab, compute.rise_of(ys, t, ab[0], ab[1], lo, hi)) for ab in eps)
          if r is not None]
    if not rs:
        return None, None, len(eps)
    imax = int(np.argmax(ys))
    ab, r = min(rs, key=lambda x: (not (x[1]["a"] <= imax <= x[1]["b"]),
                                   -x[1]["amp"]))
    return ab, float(r["days"]), len(rs)


def one_run(ssod, k, mt, prom, lo, hi, ref_win):
    ctx, err = compute.run_context(ssod, k, "post", mt)
    if err:
        return None, err
    F, t, post = ctx["F"], ctx["t_post"], ctx["post"]
    if len(t) < 3:
        return None, "window too short"
    tx, c2l = ctx["tx_all"], ctx["c2l"]
    recs = [dict(c=i, col=c, ndef=int(ctx["ndef_of"].get(c, 0) or 0),
                 ntx=int(sum(int(tx.get(m, 0)) for m in c2l.get(c, []))),
                 founding=not isinstance(ctx["par_of"].get(c), str))
            for i, c in enumerate(ctx["cols"])]
    kept = [r for r in cl.keep_designated_clades(recs)
            if not r["founding"] and r["ndef"] > 0]

    rows = []
    for r in kept:
        y = np.nan_to_num(F[r["col"]].to_numpy(float))[post]
        if not (y.max() > 0):
            continue
        # the frozen episode: detected once, at the reference window
        ref_ab, ref_days, _ = emergence(compute.smooth(y, t, ref_win), t,
                                        prom, lo, hi)
        for w in WINS:
            ys = compute.smooth(y, t, w)
            ab, days, neps = emergence(ys, t, prom, lo, hi)
            frozen = (compute.rise_of(ys, t, ref_ab[0], ref_ab[1], lo, hi)
                      if ref_ab is not None else None)
            raw = (compute.rise_of(y, t, ref_ab[0], ref_ab[1], lo, hi)
                   if ref_ab is not None else None)
            rows.append(dict(
                clade=str(r["col"]), win=w, n_eps=neps,
                free=days,
                frozen=None if frozen is None else float(frozen["days"]),
                raw=None if raw is None else float(raw["days"]),
                ref_days=ref_days))
    return rows, None


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stamp", default="latest")
    ap.add_argument("--scenarios", nargs="+", default=SCEN)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("-k", type=int, default=7)
    ap.add_argument("--mt", type=float, default=0.01)
    ap.add_argument("--prom", type=float, default=0.10)
    ap.add_argument("--lo", type=float, default=0.10)
    ap.add_argument("--hi", type=float, default=0.90)
    ap.add_argument("--ref-win", type=float, default=15.0)
    ap.add_argument("--csv", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "_rise_window_effect.csv"))
    a = ap.parse_args()
    from _scenarios import resolve_stamp
    a.stamp = resolve_stamp(a.stamp, 'impact_long_shedders')
    sys.setrecursionlimit(300000)

    allrows, failed = [], []
    for sc in a.scenarios:
        ssods = compute.dm.get_seeded_simulation_output_dirs(
            compute.one_sod(compute.EXP_NAME, sc, a.stamp))[:a.seeds]
        for ss in ssods:
            rows, err = one_run(ss, a.k, a.mt, a.prom, a.lo, a.hi, a.ref_win)
            if err:
                failed.append((sc, os.path.basename(ss), err))
                continue
            for r in rows:
                r.update(scenario=sc, seed=os.path.basename(ss))
            allrows.extend(rows)
        print(f"[rise] {sc}: {len(allrows)} rows so far", flush=True)

    df = pd.DataFrame(allrows)
    df.to_csv(a.csv, index=False)
    print(f"\n[rise] {len(df)} rows -> {a.csv}")
    if failed:
        print(f"[rise] {len(failed)} run(s) failed: {failed[:3]}")
    if df.empty:
        sys.exit(1)

    nclades = df.groupby(["scenario", "seed", "clade"]).ngroups
    print(f"\n=== k={a.k}, m_T={a.mt * 100:g}%, prom={a.prom}, band "
          f"{a.lo}-{a.hi}, frozen at win={a.ref_win:g} d ===")
    print(f"{nclades} clades over "
          f"{df.groupby(['scenario', 'seed']).ngroups} runs\n")

    g = df.groupby("win").agg(
        clades=("free", "size"),
        with_episode=("free", lambda s: int(s.notna().sum())),
        free_median=("free", "median"),
        frozen_median=("frozen", "median"),
        raw_median=("raw", "median"),
        eps_per_clade=("n_eps", "mean")).round(2)
    print("  per smoothing window, over every clade and run:")
    print(g.to_string())

    base = g.loc[WINS[0]]
    print("\n  change from win=1:")
    for w in WINS:
        r = g.loc[w]
        print(f"    win={w:4.0f}   free {r['free_median'] - base['free_median']:+7.1f} d"
              f"   frozen {r['frozen_median'] - base['frozen_median']:+7.1f} d"
              f"   raw {r['raw_median'] - base['raw_median']:+7.1f} d"
              f"   episodes/clade {r['eps_per_clade']:.2f}")

    print("\n  by scenario (median free / frozen):")
    h = df.pivot_table(index="scenario", columns="win",
                       values=["free", "frozen"], aggfunc="median").round(1)
    print(h.to_string())
