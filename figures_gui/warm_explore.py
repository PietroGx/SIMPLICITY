#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Precompute what the Sweep tab's exploration and Figure 5's live panel read,
so a slider answers in seconds instead of re-clustering 250 runs.

The clustering is the only expensive step, and it depends on k and m_T alone.
Per run this fills two disk caches, both keyed on the code that writes them:

    daily   the run's per-clade daily record (panels._daily_series memo)
    meta    both hidden-time anchors and each clade's long-born definer count
            (explore.run_meta, <store>/cache/_explore/)

about 4 s a run cold, 2 of them each. Then the Figure 5 panel at Figure 5's
own rules is warmed, so the tab opens on a cached answer. Detection itself is
not precomputed: it costs milliseconds and its rules are what you explore.

    python -I figures_gui/warm_explore.py --stamp 2 --procs 3
    python -I figures_gui/warm_explore.py --stamp 2 --k 7 --mt 0.01 --shard 3/8

k and m_T default to the studio's operating point (web/settings.js).
--pairs all does every (k, m_T) of the build grid instead, which S2 (the
clustering sensitivity tab) measures Figure 5 at: 21 pairs, 21 times the work.
After the runs, the Figure 5 panel is warmed at each pair, and at the operating
point once per server-side level of S3 (sa_levels.json: smoothing window,
prominence, ceiling f); S3's other rules are applied in the browser.

--shard I/N takes every N-th job from I; figures_gui/slurm_explore.sh runs the
shards as an array, then one --panel-only task for the panels, and
run_precompute.sh submits both with the grid.
"""
import json
import argparse
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.chdir(ROOT)

import compute                                                    # noqa: E402
import panels                                                     # noqa: E402
import sweep                                                      # noqa: E402
import explore                                                    # noqa: E402

FIG5_Q = dict(  # the detection tabs' defaults, which Figure 5 measures with
    k=None, mt=None, win="15", prom="0.1", lo="0.1", burnin="1",
              kernel="box", space="lin", prom_abs="0", sep="0", width="0", cens="1",
              start="trough", trough=str(explore.FIG5["trough"]), anchor="emergence",
              gate="start", dip="low")
# S3's levels whose answer is a different /api/explore_panel request, and the
# query parameter each one sets (web/episodes.js exploreQuery)
SA_SERVER = {"sdWin": "win", "sdProm": "prom", "brkLo": "under"}


def sa_levels():
    with open(os.path.join(HERE, "sa_levels.json")) as fh:
        return json.load(fh)["outbreak"]


def grid_pairs():
    return [(k, float(m)) for k in sweep.K for m in compute.MT_VALUES]


def panels_to_warm(pairs, k0, mt0):
    """Every Figure 5 panel query the sensitivity tabs and Figure 5 ask for at
    the default rules: each pair at Figure 5's rules, then S3's server-side
    levels at the operating point."""
    qs = [dict(FIG5_Q, k=str(k), mt=str(mt)) for k, mt in pairs]
    for name, param in SA_SERVER.items():
        qs += [dict(FIG5_Q, k=str(k0), mt=str(mt0), **{param: str(v)})
               for v in sa_levels()[name]]
    seen, out = set(), []
    for q in qs:
        key = json.dumps(explore.read_params(q), sort_keys=True)
        if key not in seen:
            seen.add(key)
            out.append(q)
    return out


def _one(job):
    store, stamp, name, scen, seed, k, mt = job
    panels.EXP_NUM, panels.EXP_NAME = stamp, name
    t0 = time.time()
    try:
        explore.pool(store, scen, seed, k, mt)
        return scen, seed, round(time.time() - t0, 1), None
    except Exception as exc:                                     # noqa: BLE001
        return scen, seed, round(time.time() - t0, 1), f"{type(exc).__name__}: {exc}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stamp", default="latest")
    ap.add_argument("--exp-name", default=compute.EXP_NAME)
    ap.add_argument("--store", default=None, help="default: the experiment's own store")
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument("--mt", type=float, default=0.01)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--pairs", choices=("op", "all"), default="op",
                    help="op: --k and --mt only; all: every (k, m_T) of the build grid")
    ap.add_argument("--no-panel", action="store_true",
                    help="skip warming the Figure 5 panels")
    ap.add_argument("--panel-only", action="store_true",
                    help="only the Figure 5 panels; the runs were warmed by the shards")
    a = ap.parse_args()
    from _scenarios import resolve_stamp
    a.stamp = resolve_stamp(a.stamp, a.exp_name)
    panels.EXP_NUM, panels.EXP_NAME = a.stamp, a.exp_name
    store = os.path.expanduser(a.store) if a.store else \
        sweep.store_dir(sweep.DEFAULT_STORE, a.stamp, a.exp_name)
    si, sn = (int(x) for x in a.shard.split("/"))

    pairs = grid_pairs() if a.pairs == "all" else [(a.k, a.mt)]
    jobs = [] if a.panel_only else [
        (store, a.stamp, a.exp_name, sc, seed, k, mt)
        for k, mt in pairs for sc in compute.SCEN_ORDER
        for seed in explore._seed_names(sc)[:a.seeds]][si::sn]
    print(f"[warm_explore] {a.exp_name} #{a.stamp}  {len(pairs)} (k, m_T) pair(s)  "
          f"{len(jobs)} runs (shard {si}/{sn})  store={store}", flush=True)
    t0, bad = time.time(), []
    with Pool(max(1, a.procs)) as pool:
        for i, (sc, seed, s, err) in enumerate(pool.imap_unordered(_one, jobs), 1):
            if err:
                bad.append(f"{sc}/{seed}: {err}")
            if err or i % 25 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  {sc}/{seed}  {s} s{'  ' + err if err else ''}"
                      f"   elapsed {time.time() - t0:.0f} s", flush=True)
    if not a.no_panel and sn == 1 and not bad:
        qs = panels_to_warm(pairs, a.k, a.mt)
        for i, q in enumerate(qs, 1):
            t1 = time.time()
            try:
                explore.explore_panel(store, dict(q, seeds=str(a.seeds)))
                err = ""
            except Exception as exc:                             # noqa: BLE001
                err = f"  FAILED {type(exc).__name__}: {exc}"
                bad.append(f"panel k={q['k']} m_T={q['mt']}: {exc}")
            print(f"[warm_explore] Figure 5 panel {i}/{len(qs)} k={q['k']} m_T={q['mt']} "
                  f"win={q['win']} prom={q['prom']} f={q.get('under', 'default')}: "
                  f"{time.time() - t1:.1f} s{err}", flush=True)
    print(f"[warm_explore] {len(jobs) - len(bad)}/{len(jobs)} runs in "
          f"{time.time() - t0:.0f} s{'; FAILED: ' + '; '.join(bad[:5]) if bad else ''}",
          flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
