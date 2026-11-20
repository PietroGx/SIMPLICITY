#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""How k and the minimum transmissions m_T shape clade formation -- the grid the
clustering page shows.

Rule B: a lineage founds a clade when it has diverged >= k mutations from its
clade root AND its subtree has produced >= m_tx of the seed's transmission
events. m_T is a SHARE, not a count: the scenarios differ 2.5x in total
transmissions (8.6k HIV_high to 22k edge_case), so an absolute minimum is a
different filter in each.

Cumulative subtree transmissions fall monotonically down the tree, so a
refusal is hereditary: the refused subtree stays in the ancestor clade, and
brings its accumulated divergence with it. That is why max_div breaks the
k-1 bound, and by more in the scenarios where long shedders generate
divergence without transmission.

Emits per (scenario, seed, k, m) rows. --json writes the per-seed MEDIAN grid
the clustering page shows; --perseed writes the PER-SEED payload the
across-the-set view draws, because the medians hide most of what is there: at
control k=5 with no minimum the 50 repeats run from 52 to 1152 clades, an IQR
wider than the median itself. Run from the repo root.
"""
import os
import sys
import json
import argparse
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from _scenarios import EXP_NAME, one_sod

SCEN = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
KS = [0, 3, 5, 7, 10, 15, 20]
# a 1-2-5 progression from "one transmission in a 20k-event run" up to 2%, densest
# over 0.003-0.2% where the singletons actually clear out
# fractions since v2.5.15; the sweep itself is unchanged, every value is the
# old percentage / 100.
MS = [0.0, 0.00002, 0.00003, 0.00005, 0.00008, 0.00012, 0.0002, 0.0003,
      0.0005, 0.0008, 0.0012, 0.002, 0.003, 0.005, 0.01, 0.02]
# clade-size histogram, doubling bins: the median of a size distribution is a
# poor summary of it (one clade of 1150 lineages also reports "median 1150")
SZ_EDGES = [1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513]




def summarise(lin2mut, croot, order, stx):
    mem = defaultdict(list)
    for v in order:
        mem[croot[v]].append(v)
    sizes = np.array([len(v) for v in mem.values()], float)
    worst = max(max(len(lin2mut[x] - lin2mut[r]) for x in v) for r, v in mem.items())
    txs = np.array([stx.get(r, 0) for r in mem], float)
    h = np.bincount(np.searchsorted(SZ_EDGES, sizes.astype(int), side="right") - 1,
                    minlength=len(SZ_EDGES))
    return dict(clades=len(mem), singleton=100.0 * float((sizes == 1).mean()),
                med_size=float(np.median(sizes)), max_div=int(worst),
                dead=100.0 * float((txs == 0).mean()), med_tx=float(np.median(txs)),
                hist=[int(x) for x in h])


def _one(job):
    scenario, ssod = job
    sys.setrecursionlimit(300000)
    try:
        ph = om.read_phylogenetic_data(ssod)
        ind = om.read_individuals_data(ssod)
    except Exception as exc:
        return scenario, os.path.basename(ssod), None, f"read: {exc}"
    lin2mut, par, order, _sub, stx = cl.lineage_tree(ph, ind, order_by="topology")
    tot = cl.transmission_total(ind)
    out = []
    for k in KS:
        for m in MS:
            croot = cl.designate_clades(lin2mut, par, order, k, stx,
                                        m * tot)
            r = summarise(lin2mut, croot, order, stx)
            r.update(k=k, m=m, tot_tx=tot)
            out.append(r)
    return scenario, os.path.basename(ssod), out, None


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--scenarios", nargs="+", default=SCEN)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--csv", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "_clade_formation_grid.csv"))
    ap.add_argument("--json", default="")
    ap.add_argument("--perseed", default="")
    a = ap.parse_args()

    jobs = []
    for sc in a.scenarios:
        for ss in dm.get_seeded_simulation_output_dirs(one_sod(EXP_NAME, sc, a.exp_num))[:a.seeds]:
            jobs.append((sc, ss))
    print(f"[grid] {len(jobs)} seeds x {len(KS)} k x {len(MS)} m", flush=True)

    rows, failed = [], []
    with Pool(a.procs) as pool:
        for i, (sc, sd, out, err) in enumerate(pool.imap_unordered(_one, jobs), 1):
            if err:
                failed.append((sc, sd, err))
            else:
                for r in out:
                    r.update(scenario=sc, seed=sd)
                    rows.append(r)
            if i % 5 == 0 or i == len(jobs):
                print(f"[grid]   {i}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows)
    df.drop(columns=["hist"]).to_csv(a.csv, index=False)
    print(f"[grid] {len(df)} rows -> {a.csv}")
    for f in failed:
        print("  FAILED", f)

    if a.json:
        g = (df.groupby(["scenario", "k", "m"])
               .agg(clades=("clades", "median"), singleton=("singleton", "median"),
                    med_size=("med_size", "median"), max_div=("max_div", "max"),
                    dead=("dead", "median"), tot_tx=("tot_tx", "median"),
                    n_seeds=("seed", "nunique")).reset_index())
        obj = {"ks": KS, "ms": MS, "scenarios": a.scenarios,
               "n_seeds": int(g["n_seeds"].min()),
               "grid": {f"{r.scenario}|{r.k}|{r.m}": [
                   int(r.clades), round(r.singleton, 1), round(r.med_size, 1),
                   int(r.max_div), round(r.dead, 1), int(r.tot_tx)]
                   for r in g.itertuples()}}
        with open(a.json, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"))
        print(f"[grid] grid -> {a.json} ({os.path.getsize(a.json)/1000:.1f} kB)")

    if a.perseed:
        import base64
        # [h0..h10, max_div, dead per mille, med_size x2] per cell, one uint16
        # block per scenario shaped [seed][k][m][F]. clades is the bin sum and
        # singletons is bin 0, so neither is stored. x2, not x10: a median of
        # integers is always a half-integer, so x2 is exact, and x10 overflowed
        # uint16 on the biggest epidemics (SOT seed_0019 runs 66,852
        # transmissions over 6,881 lineages).
        F = len(SZ_EDGES) + 3
        seeds = {sc: sorted({r["seed"] for r in rows if r["scenario"] == sc}) for sc in a.scenarios}
        by = {(r["scenario"], r["seed"], r["k"], r["m"]): r for r in rows}
        blocks = {}
        for sc in a.scenarios:
            buf = np.zeros((len(seeds[sc]), len(KS), len(MS), F), np.uint16)
            for si, sd in enumerate(seeds[sc]):
                for ki, k in enumerate(KS):
                    for mi, m in enumerate(MS):
                        r = by.get((sc, sd, k, m))
                        if r is None:
                            continue
                        v = r["hist"] + [r["max_div"], round(r["dead"] * 10),
                                         round(r["med_size"] * 2)]
                        over = [x for x in v if x > 65535]
                        if over:
                            raise SystemExit(f"[grid] uint16 overflow {sc} {sd} "
                                             f"k={k} m={m}: {over}")
                        buf[si, ki, mi] = v
            blocks[sc] = base64.b64encode(buf.tobytes()).decode()
        obj = {"ks": KS, "ms": MS, "scenarios": a.scenarios, "edges": SZ_EDGES,
               "fields": F, "seeds": {sc: [x.replace("seed_", "") for x in seeds[sc]]
                                      for sc in a.scenarios},
               "data": blocks}
        with open(a.perseed, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"))
        print(f"[grid] per-seed -> {a.perseed} "
              f"({os.path.getsize(a.perseed)/1e6:.2f} MB, {F} fields x "
              f"{sum(len(v) for v in seeds.values())} seeds x {len(KS)}k x {len(MS)}m)")
