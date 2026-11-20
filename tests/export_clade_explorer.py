#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Export per-seed data for the interactive clade explorer.

The browser must be able to reproduce the pipeline's clustering EXACTLY, so
what gets shipped is the clustering's own decisions, not an approximation:

  parent / host / nedge / depth / tem   one row per lineage
  starts[ct]                            the nodes that START a clade at that
                                        shared-mutation threshold, taken from
                                        cluster_lin_into_clades_with_meta

From those the page rebuilds clade membership, definer sets, the majority vote
(with or without the parent-root leak), size-merging and colours -- all exact.
A depth-only shortcut was tested and rejected: it puts 11% of lineages in a
different clade on edge_case seed_0001, where back-mutations and second hits
make net divergence non-monotone.

Frequencies are downsampled in time and quantised to uint16; they feed a
stacked area plot, nothing numerical.

  python tests/export_clade_explorer.py --exp-num 1 --out <dir>
"""
import os
import sys
import json
import base64
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from _scenarios import EXP_NAME, resolve_scenarios, one_sod, scenario_meta

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
THRESHOLDS = [0, 3, 5, 7, 10, 15, 20]
N_TIMES = 120


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def _export_seed(job):
    scenario, ssod, outdir = job
    sys.setrecursionlimit(200000)
    seed = os.path.basename(ssod)
    tag = f"{scenario}__{seed}"
    try:
        ph = om.read_phylogenetic_data(ssod)
    except Exception as exc:
        return tag, 0, f"read failed: {exc}"

    # ---- lineage table, ordered parents-before-children -------------------
    ph = ph.sort_values("Time_emergence", kind="mergesort").reset_index(drop=True)
    names = ph["Lineage_name"].tolist()
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    if n < 2:
        return tag, 0, "fewer than 2 lineages"

    par = np.full(n, -1, np.int32)
    for i, p in enumerate(ph["Lineage_parent"]):
        if p is not None and not (isinstance(p, float) and np.isnan(p)):
            j = idx.get(p, -1)
            if j != i:
                par[i] = j
    # a parent must precede its child; emergence order should guarantee it
    bad = int(np.sum(par >= np.arange(n)))

    lin2mut = cl.build_lineage_to_mutation_dict(ph)
    depth = np.array([len(lin2mut[m]) for m in names], np.uint16)
    nedge = np.zeros(n, np.uint8)
    for i, m in enumerate(names):
        pm = lin2mut[names[par[i]]] if par[i] >= 0 else frozenset()
        nedge[i] = min(255, len(lin2mut[m] - pm))
    host = np.array([1 if str(h).strip().lower().startswith("long") else 0
                     for h in ph["Host_type"]], np.uint8)
    tem = np.clip(ph["Time_emergence"].to_numpy(float), 0, 65535).astype(np.uint16)

    # ---- clade boundaries, from the pipeline's own clustering -------------
    starts = {}
    for ct in THRESHOLDS:
        _c2l, _l2c, _pcm, meta = cl.cluster_lin_into_clades_with_meta(
            ph, shared_mut_threshold=ct)
        roots = [idx[r] for r in meta["root_lineage"] if r in idx]
        starts[str(ct)] = _b64(np.sort(np.array(roots, np.uint32)))

    # ---- infections attributable to each lineage --------------------------
    # Clade size in lineages counts intra-host diversification; long shedders
    # produce that by construction. The epidemiological size is how many
    # individuals carry the clade, so ship it alongside.
    try:
        ind = om.read_individuals_data(ssod)
        per = ind["inherited_lineage"].value_counts()
        hosts = np.array([int(per.get(m, 0)) for m in names], np.uint32)
    except Exception as exc:
        return tag, 0, f"individuals_data failed: {exc}"

    # ---- frequency, downsampled and quantised -----------------------------
    try:
        lf = om.read_lineage_frequency(ssod)
        F = lf.pivot(index="Time_sampling", columns="Lineage_name",
                     values="Frequency_at_t").sort_index()
        all_t = F.index.to_numpy(float)
        pick = np.unique(np.linspace(0, len(all_t) - 1,
                                     min(N_TIMES, len(all_t))).astype(int))
        F = F.iloc[pick]
        times = np.clip(all_t[pick], 0, 65535).astype(np.uint16)
        F = F.reindex(columns=names).fillna(0.0)
        M = np.clip(F.to_numpy(float), 0, 1)
        q = np.rint(M * 65535).astype(np.uint16)          # time x lineage
        nzt = q > 0
        first = np.where(nzt.any(axis=0), nzt.argmax(axis=0), -1)
        last = np.where(nzt.any(axis=0),
                        len(pick) - 1 - nzt[::-1].argmax(axis=0), -2)
        f_start = first.astype(np.int32)
        f_len = (last - first + 1).clip(0).astype(np.uint16)
        vals = np.concatenate([q[f_start[i]:f_start[i] + f_len[i], i]
                               for i in range(n) if f_len[i] > 0]) \
            if (f_len > 0).any() else np.zeros(0, np.uint16)
    except Exception as exc:
        return tag, 0, f"frequency failed: {exc}"

    obj = {
        "scenario": scenario, "seed": seed, "n": int(n),
        "thresholds": THRESHOLDS,
        "parent": _b64(par), "host": _b64(host), "nedge": _b64(nedge),
        "depth": _b64(depth), "tem": _b64(tem), "hosts": _b64(hosts),
        "starts": starts,
        "times": _b64(times),
        "f_start": _b64(f_start), "f_len": _b64(f_len),
        "f_vals": _b64(vals.astype(np.uint16)),
        "order_ok": bad == 0,
    }
    path = os.path.join(outdir, f"{tag}.json")
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    return tag, os.path.getsize(path), None


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=EXP_NAME)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    present, _ = resolve_scenarios(a.exp_name, a.exp_num, required=SCEN_ORDER)
    jobs, index = [], []
    for sc in [s for s in SCEN_ORDER if s in present]:
        seeds = dm.get_seeded_simulation_output_dirs(one_sod(a.exp_name, sc, a.exp_num))
        jobs += [(sc, s, a.out) for s in seeds[:a.seeds]]
    print(f"[export] {len(jobs)} seeds x {len(THRESHOLDS)} thresholds")

    total, failed = 0, []
    with Pool(a.procs) as pool:
        for i, (tag, size, err) in enumerate(
                pool.imap_unordered(_export_seed, jobs, chunksize=1), 1):
            if err:
                failed.append((tag, err))
            else:
                total += size
                index.append({"tag": tag, "scenario": tag.split("__")[0],
                              "seed": tag.split("__")[1], "bytes": size})
            if i % 25 == 0 or i == len(jobs):
                print(f"[export]   {i}/{len(jobs)}  {total/1e6:.1f} MB", flush=True)

    meta = scenario_meta(a.exp_name)
    with open(os.path.join(a.out, "index.json"), "w") as f:
        json.dump({"exp_name": a.exp_name, "exp_num": a.exp_num,
                   "thresholds": THRESHOLDS, "n_times": N_TIMES,
                   "scenarios": {s: {"ratio": meta.get(s, {}).get("ratio", 0.0),
                                     "duration": meta.get(s, {}).get("duration")}
                                 for s in SCEN_ORDER if s in present},
                   "seeds": sorted(index, key=lambda r: (r["scenario"], r["seed"]))},
                  f, separators=(",", ":"))
    print(f"[export] total {total/1e6:.1f} MB in {len(index)} files")
    for t, e in failed[:10]:
        print(f"[export][warn] {t}: {e}")
    if failed:
        print(f"[export][warn] {len(failed)} failed")
