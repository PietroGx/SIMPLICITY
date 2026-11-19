#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is cluster_lin_into_clades_with_meta building sensible clades?

Structural checks only -- no labelling, no outcomes. Each one is a property
the clustering should have if a clade is to stand in for a real-world lineage:

  partition   every lineage in exactly one clade, nothing lost
  size        lineages per clade; a clustering of singletons aggregates nothing
  bound       within-clade divergence from the clade root. Under rule B this
              is NOT bounded by the threshold and is not meant to be: a lineage
              refused for want of transmissions folds into its ancestor and
              brings its divergence with it (BACKLOG, "the transmission
              minimum breaks the k-1 divergence bound -- ACCEPTED, no cap"), so
              the column is reported, not asserted
  definers    is any lineage's mutation counted as a definer of more than one
              clade? (siblings crossing the threshold inherit the same
              path_since_root, clustering.py:255-257)
  frequency   do clade frequencies account for the population at each time?
  depth       nesting depth of the clade tree

Read-only. Run from the repo root.
"""
import os
import sys
import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from _scenarios import EXP_NAME, one_sod

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]


def check(ssod, ct):
    sys.setrecursionlimit(100000)
    ph = om.read_phylogenetic_data(ssod)
    # Rule B, as the pipeline designates (clustering.get_clustering_for_ssod).
    # This measured the divergence-only partition until v2.5.6, which is no
    # longer what anything produces: at k=5 on control it reported 335 clades
    # and 87.8% singletons where rule B gives 14 and 0%. A structural check of
    # a partition nothing builds is worse than no check.
    ind = om.read_individuals_data(ssod)
    support = cl.subtree_transmission_counts(ph, ind)
    min_support = cl.transmission_floor(ind, cl.MIN_TX_FRAC)
    c2l, l2c, pcm, meta = cl.cluster_lin_into_clades_with_meta(
        ph, shared_mut_threshold=ct, support=support, min_support=min_support)
    lin2mut = cl.build_lineage_to_mutation_dict(ph)
    r = {"threshold": ct, "n_lineages": len(ph), "n_clades": len(c2l)}

    # -- partition ---------------------------------------------------------
    assigned = [l for members in c2l.values() for l in members]
    r["assigned"] = len(assigned)
    r["unassigned"] = len(ph) - len(set(assigned))
    r["multi_assigned"] = len(assigned) - len(set(assigned))

    # -- size --------------------------------------------------------------
    sizes = np.array([len(v) for v in c2l.values()])
    r["size_median"] = float(np.median(sizes))
    r["size_p90"] = float(np.percentile(sizes, 90))
    r["size_max"] = int(sizes.max())
    r["singleton_%"] = float(100 * (sizes == 1).mean())

    # -- divergence bound --------------------------------------------------
    root_of = dict(zip(meta["clade"], meta["root_lineage"]))
    worst, over = 0, 0
    for clade, members in c2l.items():
        rm = lin2mut.get(root_of.get(clade), frozenset())
        d = [len(lin2mut.get(m, frozenset()) - rm) for m in members]
        if d:
            worst = max(worst, max(d))
            over += sum(1 for x in d if x >= ct)
    r["max_div_from_root"] = worst
    r["members_at_or_over_thr"] = over

    # -- definer attribution ----------------------------------------------
    # One mutation event = one lineage. Count the lineages that appear as a
    # definer source for more than one clade.
    src = defaultdict(set)
    recs = 0
    for clade, d in pcm.items():
        if d is None or d.empty:
            continue
        recs += len(d)
        for lin in d["emergence_lineage"]:
            src[lin].add(clade)
        r.setdefault("_", None)
    multi = {l: cs for l, cs in src.items() if len(cs) > 1}
    r["definer_records"] = recs
    r["definer_lineages"] = len(src)
    r["lineages_defining_>1_clade"] = len(multi)
    r["worst_lineage_in_n_clades"] = max([len(c) for c in multi.values()], default=1)
    r["definer_records_double_counted_%"] = float(
        100 * sum(len(cs) - 1 for cs in multi.values()) / recs) if recs else 0.0

    # -- frequency bookkeeping --------------------------------------------
    lf = om.read_lineage_frequency(ssod)
    F_lin = lf.pivot(index="Time_sampling", columns="Lineage_name",
                     values="Frequency_at_t").sort_index()
    series = []
    for clade, members in c2l.items():
        cols = [c for c in members if c in F_lin.columns]
        if cols:
            series.append(F_lin[cols].sum(axis=1).rename(clade))
    F_cl = pd.concat(series, axis=1) if series else pd.DataFrame()
    tot = F_cl.sum(axis=1).to_numpy(float) if not F_cl.empty else np.array([np.nan])
    r["clade_freq_sum_median"] = float(np.nanmedian(tot))
    r["clade_freq_sum_min"] = float(np.nanmin(tot))
    r["lineages_missing_from_freq_%"] = float(
        100 * (1 - F_lin.shape[1] / max(len(ph), 1)))

    # -- depth -------------------------------------------------------------
    depth = [str(c).count(".") for c in c2l]
    r["depth_median"] = float(np.median(depth))
    r["depth_max"] = int(max(depth))
    r.pop("_", None)
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=EXP_NAME)
    ap.add_argument("--scenarios", nargs="+", default=["control", "HIV_high"])
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--thresholds", type=int, nargs="+", default=[5, 10])
    a = ap.parse_args()

    rows = []
    for sc in a.scenarios:
        sod = one_sod(a.exp_name, sc, a.exp_num)
        for ss in dm.get_seeded_simulation_output_dirs(sod)[:a.seeds]:
            for ct in a.thresholds:
                r = check(ss, ct)
                r["scenario"], r["seed"] = sc, os.path.basename(ss)
                rows.append(r)
                print(f"[done] {sc} {os.path.basename(ss)} ct={ct}", flush=True)
    df = pd.DataFrame(rows)
    cols = [c for c in df.columns if c not in ("scenario", "seed")]
    pd.set_option("display.width", 200)
    print()
    print(df.set_index(["scenario", "seed", "threshold"])[
        [c for c in cols if c != "threshold"]].T.to_string())
