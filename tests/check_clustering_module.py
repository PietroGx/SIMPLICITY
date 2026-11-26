#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Do simplicity/clustering.py's primitives agree with what they replaced?

The clustering primitives were extracted so each rule has one statement
(v2.5.4). Extraction is only safe while the copies it replaced still agree, and
two of those copies are out of reach:

  resolve_clade_label was one of two statements of the origin rule until
  v2.5.6, the other living in scripts/long_paper_figures. That directory is
  archived and the module's is now the only copy, so section A pins its
  BEHAVIOUR against an expected table instead of against a second
  implementation -- including the orderings ("founder long") that a reordered
  if-chain would get wrong.

  cluster_lin_into_clades_with_meta and designate_clades are two walks over one
  rule. Both call founds_clade, so they cannot disagree by construction -- but
  "cannot by construction" is what the four bad production runs were. Section B
  measures it on real runs instead.

Read-only. Run from the repo root:

    python -I tests/check_clustering_module.py [--seeds 2]
"""
import os
import sys
import argparse

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "figures_gui", "prep"))
sys.path.insert(0, os.path.join(ROOT, "tests"))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from _scenarios import EXP_NAME, one_sod

PASS = [0]
FAIL = [0]


def ck(cond, msg):
    (PASS if cond else FAIL)[0] += 1
    print(("  ok   " if cond else "  FAIL ") + msg)


# --------------------------------------------------------------- section A
# the copies in scripts/ that the module restates

# (clade_labels, clade, expected). The expected column is the behaviour the
# archived scripts/long_paper_figures copy had, measured against it on
# 2026-11-18 over these same cases with 0 disagreements, before it was removed.
LABEL_CASES = [
    ({}, "", "standard"),
    ({}, None, "standard"),
    ({}, "Clade_1", "standard"),
    ({"Clade_1": "long"}, "Clade_1", "long"),
    ({"Clade_1": "long_shedder"}, "Clade_1", "long"),
    ({"Clade_1": "LONG"}, "Clade_1", "long"),
    ({"Clade_1": "mixed"}, "Clade_1", "mixed"),
    ({"Clade_1": "MIX"}, "Clade_1", "mixed"),
    ({"Clade_1": "founder"}, "Clade_1", "founder"),
    ({"Clade_1": "standard"}, "Clade_1", "standard"),
    ({"Clade_1": "nonsense"}, "Clade_1", "standard"),
    ({"Clade_1": None}, "Clade_1", "standard"),
    ({"Clade_1": 7}, "Clade_1", "standard"),
    # orderings that change answer if the if-chain is reordered
    ({"Clade_1": "long mixed"}, "Clade_1", "long"),
    ({"Clade_1": "mixed founder"}, "Clade_1", "mixed"),
    ({"Clade_1": "founder long"}, "Clade_1", "long"),
]

HOST_CASES = [None, "long_shedder_i", "Long", "LONG_x", "standard", "",
              " long ", 0, "nonlong", "longitudinal"]


def section_a():
    print("\nA. the module against the behaviour it inherited")
    bad = [(lab, c, want, cl.resolve_clade_label(lab, c))
           for lab, c, want in LABEL_CASES
           if cl.resolve_clade_label(lab, c) != want]
    ck(not bad,
       f"resolve_clade_label gives the inherited answer on "
       f"{len(LABEL_CASES)} cases"
       + (f" -- differs on {bad}" if bad else ""))

    # the prefix rule the four inline copies all used
    bad = [h for h in HOST_CASES
           if (cl.normalise_host_type(h) == "long")
           != str(h).strip().lower().startswith("long")]
    ck(not bad,
       f"normalise_host_type is the \"long\" prefix rule on {len(HOST_CASES)} "
       f"cases" + (f" -- differs on {bad}" if bad else ""))

    ck(cl.normalise_host_type(None) == "standard",
       "normalise_host_type(None) is standard, not a crash")

    # m_T is a share of the run's own transmissions, so the floor must scale
    class _Ind(dict):
        pass
    import pandas as pd
    ind = pd.DataFrame({"inherited_lineage": ["a", "a", "b", "c"] * 25 + ["wt"] * 5,
                        "parent": ["7"] * 100 + ["root"] * 5})
    ck(cl.transmission_total(ind) == 100,
       f"transmission_total counts transmission events, not the 5 seeded "
       f"infections (got {cl.transmission_total(ind)}, want 100)")
    ck(abs(cl.transmission_floor(ind, 0.01) - 1.0) < 1e-9,
       "transmission_floor turns the fraction 0.01 into 1.0 of 100 transmissions")
    ck(abs(cl.transmission_floor(ind, 1.0) - 100.0) < 1e-9,
       "and the fraction 1.0 into all 100 -- it is a fraction, not a percentage")
    ck(cl.transmission_floor(ind, 0.0) == 0.0,
       "transmission_floor(0) is 0, so m_T=0 keeps the divergence-only set")


# --------------------------------------------------------------- section B
# the two walks over rule B, measured on real runs

def section_b(scenarios, nseeds, k, stamp=1, exp_name=EXP_NAME):
    print(f"\nB. designate_clades against cluster_lin_into_clades_with_meta "
          f"(k={k})")
    checked = 0
    for sc in scenarios:
        try:
            ssods = dm.get_seeded_simulation_output_dirs(
                one_sod(exp_name, sc, stamp))[:nseeds]
        except Exception as exc:
            print(f"  ({sc} unavailable: {exc})")
            continue
        for ss in ssods:
            name = f"{sc}/{os.path.basename(ss)}"
            ph = om.read_phylogenetic_data(ss)
            ind = om.read_individuals_data(ss)

            c2l, l2c, _pcm, _meta = cl.cluster_lin_into_clades_with_meta(ph, k)
            lin2mut, par, order, _sub, stx = cl.lineage_tree(ph, ind)
            croot = cl.designate_clades(lin2mut, par, order, k)

            # same partition, named differently: compare the GROUPINGS
            builder = {}
            for lin, clade in l2c.items():
                builder.setdefault(clade, set()).add(lin)
            walk = {}
            for lin, root in croot.items():
                walk.setdefault(root, set()).add(lin)
            same = (sorted(map(sorted, builder.values()))
                    == sorted(map(sorted, walk.values())))
            ck(same,
               f"{name}: one partition, {len(builder)} groups from the builder "
               f"vs {len(walk)} from the walk")

            # the support term: every rule B clade root's subtree clears the floor
            n_tx = cl.transmission_total(ind)
            floor = cl.transmission_floor(ind, 0.01)
            croot_b = cl.designate_clades(lin2mut, par, order, k, stx, floor)
            roots_b = {r for r in croot_b.values() if par.get(r) is not None}
            bad = [r for r in roots_b if stx.get(r, 0) < floor]
            ck(not bad,
               f"{name}: all {len(roots_b)} rule B clade roots clear m_T=1% of "
               f"{n_tx:,} transmissions (floor {floor:.1f})")

            # order_by: topology must put every parent before its children
            _l, par_t, order_t, _s, _x = cl.lineage_tree(
                ph, ind, order_by="topology")
            pos = {l: i for i, l in enumerate(order_t)}
            bad = [l for l, p in par_t.items()
                   if p is not None and p in pos and pos[p] > pos[l]]
            ck(not bad,
               f'{name}: order_by="topology" places every parent before its '
               f"children ({len(order_t):,} lineages)")
            ck(sorted(order_t) == sorted(order),
               f"{name}: both orders cover the same lineage set")
            checked += 1
    ck(checked > 0, f"{checked} run(s) compared")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", nargs="+",
                    default=["control", "HIV_high", "edge_case"])
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("-k", type=int, default=7)
    ap.add_argument("--stamp", default="latest")
    ap.add_argument("--exp-name", default=EXP_NAME)
    a = ap.parse_args()
    from _scenarios import resolve_stamp
    a.stamp = resolve_stamp(a.stamp, a.exp_name)
    sys.setrecursionlimit(300000)

    section_a()
    section_b(a.scenarios, a.seeds, a.k, a.stamp, a.exp_name)
    print(f"\n{PASS[0]} passed, {FAIL[0]} failed")
    sys.exit(1 if FAIL[0] else 0)
