#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What is each figure-2 panel's denominator, and how many clades fall outside it?

Nearly every bug in this repo has been a measurement measuring something other
than what its name claims, and a percentage is the easiest place for that to
hide: the label says "of clades" and the denominator is something narrower.
This reports, per scenario, every clade the panels drop and every class they
keep, from one clade build.

  rule B      how the build's groups become clades: the founding clade is kept
              and labelled, a group with no defining mutations and a parent is
              not a clade, and the rest must carry >= m_T of the run's own
              transmissions (simplicity/clustering.py keep_clades_by_transmission)
  panel 1     composition classes 0..k, plus `other` for a clade whose defining
              count is not k, plus `wt` for the founding clade. These are meant
              to sum to 100% of a run, so `other` is reported in full: it is
              the class most easily mistaken for a rounding error
  panel 3     the rise denominator is kept clades WITH an emergence episode;
              the latency denominator excludes truncated clades as well
  panel 2     three different denominators in one panel: every cell is a share
              of that run's own events, while the per-host intensities divide
              by the host counts -- and a run with no long shedders contributes
              to the standard intensity and not the long one, so the two are
              not averages over the same runs

Reads a build directory, not the server. Read-only. Run from the repo root:

    python -I tests/panel_denominators.py                 # the k=7 operating point
    python -I tests/panel_denominators.py --key k5_w1_p0.1_lo0.1_hi0.9_post
"""
import os
import sys
import json
import argparse

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in ("", "tests", os.path.join("tests", "studio"),
           os.path.join("scripts", "nature_plots")):
    sys.path.insert(0, os.path.join(ROOT, _p))

import simplicity.clustering as cl
import compute
import panels as P

STORE = os.path.join(ROOT, "Data", "fig2_studio")


def pick_build(cache, want_key):
    """The build directory to audit: the current code hash, preferring want_key."""
    code = compute.code_hash()
    hits = []
    for d in sorted(os.listdir(cache)):
        if d.startswith("_"):
            continue
        f = os.path.join(cache, d, "index.json")
        if not os.path.isfile(f):
            continue
        try:
            idx = json.load(open(f))
        except Exception:
            continue
        if idx.get("code") != code:
            continue
        hits.append((d, idx))
    if want_key:
        for d, idx in hits:
            if d.startswith(want_key):
                return d, idx
    for d, idx in hits:
        if idx.get("params", {}).get("k") == 7:
            return d, idx
    return hits[0] if hits else (None, None)


def audit(cache, key, idx, mt):
    k = int(idx["params"]["k"])
    print(f"\nbuild {key}")
    print(f"  params {idx['params']}   code {idx['code']}   "
          f"m_T {mt} ({mt * 100:g}%)\n")

    hdr = (f"{'scenario':10} {'runs':>5} {'groups':>8} {'kept':>7} "
           f"{'founding':>9} {'no def':>8} {'under m_T':>10} {'no phi':>7}")
    print("  rule B: which groups become clades")
    print("  " + hdr)
    rows = {}
    for sc in P.SCEN_ORDER:
        f = os.path.join(cache, key, f"clades_{sc}.json")
        if not os.path.isfile(f):
            continue
        runs = json.load(open(f))["runs"]
        rows[sc] = runs
        g = kept = found = nodef = under = nophi = 0
        for r in runs:
            g += len(r["clades"])
            ks = cl.keep_clades_by_transmission(r["clades"], r["n_tx_run"], mt)
            kept += len(ks)
            floor = mt * r["n_tx_run"]
            for c in r["clades"]:
                if c["founding"]:
                    found += 1
                elif int(c["ndef"]) <= 0:
                    nodef += 1
                elif c["ntx"] <= 0 or c["ntx"] < floor:
                    under += 1
            if r.get("phi") is None:
                nophi += 1
        print(f"  {sc:10} {len(runs):5} {g:8,} {kept:7,} {found:9,} "
              f"{nodef:8,} {under:10,} {nophi:7}")
        if g:
            print(f"  {'':10} {'':5} {'':8} {100*kept/g:6.1f}% "
                  f"{100*found/g:8.1f}% {100*nodef/g:7.1f}% {100*under/g:9.1f}%")

    print("\n  panel 1: composition classes, as % of KEPT clades and of burden")
    classes = P.classes_for(k)
    print(f"  {'scenario':10} " + " ".join(f"{c:>7}" for c in classes))
    for sc, runs in rows.items():
        n = {c: 0 for c in classes}
        b = {c: 0.0 for c in classes}
        tot = totb = 0
        for r in runs:
            for c in cl.keep_clades_by_transmission(r["clades"], r["n_tx_run"], mt):
                cls = P.class_of(c, k)
                n[cls] += 1
                b[cls] += float(c["burden"])
                tot += 1
                totb += float(c["burden"])
        print(f"  {sc:10} " + " ".join(
            f"{100*n[c]/tot:6.1f}%" if tot else "      -" for c in classes))
        print(f"  {'  (burden)':10} " + " ".join(
            f"{100*b[c]/totb:6.1f}%" if totb else "      -" for c in classes))

    print("\n  panel 3: the rise and latency denominators")
    print(f"  {'scenario':10} {'kept':>7} {'with episode':>13} {'truncated':>10} "
          f"{'latency':>8}   {'rise den.':>10} {'lat den.':>9}")
    for sc, runs in rows.items():
        kept = withep = trunc = lat = 0
        for r in runs:
            ks = {c["c"]: c for c in
                  cl.keep_clades_by_transmission(r["clades"], r["n_tx_run"], mt)}
            emg = {e["c"] for e in r["episodes"] if e["emg"] and e["c"] in ks}
            kept += len(ks)
            withep += len(emg)
            trunc += sum(1 for c in ks.values() if c["trunc"])
            lat += sum(1 for ci in emg if not ks[ci]["trunc"])
        print(f"  {sc:10} {kept:7,} {withep:13,} {trunc:10,} {lat:8,}   "
              f"{100*withep/kept if kept else 0:9.1f}% "
              f"{100*lat/kept if kept else 0:8.1f}%")
    print("\n  rise den. / lat den. are the share of KEPT clades each panel-3 "
          "statistic\n  is actually computed over.")


def audit_panel2(store, burnin=True):
    """Panel 2's denominators, read from the transmission table it reads.

    Mirrors `panels.panel2`'s own exclusions rather than restating them loosely:
    a run with no burn-in day is skipped when the window is on, a run with no
    events in the window is skipped, and the per-host intensities are averaged
    only over runs that HAVE that host type.
    """
    import base64
    d = os.path.join(store, "tx")
    f = os.path.join(d, "index.json")
    if not os.path.isfile(f):
        print("\n  panel 2: no transmission table at " + d)
        return
    idx = json.load(open(f))
    print(f"\n  panel 2: the transmission table (window "
          f"{'on' if burnin else 'off'})")
    print(f"  {'scenario':10} {'runs':>5} {'no phi':>7} {'no events':>10} "
          f"{'cell den.':>10} {'long int.':>10} {'std int.':>9}   "
          f"{'mean n_long':>11} {'mean n_std':>10}")
    for sc in idx["scenarios"]:
        runs = json.load(open(os.path.join(d, f"tx_{sc}.json")))["runs"]
        nophi = noev = used = nlong = nstd = 0
        sl, ss = [], []
        for r in runs:
            if burnin and r["phi"] is None:
                nophi += 1
                continue
            d0 = int(np.floor(float(r["phi"]))) if burnin else 0
            tot = {c: float(np.frombuffer(base64.b64decode(r[c]),
                                          np.uint16).astype(float)[d0:].sum())
                   for c in P.CELLS}
            if sum(tot.values()) <= 0:
                noev += 1
                continue
            used += 1
            if r["n_long"]:
                nlong += 1
                sl.append(r["n_long"])
            if r["n_standard"]:
                nstd += 1
                ss.append(r["n_standard"])
        print(f"  {sc:10} {len(runs):5} {nophi:7} {noev:10} {used:10} "
              f"{nlong:10} {nstd:9}   "
              f"{(sum(sl)/len(sl) if sl else 0):11.1f} "
              f"{(sum(ss)/len(ss) if ss else 0):10.1f}")
    print("  cell den. is the runs every cell percentage is averaged over; the\n"
          "  intensity columns are averaged over FEWER runs where a scenario has\n"
          "  runs with none of that host type. control has no long shedders at\n"
          "  all, so its long intensity has no denominator by design.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=STORE)
    ap.add_argument("--key", default="")
    ap.add_argument("--mt", type=float, default=0.01)
    a = ap.parse_args()
    cache = os.path.join(a.store, "cache")
    key, idx = pick_build(cache, a.key)
    if key is None:
        print(f"no build in {cache} matches compute.code_hash() "
              f"= {compute.code_hash()}")
        sys.exit(2)
    audit(cache, key, idx, a.mt)
    audit_panel2(a.store, burnin=True)
