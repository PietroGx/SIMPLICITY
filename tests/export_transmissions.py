"""Per-day transmission-event counts by host type, for panel 2.

A transmission event is one infection: a transmitter, a recipient, a day and
the lineage passed on. Counted into a 2x2 table on the host type of the
transmitter and of the recipient, per day, per run.

The table is assembled TWICE from independent records -- once from each
individual's `parent` (who infected me) and once from each individual's
`new_infections` list (who I infected) -- and the two must agree exactly. They
are the same events seen from opposite ends, so a divergence is a finding about
the run data, not a rounding detail.

Reads only individuals_data; no phylogeny, no clustering, no lineage
frequencies. ~0.85 s per run.

No windowed twins: with per-day counts the burn-in cut is a sum range, which
makes it exact and lets a panel show any window.
"""
import argparse, base64, json, os, sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))
sys.path.insert(0, os.path.join(ROOT, 'tests', 'studio'))
import codehash as ch

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import fig2_preprocess_data as pre
from _scenarios import EXP_NAME, resolve_scenarios, one_sod, scenario_meta

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
CELLS = ["ll", "ls", "sl", "ss"]          # transmitter -> recipient, l=long s=standard


def code_hash():
    """Stamped into index.json so a table built by older code is visible as
    such rather than silently trusted. AST, not bytes (codehash.py)."""
    return ch.files_hash([os.path.abspath(__file__)])


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def _is_long(h):
    return str(h).strip().lower().startswith("long")


def _one(job):
    scenario, ssod = job
    seed = os.path.basename(ssod)
    try:
        ind = om.read_individuals_data(ssod)
        phi = pre.phi_day(ssod)
        end = float(om.read_final_time(ssod))
    except Exception as exc:
        return scenario, seed, None, f"read: {exc}"

    nd = int(np.floor(end)) + 1
    lg = np.array([_is_long(h) for h in ind["type"]], bool)
    par = ind["parent"].astype(str).to_numpy()
    tinf = ind["t_infection"].astype(float).to_numpy()
    # the index is the individual id that `parent` and `individual_infected` refer to
    pos = {str(i): j for j, i in enumerate(ind.index)}

    A = {c: np.zeros(nd, np.int64) for c in CELLS}
    B = {c: np.zeros(nd, np.int64) for c in CELLS}

    def cell(tx_long, rx_long):
        return ("l" if tx_long else "s") + ("l" if rx_long else "s")

    # route A: the recipient's own record of who infected it
    seeded = 0
    for j in range(len(ind)):
        p = par[j]
        if p == "root" or p not in pos:
            seeded += 1
            continue
        d = int(np.floor(tinf[j]))
        if 0 <= d < nd:
            A[cell(lg[pos[p]], lg[j])][d] += 1

    # route B: the transmitter's own record of whom it infected
    for j, evs in enumerate(ind["new_infections"]):
        if not isinstance(evs, (list, tuple)):
            continue
        for e in evs:
            r = pos.get(str(e.get("individual_infected")))
            if r is None:
                continue
            d = int(np.floor(float(e.get("time_infection"))))
            if 0 <= d < nd:
                B[cell(lg[j], lg[r])][d] += 1

    bad = [c for c in CELLS if not np.array_equal(A[c], B[c])]
    if bad:
        gap = {c: int(np.abs(A[c] - B[c]).sum()) for c in bad}
        return scenario, seed, None, f"routes disagree on {bad}: events differing {gap}"

    tot = int(sum(A[c].sum() for c in CELLS))
    if tot + seeded != len(ind):
        return scenario, seed, None, f"events {tot} + seeded {seeded} != infected {len(ind)}"
    mx = max(int(A[c].max()) for c in CELLS)
    if mx > 65534:
        return scenario, seed, None, f"daily count {mx} overflows uint16"

    rec = dict(seed=seed, phi=None if phi is None else round(phi, 2),
               end=round(end, 2), n_days=nd,
               n_infected=int(len(ind)), n_long=int(lg.sum()),
               n_standard=int((~lg).sum()), n_seeded=seeded, n_events=tot,
               **{c: _b64(A[c].astype(np.uint16)) for c in CELLS})
    return scenario, seed, rec, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=EXP_NAME)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()

    present, _ = resolve_scenarios(a.exp_name, a.exp_num, required=SCEN_ORDER)
    jobs = []
    for sc in [s for s in SCEN_ORDER if s in present]:
        for s in dm.get_seeded_simulation_output_dirs(one_sod(a.exp_name, sc, a.exp_num))[:a.seeds]:
            jobs.append((sc, s))
    print(f"[tx] {len(jobs)} runs", flush=True)

    os.makedirs(a.out, exist_ok=True)
    per_scen, failed = {}, []
    with Pool(a.procs) as pool:
        for i, (sc, sd, rec, err) in enumerate(pool.imap_unordered(_one, jobs, chunksize=1), 1):
            if err:
                failed.append((sc, sd, err))
            else:
                per_scen.setdefault(sc, []).append(rec)
            if i % 25 == 0 or i == len(jobs):
                print(f"[tx] {i}/{len(jobs)}", flush=True)

    try:
        meta = scenario_meta(a.exp_name)
    except Exception as exc:
        print(f"[tx] scenario_meta unavailable: {exc}")
        meta = {}
    for sc in per_scen:
        per_scen[sc].sort(key=lambda r: r["seed"])
        with open(os.path.join(a.out, f"tx_{sc}.json"), "w") as fh:
            json.dump({"scenario": sc, "cells": CELLS, "runs": per_scen[sc]}, fh)

    idx = {"scenarios": [s for s in SCEN_ORDER if s in per_scen], "cells": CELLS,
           "code": code_hash(),
           "counts": {s: len(v) for s, v in per_scen.items()}, "meta": meta,
           "failed": [{"scenario": s, "seed": d, "why": w} for s, d, w in failed]}
    with open(os.path.join(a.out, "index.json"), "w") as fh:
        json.dump(idx, fh, indent=1)

    print(f"\n[tx] wrote {sum(len(v) for v in per_scen.values())} runs to {a.out}")
    for s in idx["scenarios"]:
        rs = per_scen[s]
        ev = sum(r["n_events"] for r in rs)
        print(f"  {s:12s} {len(rs):3d} runs  {ev:9,d} events  "
              f"{sum(r['n_seeded'] for r in rs):5,d} seeded")
    if failed:
        print(f"\n[tx] {len(failed)} FAILED")
        for s, d, w in failed[:10]:
            print(f"  {s} {d}: {w}")
    else:
        print("\n[tx] both routes agree on every run, every cell, every day")


if __name__ == "__main__":
    main()
