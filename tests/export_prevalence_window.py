#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Long-shedder prevalence per run, whole-run and post-burn-in.

The figure-2 builder export carries `n_infected` and `n_long_hosts` over the
WHOLE run only, so panel 1 compares a clade quantity against a denominator
that still counts the burn-in. The clades and their defining mutations are not
time series and stay as they are; the infection counts are what the burn-in
cut applies to.

Reads only the individuals table and phi, never the phylogeny or the lineage
frequencies, so this costs ~0.8 s per seed instead of the clustering export's
hour. `tests/export_fig2_builder.py` emits the same two windowed counts
natively now; this script exists so the current export can be used without
waiting for a re-export, and both paths are asserted equal on the whole-run
columns.
"""
import os
import sys
import json
import argparse
from multiprocessing import Pool

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import fig2_preprocess_data as pre
from _scenarios import EXP_NAME, resolve_scenarios, one_sod

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]


def _is_long(s):
    return s.astype(str).str.strip().str.lower().str.startswith("long")


def _one(job):
    scenario, ssod = job
    seed = os.path.basename(ssod)
    try:
        ind = om.read_individuals_data(ssod)
        phi = pre.phi_day(ssod)
        end = float(om.read_final_time(ssod))
    except Exception as exc:
        return {"scenario": scenario, "seed": seed, "error": f"read: {exc}"}

    lg = _is_long(ind["type"])
    ti = ind["t_infection"].astype(float)
    out = {"scenario": scenario, "seed": seed,
           "phi": None if phi is None else round(float(phi), 2),
           "end": round(end, 2),
           "n_infected": int(len(ind)),
           "n_long_hosts": int(lg.sum())}

    if phi is None:
        out.update(n_infected_w=None, n_long_hosts_w=None)
        return out
    post = ti >= float(phi)
    out["n_infected_w"] = int(post.sum())
    out["n_long_hosts_w"] = int((post & lg).sum())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-name", default=EXP_NAME)
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", default=os.path.join(ROOT, "tests", "artifacts", "fig2b", "prevalence.json"))
    ap.add_argument("--csv", default=os.path.join(ROOT, "tests", "_prevalence_window.csv"))
    a = ap.parse_args()

    present, missing = resolve_scenarios(a.exp_name, a.exp_num,
                                         required=SCEN_ORDER, verbose=True)
    jobs = []
    for sc in SCEN_ORDER:
        if sc not in present:
            continue
        sod = one_sod(a.exp_name, sc, a.exp_num)
        for s in dm.get_seeded_simulation_output_dirs(sod)[:a.seeds]:
            jobs.append((sc, s))
    print(f"[prev] {len(jobs)} seeds", flush=True)

    rows = []
    with Pool(a.procs) as p:
        for i, r in enumerate(p.imap_unordered(_one, jobs, chunksize=2), 1):
            rows.append(r)
            if i % 25 == 0:
                print(f"[prev]   {i}/{len(jobs)}", flush=True)

    bad = [r for r in rows if r.get("error")]
    for r in bad:
        print(f"[prev] !! {r['scenario']}/{r['seed']}: {r['error']}")
    rows = [r for r in rows if not r.get("error")]
    rows.sort(key=lambda r: (SCEN_ORDER.index(r["scenario"]), r["seed"]))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump({"exp_name": a.exp_name, "exp_num": a.exp_num,
                   "runs": rows}, f, separators=(",", ":"))
    pd.DataFrame(rows).to_csv(a.csv, index=False)
    print(f"[prev] {len(rows)} runs -> {a.out}")
    print(f"[prev] {len(rows)} runs -> {a.csv}")

    print("\nscenario     runs   long-shedder prevalence")
    print("                     whole run   post burn-in")
    df = pd.DataFrame(rows)
    for sc in SCEN_ORDER:
        d = df[df["scenario"] == sc]
        if not len(d):
            continue
        w = d.dropna(subset=["n_infected_w"])

        def pct(num, den, fr):
            if not len(fr):
                return "   --   "
            t = fr[den].sum()
            return f"{100.0 * fr[num].sum() / t:7.3f}%" if t else "   --   "
        print(f"{sc:<12}{len(d):>5}   {pct('n_long_hosts','n_infected',d)}   {pct('n_long_hosts_w','n_infected_w',w)}")


if __name__ == "__main__":
    main()
