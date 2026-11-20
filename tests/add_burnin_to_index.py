#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Add each seed's burn-in day to the clade explorer's index.json.

The cutoff is read from fig2_preprocess_data.phi_day, not re-derived here, so
the explorer's line and the figure's window cannot disagree about what burn-in
means. Only index.json is rewritten; the 250 per-seed files are untouched.
"""
import os
import sys
import json
import argparse
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import fig2_preprocess_data as pre
from _scenarios import one_sod


def _one(job):
    exp_name, exp_num, scenario, seed = job
    try:
        sods = dm.get_seeded_simulation_output_dirs(one_sod(exp_name, scenario, exp_num))
        ssod = next(s for s in sods if os.path.basename(s) == seed)
    except Exception as exc:
        return scenario, seed, None, None, f"no dir: {exc}"
    try:
        phi = pre.phi_day(ssod)
    except Exception as exc:
        return scenario, seed, None, None, f"phi: {exc}"
    try:
        end = float(om.read_final_time(ssod))
    except Exception as exc:
        return scenario, seed, phi, None, f"end: {exc}"
    return scenario, seed, phi, end, None


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()

    with open(a.index) as f:
        idx = json.load(f)
    rows = idx["seeds"]
    jobs = [(idx["exp_name"], idx["exp_num"], r["scenario"], r["seed"]) for r in rows]
    print(f"[burnin] {len(jobs)} seeds, window={pre.BURNIN_WINDOW}, "
          f"MIN_FINAL_TIME={pre.MIN_FINAL_TIME}, MIN_WINDOW_DAYS={pre.MIN_WINDOW_DAYS}")

    got = {}
    with Pool(a.procs) as pool:
        for i, (sc, sd, phi, end, err) in enumerate(pool.imap_unordered(_one, jobs, chunksize=4), 1):
            got[(sc, sd)] = (phi, end, err)
            if i % 50 == 0 or i == len(jobs):
                print(f"[burnin]   {i}/{len(jobs)}", flush=True)

    nphi = nin = 0
    for r in rows:
        phi, end, err = got.get((r["scenario"], r["seed"]), (None, None, "missing"))
        r["phi"] = None if phi is None else round(float(phi), 2)
        r["end"] = None if end is None else round(float(end), 2)
        inwin = (phi is not None and end is not None
                 and end >= pre.MIN_FINAL_TIME and (end - phi) >= pre.MIN_WINDOW_DAYS)
        r["fig2"] = bool(inwin)
        nphi += phi is not None
        nin += inwin
        if err:
            print(f"[burnin][warn] {r['scenario']}/{r['seed']}: {err}")

    idx["burnin"] = {"window": pre.BURNIN_WINDOW,
                     "fixed_days": pre.BURNIN_CUTOFF_DAYS,
                     "min_final_time": pre.MIN_FINAL_TIME,
                     "min_window_days": pre.MIN_WINDOW_DAYS}
    with open(a.index, "w") as f:
        json.dump(idx, f, separators=(",", ":"))
    print(f"[burnin] {nphi}/{len(rows)} reach phi=1; {nin} are inside fig2's window")
    print(f"[burnin] wrote {a.index} ({os.path.getsize(a.index)/1e3:.1f} kB)")
