#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Precompute everything the studio memoises, so every slider answers from disk.

The clade grid is the expensive precomputation but not the whole of it. Each
endpoint memoises on its own key, and a key the page can ask for and that is
not on disk is a wait in the middle of dragging a slider.

Measured cold, run #2, 168 builds:

    panel 1     0.12 s   key [build, m_T, burn-in]        336 keys
    panel 1c  not timed  key [build, m_T, burn-in]        336 keys  (Figure 4 points)
    panel 2        ~0 s  key [burn-in]                      2 keys
    panel 3     0.18 s   key [build, m_T, m_F]           1,344 keys
    episodes    0.04 s   key [build, m_T, m_F]           1,344 keys
    across      1.0 s    key [build, present]              168 keys
    daily       0.8-3.3s key [k, scenario, seed, m_T]    5,250 keys
    panel 4     59 s     key [build, m_T, burn-in, m_F]    112 passes
                         -- one row per build, so 112 x 168 = 18,816 rows
    dist 1      0.04 s   key [build, burn-in, dist hash]    336 keys
    dist 2      0.02 s   key [burn-in, dist hash]             2 keys
    dist 3      0.02 s   key [build, m_F, dist hash]      1,344 keys

The d* kinds are the per-run arrays behind the same CIs (dist.py), keyed as
their panels are and cheaper, since they bootstrap only the derived keys.
/api/sens assembles from d1 and d3 and is not memoised itself, so warming those
warms it. They add about 150 MB to _panels.

m_T is both a build parameter and an argument, and the two always move
together: the build is NAMED by m_T. So panels 1 and 3 are warmed at each
build's own m_T, not at the cross product. Panel 4 is the exception -- it walks
every build at the slider's m_T -- which is why it is the long pole.

Requests go through panels.panel() and the other server entry points, so this
cannot warm a key the page does not ask for. A key is skipped only if its file
exists, and the file name hashes the code that computes it (panels.scoped_hash),
so an edit to panel1 re-warms p1, p4 rows and d1 without being told to.

    python tests/studio/warm_panels.py --exp-num 2 --procs 4
    python tests/studio/warm_panels.py --exp-num 2 --what p1,p3,daily
    python tests/studio/warm_panels.py --exp-num 2 --dry-run
    python tests/studio/warm_panels.py --exp-num 2 --shard 3/16   # 4th of 16
    python tests/studio/warm_panels.py --exp-num 2 --shard 0/32 --list  # jobs per shard

--shard I/N takes every N-th job from I, so each shard gets a mix of kinds.
tests/studio/slurm_warm.sh runs the N shards as an array; a key two shards
both reach (panel 4 fills p1/p3 rows) is just computed twice.
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compute
import dist
import panels
import sweep

# The read parameters the rail offers. The build axes are not listed here: they
# are read off the cache, so a build that is missing is not warmed and a build
# that is there is, whatever the grid intended.
MF_VALUES = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]
BURNIN = (True, False)

# Panel 4 ignores these for row selection -- each row uses its own build's
# params -- but panel() validates them, so send a point that exists.
BASE = {"k": "7", "win": "15", "prom": "0.1", "lo": "0.1", "hi": "0.9",
        "window": "post"}

KINDS = ("p2", "p1", "p1c", "p3", "d2", "d1", "d3", "episodes", "across",
         "daily", "p4")


def builds_on_disk(store):
    """Each build's own params, for the builds this code can read. A build left
    by older code is skipped: _load_builds resolves by params_key, which
    carries the CURRENT hash, so it would look up a directory that is not
    there."""
    cache = os.path.join(store, "cache")
    out = []
    for d in sorted(os.listdir(cache)) if os.path.isdir(cache) else []:
        f = os.path.join(cache, d, "index.json")
        if d.startswith("_") or not os.path.exists(f):
            continue
        try:
            ix = json.load(open(f))
        except Exception:
            continue
        if ix.get("code") == compute.code_hash():
            out.append(ix["params"])
    return out


def _q(p, **extra):
    q = {k: str(v) for k, v in p.items()}
    q.update({k: str(v) for k, v in extra.items()})
    return q


def plan(store, kinds):
    """Every job, cheapest kind first, so the studio becomes usable early."""
    ps = builds_on_disk(store)
    mts = sorted({float(p["mt"]) for p in ps})
    ks = sorted({int(p["k"]) for p in ps})
    jobs = []
    if "p2" in kinds:
        jobs += [("p2", dict(n="2", burnin=int(b), **BASE)) for b in BURNIN]
    if "p1" in kinds:
        jobs += [("p1", _q(p, n=1, burnin=int(b))) for p in ps for b in BURNIN]
    if "p1c" in kinds:
        jobs += [("p1c", _q(p, n="1c", burnin=int(b))) for p in ps for b in BURNIN]
    if "p3" in kinds:
        jobs += [("p3", _q(p, n=3, mf=mf)) for p in ps for mf in MF_VALUES]
    # the per-run arrays behind the same CIs, keyed as their panels are. /api/sens
    # assembles from d1 and d3, so warming those warms it too.
    if "d2" in kinds:
        jobs += [("d2", dict(n="2", burnin=int(b), **BASE)) for b in BURNIN]
    if "d1" in kinds:
        jobs += [("d1", _q(p, n=1, burnin=int(b))) for p in ps for b in BURNIN]
    if "d3" in kinds:
        jobs += [("d3", _q(p, n=3, mf=mf)) for p in ps for mf in MF_VALUES]
    if "episodes" in kinds:
        jobs += [("episodes", _q(p, mf=mf)) for p in ps for mf in MF_VALUES]
    if "across" in kinds:
        jobs += [("across", _q(p)) for p in ps]
    if "daily" in kinds:
        jobs += [("daily", dict(BASE, k=str(k), scenario=sc, seed=str(i),
                                mt=str(mt)))
                 for k in ks for mt in mts
                 for sc in panels.SCEN_ORDER for i in range(50)]
    if "p4" in kinds:
        jobs += [("p4", dict(BASE, n="4", mt=str(mt), burnin=int(b), mf=str(mf)))
                 for mt in mts for b in BURNIN for mf in MF_VALUES]
    return ps, jobs


STORE = None


def _one(job):
    kind, q = job
    t0 = time.time()
    try:
        if kind in ("p1", "p1c", "p2", "p3", "p4"):
            panels.panel(STORE, q)
        elif kind in ("d1", "d2", "d3"):
            dist.dist(STORE, q)
        elif kind == "episodes":
            panels.episode_view(STORE, q)
        elif kind == "across":
            panels.across(STORE, q)
        elif kind == "daily":
            panels.run_daily(STORE, q)
        return kind, q, time.time() - t0, None
    except Exception as exc:
        return kind, q, time.time() - t0, f"{type(exc).__name__}: {exc}"


def _init(store, exp_num, exp_name):
    global STORE
    STORE = store
    # panels reads the RAW run output for the daily record, and defaults to
    # experiment 1; server.main() rebinds these per process and so must this.
    panels.EXP_NUM = exp_num
    panels.EXP_NAME = exp_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=compute.EXP_NAME)
    ap.add_argument("--store", default=None,
                    help="override; by default the experiment's own store")
    ap.add_argument("--base-store", default=os.path.join(ROOT, "Data", "figures_studio"))
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--what", default=",".join(KINDS),
                    help="comma-separated: " + ", ".join(KINDS))
    ap.add_argument("--shard", default="0/1",
                    help="I/N: warm only every N-th job from I")
    ap.add_argument("--dry-run", action="store_true",
                    help="say what would be warmed, warm nothing")
    ap.add_argument("--list", action="store_true",
                    help="jobs per kind in each of the N shards, warm nothing")
    a = ap.parse_args()

    store = (os.path.expanduser(a.store) if a.store else
             sweep.store_dir(a.base_store, a.exp_num, a.exp_name))
    if not os.path.isdir(os.path.join(store, "cache")):
        print(f"[warm] no cache at {store}/cache", file=sys.stderr)
        return 1
    kinds = [k.strip() for k in a.what.split(",") if k.strip()]
    bad = [k for k in kinds if k not in KINDS]
    if bad:
        print(f"[warm] unknown kind(s): {', '.join(bad)}", file=sys.stderr)
        return 1

    try:
        si, sn = (int(x) for x in a.shard.split("/"))
        assert 0 <= si < sn
    except (ValueError, AssertionError):
        print(f"[warm] --shard wants I/N with 0 <= I < N, got {a.shard}",
              file=sys.stderr)
        return 1

    ps, jobs = plan(store, kinds)
    if a.list:
        ks = [k for k in KINDS if any(j[0] == k for j in jobs)]
        print(f"[warm] {len(ps)} builds, {len(jobs)} jobs over {sn} shard(s)")
        print("shard " + " ".join(f"{k:>8s}" for k in ks) + "    total")
        for i in range(sn):
            sh = jobs[i::sn]
            print(f"{i:5d} " + " ".join(f"{sum(j[0] == k for j in sh):8d}" for k in ks)
                  + f" {len(sh):8d}")
        print("total " + " ".join(f"{sum(j[0] == k for j in jobs):8d}" for k in ks)
              + f" {len(jobs):8d}")
        return 0
    jobs = jobs[si::sn]
    print(f"[warm] store  {store}")
    print(f"[warm] builds {len(ps)} readable at code {compute.code_hash()}")
    for k in KINDS:
        n = sum(1 for j in jobs if j[0] == k)
        if n:
            print(f"[warm]   {k:9s} {n:6d} jobs")
    print(f"[warm] {len(jobs)} jobs, {a.procs} processes, shard {si}/{sn}")
    if a.dry_run:
        return 0

    t0 = time.time()
    nbad = 0
    seen = {}
    with Pool(a.procs, initializer=_init,
              initargs=(store, a.exp_num, a.exp_name)) as pool:
        for i, (kind, q, dt, err) in enumerate(
                pool.imap_unordered(_one, jobs, chunksize=4), 1):
            seen[kind] = seen.get(kind, 0) + dt
            if err:
                nbad += 1
                if nbad <= 20:
                    print(f"[warm] {i}/{len(jobs)} {kind} {q}  FAILED  {err}",
                          flush=True)
            elif i % 50 == 0 or i == len(jobs):
                el = time.time() - t0
                eta = el * (len(jobs) - i) / max(i, 1)
                print(f"[warm] {i}/{len(jobs)} {kind}  {el/60:5.1f} min in, "
                      f"~{eta/60:.0f} min left, {nbad} failed", flush=True)
    pdir = os.path.join(store, "cache", "_panels")
    size = sum(os.path.getsize(os.path.join(pdir, f))
               for f in os.listdir(pdir)) if os.path.isdir(pdir) else 0
    print(f"[warm] core-seconds by kind: "
          + ", ".join(f"{k}={v:.0f}" for k, v in sorted(seen.items())))
    print(f"[warm] done in {(time.time()-t0)/60:.1f} min, "
          f"_panels is {size/1e6:.1f} MB, {nbad} failed")
    return 1 if nbad else 0


if __name__ == "__main__":
    sys.exit(main())
