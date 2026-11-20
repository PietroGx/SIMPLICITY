#!/usr/bin/env python3
"""Cold and warm cost per key for the d1/d2/d3 memos, and the bytes on disk.

Measures only keys whose memo file is absent, so "cold" is cold.

    python tests/studio/test_dist_cold_cost.py --exp-num 2 --n 5
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import compute
import dist
import panels
import sweep


def warmed(store, tag, parts):
    return os.path.exists(panels._panel_cache_path(store, tag, parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=2)
    ap.add_argument("--exp-name", default=compute.EXP_NAME)
    ap.add_argument("--base-store",
                    default=os.path.join(os.path.dirname(os.path.dirname(HERE)),
                                         "Data", "fig2_studio"))
    ap.add_argument("--n", type=int, default=5, help="cold keys per kind")
    a = ap.parse_args()
    store = sweep.store_dir(a.base_store, a.exp_num, a.exp_name)
    ps = dist.builds_on_disk(store)
    print(f"store  {store}\nbuilds {len(ps)}")

    for kind, tag in (("d1", "d1"), ("d3", "d3")):
        cold, warm, size = [], [], []
        for p in ps:
            if len(cold) >= a.n:
                break
            # dist.py's own hash is part of the key (dist.py:51-57)
            parts = ([p, True, dist.code_hash()] if kind == "d1"
                     else [p, 0.0, dist.code_hash()])
            if warmed(store, tag, parts):
                continue
            t0 = time.time()
            if kind == "d1":
                dist.dist1_cached(store, p, True)
            else:
                dist.dist3_cached(store, p, 0.0)
            cold.append(time.time() - t0)
            t0 = time.time()
            if kind == "d1":
                d = dist.dist1_cached(store, p, True)
            else:
                d = dist.dist3_cached(store, p, 0.0)
            warm.append(time.time() - t0)
            assert d["cached"], "second call was not a cache hit"
            size.append(os.path.getsize(
                panels._panel_cache_path(store, tag, parts)))
        if not cold:
            print(f"{kind}: every key already warm, nothing to measure")
            continue
        print(f"{kind}: cold {sum(cold)/len(cold):.3f} s/key "
              f"(min {min(cold):.3f} max {max(cold):.3f}, n={len(cold)})   "
              f"warm {sum(warm)/len(warm)*1000:.1f} ms/key   "
              f"{sum(size)/len(size)/1e6:.3f} MB/key")

    parts = [True, dist.code_hash()]
    if warmed(store, "d2", parts):
        print("d2: already warm")
    else:
        t0 = time.time()
        panels._memo(store, "d2", parts, lambda: dist.dist2(store, True))
        dt = time.time() - t0
        print(f"d2: cold {dt:.3f} s/key   "
              f"{os.path.getsize(panels._panel_cache_path(store,'d2',parts))/1e6:.3f} MB/key")

    pdir = os.path.join(store, "cache", "_panels")
    per = {}
    for f in os.listdir(pdir):
        per.setdefault(f.split("_")[0], [0, 0])
        per[f.split("_")[0]][0] += 1
        per[f.split("_")[0]][1] += os.path.getsize(os.path.join(pdir, f))
    print("\n_panels on disk, by tag:")
    for t, (n, b) in sorted(per.items()):
        print(f"  {t:8s} {n:6d} files  {b/1e9:7.3f} GB")
    tot = sum(b for _, b in per.values())
    print(f"  {'total':8s} {sum(n for n,_ in per.values()):6d} files  {tot/1e9:7.3f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
