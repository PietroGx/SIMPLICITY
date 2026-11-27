#!/usr/bin/env python3
"""Does dist.py measure what panels.py measures?

For every metric key, every scenario: mean(per-run array, skipping None) against
the panel's own `mean`, and the bootstrap n against the panel's `n`. The n check
is not redundant -- a run dropped instead of held as None moves no mean but
breaks the alignment with `seeds`.

    python figures_gui/test_dist_mirrors_panels.py --stamp 2
"""
import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import compute
import dist
import panels
import sweep

TOL = 1e-9
FAIL = []
WORST = [0.0, 0.0, ""]       # abs, rel, where
NCMP = [0]


def mean_of(v):
    """Plain-Python mean over the values _boot would keep, so the comparison is
    not numpy's own summation compared with itself."""
    w = [float(x) for x in v if x is not None and math.isfinite(float(x))]
    return (sum(w) / len(w)) if w else None, len(w)


def cmp(where, arr, agg):
    NCMP[0] += 1
    got, n = mean_of(arr)
    want = (agg or {}).get("mean")
    wn = (agg or {}).get("n")
    if n != wn:
        FAIL.append(f"{where}: n {n} != panel n {wn}")
    if got is None or want is None:
        if (got is None) != (want is None):
            FAIL.append(f"{where}: mean {got} vs panel {want}")
        return
    a = abs(got - want)
    r = a / max(abs(want), 1e-300)
    if a > WORST[0]:
        WORST[:] = [a, r, where]
    if a > TOL and r > 1e-12:
        FAIL.append(f"{where}: mean {got!r} vs panel {want!r} (abs {a:.3g})")


def check_ci(d, label):
    """The `ci` block: the derived keys only, and the same estimator. A chart
    must not form its own interval off the raw runs, so what it is handed has to
    be the bootstrap of the mean panels._boot gives, not a spread."""
    want = [m for m in dist.DERIVED if m in d["series"]]
    if sorted(d["ci"]) != sorted(want):
        FAIL.append(f"{label} ci: keys {sorted(d['ci'])} != {sorted(want)}")
    for m in d["ci"]:
        for sc, agg in d["ci"][m].items():
            if set(agg) != {"n", "mean", "lo", "hi"}:
                FAIL.append(f"{label} ci {m}/{sc}: keys {sorted(agg)}")
            cmp(f"{label} ci {m}/{sc}", d["series"][m][sc], agg)
            if agg["n"] >= 3 and not (agg["lo"] <= agg["mean"] <= agg["hi"]):
                FAIL.append(f"{label} ci {m}/{sc}: mean outside the band")
            if agg["n"] < 3 and agg["lo"] is not None:
                FAIL.append(f"{label} ci {m}/{sc}: band under 3 runs")


def check_alignment(d, label):
    for key, per in d["series"].items():
        if key in d["not_seed_aligned"]:
            continue
        for sc, arr in per.items():
            if len(arr) != len(d["seeds"][sc]):
                FAIL.append(f"{label} {key}/{sc}: {len(arr)} values for "
                            f"{len(d['seeds'][sc])} seeds")


def check_p1(store, p, burnin):
    lab = f"p1 {compute.params_key(p)} burnin={int(burnin)}"
    d = dist.dist1(store, p, burnin)
    a = panels.panel1(store, p, p["mt"], burnin)
    check_alignment(d, lab)
    check_ci(d, lab)
    for sc in panels.SCEN_ORDER:
        s = a["scenarios"][sc]
        for m in ("total_burden", "total_clades"):
            cmp(f"{lab} {m}/{sc}", d["series"][m][sc], s[m])
        for c in a["classes"]:
            for m in ("frac_burden", "frac_clades", "per_clade"):
                cmp(f"{lab} {m}.{c}/{sc}", d["series"][f"{m}.{c}"][sc],
                    s["classes"][c][m])
        # the summed keys have no panel aggregate: the sum of the class means is
        # what they must land on, formed per run instead of over the CIs
        for key, grp in (("frac_burden_long", d["long_classes"]),
                         ("frac_burden_majority_long", d["majority_classes"])):
            vals = [s["classes"].get(c, {}).get("frac_burden", {}).get("mean")
                    for c in grp]
            vals = [v for v in vals if v is not None]
            want = sum(vals) if vals else None
            cmp(f"{lab} {key}/{sc}", d["series"][key][sc],
                dict(mean=want, n=mean_of(d["series"][key][sc])[1]))
    return d


def check_p2(store, burnin):
    lab = f"p2 burnin={int(burnin)}"
    d = dist.dist2(store, burnin)
    a = panels.panel2(store, burnin)
    check_alignment(d, lab)
    check_ci(d, lab)
    for sc in a["scenarios"]:
        s = a["scenarios"][sc]
        cmp(f"{lab} who_transmits/{sc}", d["series"]["who_transmits"][sc],
            s["who_transmits"])
        for c in panels.CELLS:
            cmp(f"{lab} where.{c}/{sc}", d["series"][f"where.{c}"][sc],
                s["where"][c])
        for h in ("long", "standard"):
            cmp(f"{lab} intensity.{h}/{sc}", d["series"][f"intensity.{h}"][sc],
                s["intensity"][h])
        # intensity_ratio has no panel counterpart by design: panel 2 would have
        # to divide one CI by another. Report it, and check it is a per-run ratio
        rr = d["series"]["intensity_ratio"][sc]
        il = d["series"]["intensity.long"][sc]
        isd = d["series"]["intensity.standard"][sc]
        for i, v in enumerate(rr):
            ok = (il[i] is not None and isd[i] not in (None, 0.0))
            if ok and abs(v - il[i] / isd[i]) > 1e-12:
                FAIL.append(f"{lab} intensity_ratio/{sc}[{i}] not the run ratio")
            if not ok and v is not None:
                FAIL.append(f"{lab} intensity_ratio/{sc}[{i}] has no denominator")
        m, n = mean_of(rr)
        pl = s["intensity"]["long"]["mean"]
        ps = s["intensity"]["standard"]["mean"]
        ofm = None if not (pl and ps) else pl / ps
        print(f"    intensity_ratio {sc:9s} n={n:3d} per-run mean="
              f"{'-' if m is None else round(m, 4)}"
              f"   panel long/standard means {'-' if pl is None else round(pl,5)}"
              f" / {'-' if ps is None else round(ps,5)}"
              f" = {'-' if ofm is None else round(ofm, 4)}")
    return d


def check_p3(store, p, mf):
    lab = f"p3 {compute.params_key(p)} mf={mf}"
    d = dist.dist3(store, p, mf)
    a = panels.panel3(store, p, p["mt"], mf)
    check_alignment(d, lab)
    check_ci(d, lab)
    for sc in panels.SCEN_ORDER:
        s = a["scenarios"][sc]
        for m in ("rise_days", "latent_days"):
            cmp(f"{lab} {m}/{sc}", d["series"][m][sc], s[m])
        for pre, field in (("by_class_rise", "by_class_rise"),
                           ("by_class_latent", "by_class_latent")):
            for j, agg in s[field].items():
                cmp(f"{lab} {pre}.{j}/{sc}",
                    d["series"][f"{pre}.{j}"][sc], agg)
            extra = {k.split(".", 1)[1] for k in d["series"]
                     if k.startswith(pre + ".")
                     and any(v is not None for v in d["series"][k][sc])}
            if extra - set(s[field]):
                FAIL.append(f"{lab} {pre}/{sc}: classes {extra - set(s[field])} "
                            f"not in the panel")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stamp", default="2")
    ap.add_argument("--exp-name", default=compute.EXP_NAME)
    ap.add_argument("--base-store",
                    default=os.path.join(os.path.dirname(HERE),
                                         "Data", "figures"))
    a = ap.parse_args()
    store = sweep.store_dir(a.base_store, a.stamp, a.exp_name)
    print(f"store  {store}")
    print(f"hashes compute {compute.code_hash()}  panels {panels.code_hash()}")

    on_disk = {compute.params_key(p): p for p in dist.builds_on_disk(store)}
    print(f"builds {len(on_disk)} readable")
    op = dict(compute.DEFAULTS, window="post")
    cases = [op]
    for want in ("k5_m0_w1_p0.2", "k10_m0.1_w21_p0.2", "k5_m0.001_w21_p0.1"):
        hit = next((p for kk, p in on_disk.items() if kk.startswith(want + "_")),
                   None)
        if hit:
            cases.append(hit)

    t0 = time.time()
    for p in cases:
        for burnin in (True, False):
            print(f"  panel1  {compute.params_key(p)} burnin={int(burnin)}")
            check_p1(store, p, burnin)
        for mf in (0.0, 0.1):
            print(f"  panel3  {compute.params_key(p)} mf={mf}")
            check_p3(store, p, mf)
    for burnin in (True, False):
        print(f"  panel2  burnin={int(burnin)}")
        check_p2(store, burnin)

    print(f"\n{NCMP[0]} mean/n comparisons")
    print(f"worst |diff| {WORST[0]:.3g} (rel {WORST[1]:.3g}) at {WORST[2]}")
    print(f"{len(FAIL)} failures, {time.time()-t0:.1f}s")
    for f in FAIL[:40]:
        print("  FAIL " + f)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
