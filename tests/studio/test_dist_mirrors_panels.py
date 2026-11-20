#!/usr/bin/env python3
"""Does dist.py measure what panels.py measures?

For every metric key, every scenario: mean(per-run array, skipping None) against
the panel's own `mean`, and the bootstrap n against the panel's `n`. The n check
is not redundant -- a run dropped instead of held as None moves no mean but
breaks the alignment with `seeds`.

    python tests/studio/test_dist_mirrors_panels.py --exp-num 2
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
from dist import _same

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


def check_sens(store, q, axis, expect):
    s = dist.sens(store, dict(q, axis=axis))
    lab = f"sens axis={axis}"
    lv = s["levels"]
    if lv != sorted(lv, key=float):
        FAIL.append(f"{lab}: levels not sorted numerically: {lv}")
    if len(lv) != len(set(map(float, lv))):
        FAIL.append(f"{lab}: {len(lv)} levels for {len(set(map(float,lv)))} values")
    if sorted(map(float, lv)) != sorted(map(float, expect)):
        FAIL.append(f"{lab}: levels {lv} != grid {expect}")
    for o in dist.BUILD_AXES:
        if o != axis and o not in s["pinned"]:
            FAIL.append(f"{lab}: {o} is not pinned")
    op = dist.operating_point()
    if not _same(s["operating"], op[axis]):
        FAIL.append(f"{lab}: operating {s['operating']} != {op[axis]}")
    if s["operating_key"] != compute.params_key(op):
        FAIL.append(f"{lab}: operating_key {s['operating_key']}")
    if s["operating_key"] not in s["keys"]:
        FAIL.append(f"{lab}: operating_key is not one of the levels")
    # every level's CI must come from that level's own per-run array
    for i, b in enumerate([dict(s["pinned"], **{axis: v}) for v in lv]):
        b = {k: b[k] for k in dist.BUILD_AXES}
        d1 = dist.dist1_cached(store, b, s["burnin"])
        d3 = dist.dist3_cached(store, b, s["mf"])
        for sc in panels.SCEN_ORDER:
            cmp(f"{lab}[{lv[i]}] frac_burden_majority_long/{sc}",
                d1["series"]["frac_burden_majority_long"][sc],
                s["scenarios"][sc]["frac_burden_majority_long"][i])
            for m in ("rise_days", "latent_days"):
                cmp(f"{lab}[{lv[i]}] {m}/{sc}", d3["series"][m][sc],
                    s["scenarios"][sc][m][i])
            for m in dist.SENS_METRICS:
                row = s["scenarios"][sc][m][i]
                if set(row) != {"n", "mean", "lo", "hi"}:
                    FAIL.append(f"{lab} {m}/{sc}[{i}]: keys {sorted(row)}")
    return s


def check_sens_vs_panel4(store, q):
    """At the operating point, panel 4's row and /api/sens must agree on the
    headline -- panel 4 sums the class means, /api/sens means the per-run sums."""
    p4 = panels.panel4(store, panels.clean_params(q), float(q["mt"]),
                       str(q["burnin"]) not in ("0", "false", ""),
                       panels.read_mf(q))
    # NOT p4["operating_point"]: _p4_operating_key (panels.py:457-463) matches on
    # P4_AXES only, so with m_T out of that list it names k7_m0.001, not the
    # m_T=0.01 operating point. Pick the row whose params ARE the build.
    want = panels.clean_params(q)
    row = next(r for r in p4["rows"]
               if all(_same(r["params"].get(k2), want[k2])
                      for k2 in dist.BUILD_AXES))
    s = dist.sens(store, dict(q, axis="k"))
    i = [float(v) for v in s["levels"]].index(float(row["params"]["k"]))
    print(f"    panel4 row {row['key']}")
    print(f"    panel4 operating_point says {p4['operating_point']}")
    for sc in panels.SCEN_ORDER:
        for m in dist.SENS_METRICS:
            want = row["scenarios"][sc][m]
            got = s["scenarios"][sc][m][i]["mean"]
            lo, hi = s["scenarios"][sc][m][i]["lo"], s["scenarios"][sc][m][i]["hi"]
            if (want is None) != (got is None):
                FAIL.append(f"sens vs p4 {m}/{sc}: {got} vs {want}")
            elif want is not None and abs(got - want) > TOL \
                    and abs(got - want) / max(abs(want), 1e-300) > 1e-12:
                FAIL.append(f"sens vs p4 {m}/{sc}: {got!r} vs {want!r}")
            print(f"      {sc:9s} {m:26s} p4 {want if want is None else round(want,6):>10} "
                  f"  sens {got if got is None else round(got,6):>10} "
                  f"[{'' if lo is None else round(lo,4)}, "
                  f"{'' if hi is None else round(hi,4)}]")
    return p4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=2)
    ap.add_argument("--exp-name", default=compute.EXP_NAME)
    ap.add_argument("--base-store",
                    default=os.path.join(os.path.dirname(os.path.dirname(HERE)),
                                         "Data", "figures_studio"))
    a = ap.parse_args()
    store = sweep.store_dir(a.base_store, a.exp_num, a.exp_name)
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

    q = {k: str(v) for k, v in op.items()}
    q.update(burnin="1", mf="0")
    grid = dict(k=sweep.K, mt=[float(v) for v in sweep.MT],
                win=sweep.WIN, prom=sweep.PROM)
    for axis in dist.SENS_AXES:
        print(f"  sens    axis={axis}")
        s = check_sens(store, q, axis, grid[axis])
        print(f"    levels {s['levels']}  pinned {s['pinned']}")
    print(f"  operating point {compute.params_key(dist.operating_point())}")
    print("  sens vs panel4")
    check_sens_vs_panel4(store, q)

    print(f"\n{NCMP[0]} mean/n comparisons")
    print(f"worst |diff| {WORST[0]:.3g} (rel {WORST[1]:.3g}) at {WORST[2]}")
    print(f"{len(FAIL)} failures, {time.time()-t0:.1f}s")
    for f in FAIL[:40]:
        print("  FAIL " + f)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
