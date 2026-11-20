"""Per-run values behind the panel CIs, and the sensitivity table with bands.

Every number panels.py sends is already a mean with a bootstrap CI over runs
(panels.py:146-156). The per-run values exist on the way there -- panel 1 builds
per[cls]["frac_burden"] as a list over runs at panels.py:207 and then collapses it
-- but nothing in the payload carries them.

They cannot be added by editing panels.py: the panel memo key is
[tag, parts, compute.code_hash(), panels.code_hash()] (panels.py:43-55) and
panels.code_hash() is the sha256 of panels.py's own bytes (:27-32). This file
imports panels and REUSES its functions, memoising under its own tags through
panels._memo, so panels.code_hash() never moves.

Each array here is the very list the panel builds before _boot collapses it, so
mean(array, skipping None) is the panel's own mean. Runs the panel skips hold
None at their index instead of being dropped, so every point traces back to its
seed, by_class_rise/by_class_latent included: panel 3 forms those per run
too, so `not_seed_aligned` is empty.

The derived keys -- the sums and ratios formed inside a run -- carry a `ci`
block built with panels._boot, so a chart draws the same estimator panels.py
draws everywhere else instead of forming a spread over 50 runs under the same
whisker.

panels._memo's key does not include THIS file's bytes, so code_hash() goes into
every `parts` list: an edit here discards the d* memos and nothing else.

The contract is written out in web/dist.js.
"""
import base64
import hashlib
import json
import os

import numpy as np

import compute
import codehash as ch
import panels
import simplicity.clustering as cl

# Everything params_key names, i.e. what makes a build a different build
# (compute.py:74-86).
BUILD_AXES = ["k", "mt", "win", "prom", "lo", "hi", "window"]
SENS_AXES = ["k", "mt", "win", "prom"]
SENS_METRICS = ["frac_burden_majority_long", "rise_days", "latent_days"]
# formed inside a run, so no panel field holds their interval
DERIVED = ["frac_burden_long", "frac_burden_majority_long", "intensity_ratio"]


def code_hash():
    """This file's behaviour, for the memo key.

    panels.scoped_hash names what a tag reuses from panels.py, and cannot hash
    a function defined here, so this goes into `parts` instead. Over the AST
    (codehash.py), so a comment or docstring here costs nothing either."""
    return ch.files_hash([os.path.abspath(__file__)])


def _ci(series, keys):
    """The bootstrap CI of the mean for the derived keys, panels._boot over the
    very per-run list, so it is the estimator every other CI in the studio is."""
    return {m: {sc: panels._boot(v) for sc, v in series[m].items()}
            for m in keys if m in series}


def _out(sc_keys):
    return {sc: [] for sc in sc_keys}


# ------------------------------------------------------------------ panel 1

def dist1(store, p, burnin):
    """Mirror of panels.panel1 (panels.py:183-234), uncollapsed."""
    idx, data = panels._load_builds(store, p)
    bkey = "burden" if burnin else "burden_all"
    k = int(p["k"])
    classes = panels.classes_for(k)
    jlong = [str(j) for j in range(1, k + 1)]
    maj = [str(j) for j in range(k // 2 + 1, k + 1)]      # panels.py:387
    seeds, series = {}, {}
    for m in ("frac_burden", "frac_clades", "per_clade"):
        for c in classes:
            series[f"{m}.{c}"] = _out(panels.SCEN_ORDER)
    for m in ("total_burden", "total_clades",
              "frac_burden_long", "frac_burden_majority_long"):
        series[m] = _out(panels.SCEN_ORDER)

    for sc in panels.SCEN_ORDER:
        runs = data.get(sc) or []
        seeds[sc] = [r["seed"] for r in runs]
        for r in runs:
            ks = cl.keep_designated_clades(r["clades"])
            if not ks:
                for key in series:
                    series[key][sc].append(None)
                continue
            tb = sum(c[bkey] for c in ks)
            tn = len(ks)
            series["total_burden"][sc].append(tb)
            series["total_clades"][sc].append(tn)
            byc = {}
            for c in ks:
                byc.setdefault(panels.class_of(c, k), []).append(c)
            bsum = {}
            for cls in classes:
                sel = byc.get(cls, [])
                b = sum(c[bkey] for c in sel)
                bsum[cls] = b
                series[f"frac_burden.{cls}"][sc].append(
                    b / tb if tb > 0 else None)
                series[f"frac_clades.{cls}"][sc].append(
                    len(sel) / tn if tn else None)
                series[f"per_clade.{cls}"][sc].append(
                    b / len(sel) if sel else None)
            # summed INSIDE the run: a class CI cannot be added to another's
            for key, grp in (("frac_burden_long", jlong),
                             ("frac_burden_majority_long", maj)):
                v = None if (tb <= 0 or not grp) else \
                    sum(bsum[c] for c in grp) / tb
                series[key][sc].append(v)
    return dict(kind="dist", n=1, params=p, mt=p["mt"], burnin=burnin,
                classes=classes, long_classes=jlong, majority_classes=maj,
                seeds=seeds, series=series, ci=_ci(series, DERIVED),
                not_seed_aligned=[])


# ------------------------------------------------------------------ panel 2

def dist2(store, burnin):
    """Mirror of panels.panel2 (panels.py:236-293), uncollapsed."""
    d = os.path.join(store, "tx")
    idx_f = os.path.join(d, "index.json")
    if not os.path.exists(idx_f):
        raise FileNotFoundError("transmissions")
    idx = json.load(open(idx_f))
    scens = idx["scenarios"]
    keys = (["who_transmits", "intensity_ratio"]
            + [f"where.{c}" for c in panels.CELLS]
            + ["intensity.long", "intensity.standard"])
    series = {key: _out(scens) for key in keys}
    seeds = {}
    for sc in scens:
        runs = json.load(open(os.path.join(d, f"tx_{sc}.json")))["runs"]
        seeds[sc] = [r["seed"] for r in runs]
        for r in runs:
            ser = {c: np.frombuffer(base64.b64decode(r[c]), np.uint16).astype(float)
                   for c in panels.CELLS}
            d0 = 0
            skip = False
            if burnin:
                if r["phi"] is None:
                    skip = True
                else:
                    d0 = int(np.floor(float(r["phi"])))
            tot = {c: float(ser[c][d0:].sum()) for c in panels.CELLS}
            all_ev = sum(tot.values())
            if skip or all_ev <= 0:
                for key in keys:
                    series[key][sc].append(None)
                continue
            lrow = tot["ll"] + tot["ls"]
            srow = tot["sl"] + tot["ss"]
            series["who_transmits"][sc].append(lrow / all_ev)
            for c in panels.CELLS:
                series[f"where.{c}"][sc].append(tot[c] / all_ev)
            il = lrow / r["n_long"] if r["n_long"] else None
            isd = srow / r["n_standard"] if r["n_standard"] else None
            series["intensity.long"][sc].append(il)
            series["intensity.standard"][sc].append(isd)
            # per-run ratio, not one mean over another; no ratio without both
            series["intensity_ratio"][sc].append(
                il / isd if il is not None and isd not in (None, 0.0) else None)
    return dict(kind="dist", n=2, params=None, burnin=burnin,
                cells=panels.CELLS, seeds=seeds, series=series,
                ci=_ci(series, DERIVED), not_seed_aligned=[])


# ------------------------------------------------------------------ panel 3

def dist3(store, p, mf=0.0):
    """Mirror of panels.panel3, uncollapsed.

    by_class_rise/by_class_latent are a class's mean over ONE run's clades, as
    panels.panel3 forms them, so they are seed-aligned like every other key: a
    run with no clade in that class holds None.
    """
    idx, data = panels._load_builds(store, p)
    k = int(p["k"])
    seeds, series = {}, {"rise_days": _out(panels.SCEN_ORDER),
                         "latent_days": _out(panels.SCEN_ORDER)}
    by_cls = {}
    for sc in panels.SCEN_ORDER:
        runs = data.get(sc) or []
        seeds[sc] = [r["seed"] for r in runs]
        per_run = []
        for r in runs:
            ks = {c["c"]: c for c in cl.keep_designated_clades(r["clades"])}
            emg = panels.emergence_of(r, ks, mf)
            rr, ll, rr_cls, ll_cls = [], [], {}, {}
            for ci, e in emg.items():
                c = ks[ci]
                rr.append(e["days"])
                j = panels.class_of(c, k)
                rr_cls.setdefault(j, []).append(e["days"])
                if not c["trunc"]:
                    v = e["a_day"] - c["tem"]
                    if v >= 0:
                        ll.append(v)
                        ll_cls.setdefault(j, []).append(v)
            series["rise_days"][sc].append(float(np.mean(rr)) if rr else None)
            series["latent_days"][sc].append(float(np.mean(ll)) if ll else None)
            per_run.append((rr_cls, ll_cls))
        for pre, i in (("by_class_rise", 0), ("by_class_latent", 1)):
            for j in sorted({j for pr in per_run for j in pr[i]}, key=str):
                by_cls.setdefault(f"{pre}.{j}", _out(panels.SCEN_ORDER))[sc] = [
                    float(np.mean(pr[i][j])) if j in pr[i] else None
                    for pr in per_run]
    for per in by_cls.values():
        for sc in panels.SCEN_ORDER:
            if not per[sc]:
                per[sc] = [None] * len(seeds[sc])
    series.update(by_cls)
    return dict(kind="dist", n=3, params=p, mt=p["mt"], mf=mf,
                seeds=seeds, series=series, ci=_ci(series, DERIVED),
                not_seed_aligned=[])


# ------------------------------------------------------------------ dispatch

def dist1_cached(store, p, burnin):
    return panels._memo(store, "d1", [p, burnin, code_hash()],
                        lambda: dist1(store, p, burnin))


def dist3_cached(store, p, mf=0.0):
    return panels._memo(store, "d3", [p, mf, code_hash()],
                        lambda: dist3(store, p, mf))


def dist(store, q):
    which = str(q.get("n", "1"))
    burnin = str(q.get("burnin", "1")) not in ("0", "false", "")
    if which == "2":
        return panels._memo(store, "d2", [burnin, code_hash()],
                            lambda: dist2(store, burnin))
    p = panels.clean_params(q)
    if which == "1":
        return dist1_cached(store, p, burnin)
    if which == "3":
        return dist3_cached(store, p, panels.read_mf(q))
    raise ValueError(f"no distribution for panel {which}")


# ------------------------------------------------------------------ /api/sens

def builds_on_disk(store):
    """Each build's own params, for the builds this code can read -- the filter
    warm_panels.builds_on_disk applies (warm_panels.py:61-78)."""
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


def operating_point():
    """The operating point across EVERY build parameter: compute.DEFAULTS plus
    window="post", the values a rebuild uses when nothing is overridden.

    NOT panels._p4_operating_key (panels.py:458-465), which matches on P4_AXES
    = ["k","win","prom"] only. m_T is absent from that list, so it answers with
    the first row in directory order, k7_m0.001_... -- m_T = 0.1% instead of the
    operating point's 1% -- and that key feeds _p4_summary's at_operating_point
    (panels.py:432-437).
    """
    return dict(compute.DEFAULTS, window="post")


def _same(a, b):
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return a == b


def sens(store, q):
    """Panel 4's sensitivity line, one axis at a time, with bands.

    panel4's own by_axis cannot serve it. P4_AXES is ["k","win","prom"]
    (panels.py:442), so _p4_by_axis (:503-523) leaves m_T unpinned and
    by_axis["k"] comes back with 21 levels for three values of k -- seven m_T
    rows per k, sorted lexicographically by directory name so m0 lands after
    m0.1. Here every build axis but the one asked for is pinned, m_T included,
    and levels sort numerically.

    Each metric carries a full bootstrap CI over runs. Panel 4's rows take
    .get("mean") and drop lo/hi (panels.py:388-393), and a sensitivity line
    without a band invites a 2-day wiggle to be read as a real shift.
    """
    axis = q.get("axis")
    if axis not in SENS_AXES:
        raise ValueError(f"axis must be one of {', '.join(SENS_AXES)}")
    p = panels.clean_params(q)
    burnin = str(q.get("burnin", "1")) not in ("0", "false", "")
    mf = panels.read_mf(q)
    others = [a for a in BUILD_AXES if a != axis]
    sel = [b for b in builds_on_disk(store)
           if all(_same(b.get(o), p.get(o)) for o in others)]
    sel.sort(key=lambda b: float(b[axis]))

    scens = {sc: {m: [] for m in SENS_METRICS} for sc in panels.SCEN_ORDER}
    cached = True
    for b in sel:
        # each level is read at its OWN m_T: the build is named by it
        d1 = dist1_cached(store, b, burnin)
        d3 = dist3_cached(store, b, mf)
        cached = cached and d1.get("cached") and d3.get("cached")
        for sc in panels.SCEN_ORDER:
            for m, src in ((SENS_METRICS[0], d1), (SENS_METRICS[1], d3),
                           (SENS_METRICS[2], d3)):
                scens[sc][m].append(
                    panels._boot((src["series"].get(m) or {}).get(sc) or []))
    # the level to rule. operating_key is null unless this slice passes through
    # the operating build, so a ruled level is never mistaken for it
    op = operating_point()
    okey = compute.params_key(op)
    on_slice = any(compute.params_key(b) == okey for b in sel)
    return dict(kind="sens", axis=axis, pinned={o: p[o] for o in others},
                levels=[b[axis] for b in sel],
                keys=[compute.params_key(b) for b in sel],
                operating=op[axis], operating_key=okey if on_slice else None,
                burnin=burnin, mf=mf, cached=bool(cached),
                metrics=panels.P4_METRICS, scenarios=scens)
