"""The panel cache: the aggregations each panel plots, over the clade builds.

Two rules hold everywhere in this file.

1. The run is the independent unit. Every interval is a bootstrap over RUNS,
   never over clades. A percentage is formed inside a run and then averaged
   across runs, so a run with 3000 clades does not outvote one with 300.

2. Scenarios are never pooled. The five scenarios are a designed pair of
   contrasts -- SOT against HIV_low holds prevalence and moves duration,
   HIV_low against HIV_high holds duration and moves prevalence -- so a
   statistic that averages across them destroys the experiment.
"""
import hashlib, json, os, socket

import numpy as np

import codehash as ch
import compute
import simplicity.clustering as cl

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
CELLS = ["ll", "ls", "sl", "ss"]
NBOOT = 2000
SEED = 20261116


def code_hash():
    """This whole file's behaviour. Kept for /api/index to report and as a
    tripwire recorded INSIDE each memo, not for keying one -- see scoped_hash."""
    return ch.files_hash([os.path.abspath(__file__)])


# What each memoised answer actually depends on. `@tag` means "everything that
# tag depends on", so a row built from two panels inherits both.
#
# Explicit, never derived from the call graph. A list that silently misses a
# dependency serves a stale answer under a correct-looking name, which is the
# failure this repo has paid for four times; a list that names a function that
# does not exist raises instead (codehash.names_hash), so the cheap mistake is
# the loud one. A tag with no entry has no key and raises too.
MEMO_DEPS = {
    "p1": ("panel1", "_load_builds", "_build_dir", "classes_for", "class_of",
           "_boot", "_prev", "_dur"),
    "p2": ("panel2", "_boot"),
    "p3": ("panel3", "_load_builds", "_build_dir", "classes_for", "class_of",
           "emergence_of", "_boot", "_prev", "_dur"),
    "p4row": ("_panel4_row", "panel1_cached", "panel3_cached", "@p1", "@p3"),
    "episodes": ("episode_view", "_load_builds", "_build_dir", "emergence_of",
                 "class_of", "_boot"),
    "across": ("across", "across_the_set", "values_on_disk", "_load_builds",
               "_build_dir", "_doubling_bins", "_median_per_value", "_median"),
    "daily": ("run_daily", "_daily_series", "_load_builds", "_build_dir"),
    # dist.py memoises through _memo as well. These name only what it reuses
    # from HERE; its own code goes into `parts` as dist.code_hash(), because
    # this file cannot hash a function it does not define.
    "d1": ("@p1",),
    "d2": ("@p2",),
    "d3": ("@p3",),
}

# Deliberately outside every key, and checked as such by
# tests/studio/test_memo_deps_complete.py rather than tolerated silently:
#
#   the caching machinery   decides WHERE an answer is stored, not what it is.
#                           Keying on it would mean every change to how keys
#                           are made discards every memo -- the very thing
#                           scoped_hash exists to stop.
#   the request parsers     clean_params and read_mf produce the `parts` that
#                           are ALREADY in the key, so a change to them shows
#                           up as a different key by itself.
MEMO_INFRA = ("_memo", "_panel_cache_path", "scoped_hash", "_deps",
              "code_hash", "clean_params", "read_mf")

# Module constants no function's source carries. NBOOT and SEED decide every
# interval, so they belong in every key.
_CONSTS = ("NBOOT", "SEED", "SCEN_ORDER", "CELLS")

_scoped = {}


def _deps(tag, seen=None):
    """`MEMO_DEPS[tag]` with every `@tag` reference resolved."""
    seen = seen or set()
    if tag in seen:
        return set()
    seen.add(tag)
    out = set()
    for n in MEMO_DEPS[tag]:
        out |= _deps(n[1:], seen) if n.startswith("@") else {n}
    return out


def scoped_hash(tag):
    """The code THIS answer depends on, as an AST digest.

    One whole-file hash meant a change to panel 2 discarded panels 1, 3 and 4,
    `across` and `daily` with it -- about 10,200 memos and an hour of warming
    for an edit that touched none of them. This names only what computes the
    answer, plus compute.py's hash because every one of them reads a clade
    build, plus the constants above.
    """
    if tag not in _scoped:
        consts = repr([(n, globals()[n]) for n in _CONSTS])
        _scoped[tag] = ch.names_hash(globals(), _deps(tag),
                                     extra=compute.code_hash() + consts)
    return _scoped[tag]


def tx_code_hash():
    """What the transmission exporter hashes to RIGHT NOW, so a table stamped by
    older code shows up as such instead of being trusted silently."""
    return ch.files_hash([os.path.join(compute.ROOT, "tests", "export_transmissions.py")])


def _panel_cache_path(store, tag, parts):
    """The panel cache on disk. Keyed by the request and by the code that
    computes THIS answer -- scoped_hash, not the whole file, so editing one
    panel does not discard the others.

    The set of clade builds on disk used to be in this key, because panel 4
    reads every build. It is not any more: panel 4 caches one ROW per build
    instead, so adding a build no longer discards every other panel's answer.
    """
    cache = os.path.join(store, "cache")
    raw = json.dumps([tag, parts, scoped_hash(tag)], sort_keys=True)
    h = hashlib.sha256(raw.encode()).hexdigest()[:16]
    d = os.path.join(cache, "_panels")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{tag}_{h}.json")


def _memo(store, tag, parts, build):
    f = _panel_cache_path(store, tag, parts)
    if os.path.exists(f):
        try:
            with open(f) as fh:
                out = json.load(fh)
            out["cached"] = True
            return out
        except Exception:
            pass
    out = build()
    out["cached"] = False
    # per writer: warm shards on different nodes can build the same key
    tmp = f"{f}.{socket.gethostname()}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w") as fh:
            json.dump(out, fh)
        os.replace(tmp, f)
    except Exception:
        pass
    return out


def clean_params(q):
    d = dict(compute.DEFAULTS)
    d["window"] = "post"
    for name, cast in (("k", int), ("mt", float), ("win", float),
                       ("prom", float), ("lo", float), ("hi", float)):
        if name in q and q[name] not in (None, ""):
            d[name] = cast(q[name])
    if q.get("window") in compute.WINDOWS:
        d["window"] = q["window"]
    if d["k"] not in (0, 3, 5, 7, 10, 15, 20):
        raise ValueError(f"k={d['k']} is not one of the k values")
    for n in ("prom", "lo", "hi"):
        if not 0.0 <= d[n] <= 1.0:
            raise ValueError(f"{n} out of range")
    if not 0.0 <= d["mt"] <= 1.0:
        raise ValueError("m_T is a fraction of the run's own transmissions")
    if not 1.0 <= d["win"] <= 120.0:
        raise ValueError("smoothing window out of range")
    if d["lo"] >= d["hi"]:
        raise ValueError("lo must be below hi")
    return d


def _build_dir(store, p):
    return os.path.join(store, "cache", compute.params_key(p))


def index(store):
    """What exists, and what each scenario's design values are."""
    cache = os.path.join(store, "cache")
    keys = sorted(os.listdir(cache)) if os.path.isdir(cache) else []
    built = []
    for kk in keys:
        f = os.path.join(cache, kk, "index.json")
        if os.path.exists(f):
            try:
                built.append(json.load(open(f)))
            except Exception:
                pass
    tx = os.path.join(store, "tx", "index.json")
    txi = json.load(open(tx)) if os.path.exists(tx) else None
    if txi is not None:
        txi["code_now"] = tx_code_hash()
    return dict(store=store, code=compute.code_hash(), panels_code=code_hash(),
                defaults=dict(compute.DEFAULTS, window="post"),
                builds=built, transmissions=txi,
                scen_order=SCEN_ORDER)


def _load_builds(store, p):
    d = _build_dir(store, p)
    idx_f = os.path.join(d, "index.json")
    if not os.path.exists(idx_f):
        raise FileNotFoundError(compute.params_key(p))
    idx = json.load(open(idx_f))
    out = {}
    for sc in idx["scenarios"]:
        f = os.path.join(d, f"clades_{sc}.json")
        if os.path.exists(f):
            out[sc] = json.load(open(f))["runs"]
    return idx, out


def _boot(vals, lo=2.5, hi=97.5):
    """Percentile bootstrap over runs. Returns None below 3 runs: with one or
    two there is no interval to quote, and quoting one would invent precision."""
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
    if len(v) < 3:
        return dict(n=int(len(v)),
                    mean=float(v.mean()) if len(v) else None, lo=None, hi=None)
    rng = np.random.default_rng(SEED)
    bs = rng.choice(v, size=(NBOOT, len(v)), replace=True).mean(axis=1)
    return dict(n=int(len(v)), mean=float(v.mean()),
                lo=float(np.percentile(bs, lo)), hi=float(np.percentile(bs, hi)))


def classes_for(k):
    """Composition classes. `j` means "j of the clade's k defining mutations
    arose in a long shedder". Two classes are not j at all and are labelled:
    `wt` is the founding clade, which has no defining mutations, and `other`
    holds clades whose defining-mutation count is not k -- the segment can
    overshoot, so a handful carry 8 at k=7. Neither is dropped; together with
    0..k they account for every clade, so the classes sum to 100% of a run."""
    return [str(j) for j in range(k + 1)] + ["other", "wt"]


def class_of(c, k):
    if c["founding"]:
        return "wt"
    if int(c["ndef"]) != int(k):
        return "other"
    return str(int(c["nlong"]))


def panel1(store, p, mt, burnin):
    """Composition class against burden.

    Classes are "j long of 7" for j = 0..7, plus `wt` for the founding clade,
    which has zero defining mutations -- its long-born ratio is 0/0 and
    undefined, so it is labelled, never merged into j = 0.
    """
    idx, data = _load_builds(store, p)
    bkey = "burden" if burnin else "burden_all"
    k = int(p["k"])
    classes = classes_for(k)
    out = {}
    for sc in SCEN_ORDER:
        runs = data.get(sc) or []
        per = {c: dict(frac_burden=[], frac_clades=[], per_clade=[])
               for c in classes}
        tot_b, tot_n = [], []
        for r in runs:
            ks = cl.keep_designated_clades(r["clades"])
            if not ks:
                continue
            tb = sum(c[bkey] for c in ks)
            tn = len(ks)
            tot_b.append(tb)
            tot_n.append(tn)
            byc = {}
            for c in ks:
                byc.setdefault(class_of(c, k), []).append(c)
            for cls in classes:
                sel = byc.get(cls, [])
                b = sum(c[bkey] for c in sel)
                per[cls]["frac_burden"].append(b / tb if tb > 0 else None)
                per[cls]["frac_clades"].append(len(sel) / tn if tn else None)
                per[cls]["per_clade"].append(b / len(sel) if sel else None)
        out[sc] = dict(
            runs=len(tot_b),
            total_burden=_boot(tot_b), total_clades=_boot(tot_n),
            classes={c: {m: _boot(v[m]) for m in
                         ("frac_burden", "frac_clades", "per_clade")}
                     for c, v in per.items()},
            prevalence=_prev(idx, sc), duration=_dur(idx, sc))
    return dict(kind="panel1", params=p, mt=mt, burnin=burnin,
                classes=classes, scenarios=out,
                note="fractions are of the run's clades as defined; the "
                     "founding clade is the `wt` class, so the classes sum "
                     "to 1")


def _prev(idx, sc):
    m = (idx.get("meta") or {}).get(sc) or {}
    v = m.get("ratio")
    # the config's own `long_shedders_ratio`, a fraction. It was scaled to a
    # percentage here; the studio now converts at the display layer only.
    return None if v is None else float(v)


def _dur(idx, sc):
    m = (idx.get("meta") or {}).get(sc) or {}
    return m.get("duration")


# ------------------------------------------------------------------ panel 2

def panel2(store, burnin):
    """Transmission events by host type. All three claims, from one table."""
    d = os.path.join(store, "tx")
    idx_f = os.path.join(d, "index.json")
    if not os.path.exists(idx_f):
        raise FileNotFoundError("transmissions")
    idx = json.load(open(idx_f))
    import base64
    out = {}
    series = {}
    for sc in idx["scenarios"]:
        runs = json.load(open(os.path.join(d, f"tx_{sc}.json")))["runs"]
        who, cells, inten = [], {c: [] for c in CELLS}, dict(long=[], standard=[])
        acc = None
        for r in runs:
            nd = r["n_days"]
            ser = {c: np.frombuffer(base64.b64decode(r[c]), np.uint16).astype(float)
                   for c in CELLS}
            d0 = 0
            if burnin:
                if r["phi"] is None:
                    continue
                d0 = int(np.floor(float(r["phi"])))
            tot = {c: float(ser[c][d0:].sum()) for c in CELLS}
            all_ev = sum(tot.values())
            if all_ev <= 0:
                continue
            lrow = tot["ll"] + tot["ls"]
            who.append(lrow / all_ev)
            for c in CELLS:
                cells[c].append(tot[c] / all_ev)
            if r["n_long"]:
                inten["long"].append(lrow / r["n_long"])
            if r["n_standard"]:
                inten["standard"].append((tot["sl"] + tot["ss"]) / r["n_standard"])
            # mean per-day curve across runs, on a common 0..nd axis
            cur = np.vstack([ser[c] for c in CELLS])
            if acc is None:
                acc = np.zeros((4, nd))
                cnt = np.zeros(nd)
            m = min(nd, acc.shape[1])
            acc[:, :m] += cur[:, :m]
            cnt[:m] += 1
        out[sc] = dict(
            runs=len(who),
            who_transmits=_boot(who),
            where=({c: _boot(cells[c]) for c in CELLS}),
            intensity={h: _boot(v) for h, v in inten.items()})
        if acc is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
            series[sc] = dict(days=int(acc.shape[1]),
                              mean={c: np.round(mean[i], 4).tolist()
                                    for i, c in enumerate(CELLS)})
    return dict(kind="panel2", burnin=burnin, cells=CELLS, scenarios=out,
                series=series,
                note="seeded infections have no transmitter and are excluded; "
                     "they are counted separately per run")


# ------------------------------------------------------------------ panel 3

def emergence_of(r, ks, mf):
    """{clade: its emergence episode}, after dropping episodes whose peak is
    under `mf` of the viral population.

    The stored `emg` flag is kept whenever it survives, so mf=0 is exactly what
    the build recorded. Only when it is dropped is a replacement chosen, by
    largest amplitude -- compute.py:371's own tie-break, since the clade's
    series is not here to test which episode holds its maximum.
    """
    out = {}
    for e in r["episodes"]:
        if e["c"] not in ks or e["top"] < mf:
            continue
        cur = out.get(e["c"])
        if cur is None or (e["emg"] and not cur["emg"]) or \
           (e["emg"] == cur["emg"] and e["amp"] > cur["amp"]):
            out[e["c"]] = e
    return out


def panel3(store, p, mt, mf=0.0):
    """Emergence speed and latency.

    The emergence episode is the one containing the clade's highest frequency.
    A clade already rising at the burn-in day keeps its rise time and loses its
    latent interval, so the two metrics have DIFFERENT denominators and both
    counts are reported.
    """
    idx, data = _load_builds(store, p)
    out = {}
    for sc in SCEN_ORDER:
        runs = data.get(sc) or []
        rise, lat, nrise, nlat, ntrunc, nbefore, ncens = [], [], [], [], [], [], []
        by_cls_rise, by_cls_lat = {}, {}
        for r in runs:
            ks = {c["c"]: c for c in cl.keep_designated_clades(r["clades"])}
            emg = emergence_of(r, ks, mf)
            rr, ll = [], []
            # per class, WITHIN this run: the run is the independent unit, so a
            # class's value is its mean over this run's clades and the interval
            # bootstraps those. Pooling clades across runs instead treated 329
            # clades as 329 independent draws against 50 runs, narrowing the
            # band by about 2.6x and letting a run with many clades outvote one
            # with few.
            rr_cls, ll_cls = {}, {}
            ntr = nbef = 0
            for ci, e in emg.items():
                c = ks[ci]
                rr.append(e["days"])
                j = class_of(c, p["k"])
                rr_cls.setdefault(j, []).append(e["days"])
                if c["trunc"]:
                    # already circulating when the window opened, so its first
                    # appearance is outside it: the rise stands, the latency
                    # cannot be measured
                    ntr += 1
                    continue
                v = e["a_day"] - c["tem"]
                if v >= 0:
                    ll.append(v)
                    ll_cls.setdefault(j, []).append(v)
                else:
                    # the emergence episode appears to start before the clade
                    # does. Counted, because it is the third reason a latency
                    # is missing and nothing reported it.
                    nbef += 1
            if rr:
                rise.append(float(np.mean(rr)))
                nrise.append(len(rr))
                # the three counts over the SAME runs as `runs`, so they
                # reconcile: n_with_rise = n_with_latency + n_truncated +
                # n_latency_before_start. ntrunc used to be appended for every
                # run, including runs contributing to neither estimate, which
                # made HIV_high's counts sum past its own total.
                ntrunc.append(ntr)
                nbefore.append(nbef)
                # still high when the record ends: the rise is observed, the
                # peak is a lower bound (clade_preprocess.episodes_of)
                ncens.append(sum(1 for e in emg.values() if e.get("cens")))
            if ll:
                lat.append(float(np.mean(ll)))
                nlat.append(len(ll))
            for j, v in rr_cls.items():
                by_cls_rise.setdefault(j, []).append(float(np.mean(v)))
            for j, v in ll_cls.items():
                by_cls_lat.setdefault(j, []).append(float(np.mean(v)))
        out[sc] = dict(
            runs=len(rise),
            rise_days=_boot(rise), latent_days=_boot(lat),
            n_with_rise=int(sum(nrise)), n_with_latency=int(sum(nlat)),
            n_truncated=int(sum(ntrunc)),
            n_latency_before_start=int(sum(nbefore)),
            n_censored=int(sum(ncens)),
            by_class_rise={j: _boot(v) for j, v in sorted(by_cls_rise.items())},
            by_class_latent={j: _boot(v) for j, v in sorted(by_cls_lat.items())},
            prevalence=_prev(idx, sc), duration=_dur(idx, sc))
    return dict(kind="panel3", params=p, mt=mt, scenarios=out,
                note="rise and latency have different denominators: a clade "
                     "already rising at the burn-in day keeps its rise time "
                     "and loses its latency. The counts are of CLADES, "
                     "over the runs that produced a rise, and they reconcile: "
                     "n_with_rise = n_with_latency + n_truncated + "
                     "n_latency_before_start. The intervals bootstrap over "
                     "RUNS, so their n is a run count, not a clade count.")


# ------------------------------------------------------------------ panel 4

def _panel4_row(store, q, key, mt, burnin, mf=0.0):
    """One sensitivity row: the headline of panels 1 and 3 at one clade build.

    Cached per build, so adding a build costs one row rather than the whole
    table, and moving m_T costs one row per build with each cached on its own.
    Before this, panel 4 recomputed panels 1 and 3 from scratch for every build
    on disk -- 30.4 s once the sweep put 24 builds there, against 2.8 s with 2.
    """
    p1 = panel1_cached(store, q, mt, burnin)
    p3 = panel3_cached(store, q, mt, mf)
    row = dict(key=key, params=q, scenarios={})
    for sc in SCEN_ORDER:
        a = p1["scenarios"].get(sc) or {}
        b = p3["scenarios"].get(sc) or {}
        cls = a.get("classes") or {}
        # one headline per panel: burden percentage in the majority-long
        # classes, and the mean emergence rise
        k = int(q["k"])
        maj = [str(j) for j in range(k // 2 + 1, k + 1)]   # j classes only
        vals = [cls.get(c, {}).get("frac_burden", {}).get("mean") for c in maj]
        vals = [v for v in vals if v is not None]
        row["scenarios"][sc] = dict(
            frac_burden_majority_long=sum(vals) if vals else None,
            rise_days=(b.get("rise_days") or {}).get("mean"),
            latent_days=(b.get("latent_days") or {}).get("mean"),
            runs=a.get("runs"))
    return row


def panel4(store, p, mt, burnin, mf=0.0):
    """Sensitivity of panels 1 and 3 to the clustering threshold and to the
    episode rules -- whatever clade builds exist.

    The assembly is not cached as a whole: each row is, so a new build adds one
    row's work instead of invalidating the table.
    """
    cache = os.path.join(store, "cache")
    keys = sorted(os.listdir(cache)) if os.path.isdir(cache) else []
    rows = []
    for kk in keys:
        if kk.startswith("_"):
            continue
        f = os.path.join(cache, kk, "index.json")
        if not os.path.exists(f):
            continue
        try:
            idx = json.load(open(f))
        except Exception:
            continue
        # A build stamped by older code cannot be read: _load_builds looks it up
        # by params_key, which carries the CURRENT hash, so it would resolve to a
        # directory that does not exist. Skip it rather than fill the table with
        # rows whose only content is that 404.
        if idx.get("code") != compute.code_hash():
            continue
        q = idx["params"]
        try:
            rows.append(_memo(store, "p4row", [kk, mt, burnin, mf],
                              lambda q=q, kk=kk: _panel4_row(store, q, kk, mt,
                                                             burnin, mf)))
        except Exception as exc:
            rows.append(dict(key=kk, params=q, error=str(exc)))
    ok = [r for r in rows if "error" not in r]
    return dict(kind="panel4", mt=mt, burnin=burnin, rows=rows,
                cached=all(r.get("cached") for r in rows if "error" not in r),
                axes=_p4_axes(ok), operating_point=_p4_operating_key(ok),
                summary=_p4_summary(ok), by_axis=_p4_by_axis(ok),
                metrics=P4_METRICS,
                note="each row is one clade build. Only builds present in the "
                     "cache appear; rebuild at another k or another episode "
                     "setting to add a row.")


# m_T is a BUILD parameter -- the build is named by it -- so it has to be
# pinned like the rest. It was missing, and that single omission produced two
# defects: _p4_operating_key matched k/win/prom only and returned the first row
# in directory order, m_T = 0.1% rather than the operating point's 1% (control
# latency 137.88 d against 183.53 d); and _p4_by_axis left it unpinned, so
# by_axis["k"] came back with 21 levels for three values of k.
P4_AXES = ["k", "mt", "win", "prom"]
P4_METRICS = {
    "frac_burden_majority_long": "burden in the majority-long classes "
                                 "(fraction)",
    "rise_days": "mean emergence rise (days)",
    "latent_days": "mean latency (days)",
}


def _p4_axes(rows):
    """The grid the builds on disk actually span, one sorted list per axis."""
    return {a: sorted({r["params"][a] for r in rows if a in r["params"]})
            for a in P4_AXES}


def _p4_operating_key(rows):
    """The build at the operating point: compute.DEFAULTS, the values a
    rebuild uses when nothing is overridden."""
    for r in rows:
        q = r["params"]
        if all(q.get(a) == compute.DEFAULTS[a] for a in P4_AXES):
            return r["key"]
    return None


def _median(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    m = len(v) // 2
    return v[m] if len(v) % 2 else (v[m - 1] + v[m]) / 2.0


def _p4_summary(rows):
    """Per scenario and metric, what the WHOLE grid does to the number.

    `span` is max - min across every build on disk: the honest width of the
    claim, as opposed to its value at one chosen setting. `at_operating_point`
    is that chosen setting, so the two can be read against each other.
    """
    opkey = _p4_operating_key(rows)
    out = {}
    for sc in SCEN_ORDER:
        per = {}
        for m in P4_METRICS:
            vals = [(r["scenarios"].get(sc) or {}).get(m) for r in rows]
            vals = [v for v in vals if v is not None]
            at = None
            if opkey:
                hit = next((r for r in rows if r["key"] == opkey), None)
                at = (hit["scenarios"].get(sc) or {}).get(m) if hit else None
            per[m] = dict(
                n=len(vals),
                min=min(vals) if vals else None,
                max=max(vals) if vals else None,
                median=_median(vals),
                at_operating_point=at,
                span=(max(vals) - min(vals)) if vals else None)
        out[sc] = per
    return out


def _p4_by_axis(rows):
    """One-at-a-time sensitivity: move one axis, pin the others at the
    operating point.

    The whole-grid span mixes the axes together, so a number can look unstable
    because of a setting nobody would choose. This isolates each axis, which is
    the form the question is usually asked in -- "does k change the answer?"
    """
    out = {}
    for axis in P4_AXES:
        others = [a for a in P4_AXES if a != axis]
        sel = [r for r in rows
               if all(r["params"].get(o) == compute.DEFAULTS[o] for o in others)]
        sel.sort(key=lambda r: r["params"].get(axis))
        out[axis] = dict(
            levels=[r["params"].get(axis) for r in sel],
            keys=[r["key"] for r in sel],
            scenarios={sc: {m: [(r["scenarios"].get(sc) or {}).get(m)
                                for r in sel]
                            for m in P4_METRICS} for sc in SCEN_ORDER},
            pinned={o: compute.DEFAULTS[o] for o in others})
    return out


# ------------------------------------------------------------------ dispatch

def panel1_cached(store, p, mt, burnin):
    """Panel 1 through the cache. Panel 4 builds a row per clade build from
    panels 1 and 3, and used to call them directly -- recomputing both from
    scratch for every build on disk, 24 times over after the parameter sweep.
    Going through the cache means a row shares the very answer the Panel 1 tab
    shows, so whichever is asked for first pays and the other is free."""
    return _memo(store, "p1", [p, mt, burnin],
                 lambda: panel1(store, p, mt, burnin))


def panel3_cached(store, p, mt, mf=0.0):
    return _memo(store, "p3", [p, mt, mf], lambda: panel3(store, p, mt, mf))


def read_mf(q):
    """m_F: the smallest peak frequency an episode must reach to be counted.

    A READ parameter, like the burn-in toggle: episodes are stored, so moving
    it filters an existing build instead of naming a new one. 0 counts every
    episode, which is what the build recorded.
    """
    mf = float(q.get("mf", 0.0))
    if not 0.0 <= mf <= 1.0:
        raise ValueError("m_F is a fraction of the viral population")
    return mf


def panel(store, q):
    which = q.get("n", "1")
    mt = float(q.get("mt", compute.DEFAULTS["mt"]))
    if not 0.0 <= mt <= 1.0:
        raise ValueError("m_T is a fraction of the run's transmissions")
    burnin = str(q.get("burnin", "1")) not in ("0", "false", "")
    mf = read_mf(q)
    if which == "2":
        return _memo(store, "p2", [burnin], lambda: panel2(store, burnin))
    p = clean_params(q)
    if which == "1":
        return panel1_cached(store, p, mt, burnin)
    if which == "3":
        return panel3_cached(store, p, mt, mf)
    if which == "4":
        return panel4(store, p, mt, burnin, mf)
    raise ValueError(f"no panel {which}")


def episode_view(store, q):
    """Page 2: what the current episode rules produce, per scenario, plus a
    sample of clade trajectories' episode structure for eyeballing."""
    p = clean_params(q)
    mt = float(q.get("mt", compute.DEFAULTS["mt"]))
    mf = read_mf(q)
    idx, data = _load_builds(store, p)
    out = {}
    for sc in idx["scenarios"]:
        runs = data.get(sc) or []
        per_clade, rates, amps = [], [], []
        for r in runs:
            ks = {c["c"]: c for c in cl.keep_designated_clades(r["clades"])}
            eps = [e for e in r["episodes"]
                   if e["c"] in ks and e["top"] >= mf]
            if ks:
                rates.append(len(eps) / len(ks))
            per_clade += [e["days"] for e in eps if e["emg"]]
            amps += [e["amp"] for e in eps]
        out[sc] = dict(runs=len(runs),
                       episodes_per_clade=_boot(rates),
                       emergence_rise_days=_boot(per_clade),
                       n_episodes=len(amps))
    return dict(kind="episodes", params=p, mt=mt, scenarios=out,
                counts=idx.get("counts"), failed=idx.get("failed"))


# ------------------------------------------------- the inspection endpoints

RUNS_DIR = "runs"

# The m_T axis of the build grid, from compute so the sweep driver and the
# panels cannot disagree about which values exist. Clade sizes stay in doubling
# bins, as the artifact drew them.
MT_VALUES = compute.MT_VALUES
ACROSS_METRICS = ["clades", "singleton_frac", "median_size", "max_div",
                  "never_tx_frac"]


def run_index(store):
    """Which runs the per-seed export holds, and each one's burn-in day.

    The burn-in comes from this index rather than being re-derived, because the
    exporter took it from `fig2_preprocess_data.phi_day` -- the same function
    the clade pass and the figure use. Re-deriving it here would be a second
    answer to a question that already has one.
    """
    f = os.path.join(store, RUNS_DIR, "index.json")
    if not os.path.exists(f):
        raise FileNotFoundError("the per-seed export (runs/) is not present")
    idx = json.load(open(f))
    by_scen = {}
    for s in idx.get("seeds", []):
        by_scen.setdefault(s["scenario"], []).append(
            dict(seed=s["seed"], phi=s.get("phi"), end=s.get("end")))
    for v in by_scen.values():
        v.sort(key=lambda r: r["seed"])
    return dict(scenarios=by_scen,
                thresholds=idx.get("thresholds", []),
                n_times=idx.get("n_times"),
                counts={k: len(v) for k, v in by_scen.items()})


def run_file(store, q):
    """One run's packed arrays, with its burn-in day attached.

    The arrays are passed through untouched: the browser decodes them and
    rebuilds clade membership from `starts`, which is the clustering's own
    decision rather than an approximation of it.
    """
    scen = q.get("scenario") or SCEN_ORDER[0]
    if scen not in SCEN_ORDER:
        raise ValueError(f"unknown scenario {scen!r}")
    ri = run_index(store)
    runs = ri["scenarios"].get(scen) or []
    if not runs:
        raise FileNotFoundError(f"no exported runs for {scen}")
    i = int(q.get("seed", 0))
    if not 0 <= i < len(runs):
        raise ValueError(f"seed index {i} outside 0..{len(runs) - 1} for {scen}")
    rec = runs[i]
    f = os.path.join(store, RUNS_DIR, f"{scen}__{rec['seed']}.json")
    if not os.path.exists(f):
        raise FileNotFoundError(os.path.basename(f))
    obj = json.load(open(f))
    obj["phi"] = rec["phi"]
    obj["end"] = rec["end"]
    obj["seed_index"] = i
    return obj


def _doubling_bins(sizes):
    """Clade sizes as a histogram over doubling bins: 1, 2-3, 4-7, 8-15, ...
    A linear axis is useless here -- most clades hold one lineage and the
    largest hold thousands."""
    out = {}
    for n in sizes:
        b = 0 if n <= 1 else int(n).bit_length() - 1
        out[b] = out.get(b, 0) + 1
    return [dict(bin=b, lo=1 << b, hi=(1 << (b + 1)) - 1, n=out[b])
            for b in sorted(out)]


def across_the_set(store, p):
    """The clustering's shape along the m_T axis of the build grid.

    This is the across-the-set tab's data. The artifact drew it from a
    precomputed (seed x k x m) grid, `runs/formation.json`; it reads clade
    builds instead so that one definition answers everywhere, and it keeps that
    grid's values so the axes are unchanged.

    Until v2.5.5 it swept them by FILTERING one build, which was right
    while the studio designated on divergence and dropped clades under m_T.
    Under rule B, m_T is part of designation: each value is a different
    partition, not the same one filtered, so each needs its own build. A
    value with no build on disk is reported as missing rather than guessed --
    `tests/studio/slurm_precompute.sh` is what fills the grid.

    m_T = 0 means no transmission minimum, so it is the value where the
    definition is divergence alone. Above it a clade had to have transmitted,
    though a designated clade can still carry no transmissions OF ITS OWN once
    its descendants have split off, which is why `never_tx_pct` is reported at
    every value.
    """
    out = {sc: dict(runs=0, seeds={}, size_bins=None, median=[])
           for sc in SCEN_ORDER}
    values_present, values_missing = [], []
    sizes_all = {sc: [] for sc in SCEN_ORDER}

    for m in MT_VALUES:
        try:
            idx, data = _load_builds(store, dict(p, mt=m))
        except Exception:
            values_missing.append(m)
            continue
        values_present.append(m)
        for sc in SCEN_ORDER:
            for r in (data.get(sc) or []):
                kept = cl.keep_designated_clades(r["clades"])
                n = len(kept)
                nlin = np.array([c["nlin"] for c in kept], float)
                ntx = np.array([c["ntx"] for c in kept], float)
                mdiv = np.array([c.get("max_div", 0) for c in kept], float)
                sizes_all[sc].extend(nlin.tolist())
                out[sc]["seeds"].setdefault(r["seed"], dict(
                    seed=r["seed"], per_value=[], n_tx_run=r["n_tx_run"],
                    run_max_div=int(r.get("max_div", 0))))
                out[sc]["seeds"][r["seed"]]["per_value"].append(dict(
                    m=m, clades=n,
                    singleton_frac=float((nlin <= 1).sum() / max(n, 1)),
                    median_size=float(np.median(nlin)) if n else 0.0,
                    max_div=float(mdiv.max()) if n else 0.0,
                    never_tx_frac=float((ntx <= 0).sum() / max(n, 1)),
                ))

    for sc in SCEN_ORDER:
        series = list(out[sc]["seeds"].values())
        out[sc] = dict(runs=len(series), seeds=series,
                       size_bins=_doubling_bins(sizes_all[sc]),
                       median=_median_per_value(series, values_present))
    return dict(kind="across", params=p, values=values_present,
                values_missing=values_missing, values_full=MT_VALUES,
                metrics=ACROSS_METRICS, scenarios=out,
                note=("each m_T value is its own clade build, because rule B makes "
                      "m_T part of designation"
                      + (f" -- {len(values_missing)} of {len(MT_VALUES)} values "
                         f"are not built yet" if values_missing else "")))


def _median_per_value(series, values=None):
    """The median across runs at each m_T value, for the bold line the artifact drew
    over the faint per-seed ones. Median over RUNS, never over clades."""
    if not series:
        return []
    out = []
    for j, m in enumerate(values if values is not None else MT_VALUES):
        row = dict(m=m)
        for key in ACROSS_METRICS:
            vals = [s["per_value"][j][key] for s in series
                    if j < len(s["per_value"])]
            row[key] = float(np.median(vals)) if vals else 0.0
        row["n"] = len(series)
        out.append(row)
    return out


def values_on_disk(store, p):
    """Which m_T values of this (k, win, prom) have a build, cheaply.

    The answer is part of the across-the-set memo key: the tab now reads one
    build per value, so a newly built one has to invalidate it. Panel 4 had the
    same problem and solved it by caching per row; there is one answer here, so
    the key carries the set instead.
    """
    return [m for m in MT_VALUES
            if os.path.isdir(_build_dir(store, dict(p, mt=m)))]


def across(store, q):
    """The across-the-set answer, memoised like a panel: it reads one clade
    build per m_T value, tens of MB of JSON each, so it must not be recomputed
    per redraw."""
    p = clean_params(q)
    present = values_on_disk(store, p)
    return _memo(store, "across", [p, present],
                 lambda: across_the_set(store, p))


EXP_NUM = 1
EXP_NAME = compute.EXP_NAME


def _daily_series(store, k, scen, seed, mt):
    """One run's per-clade DAILY trajectories, memoised.

    Keyed on only what the clustering depends on -- k, the run, m_T -- so moving
    the smoothing window or the prominence bar does not re-cluster. The whole
    record is served and the caller slices it, which is why `window` is not in
    the key.
    """
    def build():
        ssod = compute.ssod_for(scen, seed, EXP_NUM, EXP_NAME)
        if ssod is None:
            raise FileNotFoundError(f"no simulation output for {scen}/{seed}")
        # m_T must be in the params, not only in the argument: run_context
        # designates under rule B from p["mt"], and leaving it out silently
        # designated on divergence alone -- 199 groups where the build has 10
        payload, err = compute.run_daily(
            ssod, dict(k=k, mt=mt, window="all"), mt)
        if err:
            raise RuntimeError(err)
        return payload

    return _memo(store, "daily", [k, scen, seed, mt], build)


def run_daily(store, q):
    """One run's per-clade daily trajectories AND its episodes, for the Sweep
    time tab.

    Two things used to go wrong here. The tab drew its episodes from runs/,
    which samples about every 9.2 days, while panel 3 measured on the daily
    record -- two numbers for one definition. And the browser carried its own
    reimplementation of scipy's prominence test, which chose a different
    trough: on HIV_high seed 0 Clade_0.366.5.11 it started the episode at index
    496 where compute.py starts it at 600, a 71-day difference in the rise.

    So the episodes are measured HERE, by the same `compute.episodes_of` and
    `compute.rise_of` the clade build calls, and the browser draws what it is
    given. The series is memoised on the clustering alone, so a slider move
    costs a peak-finding pass over a few hundred points, not a re-cluster.
    """
    p = clean_params(q)
    scen = q.get("scenario") or SCEN_ORDER[0]
    if scen not in SCEN_ORDER:
        raise ValueError(f"unknown scenario {scen!r}")
    idx, data = _load_builds(store, p)
    runs = data.get(scen) or []
    if not runs:
        raise FileNotFoundError(f"{scen} is not in this build")
    i = int(q.get("seed", 0))
    if not 0 <= i < len(runs):
        raise ValueError(f"seed index {i} outside 0..{len(runs) - 1}")
    seed = runs[i]["seed"]
    mt = float(q.get("mt", compute.DEFAULTS["mt"]))
    mf = read_mf(q)
    burnin = str(q.get("burnin", "1")) not in ("0", "", "false", "False")

    ser = _daily_series(store, int(p["k"]), scen, seed, mt)
    t_all = np.asarray(ser["t"], float)
    phi = ser.get("phi")
    j0 = 0
    if burnin and phi is not None:
        j0 = int(np.searchsorted(t_all, float(phi), side="left"))
    t = t_all[j0:]

    out = []
    for c in ser["clades"]:
        y = np.asarray(c["y"], float)[j0:]
        rec = dict(c, y=c["y"], j0=j0, eps=[], main=None, h=0.0,
                   peak=float(y.max()) if len(y) else 0.0)
        if len(t) >= 3:
            ys = compute.smooth(y, t, float(p["win"]))
            h, eps = compute.episodes_of(ys, float(p["prom"]))
            rs = [dict(r, cens=c) for r, c in
                  ((compute.rise_of(ys, t, a, b, float(p["lo"]), float(p["hi"])), c)
                   for a, b, c in eps) if r is not None]
            # before m_F, so the tab can say WHY a clade has no episode
            rec["n_found"] = len(rs)
            # m_F drops an episode that never reaches that share of the viral
            # population, before the emergence pick, so the headline episode is
            # chosen among the ones that count.
            rs = [r for r in rs if r["top"] >= mf]
            rec["h"] = round(float(h), 6)
            if rs:
                imax = int(np.argmax(ys))
                emg = min(rs, key=lambda r: (not (r["a"] <= imax <= r["b"]),
                                             -r["amp"]))
                # indices are relative to the slice, and `j0` travels with
                # them, so the browser slices where this did rather than
                # recomputing the burn-in index and landing one sample off
                rec["eps"] = [dict(a=int(r["a"]), b=int(r["b"]),
                                   ia=int(r["ia"]), ib=int(r["ib"]),
                                   days=round(float(r["days"]), 4),
                                   amp=round(float(r["amp"]), 6),
                                   top=round(float(r["top"]), 6),
                                   emg=bool(r is emg),
                                   cens=bool(r["cens"])) for r in rs]
                rec["main"] = next(j for j, r in enumerate(rs) if r is emg)
        out.append(rec)

    return dict(kind="run_daily", params=p, scenario=scen, seed_index=i,
                seed=seed, mt=mt, burnin=bool(burnin), j0=j0,
                t=ser["t"], phi=phi, end=ser.get("end"),
                n_tx_run=ser.get("n_tx_run"), n_groups=ser.get("n_groups"),
                n_kept=ser.get("n_kept"), capped=ser.get("capped"),
                clades=out)


def build_run(store, q):
    """The clade build's own record for ONE run: its clades and its episodes,
    exactly as the pass wrote them.

    Added for tests/studio/check_episodes.js, which pins the browser's episode
    detection to compute.py's. Every other endpoint aggregates, so there was no
    way to compare the two implementations run by run -- and an aggregate can
    agree while the per-clade numbers do not.
    """
    p = clean_params(q)
    scen = q.get("scenario") or SCEN_ORDER[0]
    if scen not in SCEN_ORDER:
        raise ValueError(f"unknown scenario {scen!r}")
    idx, data = _load_builds(store, p)
    runs = data.get(scen) or []
    if not runs:
        raise FileNotFoundError(f"{scen} is not in this build")
    i = int(q.get("seed", 0))
    if not 0 <= i < len(runs):
        raise ValueError(f"seed index {i} outside 0..{len(runs) - 1}")
    r = runs[i]
    return dict(kind="build_run", params=p, scenario=scen, seed_index=i, run=r)
