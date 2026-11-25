"""Episode, peak and hidden-time detection with every rule exposed, for the
Sweep time tab's exploration controls (/api/explore).

Figure 5 reads the clade build, which measured its episodes at ONE set of
rules (compute.py, clade_preprocess.episodes_of). /api/run_daily re-measures on
the daily record but still needs a build at the requested window and
prominence, so the Sweep sliders could only reach the grid's own values. This
module measures on the same daily record with no build at all, and lets every
rule move:

    smoothing   kernel (box | gauss | savgol | median | none), window in days,
                and whether to smooth the frequency or its log10
    peaks       prominence relative to the clade's height AND/OR absolute, a
                minimum separation and a minimum width in days, censored peaks
    start       where an episode begins: the last trough deeper than a bar
                (Figure 5), the last day below a fraction of the peak, or
                scipy's prominence base
    anchor      where hidden time starts: the clade's intra-host emergence
                (Time_emergence of its root lineage -- what Figure 5 uses) or
                the first transmission of any of its lineages

At the build's values every rule reduces to the build's own code, called, not
copied: `cp.smooth`, scipy's `find_peaks` with the same arguments, the same
trough bar and censored test, `cp.rise_of`, the same m_F filter and emergence
pick -- with one deliberate difference: a clade still rising on the last day
counts as a censored peak (see `episodes`). figures_gui/check_explore.py holds
the two to equality against /api/run_daily on every other clade.

Not hashed into any cache key, on purpose: nothing here is stored. compute.py,
panels.py and clade_preprocess.py are untouched -- their AST names the builds
and panel memos, and an edit there orphans them.
"""
import functools
import hashlib
import json
import os
import socket

import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks, peak_widths, savgol_filter

import codehash as ch
import compute
import panels
import clade_preprocess as cp
import simplicity.clustering as cl
import simplicity.dir_manager as dm
from _scenarios import one_sod

KERNELS = ("box", "gauss", "savgol", "median", "none")
SPACES = ("lin", "log")
STARTS = ("trough", "frac", "base")
ANCHORS = ("emergence", "transmission")
# which clades the burn-in day lets through. trunc (Figure 5): measure from the
# burn-in day and drop a clade already circulating then. start: measure on the
# whole record and keep a clade whose breakout STARTS after the burn-in day --
# the browser applies that test, since it picks the breakout
GATES = ("trunc", "start")
# which dips may start an episode under start=trough. any (the build's): every
# dip clearing the trough bar, so a wobble at the top of a plateau can start an
# episode a day before its peak. low: only a dip at least the prominence bar
# below the peak -- one the episode can be said to rise from
DIPS = ("any", "low")
LOG_FLOOR = 1e-6          # added before log10, so an absent clade is -6, not -inf

# the rules Figure 5 measured with; win, prom and lo come from the build axes
FIG5 = dict(kernel="box", space="lin", prom_abs=0.0, sep=0.0, width=0.0,
            cens=True, start="trough", trough=cp.TROUGH_PROM, anchor="emergence",
            gate="trunc", dip="any")


def read_params(q):
    """The detection rules from a query, validated. k and m_T name the
    partition, as everywhere else; the rest is this module's own."""
    def num(name, default, lo, hi):
        v = float(q.get(name, default))
        if not lo <= v <= hi:
            raise ValueError(f"{name}={v} outside {lo}..{hi}")
        return v

    def choice(name, options, default):
        v = q.get(name, default)
        if v not in options:
            raise ValueError(f"{name} must be one of {options}, got {v!r}")
        return v

    k = int(q.get("k", compute.DEFAULTS["k"]))
    if k not in (0, 3, 5, 7, 10, 15, 20):
        raise ValueError(f"k={k} is not one of the k values")
    lo = num("lo", 0.10, 0.0, 0.49)
    return dict(
        k=k, mt=num("mt", compute.DEFAULTS["mt"], 0.0, 1.0),
        kernel=choice("kernel", KERNELS, FIG5["kernel"]),
        space=choice("space", SPACES, FIG5["space"]),
        win=num("win", 15.0, 1.0, 365.0),
        prom=num("prom", 0.10, 0.0, 1.0),
        prom_abs=num("prom_abs", 0.0, 0.0, 1.0),
        sep=num("sep", 0.0, 0.0, 1000.0),
        width=num("width", 0.0, 0.0, 1000.0),
        cens=str(q.get("cens", "1")) not in ("0", "", "false", "False"),
        start=choice("start", STARTS, FIG5["start"]),
        trough=num("trough", FIG5["trough"], 0.0, 1.0),
        frac=num("frac", 0.10, 0.0, 1.0),
        lo=lo, hi=round(1.0 - lo, 6),
        anchor=choice("anchor", ANCHORS, FIG5["anchor"]),
        gate=choice("gate", GATES, FIG5["gate"]),
        dip=choice("dip", DIPS, FIG5["dip"]),
        mf=panels.read_mf(q),
        burnin=str(q.get("burnin", "1")) not in ("0", "", "false", "False"))


def smooth(y, t, kernel, win, space):
    """The series the peaks are found on. `win` is in days for every kernel:
    the box's full width (as Figure 5), the Gaussian's full width at half
    maximum, the Savitzky-Golay and median windows' width. log space smooths
    log10(y + 1e-6), so a clade's early climb through 1e-4..1e-2 weighs as
    much as its last doubling."""
    y = np.asarray(y, float)
    if kernel == "none" or len(y) < 3:
        return y.copy()
    z = np.log10(y + LOG_FLOOR) if space == "log" else y
    if kernel == "box":
        z = cp.smooth(z, t, win)
    else:
        step = (t[-1] - t[0]) / (len(t) - 1)
        w = max(1, int(round(win / max(step, 1e-9))))
        if kernel == "gauss":
            z = gaussian_filter1d(z, sigma=w / 2.3548, mode="nearest")
        elif kernel == "median":
            z = median_filter(z, size=w | 1, mode="nearest")
        elif kernel == "savgol":
            wl = min(max(w | 1, 5), len(z) if len(z) % 2 else len(z) - 1)
            z = savgol_filter(z, wl, 2, mode="interp") if wl >= 5 else z
    out = np.power(10.0, z) - LOG_FLOOR if space == "log" else z
    return np.clip(out, 0.0, None) if kernel != "box" or space == "log" else out


def episodes(ys, t, r):
    """(height, [(start, peak, censored)]) under rules `r`: clade_preprocess.
    episodes_of, except that a clade still rising on the last day counts as a
    censored peak there. The build leaves it out (no observed peak); for a
    breakout only its start matters, and that is observed."""
    ys = np.asarray(ys, float)
    n = len(ys)
    h = float(np.max(ys)) if n else 0.0
    if not (h > 0):
        return 0.0, []
    bar = max(r["prom"] * h, r["prom_abs"])
    step = (t[-1] - t[0]) / (n - 1)
    kw = dict(prominence=bar)
    if r["sep"] > 0:
        kw["distance"] = max(1, int(round(r["sep"] / step)))
    if r["width"] > 0:
        kw["width"] = r["width"] / step
    pk, props = find_peaks(ys, **kw)
    peaks = [(int(i), False, int(lb)) for i, lb in zip(pk, props["left_bases"])]
    if r["cens"]:
        m = int(np.argmax(ys >= h * (1 - cp.CENSOR_TOL)))
        if h - ys[:m + 1].min() >= bar and h - ys[m:].min() < bar:
            peaks = sorted(peaks + [(m, True, int(np.argmin(ys[:m + 1])))])
    if not peaks:
        return h, []
    if r["start"] == "trough":
        mins, _ = find_peaks(-ys, prominence=r["trough"] * h)
        mins = np.concatenate([[0], mins])
    out = []
    for i, c, base in peaks:
        if r["start"] == "trough":
            before = mins[mins < i]
            if r["dip"] == "low":
                before = before[ys[before] <= ys[i] - bar]
            a = int(before[-1]) if len(before) else 0
        elif r["start"] == "frac":
            below = np.nonzero(ys[:i] <= r["frac"] * ys[i])[0]
            a = int(below[-1]) if len(below) else 0
        else:
            a = base
        out.append((a, i, c))
    return h, out


def widths(ys, t, r):
    """Every peak the prominence (and separation) rules pass, with the width
    Min width is tested on: scipy's at half the peak's prominence. Per peak
    [peak, left, right, height, days, kept], indices into `ys`, left and right
    interpolated; kept is False where Min width dropped it. The censored peak is
    not here: it never comes down, so it has no width and is never tested."""
    ys = np.asarray(ys, float)
    n = len(ys)
    h = float(np.max(ys)) if n else 0.0
    if n < 3 or not (h > 0):
        return []
    bar = max(r["prom"] * h, r["prom_abs"])
    step = (t[-1] - t[0]) / (n - 1)
    kw = dict(prominence=bar)
    if r["sep"] > 0:
        kw["distance"] = max(1, int(round(r["sep"] / step)))
    pk, props = find_peaks(ys, **kw)
    if not len(pk):
        return []
    w, wh, li, ri = peak_widths(ys, pk, rel_height=0.5, prominence_data=(
        props["prominences"], props["left_bases"], props["right_bases"]))
    wmin = r["width"] / step if r["width"] > 0 else 0.0
    return [[int(p), round(float(a), 2), round(float(b), 2), round(float(y), 7),
             round(float(x * step), 2), bool(x >= wmin)]
            for p, x, y, a, b in zip(pk, w, wh, li, ri)]


# ------------------------------------------------------------------ caches
# On disk under <store>/cache/_explore/, keyed on the request and on the code
# that computes it: the named functions here (codehash.names_hash, as
# panels.scoped_hash does) plus compute.py's build hash, which covers the
# clustering. Editing a function re-keys only what it computes.

def _xmemo(store, tag, parts, deps, build):
    key = json.dumps([tag, parts, ch.names_hash(globals(), deps, extra=compute.code_hash())],
                     sort_keys=True, default=str)
    d = os.path.join(store, "cache", "_explore")
    os.makedirs(d, exist_ok=True)
    f = os.path.join(d, f"{tag}_{hashlib.sha256(key.encode()).hexdigest()[:16]}.json")
    if os.path.exists(f):
        try:
            with open(f) as fh:
                return json.load(fh)
        except Exception:                                        # noqa: BLE001
            pass
    out = build()
    tmp = f"{f}.{socket.gethostname()}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w") as fh:
            json.dump(out, fh)
        os.replace(tmp, f)
    except Exception:                                            # noqa: BLE001
        pass
    return out


def _run_meta_build(scen, seed, k, mt):
    """Per clade name: [t_emerge, t_tx1, nlong, ndef, founding].

    t_emerge is intra-host emergence, its root lineage's Time_emergence (the
    build's `tem`); t_tx1 the earliest infection inheriting any of its
    lineages, None if none did; nlong how many of its defining mutations arose
    in a long shedder, as compute._one counts them. From the run's raw output,
    clustered as the build clusters it."""
    ssod = compute.ssod_for(scen, seed, panels.EXP_NUM, panels.EXP_NAME)
    if ssod is None:
        raise FileNotFoundError(f"no simulation output for {scen}/{seed}: the "
                                f"exploration reads the run itself")
    ctx, err = compute.run_context(ssod, k, "all", mt)
    if err:
        raise RuntimeError(err)
    ind = ctx["ind"]
    lin = cl.transmitted_lineages(ind)
    first = (ind.loc[lin.index, "t_infection"].astype(float)
             .groupby(lin.astype(str)).min().to_dict())
    out = {}
    for c in ctx["cols"]:
        tx = [first[m] for m in ctx["c2l"].get(c, []) if m in first]
        d = ctx["pcm"].get(c)
        nlong = int((d["host_type"].astype(str).str.lower()
                     .str.startswith("long")).sum()) if d is not None and len(d) else 0
        out[str(c)] = [float(ctx["tem_of"].get(c, 0) or 0),
                       float(min(tx)) if tx else None, nlong,
                       int(ctx["ndef_of"].get(c, 0) or 0),
                       not isinstance(ctx["par_of"].get(c), str)]
    return out


_META = {}


def run_meta(store, scen, seed, k, mt):
    key = (store, scen, seed, k, mt)
    if key not in _META:
        if len(_META) > 2000:
            _META.clear()
        _META[key] = _xmemo(store, "meta", [k, scen, seed, mt], ("_run_meta_build",),
                            lambda: _run_meta_build(scen, seed, k, mt))
    return _META[key]


@functools.lru_cache(maxsize=8)
def _seed_names(scen):
    sod = one_sod(panels.EXP_NAME, scen, panels.EXP_NUM)
    if sod is None:
        raise FileNotFoundError(f"{scen} has no simulation output here")
    return [os.path.basename(d) for d in dm.get_seeded_simulation_output_dirs(sod)]


def seed_name(scen, i):
    names = _seed_names(scen)
    if not 0 <= int(i) < len(names):
        raise ValueError(f"seed index {i} outside 0..{len(names) - 1}")
    return names[int(i)]


_POOL = {}


def pool(store, scen, seed, k, mt):
    """Figure 5's candidates from one run, compact: the daily axis, phi, end,
    and each non-founding clade with a defining mutation, its series as an
    array. Held in memory, about 0.2 MB a run, so a panel over 250 runs does
    not re-read 250 daily records on every slider move."""
    key = (store, scen, seed, k, mt)
    if key not in _POOL:
        if len(_POOL) > 2000:
            _POOL.clear()
        ser = panels._daily_series(store, k, scen, seed, mt)
        meta = run_meta(store, scen, seed, k, mt)
        cl = []
        for c in ser["clades"]:
            m = meta.get(str(c.get("name", "")))
            if m is None or c["founding"] or int(c["ndef"]) <= 0:
                continue
            cl.append(dict(name=str(c["name"]), y=np.asarray(c["y"], float),
                           t_em=m[0], t_tx=m[1], nlong=m[2], ndef=int(c["ndef"])))
        _POOL[key] = dict(t=np.asarray(ser["t"], float), phi=ser.get("phi"),
                          end=ser.get("end"), clades=cl)
    return _POOL[key]


def detect(y, t, r):
    """Every episode of one clade on the axis `t` under rules `r`, unfiltered
    by m_F: (smoothed series, height, [rise dicts with `cens`])."""
    ys = smooth(y, t, r["kernel"], r["win"], r["space"])
    h, eps = episodes(ys, t, r)
    rs = [dict(e, cens=cens) for e, cens in
          ((cp.rise_of(ys, t, a, b, r["lo"], r["hi"]), cens) for a, b, cens in eps)
          if e is not None]
    return ys, h, rs


def _pre(ys, t, tem, a):
    """The smoothed curve's highest value from the anchor up to (not
    including) an episode's start: what the clade reached during the hidden
    time that episode would end. None when the start is at or before the anchor."""
    if tem is None:
        return None
    ja = int(np.searchsorted(t, float(tem)))
    return round(float(np.max(ys[ja:a])), 7) if a > ja else None


def _slice(t_all, phi, r):
    return int(np.searchsorted(t_all, float(phi), side="left")) \
        if r["burnin"] and r["gate"] == "trunc" and phi is not None else 0


def _trunc(r, phi, tem, y):
    # as compute._one: circulating above a trace at the burn-in day, with its
    # anchor before it, so its hidden time began outside the window
    return bool(r["burnin"] and r["gate"] == "trunc" and phi is not None and tem is not None
                and tem < float(phi) and len(y) and y[0] > 1e-4)


_RUNS = {}


def explore(store, q):
    """One run's clades, measured under the rules in `q`, in /api/run_daily's
    shape so the Sweep tab draws either. Each clade also carries both anchors,
    t_emerge and t_tx1; `tem` is whichever `anchor` picks. The last answers are
    held in memory: going back to a slider stop costs nothing."""
    r = read_params(q)
    scen = q.get("scenario") or compute.SCEN_ORDER[0]
    if scen not in compute.SCEN_ORDER:
        raise ValueError(f"unknown scenario {scen!r}")
    i = int(q.get("seed", 0))
    seed = seed_name(scen, i)
    key = json.dumps([store, scen, seed, r], sort_keys=True)
    if key in _RUNS:
        return _RUNS[key]
    ser = panels._daily_series(store, r["k"], scen, seed, r["mt"])
    meta = run_meta(store, scen, seed, r["k"], r["mt"])
    t_all = np.asarray(ser["t"], float)
    phi = ser.get("phi")
    j0 = _slice(t_all, phi, r)
    t = t_all[j0:]

    out = []
    for c in ser["clades"]:
        y = np.asarray(c["y"], float)[j0:]
        m = meta.get(str(c.get("name", "")))
        t_em, t_tx = (m[0], m[1]) if m else (None, None)
        tem = t_em if r["anchor"] == "emergence" else t_tx
        rec = dict(c, y=c["y"], j0=j0, eps=[], main=None, h=0.0,
                   peak=float(y.max()) if len(y) else 0.0,
                   tem=tem, t_emerge=t_em, t_tx1=t_tx, trunc=_trunc(r, phi, tem, y))
        if len(t) >= 3:
            ys, h, rs = detect(y, t, r)
            rec["n_found"] = len(rs)
            rs = [e for e in rs if e["top"] >= r["mf"]]
            rec["h"] = round(float(h), 6)
            if rs:
                imax = int(np.argmax(ys))
                emg = min(rs, key=lambda e: (not (e["a"] <= imax <= e["b"]), -e["amp"]))
                rec["eps"] = [dict(a=int(e["a"]), b=int(e["b"]),
                                   ia=int(e["ia"]), ib=int(e["ib"]),
                                   days=round(float(e["days"]), 4),
                                   amp=round(float(e["amp"]), 6),
                                   top=round(float(e["top"]), 6),
                                   emg=bool(e is emg), cens=bool(e["cens"]),
                                   pre=_pre(ys, t, tem, int(e["a"])))
                              for e in rs]
                rec["main"] = next(j for j, e in enumerate(rs) if e is emg)
            # the smoothed series the peaks were found on, so the tab draws
            # what was measured rather than re-smoothing with a box
            rec["ys"] = [round(float(v), 7) for v in ys]
            rec["wid"] = widths(ys, t, r)
        out.append(rec)

    res = dict(kind="explore", rules=r, scenario=scen, seed_index=i, seed=seed,
               mt=r["mt"], burnin=r["burnin"], j0=j0, t=ser["t"], phi=phi,
               end=ser.get("end"), n_tx_run=ser.get("n_tx_run"),
               n_groups=ser.get("n_groups"), n_kept=ser.get("n_kept"),
               capped=ser.get("capped"), clades=out)
    if len(_RUNS) >= 48:
        _RUNS.pop(next(iter(_RUNS)))
    _RUNS[key] = res
    return res


def _panel_build(store, r, seeds):
    out = {}
    for sc in compute.SCEN_ORDER:
        runs = []
        for seed in _seed_names(sc)[:seeds]:
            P = pool(store, sc, seed, r["k"], r["mt"])
            phi = P["phi"]
            # the build measures on its post-burn-in window, which is empty for
            # a run that never reaches phi: its clades stay, with no episodes
            empty = r["burnin"] and phi is None
            j0 = _slice(P["t"], phi, r)
            t = P["t"][j0:]
            clades = []
            for c in P["clades"]:
                y = c["y"][j0:] if not empty else c["y"][:0]
                tem = c["t_em"] if r["anchor"] == "emergence" else c["t_tx"]
                if tem is None or _trunc(r, phi, tem, y):
                    continue
                eps = []
                if len(t) >= 3 and not empty:
                    ys, _, rs = detect(y, t, r)
                    eps = sorted([round(float(t[e["a"]]), 2), round(float(e["top"]), 5),
                                  round(float(e["amp"]), 5), round(float(e["days"]), 2),
                                  _pre(ys, t, tem, int(e["a"]))]
                                 for e in rs)
                clades.append([int(c["nlong"]), int(c["ndef"]), round(float(tem), 2),
                               eps, c["name"]])
            runs.append(dict(seed=seed, phi=phi, end=P["end"], clades=clades))
        out[sc] = dict(runs=runs)
    return out


def explore_panel(store, q):
    """Figure 5 under the rules in `q`: /api/panel?n=3b's shape -- per
    scenario, per run, each candidate clade as [nlong, ndef, anchor, [[a_day,
    top, amp, days], ...], name] -- so the Figure 5 tab plots it unchanged and
    applies F, H, f and m_F itself. Cached on disk per rule set; the first
    answer for new rules re-measures every run (seconds once the runs' daily
    records and anchors are on disk, see warm_explore.py)."""
    r = read_params(q)
    seeds = int(q.get("seeds", 50))
    rules = {k: v for k, v in r.items() if k != "mf"}
    scen = _xmemo(store, "panel", [rules, seeds],
                  ("_panel_build", "detect", "smooth", "episodes", "_slice", "_trunc", "_pre",
                   "pool", "_run_meta_build"),
                  lambda: _panel_build(store, r, seeds))
    return dict(kind="explore_panel", rules=r, seeds=seeds, scenarios=scen)
