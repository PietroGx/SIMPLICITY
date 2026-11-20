#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-clade and per-seed metric tables for the figure-2 plot builder.

Every metric is computed by the function the figure itself uses --
_peak_values, _survival_durations, _growth_crossing_times from
long_shedders_preprocess, and phi_day from fig2_preprocess_data. Nothing is
re-derived here; this module only assembles and encodes.

The pool is fig2's pool: non-root clades with at least one defining mutation.
Each metric is emitted twice, over the whole run and over the post-burn-in
window, so the page can switch without a second export.
"""
import os
import sys
import json
import base64
import shutil
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
import fig2_preprocess_data as pre
from clade_preprocess import (_peak_values, _survival_durations,
                              _growth_crossing_times,
                              smooth, episodes_of, rise_of)
from _scenarios import EXP_NAME, resolve_scenarios, one_sod, scenario_meta

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
THRESHOLDS = [0, 3, 5, 7, 10, 15, 20]
LAB = {"standard": 0, "long": 1, "mixed": 2, "founder": 3}
NEVER = 65535                      # packed uint16: "not measurable"
RISE_WIN = 15                      # days, centred, for the rise-time smoothing
RISE_PROM = 0.10                   # episode prominence bar, as the studio's


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def _u16(x, scale=1.0):
    """Pack to uint16, keeping "not measurable" distinct from a measured zero.

    Real values clip into 0..65534; a non-finite one -- which is what a metric
    yields when its own precondition fails, e.g. a clade that never crosses 1%
    has no crossing time -- packs to NEVER.

    Both facts used to land on the same number twice over: nan_to_num turned an
    unmeasurable metric into 0, indistinguishable from a measured 0, and NEVER
    itself was clipped to 65534 on the way out, so the page could only
    recognise it by a >= 10000 test. Measured on 2,657 clades at k=7: 2,643
    (99.5%) of grow values are this sentinel.
    """
    a = np.asarray(x, float) * scale
    out = np.clip(np.nan_to_num(a, nan=0.0, posinf=65534, neginf=0.0),
                  0, 65534).astype(np.uint16)
    out[~np.isfinite(a)] = NEVER
    return out





def _metrics(F):
    """Per-clade trajectory metrics for one frequency frame.

    peak / survival / growth come from the figure's own functions. The rest are
    assembled here: area under the trajectory, where it ends, the steepest
    logit slope (a selection coefficient), and the width at half its own peak.
    """
    empty = ({},) * 9
    if F is None or F.empty or F.shape[1] == 0:
        return empty
    peak = _peak_values(F)
    surv = _survival_durations(F)
    grow = _growth_crossing_times(F)
    t = F.index.to_numpy(float)
    tpk, auc, fin, rmax, whm = {}, {}, {}, {}, {}
    risef = {}
    for c in F.columns:
        y = np.nan_to_num(F[c].to_numpy(float))
        # first day AT the top, not the first day argmax happens to pick: a
        # clade at fixation holds its maximum for weeks (up to 181 days here),
        # and plain argmax resolved that tie on floating-point noise -- a
        # change of summation order alone moved it by up to 27 days.
        ymax = y.max()
        tpk[c] = float(t[(y >= ymax * (1 - 1e-9)).argmax()])
        auc[c] = float(np.trapezoid(y, t)) if len(t) > 1 else 0.0
        fin[c] = float(y[-1])
        # steepest rise in logit space, per day: the selection-coefficient scale
        if len(t) > 1:
            e = 1e-4
            z = np.log(np.clip(y, e, 1 - e) / (1 - np.clip(y, e, 1 - e)))
            dt = np.diff(t)
            sl = np.divide(np.diff(z), dt, out=np.zeros_like(dt), where=dt > 0)
            rmax[c] = float(np.max(sl)) if sl.size else float("nan")
        else:
            rmax[c] = float("nan")
        # ONE rise: the emergence episode's, the same definition the studio
        # reports and `compute._one` writes into a clade build. This used to be
        # two columns from two local functions, neither of which was the
        # studio's -- see CHANGELOG v2.5.7.
        if ymax > 0 and len(t) > 2:
            ys = smooth(y, t, RISE_WIN)
            h, eps = episodes_of(ys, RISE_PROM)
            rs = [r for r in (rise_of(ys, t, a, b, 0.10, 0.90) for a, b, _ in eps)
                  if r is not None]
            if rs:
                imax = int(np.argmax(ys))
                emg = min(rs, key=lambda r: (not (r["a"] <= imax <= r["b"]),
                                             -r["amp"]))
                risef[c] = float(emg["days"])
        half = 0.5 * ymax
        m = y >= half
        whm[c] = float(t[len(m) - 1 - m[::-1].argmax()] - t[m.argmax()]) if m.any() else 0.0
    return peak, surv, grow, tpk, auc, fin, rmax, whm, risef


def _one(job):
    scenario, ssod, thresholds, min_tx_frac = job
    seed = os.path.basename(ssod)
    sys.setrecursionlimit(300000)
    try:
        ph = om.read_phylogenetic_data(ssod)
        lf = om.read_lineage_frequency(ssod)
        ind = om.read_individuals_data(ssod)
    except Exception as exc:
        return scenario, seed, None, None, f"read: {exc}"

    try:
        phi = pre.phi_day(ssod)
        end = float(om.read_final_time(ssod))
    except Exception as exc:
        return scenario, seed, None, None, f"window: {exc}"
    cutoff = phi if phi is not None else float("inf")
    fig2 = bool(phi is not None and end >= pre.MIN_FINAL_TIME
                and (end - phi) >= pre.MIN_WINDOW_DAYS)

    # The old code pivoted to a dense times x LINEAGES frame and then summed
    # column groups into clades. That frame is ~97% empty (HIV_high: 1096 x
    # 12,333 float64 holding 2.9% non-null) and at k=0 the concat that followed
    # built one pandas Series per lineage. Factorise once here and let bincount
    # accumulate straight into times x CLADES, per threshold, below.
    lin_code, lin_names = pd.factorize(lf["Lineage_name"], sort=False)
    t_code, t_vals = pd.factorize(lf["Time_sampling"], sort=True)
    freq = lf["Frequency_at_t"].to_numpy(float)
    n_t = len(t_vals)
    del lf
    lin2mut = cl.build_lineage_to_mutation_dict(ph)
    lin2host = dict(zip(ph["Lineage_name"], ph["Host_type"]))

    tx_all = ind["inherited_lineage"].value_counts().to_dict()
    post = ind[ind["t_infection"].astype(float) >= cutoff]
    tx_win = post["inherited_lineage"].value_counts().to_dict()

    # rule B: a lineage founds a clade only if its subtree also carries at
    # least min_tx_frac of this run's own transmission events. Designation, not
    # a filter afterwards -- a refused lineage FOLDS into its surviving
    # ancestor and brings its mutations with it, where dropping would discard
    # it. Measured at k=7 and 1%: the same clade counts, 1.12-1.51x more
    # lineages covered. min_tx_frac = 0 is divergence alone.
    support = min_support = None
    if float(min_tx_frac) > 0:
        support = cl.subtree_transmission_counts(ph, ind)
        min_support = cl.transmission_floor(ind, float(min_tx_frac))

    rows, per_thr = [], {}
    for k in thresholds:
        c2l, _l2c, pcm, meta = cl.cluster_lin_into_clades_with_meta(
            ph, shared_mut_threshold=k, support=support,
            min_support=min_support or 0.0)
        # k=0 gives one clade per lineage, so a lineage whose edge adds no
        # mutation under a parent whose edge adds none either produces a clade
        # with no definers, which label_clades_from_definers refuses (it only
        # tolerates that for the root). Those clades are outside fig2's pool
        # anyway -- it takes n_defining > 0 -- so they are withheld from the
        # labelling. The function is per-clade, so every clade that survives
        # gets exactly the label it would have got.
        par_of = dict(zip(meta["clade"], meta["parent_clade"]))
        votable = {c: d for c, d in pcm.items()
                   if not ((d is None or len(d) == 0) and isinstance(par_of.get(c), str))}
        labels_s, _ = cl.label_clades_from_definers(votable, meta)
        labels = labels_s.to_dict() if labels_s is not None else {}

        nonroot = meta.loc[meta["parent_clade"].notna() & (meta["n_defining"] > 0), "clade"].tolist()
        nchild = meta["parent_clade"].value_counts().to_dict()
        root_of = dict(zip(meta["clade"], meta["root_lineage"]))
        ndef_of = dict(zip(meta["clade"], meta["n_defining"]))
        tstart = dict(zip(meta["clade"], meta["start_time"]))

        clades = list(c2l)
        cpos = {c: j for j, c in enumerate(clades)}
        lin2cl = np.full(len(lin_names), -1, np.int64)
        for j, nm in enumerate(lin_names):
            c = _l2c.get(nm)
            if c is not None:
                lin2cl[j] = cpos[c]
        row_cl = lin2cl[lin_code]
        keep = row_cl >= 0
        if keep.any():
            flat = t_code[keep] * len(clades) + row_cl[keep]
            M = np.bincount(flat, weights=freq[keep],
                            minlength=n_t * len(clades)).reshape(n_t, len(clades))
            # a clade with no frequency record at all was absent from the old
            # concat, so it must stay absent here for reindex/dropna to agree
            seen = np.zeros(len(clades), bool)
            seen[row_cl[keep]] = True
            FC = pd.DataFrame(M[:, seen], index=t_vals,
                              columns=[c for c, f in zip(clades, seen) if f])
        else:
            FC = pd.DataFrame()
        F = FC.reindex(columns=[c for c in nonroot if c in FC.columns]).dropna(axis=1, how="all")
        FW = F.loc[F.index.astype(float) >= cutoff] if not F.empty else F

        pk, sv, gw, tp, ac, fn, rm, wh, rf = _metrics(F)
        pkw, svw, gww, tpw, acw, fnw, rmw, whw, rfw = _metrics(FW)

        # how deep the clade's own descendants reach in the clade tree
        kids = {}
        for c0, p0 in zip(meta["clade"], meta["parent_clade"]):
            if isinstance(p0, str):
                kids.setdefault(p0, []).append(c0)
        def _sdepth(c0, seen=None):
            seen = seen or set()
            if c0 in seen:
                return 0
            seen.add(c0)
            return 1 + max((_sdepth(x, seen) for x in kids.get(c0, ())), default=-1)

        rec = []
        for c in nonroot:
            mem = c2l.get(c, ())
            rootlin = root_of.get(c)
            d = pcm.get(c)
            nl = int((d["host_type"].astype(str).str.lower().str.startswith("long")).sum()) if d is not None and len(d) else 0
            # the denominator counts what the numerator counts: ROWS of the
            # definer frame. `n_defining` counts DISTINCT mutations, and the
            # frame carries one row per (node, mutation), so a mutation
            # recurring on two edges of one segment made nl/nd exceed 1 --
            # which then wrapped in the uint8 cast rather than clipping, and a
            # fully long-origin clade would have read 13.7% instead of 100%.
            # The two are equal whenever no mutation recurs, which is every
            # clade measured so far (0 of 2,283,098 records), so this changes
            # nothing on the current data and is correct if one ever does.
            nd = (len(d) if d is not None else 0) or 1
            rmut = lin2mut.get(rootlin, frozenset())
            mdiv = max((len(lin2mut.get(m, frozenset()) - rmut) for m in mem), default=0)
            rec.append((
                LAB.get(cl.resolve_clade_label(labels, c), 0),
                LAB["long"] if cl.normalise_host_type(lin2host.get(rootlin)) == "long" else LAB["standard"],
                min(255, max(0, int(round(255 * nl / nd)))),
                int(ndef_of.get(c, 0)), len(mem),
                sum(int(tx_all.get(m, 0)) for m in mem),
                sum(int(tx_win.get(m, 0)) for m in mem),
                float(tstart.get(c, 0) or 0),
                pk.get(c, 0.0), pkw.get(c, 0.0),
                sv.get(c, 0.0), svw.get(c, 0.0),
                gw.get(c, float("nan")), gww.get(c, float("nan")),
                tp.get(c, 0.0), tpw.get(c, 0.0),
                int(nchild.get(c, 0)), int(mdiv),
                ac.get(c, 0.0), acw.get(c, 0.0),
                fn.get(c, 0.0), fnw.get(c, 0.0),
                rm.get(c, 0.0), rmw.get(c, 0.0),
                _sdepth(c),
                rf.get(c, float("nan")), rfw.get(c, float("nan")),
                wh.get(c, 0.0), whw.get(c, 0.0),
            ))
        rows.append((k, rec))
        per_thr[k] = dict(n_clades=len(rec),
                          n_long=sum(1 for r in rec if r[0] == LAB["long"]),
                          n_std=sum(1 for r in rec if r[0] == LAB["standard"]),
                          n_mixed=sum(1 for r in rec if r[0] == LAB["mixed"]))

    try:
        subs = pre._substitutions_by_type(ind, cutoff=None)
        subs_w = pre._substitutions_by_type(ind, cutoff=cutoff if phi is not None else None)
    except Exception:
        subs, subs_w = {"long": 0, "standard": 0}, {"long": 0, "standard": 0}

    # Clades and their defining mutations are not time series and are never
    # windowed. The host counts are: n_infected/n_long_hosts span the whole run,
    # the _w pair only individuals infected on or after phi. Both are a check on
    # the design value, which is the number a panel should quote -- the
    # long-shedder ratio is a fixed per-individual probability, so the realised
    # fraction is the same before and after phi to three decimals.
    _lg = ind["type"].astype(str) == "long_shedder"
    _post = ind["t_infection"].astype(float) >= cutoff

    seedrec = dict(scenario=scenario, seed=seed, phi=None if phi is None else round(phi, 2),
                   end=round(end, 2), fig2=fig2,
                   n_lineages=int(len(ph)), n_infected=int(len(ind)),
                   n_long_hosts=int(_lg.sum()),
                   n_infected_w=None if phi is None else int(_post.sum()),
                   n_long_hosts_w=None if phi is None else int((_post & _lg).sum()),
                   subs_long=int(subs.get("long", 0)), subs_std=int(subs.get("standard", 0)),
                   subs_long_w=int(subs_w.get("long", 0)), subs_std_w=int(subs_w.get("standard", 0)),
                   per_thr={str(k): v for k, v in per_thr.items()})
    return scenario, seed, rows, seedrec, None


COLS = ["lab_vote", "lab_root", "plong", "ndef", "nlin", "ntx", "ntx_w", "tem",
        "peak", "peak_w", "surv", "surv_w", "grow", "grow_w", "tpeak", "tpeak_w",
        "nchild", "mdiv", "auc", "auc_w", "ffin", "ffin_w", "rmax", "rmax_w", "sdepth",
        "risef", "risef_w", "whm", "whm_w"]


# The artifact host refuses a published text file over 16 MB, and the k=0
# baseline is one clade per lineage — HIV_high lands near 39 MB. Split such a
# table into row ranges; the page concatenates the parts back.
PART_MAX = 12_000_000
_U8SET = {"seed", "lab_vote", "lab_root", "plong"}


def _shard(obj, cap=PART_MAX):
    """One encoded table -> [(part index or None, json text)] small enough to publish."""
    blob = json.dumps(obj, separators=(",", ":"))
    if len(blob) <= cap:
        return [(None, blob)]
    nrow = int(obj["n"])
    nparts = -(-len(blob) // cap)
    step = -(-nrow // nparts)
    cols = [c for c in obj if isinstance(obj[c], str) and c != "scenario"]
    dec = {c: np.frombuffer(base64.b64decode(obj[c]),
                            np.uint8 if c in _U8SET else np.uint16) for c in cols}
    out = []
    for i, lo in enumerate(range(0, nrow, step)):
        hi = min(lo + step, nrow)
        o = {k: v for k, v in obj.items() if not isinstance(v, str) or k == "scenario"}
        o["n"] = hi - lo
        for c in cols:
            o[c] = _b64(dec[c][lo:hi])
        out.append((i, json.dumps(o, separators=(",", ":"))))
    return out


def _write_table(outdir, sc, k, obj):
    """Write one table, split if need be, and return its index entry."""
    names, sz = [], 0
    for i, text in _shard(obj):
        name = f"clades_{sc}_k{k}.json" if i is None else f"clades_{sc}_k{k}.p{i}.json"
        with open(os.path.join(outdir, name), "w") as f:
            f.write(text)
        names.append(name)
        sz += os.path.getsize(os.path.join(outdir, name))
    return {"file": names[0] if len(names) == 1 else names, "n": obj["n"], "bytes": sz}


def _encode(a, sid):
    """Column-major binary for one (scenario, threshold).

    Takes the already-assembled (n, len(COLS)) float64 block and its per-row
    seed index, so the caller can stream it off disk one table at a time
    rather than holding every seed's records until the pool drains.
    """
    if a is None or a.shape[0] == 0:
        return None
    out = {"n": int(a.shape[0]), "seed": _b64(sid)}
    out["lab_vote"] = _b64(a[:, 0].astype(np.uint8))
    out["lab_root"] = _b64(a[:, 1].astype(np.uint8))
    # clipped, not merely cast: uint8 WRAPS, so an out-of-range share would
    # read as a small one rather than as an obvious error
    out["plong"] = _b64(np.clip(a[:, 2], 0, 255).astype(np.uint8))
    for j, name, scale in [(3, "ndef", 1), (4, "nlin", 1), (5, "ntx", 1), (6, "ntx_w", 1),
                           (7, "tem", 1), (8, "peak", 65534), (9, "peak_w", 65534),
                           (10, "surv", 1), (11, "surv_w", 1), (12, "grow", 1), (13, "grow_w", 1),
                           (14, "tpeak", 1), (15, "tpeak_w", 1), (16, "nchild", 1), (17, "mdiv", 1),
                           (18, "auc", 10), (19, "auc_w", 10),
                           (20, "ffin", 65534), (21, "ffin_w", 65534),
                           (22, "rmax", 1000), (23, "rmax_w", 1000), (24, "sdepth", 1),
                           (25, "risef", 1), (26, "risef_w", 1),
                           (27, "whm", 1), (28, "whm_w", 1)]:
        out[name] = _b64(_u16(a[:, j], scale))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=EXP_NAME)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=8)
    # Each seed spools to disk as it finishes and tables are encoded one at a
    # time, so the parent stays near 150 MB whatever the run length. A worker
    # still peaks around 0.9 GB on an HIV_high seed at every m_T value, so the
    # run costs roughly 0.15 + procs x 0.9 GB. Scenarios can be run one at a
    # time and merged into the same --out.
    ap.add_argument("--scenarios", type=str, default="",
                    help="comma-separated subset; results merge into an existing --out")
    ap.add_argument("--mt", type=float, default=0.01,
                    help="m_T for rule B designation, as a fraction of the "
                         "run's own transmission events. 0 designates on "
                         "divergence alone. One export is one m_T: it changes "
                         "the partition, so it is not a dimension of the "
                         "tables.")
    ap.add_argument("--reshard", action="store_true",
                    help="split oversized tables in an existing --out and exit")
    a = ap.parse_args()
    want = [x for x in a.scenarios.split(",") if x] or SCEN_ORDER
    os.makedirs(a.out, exist_ok=True)

    if a.reshard:
        ipath = os.path.join(a.out, "index.json")
        with open(ipath) as f:
            idx = json.load(f)
        done = 0
        for key, rec in idx["files"].items():
            if isinstance(rec["file"], list):
                continue
            src = os.path.join(a.out, rec["file"])
            if os.path.getsize(src) <= PART_MAX:
                continue
            with open(src) as f:
                obj = json.load(f)
            sc, k = key.split("|")
            ent = _write_table(a.out, sc, k, obj)
            if isinstance(ent["file"], list):
                os.remove(src)
            idx["files"][key] = ent
            done += 1
            print(f"[reshard] {key}: {rec['bytes']/1e6:.1f} MB -> "
                  f"{len(ent['file'])} parts, {ent['bytes']/1e6:.1f} MB")
        with open(ipath, "w") as f:
            json.dump(idx, f)
        print(f"[reshard] {done} table(s) split")
        sys.exit(0)

    present, _ = resolve_scenarios(a.exp_name, a.exp_num, required=SCEN_ORDER)
    jobs = []
    for sc in [s for s in SCEN_ORDER if s in present and s in want]:
        for s in dm.get_seeded_simulation_output_dirs(one_sod(a.exp_name, sc, a.exp_num))[:a.seeds]:
            jobs.append((sc, s, THRESHOLDS, a.mt))
    print(f"[fig2b] {len(jobs)} seeds x {len(THRESHOLDS)} thresholds "
          f"at m_T={a.mt}%", flush=True)

    # Each seed's records go straight to a spool file and are dropped from
    # memory. Holding them all until the pool drained cost ~1.2 GB on HIV_high,
    # as Python tuples at ~650 bytes per clade, and _encode then doubled it.
    # float64 on disk, not float32: _u16 truncates, so a rounding difference
    # would move a value across an integer boundary and change the output.
    spool = os.path.join(a.out, "_spool")
    os.makedirs(spool, exist_ok=True)
    nf = len(COLS)
    counts = {}                                     # (scen, k) -> {seed: nrows}
    seeds_meta, failed = [], []
    with Pool(a.procs) as pool:
        for i, (sc, sd, rows, srec, err) in enumerate(pool.imap_unordered(_one, jobs, chunksize=1), 1):
            if err:
                failed.append((sc, sd, err))
            else:
                seeds_meta.append(srec)
                for k, rec in rows:
                    if not rec:
                        continue
                    arr = np.asarray(rec, float)
                    assert arr.shape[1] == nf, f"{sc} {sd} k={k}: {arr.shape[1]} fields, want {nf}"
                    arr.tofile(os.path.join(spool, f"{sc}|{k}|{sd}.f64"))
                    counts.setdefault((sc, k), {})[sd] = int(arr.shape[0])
            rows = None
            if i % 10 == 0 or i == len(jobs):
                print(f"[fig2b]   {i}/{len(jobs)}", flush=True)

    seeds_meta.sort(key=lambda r: (r["scenario"], r["seed"]))
    seedidx = {}
    for sc in SCEN_ORDER:
        for j, r in enumerate([x for x in seeds_meta if x["scenario"] == sc]):
            seedidx[r["seed"]] = j

    total = 0
    files = {}
    for (sc, k), per_seed in sorted(counts.items()):
        # read back in seed order, so the row order matches what the in-memory
        # version produced and the encoded tables stay byte-identical
        order = sorted(per_seed)
        ntot = sum(per_seed[sd] for sd in order)
        a_blk = np.empty((ntot, nf), float)
        sid = np.empty(ntot, np.uint8)
        o = 0
        for sd in order:
            m = per_seed[sd]
            p = os.path.join(spool, f"{sc}|{k}|{sd}.f64")
            a_blk[o:o + m] = np.fromfile(p, float).reshape(m, nf)
            sid[o:o + m] = seedidx[sd]
            o += m
        obj = _encode(a_blk, sid)
        del a_blk, sid
        if obj is None:
            continue
        obj.update(scenario=sc, threshold=k, cols=COLS)
        ent = _write_table(a.out, sc, k, obj)
        total += ent["bytes"]
        files[f"{sc}|{k}"] = ent

    # merge with whatever a previous scenario run left here
    ipath = os.path.join(a.out, "index.json")
    if os.path.exists(ipath):
        with open(ipath) as f:
            old = json.load(f)
        files = {**old.get("files", {}), **files}
        # seed names repeat across scenarios, so the key must carry both
        keep = {(r["scenario"], r["seed"]): r for r in old.get("seeds", [])}
        keep.update({(r["scenario"], r["seed"]): r for r in seeds_meta})
        seeds_meta = sorted(keep.values(), key=lambda r: (r["scenario"], r["seed"]))

    meta = scenario_meta(a.exp_name)
    with open(ipath, "w") as f:
        json.dump({"exp_name": a.exp_name, "exp_num": a.exp_num,
                   "thresholds": THRESHOLDS, "cols": COLS, "labels": LAB, "never": NEVER,
                   "min_tx_frac": a.mt,
                   "window": {"mode": pre.BURNIN_WINDOW, "fixed_days": pre.BURNIN_CUTOFF_DAYS,
                              "min_final_time": pre.MIN_FINAL_TIME,
                              "min_window_days": pre.MIN_WINDOW_DAYS},
                   "scenarios": {s: {"ratio": meta.get(s, {}).get("ratio", 0.0),
                                     "duration": meta.get(s, {}).get("duration")}
                                 for s in SCEN_ORDER if s in present
                                 and any(k.startswith(s + "|") for k in files)},
                   "files": files, "seeds": seeds_meta}, f, separators=(",", ":"))
    shutil.rmtree(spool, ignore_errors=True)
    print(f"[fig2b] {total/1e6:.1f} MB in {len(files)} files, {len(seeds_meta)} seeds")
    for sc, sd, e in failed[:10]:
        print(f"[fig2b][warn] {sc}/{sd}: {e}")
