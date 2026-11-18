"""The clade pass for the figure-2 studio: clades, composition, burden and
episodes at one clustering threshold, at DAILY resolution.

One run of this writes one CLADE BUILD: a directory under `cache/`, named for
its parameters and a hash of the code that produced it.

Why this reads the raw run output rather than the clade-explorer export: that
export downsamples to 120 time points, about 9 days apart. Episode rise times
land at 30-39 days, so on that axis a whole rise is 3-4 samples and every
crossing day is quantised to ~9 days. Episodes need the daily record.

Clade formation is NOT reimplemented here: it calls
`simplicity.clustering.cluster_lin_into_clades_with_meta`, the same function the
pipeline and `tests/export_fig2_builder.py` use. The frame construction mirrors
`tests/export_fig2_builder.py:190-209`.

Episode detection mirrors `tests/artifacts/rise_lab.html` -- `smooth`,
`episodesOf`, `riseOf` -- whose own comment states it mirrors
`scipy.signal.find_peaks(ys, prominence=...)`, which is what is used here.
"""
import argparse, base64, hashlib, json, os, sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
import fig2_preprocess_data as pre
# the episode definition lives with the other clade metrics; `compute.smooth`
# and friends stay valid names here because they are imported, and panels.py
# reaches them that way
import clade_preprocess as cp
from clade_preprocess import smooth, episodes_of, rise_of
from _scenarios import EXP_NAME, resolve_scenarios, one_sod, scenario_meta

SCEN_ORDER = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]

DEFAULTS = dict(k=7, mt=1.0, win=15.0, prom=0.10, lo=0.10, hi=0.90)
# origin-label codes, same as tests/export_fig2_builder.py:39
LAB = {"standard": 0, "long": 1, "mixed": 2, "founder": 3}
WINDOWS = ("post", "all")   # "post" = from the burn-in day; "all" = whole record

# The m_T ladder, kept identical to the artifact's runs/formation.json so the
# across-the-set tab's axis is unchanged. It lives here rather than in
# panels.py because m_T became a BUILD parameter in v2.5.5: this is the grid's
# m_T axis, and the sweep driver and the panels must agree on it.
MT_LADDER = [0.0, 0.002, 0.003, 0.005, 0.008, 0.012, 0.02, 0.03, 0.05, 0.08,
             0.12, 0.2, 0.3, 0.5, 1.0, 2.0]


def code_hash():
    """Hash of this file, the clustering it calls, and the episode definition.

    A clade build is the output of all three, so all three belong in the key:
    this file, the clustering it calls, and the episode definition it measures
    with (`clade_preprocess`, which owns `smooth`/`episodes_of`/`rise_of` since
    v2.5.7). The reason is on the record: an uncommitted one-line change to the
    clade walk (clustering.py, v2.5.3) altered what a clade is while this file
    was untouched, so 24 builds carried the new maths under a key that claimed
    the old. A build now names the code that computed it.
    """
    h = hashlib.sha256()
    for f in (os.path.abspath(__file__), os.path.abspath(cl.__file__),
              os.path.abspath(cp.__file__)):
        with open(f, "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:12]


def params_key(p):
    """The build's directory name.

    m_T is in here because rule B designation uses it: a lineage founds a clade
    only if its subtree carries at least m_T of the run's transmissions, so
    moving m_T moves the PARTITION and cannot be applied over an existing
    build. It was a read parameter until v2.5.5, when designation changed from
    divergence-alone-then-drop to rule B.
    """
    s = (f"k{p['k']}_m{p.get('mt', 0.0):g}"
         f"_w{p['win']:g}_p{p['prom']:g}_lo{p['lo']:g}_hi{p['hi']:g}"
         f"_{p.get('window', 'post')}")
    return f"{s}__{code_hash()}"







def run_context(ssod, k, window, mt=0.0):
    """Everything one run needs before the per-clade loop, on the DAILY record.

    `_one` and `run_daily` both start here, because the episode numbers only
    agree while both measure on the same axis. That axis is the daily
    lineage_frequency record: the 120-point export in runs/ samples about every
    9.2 days, which quantises every crossing day and turns a 30-day rise into
    three samples (BACKLOG, "Sweep time 9.2-day gap").

    Clades are designated under RULE B: a lineage founds a clade when it has
    diverged >= k mutations from its clade root AND its subtree carries at
    least `mt` percent of the run's own transmission events. `mt` = 0 falls
    back to divergence alone.

    Rule B FOLDS a lineage that fails the support term into its surviving
    ancestor. The studio used to designate on divergence and then DROP clades
    under m_T, which discarded those lineages instead -- 94-98% of groups at
    the operating point -- so the two are different partitions, not the same
    one filtered.

    Returns (ctx, None), or (None, "stage: reason") for the caller to report.
    """
    k = int(k)
    try:
        ph = om.read_phylogenetic_data(ssod)
        lf = om.read_lineage_frequency(ssod)
        ind = om.read_individuals_data(ssod)
        phi = pre.phi_day(ssod)
        end = float(om.read_final_time(ssod))
    except Exception as exc:
        return None, f"read: {exc}"

    cutoff = (phi if phi is not None else float("inf")) \
        if window == "post" else float("-inf")

    # frame: times x clades, mirroring export_fig2_builder.py:190-209
    lin_code, lin_names = pd.factorize(lf["Lineage_name"], sort=False)
    t_code, t_vals = pd.factorize(lf["Time_sampling"], sort=True)
    freq = lf["Frequency_at_t"].to_numpy(float)
    n_t = len(t_vals)
    del lf
    try:
        support = min_support = None
        if float(mt) > 0:
            support = cl.subtree_transmission_counts(ph, ind)
            min_support = cl.transmission_floor(ind, float(mt))
        c2l, _l2c, pcm, meta = cl.cluster_lin_into_clades_with_meta(
            ph, shared_mut_threshold=k, support=support,
            min_support=min_support or 0.0)
    except Exception as exc:
        return None, f"cluster: {exc}"

    lin2host = dict(zip(ph["Lineage_name"], ph["Host_type"]))
    # the clade-explorer export numbers lineages by ph row order
    # (export_clade_explorer.py:63), so that order is how a daily
    # series is matched to the browser's own clade rows
    lin_pos = {n: i for i, n in enumerate(ph["Lineage_name"])}
    tx_all = ind["inherited_lineage"].value_counts().to_dict()
    n_tx_run = int(sum(tx_all.values()))

    clades = list(c2l)
    cpos = {c: j for j, c in enumerate(clades)}
    lin2cl = np.full(len(lin_names), -1, np.int64)
    for j, nm in enumerate(lin_names):
        c = _l2c.get(nm)
        if c is not None:
            lin2cl[j] = cpos[c]
    row_cl = lin2cl[lin_code]
    keep = row_cl >= 0
    if not keep.any():
        return None, "no lineage mapped to a clade"
    flat = t_code[keep] * len(clades) + row_cl[keep]
    M = np.bincount(flat, weights=freq[keep],
                    minlength=n_t * len(clades)).reshape(n_t, len(clades))
    seen = np.zeros(len(clades), bool)
    seen[row_cl[keep]] = True
    cols = [c for c, f in zip(clades, seen) if f]
    F = pd.DataFrame(M[:, seen], index=t_vals, columns=cols)

    md = meta.set_index("clade")
    par_of = md["parent_clade"].to_dict()
    ndef_of = md["n_defining"].to_dict()
    tem_of = md["start_time"].to_dict()
    root_of = md["root_lineage"].to_dict()

    # Origin by majority vote over the defining mutations. Definer-less
    # non-root clades are withheld from the vote exactly as
    # tests/export_fig2_builder.py:176-182 does: the function refuses them,
    # and they carry n_defining = 0 so no panel reads them anyway.
    lin2mut = cl.build_lineage_to_mutation_dict(ph)
    votable = {c: d for c, d in pcm.items()
               if not ((d is None or len(d) == 0) and isinstance(par_of.get(c), str))}
    try:
        labels_s, _ = cl.label_clades_from_definers(votable, meta)
        vote_of = labels_s.to_dict() if labels_s is not None else {}
    except Exception as exc:
        return None, f"label: {exc}"

    t_all = F.index.to_numpy(float)
    post = t_all >= cutoff
    t_post = t_all[post]
    return dict(
        F=F,
        cols=cols,
        t_all=t_all,
        t_post=t_post,
        post=post,
        cutoff=cutoff,
        phi=phi,
        end=end,
        ind=ind,
        n_tx_run=n_tx_run,
        c2l=c2l,
        pcm=pcm,
        par_of=par_of,
        ndef_of=ndef_of,
        tem_of=tem_of,
        root_of=root_of,
        lin2mut=lin2mut,
        lin2host=lin2host,
        lin_pos=lin_pos,
        tx_all=tx_all,
        vote_of=vote_of,
    ), None


MAX_DAILY_CLADES = 200


def ssod_for(scenario, seed, exp_num=1, exp_name=EXP_NAME):
    """The simulation output directory for one scenario and seed name."""
    for d in dm.get_seeded_simulation_output_dirs(
            one_sod(exp_name, scenario, exp_num)):
        if os.path.basename(d) == seed:
            return d
    return None


def run_daily(ssod, p, mt=None, max_clades=MAX_DAILY_CLADES):
    """Per-clade DAILY trajectories for one run, on the axis the panels use.

    The browser used to measure episodes on runs/, the clade-explorer export,
    which samples 120 time points about 9.2 days apart -- so a 30-day rise was
    three samples and every crossing day was quantised. This serves the same
    daily record `_one` measures on, which is what lets web/episodes.js and
    `_one` report one number instead of two.

    Each clade carries `root`, the index of its root lineage in ph row order,
    which is how the export numbers lineages -- so the browser can match a
    daily series to the clade row it already drew, rather than by name.

    Only the clades rule B keeps at `mt` are packed, because those are the ones
    the tab draws, and `max_clades` caps the payload at the largest of them.
    Plain JSON numbers, not a packed buffer: this is one run over loopback, and
    a second packing path is a second thing to keep in step. They are NOT
    rounded: an episode starts at the local minimum directly before its peak,
    and rounding to 5 decimals reshaped those minima -- measured over 133
    clades in 15 runs, one of them (HIV_high seed 0 Clade_0.366.5.11) moved
    from a 35-day rise to a 106-day one. The median was unaffected at 34.0 d,
    so this is a knife edge rather than a bias, and BACKLOG carries it as a
    question about the episode-start definition.

    Returns (payload, None) or (None, "stage: reason").
    """
    ctx, err = run_context(ssod, p["k"], p.get("window", "post"),
                           p.get("mt", 0.0))
    if err:
        return None, err

    F, cols, post, t_post = ctx["F"], ctx["cols"], ctx["post"], ctx["t_post"]
    tx_all, c2l = ctx["tx_all"], ctx["c2l"]
    lin_pos, root_of = ctx["lin_pos"], ctx["root_of"]
    recs = []
    for ci, c in enumerate(cols):
        mem = c2l.get(c, [])
        recs.append(dict(c=ci, name=str(c), col=c,
                         ndef=int(ctx["ndef_of"].get(c, 0) or 0),
                         ntx=int(sum(int(tx_all.get(m, 0)) for m in mem)),
                         founding=not isinstance(ctx["par_of"].get(c), str),
                         root=int(lin_pos.get(root_of.get(c), -1))))
    # rule B designation already applied the support term, to the founding
    # lineage's SUBTREE; a clade's own ntx counts only the lineages that stayed,
    # so re-applying the floor here would drop legitimate clades
    kept = cl.keep_designated_clades(recs)

    ys = {r["c"]: np.nan_to_num(F[r["col"]].to_numpy(float))[post] for r in kept}
    kept.sort(key=lambda r: -(float(ys[r["c"]].max()) if len(ys[r["c"]]) else 0.0))
    kept = kept[:max_clades]

    return dict(
        t=[round(float(v), 2) for v in t_post],
        phi=None if ctx["phi"] is None else round(ctx["phi"], 2),
        end=round(ctx["end"], 2),
        n_tx_run=ctx["n_tx_run"],
        n_groups=len(cols),
        n_kept=len(kept),
        capped=bool(len(kept) == max_clades),
        clades=[dict(c=r["c"], name=r["name"], ndef=r["ndef"], ntx=r["ntx"],
                     founding=bool(r["founding"]), root=r["root"],
                     y=[float(v) for v in ys[r["c"]]])
                for r in kept]), None


def _one(job):
    scenario, ssod, p = job
    seed = os.path.basename(ssod)
    sys.setrecursionlimit(300000)
    ctx, err = run_context(ssod, p["k"], p.get("window", "post"),
                           p.get("mt", 0.0))
    if err:
        return scenario, seed, None, err
    F = ctx['F']
    cols = ctx['cols']
    t_all = ctx['t_all']
    t_post = ctx['t_post']
    post = ctx['post']
    cutoff = ctx['cutoff']
    phi = ctx['phi']
    end = ctx['end']
    ind = ctx['ind']
    n_tx_run = ctx['n_tx_run']
    c2l = ctx['c2l']
    pcm = ctx['pcm']
    par_of = ctx['par_of']
    ndef_of = ctx['ndef_of']
    tem_of = ctx['tem_of']
    root_of = ctx['root_of']
    lin2mut = ctx['lin2mut']
    lin2host = ctx['lin2host']
    tx_all = ctx['tx_all']
    vote_of = ctx['vote_of']

    crec, erec = [], []
    for ci, c in enumerate(cols):
        nd = int(ndef_of.get(c, 0) or 0)
        parent = par_of.get(c)
        is_founding = not isinstance(parent, str)
        # composition: long-born share of the DEFINING mutations only
        # same expression as tests/export_fig2_builder.py:233
        d = pcm.get(c)
        nlong = int((d["host_type"].astype(str).str.lower()
                     .str.startswith("long")).sum()) if d is not None and len(d) else 0
        mem = c2l.get(c, [])
        ntx = int(sum(int(tx_all.get(m, 0)) for m in mem))
        rootlin = root_of.get(c)
        max_div = cl.clade_max_divergence(mem, rootlin, lin2mut)
        y_all = np.nan_to_num(F[c].to_numpy(float))
        burden_all = float(np.trapezoid(y_all, t_all)) if len(t_all) > 1 else 0.0
        y = y_all[post]
        burden = float(np.trapezoid(y, t_post)) if len(t_post) > 1 else 0.0
        tem = float(tem_of.get(c, 0) or 0)
        # truncated: already above a trace frequency at the burn-in day, so its
        # first appearance is outside the window and its latency is unmeasurable
        trunc = bool(phi is not None and tem < cutoff and len(y) and y[0] > 1e-4)

        crec.append(dict(c=ci, name=str(c), ndef=nd, nlong=nlong, ntx=ntx,
                         tem=tem, founding=is_founding, trunc=trunc,
                         burden=round(burden, 4), burden_all=round(burden_all, 4),
                         nlin=len(mem), max_div=int(max_div),
                         lab_vote=LAB.get(cl.resolve_clade_label(vote_of, c), 0),
                         lab_root=LAB["long"]
                                  if cl.normalise_host_type(
                                      lin2host.get(rootlin)) == "long"
                                  else LAB["standard"]))

        if len(t_post) < 3:
            continue
        ys = smooth(y, t_post, p["win"])
        h, eps = episodes_of(ys, p["prom"])
        if not eps:
            continue
        rs = [r for r in (rise_of(ys, t_post, a, b, p["lo"], p["hi"]) for a, b in eps)
              if r is not None]
        if not rs:
            continue
        # the emergence episode is the one containing the clade's own maximum
        imax = int(np.argmax(ys))
        emg = min(rs, key=lambda r: (not (r["a"] <= imax <= r["b"]), -r["amp"]))
        for r in rs:
            erec.append(dict(c=ci,
                             a_day=round(float(t_post[r["a"]]), 2),
                             pk_day=round(float(t_post[r["b"]]), 2),
                             lo_day=round(float(t_post[r["ia"]]), 2),
                             hi_day=round(float(t_post[r["ib"]]), 2),
                             days=round(r["days"], 2),
                             amp=round(r["amp"], 5),
                             top=round(r["top"], 5),
                             emg=bool(r is emg)))

    rec = dict(seed=seed, phi=None if phi is None else round(phi, 2),
               end=round(end, 2), n_tx_run=n_tx_run,
               n_infected=int(len(ind)),
               n_long_hosts=int(sum(
                   1 for h in ind["type"]
                   if cl.normalise_host_type(h) == "long")),
               max_div=max((c["max_div"] for c in crec), default=0),
               clades=crec, episodes=erec)
    return scenario, seed, rec, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, default=1)
    ap.add_argument("--exp-name", type=str, default=EXP_NAME)
    ap.add_argument("--out", required=True, help="cache root; a keyed subdir is created")
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--scenarios", type=str, default="")
    for name, val in DEFAULTS.items():
        ap.add_argument(f"--{name}", type=type(val), default=val)
    ap.add_argument("--window", choices=WINDOWS, default="post")
    a = ap.parse_args()

    p = {n: getattr(a, n) for n in DEFAULTS}
    p["window"] = a.window
    key = params_key(p)
    outdir = os.path.join(a.out, key)
    os.makedirs(outdir, exist_ok=True)

    want = set(x for x in a.scenarios.split(",") if x) or set(SCEN_ORDER)
    present, _ = resolve_scenarios(a.exp_name, a.exp_num, required=SCEN_ORDER)
    jobs = []
    for sc in [s for s in SCEN_ORDER if s in present and s in want]:
        for s in dm.get_seeded_simulation_output_dirs(one_sod(a.exp_name, sc, a.exp_num))[:a.seeds]:
            jobs.append((sc, s, p))
    print(f"[clades] {len(jobs)} runs at {key}", flush=True)

    per, failed = {}, []
    with Pool(a.procs) as pool:
        for i, (sc, sd, rec, err) in enumerate(pool.imap_unordered(_one, jobs, chunksize=1), 1):
            if err:
                failed.append((sc, sd, err))
            else:
                per.setdefault(sc, []).append(rec)
            if i % 10 == 0 or i == len(jobs):
                print(f"[clades] {i}/{len(jobs)}", flush=True)

    for sc in per:
        per[sc].sort(key=lambda r: r["seed"])
        with open(os.path.join(outdir, f"clades_{sc}.json"), "w") as fh:
            json.dump({"scenario": sc, "runs": per[sc]}, fh)

    try:
        meta = scenario_meta(a.exp_name)
    except Exception:
        meta = {}
    idx = dict(params=p, key=key, code=code_hash(),
               scenarios=[s for s in SCEN_ORDER if s in per],
               counts={s: len(v) for s, v in per.items()}, meta=meta,
               failed=[{"scenario": s, "seed": d, "why": w} for s, d, w in failed])
    with open(os.path.join(outdir, "index.json"), "w") as fh:
        json.dump(idx, fh, indent=1)

    print(f"\n[clades] {sum(len(v) for v in per.values())} runs -> {outdir}")
    for s in idx["scenarios"]:
        rs = per[s]
        nc = sum(len(r["clades"]) for r in rs)
        ne = sum(len(r["episodes"]) for r in rs)
        # clades as DESIGNATED: rule B already applied the support term during
        # designation, so this only drops definer-less non-root groups. It is
        # not the old "at m_T=1%" number, which filtered a divergence-only
        # partition and discarded the lineages it dropped.
        oc = oe = 0
        for r in rs:
            keep = {c["c"] for c in cl.keep_designated_clades(r["clades"])}
            oc += len(keep)
            oe += sum(1 for e in r["episodes"] if e["c"] in keep)
        print(f"  {s:12s} {len(rs):3d} runs  {nc:7,d} groups  {ne:7,d} episodes"
              f"   designated: {oc:5,d} clades  {oe:5,d} episodes"
              f"  {oe/max(oc,1):.2f} per clade")
    if failed:
        print(f"\n[clades] {len(failed)} FAILED")
        for s, d, w in failed[:10]:
            print(f"  {s} {d}: {w}")


if __name__ == "__main__":
    main()
