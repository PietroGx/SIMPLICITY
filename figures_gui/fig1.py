"""Figure 1 for the studio, read back off the pipeline's own axes.

    python figures_gui/fig1.py build --stamp 2 [--exp-name NAME] --store <experiment store>

Runs figures_gui/prep/fig1_save.build_figure_1 with savefig suppressed and
records what each fig1_plots painter was handed and drew: the frames, the axis
limits, the violin outlines, the box statistics, the density curves. Nothing
here re-derives a clock, a fit or a density. The clock fits are recomputed with
fig1_plots._through_origin only to keep full precision, and the build stops if
they differ from the lines and annotations the pipeline drew.

Writes <store>/fig1/fig1.json and <store>/fig1/index.json; `load(store)` reads
them for /api/fig1, stale when the code or the real-world inputs differ from
the build's. Built locally only: the inputs include patient data that is not on
the cluster, and fig1/ is never packed (layout.SECTIONS).
"""
import argparse, contextlib, hashlib, io, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NP = os.path.join(HERE, "prep")
sys.path.insert(0, HERE)
from codehash import files_hash

SOURCES = [os.path.abspath(__file__)] + [os.path.join(NP, f) for f in (
    "fig1_save.py", "fig1_preprocess_data.py", "fig1_plots.py", "_scenarios.py")] + [
    os.path.join(ROOT, "simplicity", "tuning", "evolutionary_rate.py")]
REAL = "Data/RealWorldData/"
INPUTS = ("literature_long_shedders_data.csv", "patient_global_clock.csv",
          "patient_intrahost_clock.csv")


def inputs_sig():
    """{file: sha256 prefix, or None when absent}: a build without the patient
    files quietly loses D's and E's patient series, so the build records them."""
    out = {}
    for f in INPUTS:
        p = os.path.join(ROOT, REAL, f)
        try:
            with open(p, "rb") as fh:
                out[f] = hashlib.sha256(fh.read()).hexdigest()[:16]
        except OSError:
            out[f] = None
    return out


def code_hash():
    return files_hash(SOURCES)


def out_dir(store):
    return os.path.join(os.path.expanduser(store), "fig1")


def load(store):
    d = out_dir(store)
    fj, ij = os.path.join(d, "fig1.json"), os.path.join(d, "index.json")
    for p in (fj, ij):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"{p} -- build it with python figures_gui/fig1.py "
                                    f"build --stamp N --store {store}")
    with open(fj) as fh:
        out = json.load(fh)
    with open(ij) as fh:
        ix = json.load(fh)
    now, sig = code_hash(), inputs_sig()
    why = (["code"] if ix.get("code_hash") != now else []) + \
          (["inputs"] if ix.get("inputs") != sig else [])
    out["index"] = dict(ix, current_hash=now, stale=bool(why), stale_why=why,
                        inputs_missing=[f for f, v in sig.items() if v is None])
    return out


# ------------------------------------------------------------------ build

def _r(v, sig=6):
    return [float(f"{x:.{sig}g}") for x in v]


def _run_pipeline(stamp, exp_name):
    """build_figure_1 as fig1_save runs it, every painter call recorded."""
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, NP)
    import matplotlib.pyplot as plt
    import fig1_save
    import fig1_plots as plots
    painters = ("plot_fig1_intra_host", "plot_fig1_violins", "plot_infection_duration",
                "plot_clock_overlay", "plot_divergence_density", "add_figure_legend")
    calls, matched, saved = {}, [], {}

    def wrap(name, f):
        def g(*a, **k):
            if name == "plot_infection_duration":
                bxp = a[0].bxp
                def rec(*aa, **kk):
                    calls.setdefault("bxp", []).append((kk.get("bxpstats", aa[0] if aa else None),
                                                        kk.get("positions")))
                    return bxp(*aa, **kk)
                a[0].bxp = rec
            calls.setdefault(name, []).append((a, k))
            return f(*a, **k)
        return g

    def tm(*a, **k):
        out = saved["_time_match"](*a, **k)
        matched.append(out[1])
        return out

    for n in painters + ("_time_match",):
        saved[n] = getattr(plots, n)
        setattr(plots, n, tm if n == "_time_match" else wrap(n, saved[n]))
    keep = (fig1_save.figure_path, plt.savefig, plots.SHOW_FIT_STATS)
    fig1_save.figure_path = lambda *a, **k: "<nothing: savefig suppressed by figures_gui/fig1.py>"
    plt.savefig = lambda *a, **k: None
    plots.SHOW_FIT_STATS = True
    log = io.StringIO()
    tee = type("Tee", (), {"write": lambda s, t: (sys.__stdout__.write(t), log.write(t))[0],
                           "flush": lambda s: sys.__stdout__.flush()})()
    try:
        with contextlib.redirect_stdout(tee):
            fig1_save.build_figure_1(stamp, exp_name, "png")
    finally:
        for n, f in saved.items():
            setattr(plots, n, f)
        fig1_save.figure_path, plt.savefig, plots.SHOW_FIT_STATS = keep
    return plots, calls, matched, [l for l in log.getvalue().splitlines() if l.startswith("[fig1]")]


def _scenario_of(exp_name):
    import fig1_preprocess_data as preproc
    return {preproc.get_clinical_label(s): s for s in preproc.duration_groups(exp_name)}


def _durations(stamp, exp_name):
    """Each cohort's designed infection length, tau_1 + tau_2 + tau_3 (tau_3_long
    for long shedders) + tau_4, read from one simulation output dir per scenario
    as panel A's curves are: for a long shedder this is inf_duration_long."""
    import fig1_preprocess_data as preproc
    out = {}
    for sc in preproc.duration_groups(exp_name):
        sod = preproc.one_sod(exp_name, sc, stamp)
        if sod is None:
            continue
        get = lambda n: float(preproc.sm.get_parameter_value_from_simulation_output_dir(sod, n))
        out[preproc.get_clinical_label(sc)] = round(
            get("tau_1") + get("tau_2") + get("tau_3" if sc == "control" else "tau_3_long") + get("tau_4"), 2)
    return out


def _lims(ax):
    return dict(xlim=list(map(float, ax.get_xlim())), ylim=list(map(float, ax.get_ylim())),
                xlabel=ax.get_xlabel(), ylabel=ax.get_ylabel())


def _panel_a(calls):
    (ax, df), k = calls["plot_fig1_intra_host"][0]
    x_max = k.get("x_max", 800)
    curves = []
    if df.empty:
        return dict(curves=[], x_max=x_max, **_lims(ax))
    for c, g in df.groupby("cohort", sort=False):
        g = g[g["time"] <= x_max].iloc[::2]
        curves.append(dict(cohort=c, t=g["time"].astype(int).tolist(), p=_r(g["p_infectious"], 4)))
    return dict(curves=curves, x_max=x_max, **_lims(ax))


def _panel_b(calls, plots):
    import numpy as np
    (ax, df), _ = calls["plot_fig1_violins"][0]
    if df.empty:
        return dict(violins=[], **_lims(ax))
    shown = plots.violin_sample(df)   # the painter's own draw: same seed, same rows
    order = [t.get_text() for t in ax.get_xticklabels()]
    polys = [c.get_paths()[0].vertices for c in ax.collections]
    qs = sorted(((float(np.mean(l.get_xdata())), float(l.get_ydata()[0])) for l in ax.lines))
    out = []
    for i, c in enumerate(order):
        v = next(p for p in polys if abs(p[:, 0].mean() - i) < 0.45)
        q = [y for x, y in qs if abs(x - i) < 0.45]
        d = shown.loc[shown["cohort"] == c, "duration"]
        out.append(dict(cohort=c, n=int(len(d)), n_all=int((df["cohort"] == c).sum()),
                        mean=float(d.mean()), q=sorted(q),
                        dx=_r(v[:, 0] - i, 4), y=_r(v[:, 1], 5)))
    return dict(violins=out, **_lims(ax))


def _role(plots, colour):
    from matplotlib.colors import to_hex
    h = to_hex(colour).lower()
    for k, v in plots.SCEN_COLORS.items():
        if v.lower() == h:
            return k
    return "unmodelled"


def _panel_c(calls, plots):
    (ax, df), _ = calls["plot_infection_duration"][0]
    if df.empty:
        return dict(categories=[], n=0, **_lims(ax))
    order = [t.get_text() for t in ax.get_xticklabels()]
    stats = {int(round(p)): s for st, ps in calls.get("bxp", []) for s, p in zip(st, ps)}
    cats = []
    for i, c in enumerate(order):
        g = df[df["clinical_category"] == c]
        s = stats[i]
        cats.append(dict(cat=c, n=int(len(g)), mean=float(g["duration"].mean()),
                         role=_role(plots, plots.SCEN_COLORS.get(
                             plots.CATEGORY_TO_SCENARIO.get(c), plots.UNMODELLED_COLOR)),
                         box={k: float(s[k]) for k in ("whislo", "q1", "med", "q3", "whishi")},
                         pts=[[float(d), _role(plots, col)] for d, col in zip(g["duration"], g["pt_color"])]))
    return dict(categories=cats, n=int(len(df)), cutoff=plots.LONG_SHEDDER_CUTOFF_DAYS, **_lims(ax))


def _clock(call, plots):
    """One plot_clock_overlay call: its inputs, its fits, checked against what it drew."""
    import numpy as np
    from matplotlib.colors import to_hex
    (ax, df_real, df_model), k = call
    real_x, real_y, x_max = k["real_x"], k["real_y"], k.get("x_max", 400)
    entries, out = [], dict(title=k.get("title"), real_x=real_x, real_y=real_y,
                            data_label=k.get("data_label", "Patient data"), x_max=x_max,
                            y_max=k.get("y_max"), model=[], **_lims(ax))
    out["data"] = None
    r = df_real.dropna(subset=[real_x, real_y]) if not df_real.empty else df_real
    if len(r) > 1:
        xr = r[real_x].astype(float).values
        yr = r[real_y].astype(float).values
        xr = xr - np.min(xr)
        s, r2 = plots._through_origin(xr, yr)
        out["data"] = dict(n=int(len(xr)), slope=float(s), rate=float(s * 365.25), r2=float(r2),
                           x_offset=float(r[real_x].astype(float).min()), x=_r(xr, 6), y=_r(yr, 5))
        entries.append((out["data_label"], s * 365.25, r2, plots.DATA_COLOR, s, 0.0))
    cohorts = set(df_model["cohort"]) if not df_model.empty else set()
    for c in [c for c in ("Control", "SOT", "HIV", "Edge Case") if c in cohorts]:
        cdf = df_model[df_model["cohort"] == c]
        xm = cdf["Sequencing_time"].astype(float).values * 365.25
        ym = cdf["Distance_from_root"].astype(float).values
        _, r2 = plots._through_origin(xm, ym)
        per = []
        for _, g in cdf.groupby("seed"):
            si, _ = plots._through_origin(g["Sequencing_time"].astype(float).values * 365.25,
                                          g["Distance_from_root"].astype(float).values)
            if np.isfinite(si):
                per.append(float(si))
        slope = float(np.mean(per))
        sd = float(np.std(per, ddof=1)) if len(per) > 1 else 0.0
        t = plots._thin(cdf)
        out["model"].append(dict(cohort=c, n=int(len(cdf)), seeds=int(cdf["seed"].nunique()),
                                 fitted_seeds=len(per), slope=slope, sd=sd, rate=slope * 365.25,
                                 sd_rate=sd * 365.25, r2=float(r2), per_seed=_r(per, 5),
                                 thin=dict(x=_r(t["Sequencing_time"].astype(float) * 365.25, 5),
                                           y=_r(t["Distance_from_root"].astype(float), 5))))
        entries.append((f"Model, {c}", slope * 365.25, r2, plots.SCEN_COLORS[c], slope, sd))
    # the check: every slope drawn, every band, and the annotation text verbatim
    drawn = {to_hex(l.get_color()).lower(): l.get_ydata()[1] / l.get_xdata()[1] for l in ax.lines}
    bands = {to_hex(p.get_facecolor()[0]).lower(): p.get_paths()[0].vertices[:, 1].max() / x_max
             for p in ax.collections if hasattr(p, "get_facecolor") and len(p.get_paths())
             and p.get_paths()[0].vertices[:, 0].max() > 0 and not len(p.get_offsets()) > 1}
    text = "\n".join(t.get_text() for t in ax.texts if "s/s/y" in t.get_text())
    want = "\n".join(f"{n}: {rate:.5f} s/s/y ($R^2$={r2:.2f})" for n, rate, r2, *_ in entries)
    bad = []
    if text != want:
        bad.append(f"annotation:\n{text}\n!=\n{want}")
    for n, _, _, col, slope, sd in entries:
        got = drawn.get(col.lower())
        if got is None or not np.isclose(got, slope, rtol=1e-9, atol=0):
            bad.append(f"{n}: drawn slope {got} != {slope}")
        if sd > 0 and not np.isclose(bands.get(col.lower(), np.nan), slope + sd, rtol=1e-9):
            bad.append(f"{n}: drawn band top {bands.get(col.lower())} != {slope + sd}")
    if bad:
        raise RuntimeError(f"{out['title']}: the studio's fit is not the figure's\n" + "\n".join(bad))
    out["annotation"] = text
    return out


def _panel_e(calls, matched, plots):
    from matplotlib.colors import to_hex
    (ax, df_real, df_model), k = calls["plot_divergence_density"][0]
    data_hex, ctrl_hex = plots.DATA_COLOR.lower(), plots.SCEN_COLORS["Control"].lower()
    series = {}
    for l in ax.lines:
        h = to_hex(l.get_color()).lower()
        name = "data" if h == data_hex else "model" if h == ctrl_hex else None
        if name is None:
            continue
        kind = "hist" if l.get_drawstyle() == "steps-mid" else "kde"
        series.setdefault(name, {})[kind] = dict(x=_r(l.get_xdata(), 5), y=_r(l.get_ydata(), 5))
    rx = k.get("real_x", "days_since_ref")
    ry = k.get("real_y", "hamming_per_overlap")
    import numpy as np
    rng = np.random.default_rng(42)

    def pts(v, n=1500):
        """the series' own values, a random 1,500 where there are more, for the studio's dots"""
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        return _r(np.sort(rng.choice(v, n, replace=False) if len(v) > n else v), 5)
    if "data" in series:
        r = df_real.dropna(subset=[rx, ry])
        series["data"].update(n=int(len(r)), pts=pts(r[ry].astype(float)),
                              label=k.get("data_label", "Patient data"))
    if "model" in series:
        series["model"].update(n=int((df_model["cohort"] == "Control").sum()),
                               matched=len(matched[0]) if matched else None,
                               pts=pts(matched[0]) if matched else [], label="Model, Control")
    return dict(series=series, x_max=k.get("x_max"), title=k.get("title"),
                n_time_bins=10, **_lims(ax))


def build(stamp, exp_name, store):
    t0 = time.time()
    plots, calls, matched, log = _run_pipeline(stamp, exp_name)
    leg = calls["add_figure_legend"][0][1]
    out = dict(
        experiment=dict(name=exp_name, num=stamp),
        cohort_scenario=_scenario_of(exp_name),
        cohort_duration=_durations(stamp, exp_name),
        sources=dict(
            A="fig1_preprocess_data.get_panel_a_data: Host(tau_1, tau_2, tau_3 or tau_3_long, tau_4)"
              ".data_plot_ih_solution, taus read from each scenario's simulation output dir",
            B="fig1_preprocess_data.get_panel_b_data: individuals_data.csv, t_not_infected - t_infection, "
              "10 seeds drawn with random.Random(42) from control's, the same seeds in every scenario",
            C=REAL + "literature_long_shedders_data.csv",
            D_data=REAL + "patient_global_clock.csv (overlap_coverage >= 0.95)",
            D_model="sequencing_data_regression.csv of control, standard individuals, "
                    "10% of hosts per seed (get_model_global_clock)",
            F_data=REAL + "patient_intrahost_clock.csv",
            F_model="evolutionary_rate.extract_ih_regression_data on individuals_data.csv + "
                    "phylogenetic_data.csv, every seed, every long shedder (get_model_intrahost_clock)"),
        A=_panel_a(calls), B=_panel_b(calls, plots), C=_panel_c(calls, plots),
        D=_clock(calls["plot_clock_overlay"][0], plots),
        E=_panel_e(calls, matched, plots),
        F=_clock(calls["plot_clock_overlay"][1], plots),
        legend=dict(has_edge_case=bool(leg.get("has_edge_case", True))),
        log=log)
    d = out_dir(store)
    os.makedirs(d, exist_ok=True)
    body = json.dumps(out, separators=(",", ":"), allow_nan=False)
    with open(os.path.join(d, "fig1.json.tmp"), "w") as fh:
        fh.write(body)
    os.replace(os.path.join(d, "fig1.json.tmp"), os.path.join(d, "fig1.json"))
    ix = dict(code_hash=code_hash(), inputs=inputs_sig(), built=time.strftime("%Y-%m-%dT%H:%M:%S"),
              stamp=stamp, exp_name=exp_name, bytes=len(body),
              secs=round(time.time() - t0, 1), sources=[os.path.relpath(p, ROOT) for p in SOURCES])
    with open(os.path.join(d, "index.json"), "w") as fh:
        json.dump(ix, fh, indent=1)
    print(f"[fig1-studio] wrote {d}/fig1.json ({len(body) / 1024:.0f} KB) in {ix['secs']} s, "
          f"code {ix['code_hash']}")
    for p in ("D", "F"):
        for e in [x for x in [out[p]["data"]] + out[p]["model"] if x]:
            print(f"[fig1-studio] {p} {e.get('cohort', 'data')}: {e['rate']:.5f} s/s/y "
                  f"R2={e['r2']:.4f} n={e['n']}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--stamp", required=True)
    b.add_argument("--exp-name", default="impact_long_shedders")
    b.add_argument("--store", required=True, help="the experiment's store, sweep.store_dir")
    a = ap.parse_args()
    build(a.stamp, a.exp_name, a.store)
    return 0


if __name__ == "__main__":
    sys.exit(main())
