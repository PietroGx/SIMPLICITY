"""The Figure 1 tab's stored JSON must be what fig1_save draws and prints.

Runs scripts/nature_plots/fig1_save.build_figure_1 afresh (savefig suppressed,
nothing written), takes its six axes off the figure -- not through
tests/studio/fig1.py's painter hooks -- and compares them with
<store>/fig1/fig1.json: the D and F annotations verbatim, every drawn slope and
band, the printed point and seed counts, the violin quartiles, the box
statistics, and the A and E curves. Reads only.

    python tests/studio/test_fig1_mirrors_pipeline.py --exp-num 2 --store <experiment store>
"""
import argparse, contextlib, io, os, re, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cbook
from matplotlib.colors import to_hex

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "nature_plots"))
import fig1
import fig1_plots as plots
import fig1_preprocess_data as preproc
import fig1_save

bad = 0


def ck(ok, msg):
    global bad
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {msg}")


def close(a, b, rtol):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=0)


def pipeline(exp_num, exp_name):
    keep = (fig1_save.figure_path, plt.savefig)
    fig1_save.figure_path = lambda *a, **k: "<not written>"
    plt.savefig = lambda *a, **k: None
    out = io.StringIO()
    try:
        with contextlib.redirect_stdout(out):
            fig1_save.build_figure_1(exp_num, exp_name, "png")
    finally:
        fig1_save.figure_path, plt.savefig = keep
    return plt.gcf().axes, out.getvalue()


def clock(name, ax, P, printed):
    text = "\n".join(t.get_text() for t in ax.texts if "s/s/y" in t.get_text())
    ck(text == P["annotation"], f"{name}: annotation drawn = annotation stored")
    rows = [(P["data_label"], P["data"])] + [(f"Model, {m['cohort']}", m) for m in P["model"]]
    want = "\n".join(f"{n}: {e['rate']:.5f} s/s/y ($R^2$={e['r2']:.2f})" for n, e in rows)
    ck(text == want, f"{name}: stored rates and R² print as the figure's:\n      " + want.replace("\n", "\n      "))
    drawn = {to_hex(l.get_color()).lower(): l.get_ydata()[1] / l.get_xdata()[1] for l in ax.lines}
    for n, e in rows:
        col = (plots.DATA_COLOR if e is P["data"] else plots.SCEN_COLORS[e["cohort"]]).lower()
        ck(np.isclose(drawn.get(col, np.nan), e["slope"], rtol=1e-12, atol=0),
           f"{name}: {n} slope drawn {drawn.get(col, np.nan):.6e} = stored {e['slope']:.6e}")
    pts = [c for c in ax.collections if len(c.get_offsets()) > 1]
    ck(len(pts) == 1 and len(pts[0].get_offsets()) == P["data"]["n"],
       f"{name}: {P['data']['n']} patient points stored, {len(pts[0].get_offsets()) if pts else 0} drawn")
    for m in P["model"]:
        n = printed.get(m["cohort"])
        ck(n == m["n"], f"{name}: {m['cohort']} {m['n']} points stored, {n} printed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-num", type=int, required=True)
    ap.add_argument("--exp-name", default="impact_long_shedders")
    ap.add_argument("--store", required=True)
    a = ap.parse_args()
    d = fig1.load(a.store)
    ck(not d["index"]["stale"], f"code hash {d['index']['code_hash']} is current ({d['index']['current_hash']})")
    ck(d["experiment"] == dict(name=a.exp_name, num=a.exp_num), f"stored experiment {d['experiment']}")

    axs, log = pipeline(a.exp_num, a.exp_name)
    ck(len(axs) == 6, f"{len(axs)} axes")
    aA, aB, aC, aD, aE, aF = axs[:6]

    m = re.search(r"global clock: (\d+) points from \S+ \((\d+) seeds", log)
    clock("D", aD, d["D"], {"Control": int(m.group(1))} if m else {})
    ck(m and int(m.group(2)) == d["D"]["model"][0]["seeds"], f"D: {d['D']['model'][0]['seeds']} seeds stored, printed {m and m.group(2)}")
    clock("F", aF, d["F"], {c: int(n) for n, c in re.findall(r"intra-host clock: (\d+) points for ([^(]+?) \(", log)})
    m = re.search(r"global clock cohort: (\d+) of", log)
    ck(m and int(m.group(1)) == d["D"]["data"]["n"], f"D: patient sequences printed {m and m.group(1)}, stored {d['D']['data']['n']}")

    # A: the curves seaborn drew, every second day as stored
    for c in d["A"]["curves"]:
        col = plots.SCEN_COLORS[c["cohort"]].lower()
        ln = [l for l in aA.lines if to_hex(l.get_color()).lower() == col]
        ok = bool(ln)
        if ok:
            x, y = ln[0].get_xdata(), ln[0].get_ydata()
            ok = close(np.interp(c["t"], x, y), c["p"], 1e-3)
        ck(ok, f"A: {c['cohort']} curve as drawn")
    ck(list(aA.get_xlim()) == d["A"]["xlim"], f"A: x limit {d['A']['xlim']}")

    # B: quartile lines and counts
    qs = sorted((float(np.mean(l.get_xdata())), float(l.get_ydata()[0])) for l in aB.lines)
    dfb = preproc.get_panel_b_data(exp_num=a.exp_num, exp_name=a.exp_name)
    for i, v in enumerate(d["B"]["violins"]):
        q = sorted(y for x, y in qs if abs(x - i) < 0.45)
        ck(close(q, v["q"], 1e-12), f"B: {v['cohort']} quartiles {[round(x, 1) for x in v['q']]}")
        n = int((dfb["cohort"] == v["cohort"]).sum())
        ck(n == v["n"], f"B: {v['cohort']} n stored {v['n']}, get_panel_b_data {n}")

    # C: box statistics from matplotlib's own rule, over the figure's categories
    dfc = preproc.get_panel_c_data()
    plots.plot_infection_duration(plt.figure().add_subplot(), dfc)
    order = [t.get_text() for t in aC.get_xticklabels()]
    ck(order == [k["cat"] for k in d["C"]["categories"]], f"C: category order {order}")
    m = re.search(r"panel C: (\d+) patients", log)
    ck(m and int(m.group(1)) == d["C"]["n"], f"C: {d['C']['n']} patients stored, printed {m and m.group(1)}")
    for k in d["C"]["categories"]:
        v = dfc.loc[dfc["clinical_category"] == k["cat"], "duration"].values
        s = cbook.boxplot_stats(v, whis=1.5)[0]
        ck(len(v) == k["n"] and all(np.isclose(s[f], k["box"][f]) for f in k["box"]),
           f"C: {k['cat']} n={k['n']} median {k['box']['med']:g}")

    # E: the density curves and histograms
    S = d["E"]["series"]
    for name, col in (("data", plots.DATA_COLOR), ("model", plots.SCEN_COLORS["Control"])):
        for kind in ("hist", "kde"):
            ln = [l for l in aE.lines if to_hex(l.get_color()).lower() == col.lower()
                  and (l.get_drawstyle() == "steps-mid") == (kind == "hist")]
            ok = bool(ln) and kind in S.get(name, {})
            if ok:
                ok = close(ln[0].get_xdata(), S[name][kind]["x"], 1e-4) and close(ln[0].get_ydata(), S[name][kind]["y"], 1e-4)
            ck(ok, f"E: {name} {kind} as drawn ({len(S.get(name, {}).get(kind, {}).get('x', []))} points)")

    print(f"\n{'PASS' if not bad else f'{bad} FAILED'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
