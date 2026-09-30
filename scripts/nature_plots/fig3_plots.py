import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

METRIC_COLS = ["Peak", "Burden", "Survival", "Growth"]


def format_clean_axis(ax, remove_ticks=True):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])


def _missing(ax, label):
    format_clean_axis(ax, remove_ticks=True)
    ax.text(0.5, 0.5, f"{label}\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)


def plot_fig3_metrics(axes, metrics_df, palette=None, scenario_order=None):
    """The four clade metrics, one small panel each, scenario on the x-axis.

    Replaces a PCA of the same four metrics. PC1/PC2 carried no units and no
    stated loadings, so the only readable statement was "control clusters,
    the rest spread" -- it could not say WHICH metric separated them. These
    show the metrics directly.
    """
    if metrics_df is None or metrics_df.empty:
        for ax in axes:
            _missing(ax, "Panel B")
        return

    scenarios = scenario_order or metrics_df["scenario"].unique().tolist()
    scenarios = [s for s in scenarios if s in set(metrics_df["scenario"])]
    palette = palette or {}

    for ax, metric in zip(axes, METRIC_COLS):
        for i, scenario in enumerate(scenarios):
            vals = metrics_df.loc[metrics_df["scenario"] == scenario, metric].dropna()
            if vals.empty:
                continue
            colour = palette.get(scenario, "#888888")
            jitter = np.random.default_rng(42).normal(i, 0.07, size=len(vals))
            ax.scatter(jitter, vals, s=5, alpha=0.45, color=colour,
                       edgecolors="none", zorder=3)
            ax.hlines(np.median(vals), i - 0.28, i + 0.28,
                      color=colour, lw=1.6, zorder=4)
        ax.set_xticks(range(len(scenarios)))
        ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=5.4)
        ax.set_title(metric, fontsize=6.5, pad=2)
        ax.tick_params(axis='y', labelsize=5.4)
        format_clean_axis(ax, remove_ticks=False)


def _density_shade(ax, pts, colour, zorder=1, levels=4):
    """Filled KDE contours for one point cloud; silently skipped if the cloud
    is degenerate (too few points, or all identical)."""
    if pts.shape[0] < 20:
        return
    try:
        from scipy.stats import gaussian_kde
        from matplotlib.colors import LinearSegmentedColormap
        jitter = np.random.default_rng(0).normal(0, 1e-9, size=pts.shape)
        kde = gaussian_kde((pts + jitter).T)
        xs = np.linspace(pts[:, 0].min(), pts[:, 0].max(), 80)
        ys = np.linspace(pts[:, 1].min(), pts[:, 1].max(), 80)
        X, Y = np.meshgrid(xs, ys)
        Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        cmap = LinearSegmentedColormap.from_list("d", [(1, 1, 1, 0), colour])
        ax.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.30, zorder=zorder)
    except Exception:
        return


def plot_fig3_sequence_pca(ax, snp_df, labels, palette=None):
    """
    One point per sequenced genome, colored by individual_type. Features are
    already 0/1 presence-absence on a comparable scale, so only centering is
    needed (no variance standardization).
    """
    if snp_df is None or snp_df.empty or snp_df.shape[0] < 3 or snp_df.shape[1] == 0:
        _missing(ax, "Panel C (scatter)")
        return

    X = snp_df.to_numpy(dtype=float)
    Xc = X - X.mean(axis=0)
    n_components = min(2, Xc.shape[0] - 1, Xc.shape[1])
    if n_components < 2:
        _missing(ax, "Panel C (scatter)")
        return
    pcs = PCA(n_components=2).fit_transform(Xc)

    default_palette = {"standard": "#3b528b", "long_shedder": "#e41a1c"}
    palette = palette or default_palette
    labels = pd.Series(labels).reset_index(drop=True)

    # Standard genomes are far more numerous but occupy FAR fewer distinct
    # points (a standard host carries 1-4 near-identical lineages, a long
    # shedder 5-15 diverse ones), so they stack invisibly under the long
    # shedders. Draw them last and report both counts: the contrast between
    # "many genomes, few positions" and "fewer genomes, many positions" is the
    # panel's actual result, not something to hide.
    order = [l for l in ("long_shedder", "standard") if l in set(labels)]
    order += [l for l in labels.unique() if l not in order]
    for z, label in enumerate(order):
        mask = (labels == label).to_numpy()
        n_total = int(mask.sum())
        pts = pcs[mask]
        n_unique = len(np.unique(pts, axis=0))
        colour = palette.get(label, "black")
        # Filled density contours: the panel's point is how much sequence space
        # each type OCCUPIES, and identical genomes stack invisibly in a scatter.
        _density_shade(ax, pts, colour, zorder=1 + z)
        ax.scatter(pts[:, 0], pts[:, 1], s=9, alpha=0.40, color=colour,
                   label=f"{label} (n={n_total}, {n_unique} distinct)",
                   edgecolors="none", zorder=4 + z)

    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    ax.legend(fontsize=5.6, frameon=False, loc="best")
    format_clean_axis(ax, remove_ticks=False)


def plot_fig3_consistency(ax, consistency_df):
    """Panel D: the three mean pairwise distances, one box per comparison.

    Ordering is the result. If long shedders were a DISPLACED population the
    between-group bar would be tallest; it is not -- long-long is, which means
    a broader cloud around the same centre.
    """
    if consistency_df is None or consistency_df.empty:
        _missing(ax, "Panel D")
        return

    order = ["standard-standard", "long-long", "long-standard"]
    order = [o for o in order if o in set(consistency_df["comparison"])]
    colours = {"standard-standard": "#3b528b",
               "long-long": "#e41a1c",
               "long-standard": "#7a7a7a"}

    data = [consistency_df.loc[consistency_df["comparison"] == o, "distance"].dropna().to_numpy()
            for o in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.55,
                    showfliers=False, patch_artist=True)
    for patch, o in zip(bp["boxes"], order):
        patch.set_facecolor(colours.get(o, "#999999"))
        patch.set_alpha(0.30)
        patch.set_linewidth(0.8)
    for part in ("medians", "whiskers", "caps"):
        for line in bp[part]:
            line.set_color("0.25")
            line.set_linewidth(0.9)

    rng = np.random.default_rng(42)
    for i, (o, vals) in enumerate(zip(order, data)):
        if len(vals) == 0:
            continue
        ax.scatter(rng.normal(i, 0.055, size=len(vals)), vals, s=6, alpha=0.5,
                   color=colours.get(o, "#999999"), edgecolors="none", zorder=4)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([o.replace("-", "\nvs ") for o in order], fontsize=5.8)
    ax.set_ylabel("Mean pairwise Hamming\ndistance (substitutions)", fontsize=6.5)
    format_clean_axis(ax, remove_ticks=False)


# =============================================================================
# Conversion efficiency panels
# =============================================================================
ORIGIN_COLORS = {"long": "#e41a1c", "standard": "#3b528b"}

# Panels B, C and D all show ONE quantity. Naming it identically in each title is
# what makes the figure readable at a glance; the subtitle carries the only thing
# that differs between them, which is how it was aggregated.
METRIC_NAME = "Substitutions that founded a dominant clade"
METRIC_NOTE = "dominant = reached >50% of circulating infections"


def _titled(ax, title, subtitle=None, size=6.4, pad=3):
    ax.set_title(title, fontsize=size, pad=pad + (7 if subtitle else 0))
    if subtitle:
        ax.text(0.5, 1.015, subtitle, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=4.9, color="0.45")


def plot_efficiency_bars(ax, eff, scenario_order=None):
    """Panel A: what share of a host type's substitutions found a clade that
    reached more than half of circulating infections, pooled over simulations.

    Scenarios on the x-axis, percentage on the y, the two host types joined so
    the direction of the gap reads first. SOT pointing the other way is the
    result.
    """
    if eff is None or eff.empty:
        _missing(ax, "Panel A")
        return
    order = [x for x in (scenario_order or eff.index) if x in eff.index]
    order = [x for x in order if np.isfinite(eff.loc[x, "eff_long"])]
    if not order:
        _missing(ax, "Panel A")
        return

    x = np.arange(len(order))
    for xi, sc in zip(x, order):
        lo, st = eff.loc[sc, "eff_long"] * 100, eff.loc[sc, "eff_std"] * 100
        if np.isfinite(st) and np.isfinite(lo):
            ax.plot([xi, xi], [st, lo], color="0.55", lw=1.2, zorder=2,
                    solid_capstyle="round")
            ax.scatter([xi], [st], s=36, color=ORIGIN_COLORS["standard"],
                       zorder=4, edgecolors="none")
        ax.scatter([xi], [lo], s=36, color=ORIGIN_COLORS["long"], zorder=4,
                   edgecolors="none")
        r = eff.loc[sc, "ratio"]
        if np.isfinite(r):
            ax.text(xi + 0.14, max(lo, st if np.isfinite(st) else 0), f"{r:.1f}x",
                    va="center", ha="left", fontsize=5.8, fontweight="bold",
                    color=ORIGIN_COLORS["long"] if r > 1 else "0.35")

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=5.8)
    ax.set_xlim(-0.5, len(order) - 0.1)
    top = max(max(eff.loc[v, "eff_long"] * 100,
                  (eff.loc[v, "eff_std"] * 100) if np.isfinite(eff.loc[v, "eff_std"]) else 0)
              for v in order)
    ax.set_ylim(0, top * 1.42)        # headroom for the legend, top right
    ax.set_ylabel("% substitutions", fontsize=6.5)
    _titled(ax, METRIC_NAME, "all simulations pooled")
    ax.tick_params(axis="y", labelsize=5.6)

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="none", color=ORIGIN_COLORS["long"],
               markersize=3.6, label="long shedders"),
        Line2D([0], [0], marker="o", ls="none", color=ORIGIN_COLORS["standard"],
               markersize=3.6, label="standard hosts")],
        fontsize=5.4, frameon=False, loc="upper right")
    format_clean_axis(ax, remove_ticks=False)


def plot_efficiency_vs_duration(ax, raw, meta, threshold=5, peak=0.50,
                                scenario_order=None):
    """Panel B: the same rate, once per simulation.

    Scenarios on the x-axis rather than duration -- with four scenarios, two of
    which share a duration, a numeric axis left most of its range empty and
    overplotted the pair. Duration is annotated under each label instead.
    """
    if raw is None or raw.empty:
        _missing(ax, "Panel B")
        return
    d = raw[raw["threshold"] == threshold]
    if d.empty:
        _missing(ax, "Panel B")
        return

    order = [x for x in (scenario_order or d["scenario"].unique())
             if x in set(d["scenario"])]
    rng = np.random.default_rng(3)
    drawn = False
    for i, scenario in enumerate(order):
        g = d[d["scenario"] == scenario]
        for key, colour, off in (("std", ORIGIN_COLORS["standard"], -0.17),
                                 ("long", ORIGIN_COLORS["long"], +0.17)):
            subs = g[f"subs_{key}"].to_numpy(float)
            succ = g[f"succ_{key}_{peak}"].to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                e = np.where(subs > 0, 100.0 * succ / subs, np.nan)
            e = e[np.isfinite(e)]
            if not len(e):
                continue
            ax.scatter(rng.normal(i + off, 0.045, size=len(e)), e, s=6.5,
                       alpha=0.45, color=colour, edgecolors="none", zorder=3)
            ax.hlines(np.median(e), i + off - 0.12, i + off + 0.12,
                      color=colour, lw=1.7, zorder=5)
            drawn = True
    if not drawn:
        _missing(ax, "Panel B")
        return

    labels = []
    for sc in order:
        dur = (meta.get(sc) or {}).get("duration")
        labels.append(f"{sc}\n{dur:.0f} d" if dur else sc)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontsize=5.4)
    ax.set_ylabel("% substitutions", fontsize=6.5)
    _titled(ax, METRIC_NAME, "one point per simulation")
    ax.set_ylim(0, None)          # a percentage: the axis starts at zero
    ax.margins(y=0.05)
    ax.tick_params(axis="y", labelsize=5.6)
    format_clean_axis(ax, remove_ticks=False)


# One title per subplot, no block title. Peak was dropped: its maximum
# saturates at 1.0 in 88-98% of simulations and ties are broken by dict order
# (see BACKLOG). B/C/D still use _peak_values, thresholded rather than ranked.
WIN_TITLES = {
    "Burden":   "Most infections caused",
    "Survival": "Longest time in circulation",
    "Growth":   "Spread fastest",
}
WIN_METRICS = ["Burden", "Survival", "Growth"]
WIN_YLABEL = "% of simulations where a long-shedder clade ranked first"


def plot_winrate(axes, summary, scenario_order=None, metrics=None):
    """Panel C: how often the top-scoring clade of a simulation came from a
    long shedder, against what chance alone would give.

    The measures rank clades WITHIN a simulation, so only a proportion across
    simulations is comparable -- absolute values are not, each being its own
    stochastic trajectory. The black bar is the expectation: the share of
    candidate clades that were long-origin. The gap is the result.
    """
    metrics = metrics or WIN_METRICS
    if summary is None or summary.empty:
        for ax in axes:
            _missing(ax, "Panel C")
        return

    for ax, metric in zip(axes, metrics):
        sub = summary[summary["metric"] == metric]
        order = [x for x in (scenario_order or sub["scenario"].unique())
                 if x in set(sub["scenario"])]
        if not order:
            _missing(ax, metric)
            continue
        sub = sub.set_index("scenario")
        x = np.arange(len(order))
        win = [sub.loc[v, "win"] for v in order]
        lo = [max(0.0, sub.loc[v, "win"] - sub.loc[v, "lo"]) for v in order]
        hi = [max(0.0, sub.loc[v, "hi"] - sub.loc[v, "win"]) for v in order]

        for i, v in enumerate(order):
            ax.hlines(sub.loc[v, "null"], i - 0.3, i + 0.3, color="0.1", lw=1.1,
                      zorder=3)
        ax.errorbar(x, win, yerr=[lo, hi], fmt="o", ms=2.8, lw=0,
                    elinewidth=0.9, capsize=1.5, color=ORIGIN_COLORS["long"],
                    zorder=5)

        ax.set_xticks(x)
        last = metric == metrics[-1]
        ax.set_xticklabels(order if last else [""] * len(order),
                           rotation=45 if last else 0, ha="right" if last else "center",
                           fontsize=5.4)
        ax.set_ylim(0, 1.10)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(["0", "50", "100"], fontsize=5.4)
        ax.set_title(WIN_TITLES.get(metric, metric), fontsize=6.2, pad=2.5)
        format_clean_axis(ax, remove_ticks=False)


def winrate_legend(fig, **kw):
    """Shared legend for the win-rate panels."""
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", ls="none", color=ORIGIN_COLORS["long"],
               markersize=3.0, label="observed (95% CI)"),
        Line2D([0], [0], color="0.1", lw=1.1,
               label="expected from the long-origin share of clades"),
    ]
    return fig.legend(handles=handles, frameon=False, fontsize=5.2, **kw)


def plot_efficiency_robustness(ax, rob, scenario_order=None, palette=None):
    """Panel D: does the ratio survive the two arbitrary analysis choices?

    One line per scenario against the clade threshold, with a band spanning the
    four incidence cutoffs. That replaces a strip plus four heat-maps: the
    earlier version encoded the same number three ways (position, colour, text)
    across six sub-axes, and none of them said what the number was. Here the
    reader sees directly that every scenario's band stays on one side of parity
    whatever either choice is set to.
    """
    if rob is None or rob.empty:
        _missing(ax, "Panel D")
        return None

    order = [x for x in (scenario_order or rob["scenario"].unique())
             if x in set(rob["scenario"])]
    palette = palette or {}
    taus = sorted(rob["threshold"].unique())
    ax.axhline(1.0, color="0.4", lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(0.995, 1.0, "parity", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=4.8, color="0.4")

    drawn = False
    for scenario in order:
        g = rob[rob["scenario"] == scenario]
        med, lo, hi, xs = [], [], [], []
        for t in taus:
            v = (g.loc[g.threshold == t, "ratio"]
                 .replace([np.inf, -np.inf], np.nan).dropna())
            v = v[v > 0]
            if v.empty:
                continue
            xs.append(t); med.append(float(v.median()))
            lo.append(float(v.min())); hi.append(float(v.max()))
        if not xs:
            continue
        colour = palette.get(scenario, "#888888")
        ax.fill_between(xs, lo, hi, color=colour, alpha=0.16, lw=0, zorder=3)
        ax.plot(xs, med, color=colour, lw=1.4, zorder=4)
        ax.scatter(xs, med, s=9, color=colour, zorder=5, edgecolors="none")
        ax.text(xs[-1] + 0.35, med[-1], scenario, fontsize=5.0, color=colour,
                va="center", ha="left")
        drawn = True
    if not drawn:
        _missing(ax, "Panel D")
        return None

    ax.set_yscale("log")
    ax.set_xticks(taus)
    ax.set_xticklabels([str(t) for t in taus], fontsize=5.4)
    ax.set_xlim(min(taus) - 0.6, max(taus) + 4.2)
    ax.set_xlabel("clade threshold (substitutions)", fontsize=6.0)
    ax.set_ylabel("ratio, long \u00f7 standard", fontsize=6.0)
    ax.tick_params(axis="y", labelsize=5.4)
    _titled(ax, METRIC_NAME,
            "ratio of that %, under every analysis choice "
            "(band = 10\u201375% incidence cutoffs)")
    format_clean_axis(ax, remove_ticks=False)
    return None
