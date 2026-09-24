import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
