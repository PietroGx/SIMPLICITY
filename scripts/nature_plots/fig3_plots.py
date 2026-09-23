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
        n_unique = len(np.unique(pcs[mask], axis=0))
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=11, alpha=0.45,
                   color=palette.get(label, "black"),
                   label=f"{label} (n={n_total}, {n_unique} distinct)",
                   edgecolors="none", zorder=3 + z)

    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    ax.legend(fontsize=5.6, frameon=False, loc="best")
    format_clean_axis(ax, remove_ticks=False)


def plot_fig3_consistency(ax, consistency_df):
    """
    Per-seed excess Hamming divergence (long vs. standard, minus standard's
    own within-cohort divergence). Tight and consistently positive = robust
    effect; scattered around zero = not much of one.
    """
    if consistency_df is None or consistency_df.empty:
        _missing(ax, "Panel C (consistency)")
        return

    y = consistency_df["excess_divergence"].to_numpy(dtype=float)
    jitter = np.random.RandomState(42).normal(loc=0, scale=0.04, size=len(y))

    ax.axhline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.boxplot(y, positions=[0], widths=0.5, showfliers=False,
              boxprops=dict(color="black"), medianprops=dict(color="black"))
    ax.scatter(jitter, y, s=16, alpha=0.7, color="#4daf4a", edgecolors="none", zorder=3)

    ax.set_xticks([])
    ax.set_ylabel("Excess divergence (Hamming)\nlong vs. standard", fontsize=7)
    format_clean_axis(ax, remove_ticks=False)
    ax.set_xticks([])
