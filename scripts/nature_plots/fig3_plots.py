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


def plot_fig3_metrics_pca(ax, metrics_df, palette=None):
    """
    One point per seed, all scenarios, colored by scenario. Standardized
    before PCA -- Burden is a raw count that can be orders of magnitude
    larger than Peak's 0-1 frequency, so skipping standardization would let
    Burden's scale dominate the PCA for no biological reason.
    """
    if metrics_df.empty or len(metrics_df) < 2:
        _missing(ax, "Panel B")
        return

    X = metrics_df[METRIC_COLS].fillna(0.0).to_numpy(dtype=float)
    Xs = StandardScaler().fit_transform(X)
    pcs = PCA(n_components=2).fit_transform(Xs)

    scenarios = metrics_df["scenario"].unique().tolist()
    if palette is None:
        cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])
        palette = {s: cycle[i % len(cycle)] for i, s in enumerate(scenarios)}

    for scenario in scenarios:
        mask = (metrics_df["scenario"] == scenario).to_numpy()
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=14, alpha=0.75,
                  color=palette.get(scenario, "black"), label=scenario, edgecolors="none")

    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    ax.legend(fontsize=6, frameon=False, loc="best")
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

    for label in labels.unique():
        mask = (labels == label).to_numpy()
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=10, alpha=0.5,
                  color=palette.get(label, "black"), label=label, edgecolors="none")

    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    ax.legend(fontsize=6, frameon=False, loc="best")
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
