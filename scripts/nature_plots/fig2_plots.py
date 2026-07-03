import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import simplicity.plots_manager as pm

def format_clean_axis(ax, remove_ticks=True):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

def plot_fig2_lineage_freq(ax, lf, cmap_df, global_max_t, title, freq_threshold=0.05):
    """Row 1: Standard Lineage Frequencies locked to global time axis."""
    if lf is None or lf.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, "Data Missing", ha='center', va='center', transform=ax.transAxes)
        return

    pivot = lf.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t")
    cols = [c for c in pivot.columns if (pivot[c].fillna(0) > 0).any()]
    pivot = pivot[cols]
    
    if freq_threshold is not None and freq_threshold > 0:
        max_by_col = pivot.fillna(0).max(axis=0)
        keep_cols = [c for c in pivot.columns if max_by_col.get(c, 0.0) >= freq_threshold]
        pivot = pivot[keep_cols]
        
    colors = [pm.get_lineage_color(c, cmap_df) for c in pivot.columns]
    
    pivot.plot(kind="area", stacked=False, color=colors, alpha=0.6, ax=ax, legend=False)
    
    ax.set_xlim(0, global_max_t)
    ax.set_ylim(0, 1.0)
    ax.set_title(title, fontsize=8, fontweight='bold')
    format_clean_axis(ax, remove_ticks=False)

def plot_fig2_clustered_freq(ax, lf, clade_to_lineages, clade_meta_df, cmap_df, global_max_t):
    """Row 2: Clustered Lineage Frequencies locked to global time axis."""
    if lf is None or lf.empty or not clade_to_lineages:
        format_clean_axis(ax, remove_ticks=True)
        return

    full_pivot = lf.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t")
    
    clade_series = []
    for clade, members in clade_to_lineages.items():
        cols = [c for c in members if c in full_pivot.columns]
        if cols: 
            clade_series.append(full_pivot[cols].sum(axis=1).rename(clade))
        
    if not clade_series:
        format_clean_axis(ax, remove_ticks=True)
        return
        
    clade_freq = pd.concat(clade_series, axis=1)
    
    root_by_clade = clade_meta_df.set_index("clade")["root_lineage"].to_dict() if not clade_meta_df.empty else {}
    colors = [pm.get_lineage_color(root_by_clade.get(c), cmap_df) for c in clade_freq.columns]
    
    clade_freq.plot(kind="area", stacked=True, color=colors, alpha=0.8, ax=ax, legend=False)
    
    ax.set_xlim(0, global_max_t)
    ax.set_ylim(0, 1.0)
    format_clean_axis(ax, remove_ticks=False)

def plot_fig2_divergence_violin(ax, dists, color):
    """Row 3: Divergence gained during infection."""
    if not dists:
        format_clean_axis(ax, remove_ticks=True)
        return
        
    df = pd.DataFrame({'dist': dists, 'group': '1'})
    
    sns.violinplot(data=df, x='group', y='dist', color=color, cut=0, inner="quartile", ax=ax, linewidth=0.8, density_norm='width')
    
    ax.set_xlabel("")
    ax.set_xticks([])
    ax.set_ylim(bottom=0)
    format_clean_axis(ax, remove_ticks=False)
