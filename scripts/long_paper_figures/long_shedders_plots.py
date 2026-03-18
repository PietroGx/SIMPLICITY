#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import simplicity.intra_host_model as ih
import simplicity.plots_manager as pm

PIE_COLORS = {"standard": "#4daf4a", "long": "#e41a1c", "mixed": "#377eb8", "founder": "#8172B2"}

def plot_fig1_intra_host(ax, tau_values, baseline_tau=7.5, time_grid=300):
    tvals = np.arange(0, float(time_grid), 1.0)
    color_map = {}
    
    host = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=float(baseline_tau), tau_4=8)
    y = host.data_plot_ih_solution(0, time_grid, 1)[0]
    line = ax.plot(tvals, y, lw=2.0, label=f"Baseline ({baseline_tau}d)")[0]
    color_map["Baseline"] = line.get_color()
    
    for tau in sorted(tau_values):
        host = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=float(tau), tau_4=8)
        y = host.data_plot_ih_solution(0, time_grid, 1)[0]
        line = ax.plot(tvals, y, lw=1.8, linestyle="--", label=f"Long ({float(tau):g}d)")[0]
        color_map[f"{float(tau):g}d"] = line.get_color()

    ax.set_xlabel("time [days]")
    ax.set_ylabel("P(infectious after t)")
    ax.set_xlim(0, time_grid)
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=9)
    return color_map

def plot_fig1_violins(ax, data_list, labels, color_map):
    parts = ax.violinplot(dataset=data_list, showmeans=False, showmedians=False, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        c = color_map.get(labels[i], "gray")
        body.set_facecolor(c)
        body.set_edgecolor(c)
        body.set_alpha(0.8)
        vals = data_list[i]
        if len(vals) >= 2:
            q1, q2, q3 = np.percentile(vals, [25, 50, 75])
            ax.scatter([i+1], [q2], marker="o", s=12, zorder=3, color='black')
            ax.vlines([i+1], q1, q3, lw=2, zorder=3, color='black')

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_ylabel("infection duration [days]")

def plot_fig2_timelines(ax, polished_data, colormap_df, t_final, only_long, label):
    jitter_step = 0.15
    data = polished_data[polished_data['type'] == 'long_shedder'] if only_long else polished_data
    data = data.sort_values(by='t_infection').reset_index(drop=True)

    for idx, row in data.iterrows():
        is_long = (row['type'] == 'long_shedder')
        traj = row['IH_lineages_trajectory']
        base_y = idx

        segments = [(lin, times['ih_birth'], times['ih_death']) for lin, times in traj.items() if isinstance(times, dict)]
        segments.sort(key=lambda x: x[1])
        overlaps = any(segments[i][1] < segments[i - 1][2] for i in range(1, len(segments)))

        for j, (lineage, ih_start, ih_end) in enumerate(segments):
            y_j = base_y + (j * jitter_step - jitter_step / 2) if overlaps else base_y
            color = pm.get_lineage_color(lineage, colormap_df)
            alpha = 1.0 if is_long else 0.3
            zorder = 5 if is_long else 1

            ax.hlines(y_j, ih_start, ih_end, color=color, linewidth=2, alpha=alpha, zorder=zorder)
            if j == 0 and is_long:
                ax.plot(ih_start, y_j, marker='o', color=color, markersize=4, zorder=zorder + 1)

    ax.set_ylabel("Infected number")
    ax.set_xlabel("Time (days)")
    ax.set_xlim(0, t_final)
    ax.set_ylim(0)
    ax.text(-0.05, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

def plot_fig3_lineage_freq(ax, ssod, label_char, freq_threshold=0.05):
    import simplicity.output_manager as om
    lf = om.read_lineage_frequency(ssod)
    t_final = om.read_final_time(ssod)
    cmap_df = pm.make_lineages_colormap(ssod)
    
    pivot = lf.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t")
    cols = [c for c in pivot.columns if (pivot[c].fillna(0) > 0).any()]
    pivot = pivot[cols]
    
    if freq_threshold is not None and freq_threshold > 0:
        max_by_col = pivot.fillna(0).max(axis=0)
        keep_cols = [c for c in pivot.columns if max_by_col.get(c, 0.0) >= freq_threshold]
        pivot = pivot[keep_cols]
        
    colors = [pm.get_lineage_color(c, cmap_df) for c in pivot.columns]
    pivot.plot(kind="area", stacked=False, color=colors, alpha=0.6, ax=ax, legend=False)
    
    ax.set_xlim(0, t_final)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Lineage freq.")
    ax.set_xlabel("Time (d)")
    ax.set_title(label_char, loc="left", pad=8, fontsize=16, fontweight="bold")
    return cmap_df, t_final

def plot_fig3_clustered_freq(ax, ssod, cluster_threshold, colormap_df, t_final, label_char):
    import simplicity.output_manager as om
    import simplicity.clustering as cl
    
    phylo_df = om.read_phylogenetic_data(ssod)
    lf = om.read_lineage_frequency(ssod)
    full_pivot = lf.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t")
    
    clade_to_lineages, _, per_clade_mut_df, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
        phylo_df, shared_mut_threshold=cluster_threshold
    )
    
    clade_series = []
    for clade, members in clade_to_lineages.items():
        cols = [c for c in members if c in full_pivot.columns]
        if cols: clade_series.append(full_pivot[cols].sum(axis=1).rename(clade))
        
    if not clade_series: return
    clade_freq = pd.concat(clade_series, axis=1)
    
    root_by_clade = clade_meta_df.set_index("clade")["root_lineage"].to_dict() if not clade_meta_df.empty else {}
    colors = [pm.get_lineage_color(root_by_clade.get(c), colormap_df) for c in clade_freq.columns]
    
    clade_freq.plot(kind="area", stacked=True, color=colors, alpha=0.6, ax=ax, legend=False)
    ax.set_xlim(0, t_final)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Clade freq.")
    ax.set_xlabel("Time (d)")
    ax.set_title(label_char, loc="left", pad=8, fontsize=16, fontweight="bold")

def plot_fig3_violins(ax, dfE, labelsE):
    if dfE.empty: return
    color_map = {}
    cycle = plt.rcParams.get("axes.prop_cycle").by_key()["color"]
    for i, lab in enumerate(labelsE): color_map[lab] = cycle[i % len(cycle)]
    
    sns.violinplot(data=dfE, x="group", y="dist", hue="group", order=labelsE, palette=color_map, cut=0, inner="quartile", ax=ax, legend=False)
    ax.set_xlabel("")
    ax.set_ylabel("Hamming distance")
    ax.set_title("E. Genetic Distance per Transmission Event", loc="left", pad=8, fontsize=16, fontweight="bold")

def plot_fig4_pies(axes_2x2, pie_data_dict):
    metrics = ["Peak", "Burden", "Survival", "Growth"]
    label_order = ["standard", "long", "mixed", "founder"]
    
    for idx, metric in enumerate(metrics):
        ax = axes_2x2[idx]
        counts = pie_data_dict[metric]
        
        if sum(counts.values()) == 0:
            ax.text(0.5, 0.5, "No Data", ha="center", va="center", fontsize=8)
            ax.axis("off")
            continue
            
        values = [counts.get(k, 0) for k in label_order]
        colors = [PIE_COLORS[k] for k in label_order]
        
        ax.pie(values, colors=colors, wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        ax.set_title(metric, fontsize=9)

def add_global_pie_legend(fig):
    handles = [Patch(facecolor=PIE_COLORS[k], edgecolor="white", label=k) for k in ["standard", "long", "mixed", "founder"]]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98), frameon=False, ncol=4)

def plot_fig5_stats_cell(axes_2x2, stats_df):
    if stats_df.empty:
        for ax in axes_2x2: ax.axis("off")
        axes_2x2[0].text(0.5, 0.5, "No Data", ha="center", va="center", fontsize=8)
        return
        
    ax_dur, ax_inf, ax_diag, ax_pie = axes_2x2
    
    # 1. Duration Histogram
    sns.histplot(stats_df['duration'], ax=ax_dur, color='gray', bins=10, stat="probability")
    ax_dur.set_title("Duration (d)", fontsize=9)
    ax_dur.set_ylabel("Probability", fontsize=8)
    ax_dur.set_xlabel("Final Time", fontsize=8)
    
    # 2. Total Infected (Twin Axes)
    sns.histplot(stats_df['inf_std'], ax=ax_inf, color=PIE_COLORS['standard'], alpha=0.5, bins=10, stat="probability")
    ax_inf.set_ylabel("Prob. Normal", color=PIE_COLORS['standard'], fontsize=8)
    ax_inf.tick_params(axis='y', labelcolor=PIE_COLORS['standard'])
    ax_inf.set_xlabel("Infected Count", fontsize=8)
    
    ax_inf_long = ax_inf.twinx()
    sns.histplot(stats_df['inf_long'], ax=ax_inf_long, color=PIE_COLORS['long'], alpha=0.5, bins=10, stat="probability")
    ax_inf_long.set_ylabel("Prob. Long", color=PIE_COLORS['long'], fontsize=8)
    ax_inf_long.tick_params(axis='y', labelcolor=PIE_COLORS['long'])
    ax_inf.set_title("Total Infected", fontsize=9)
    
    # Bring twinx to foreground
    ax_inf.set_zorder(1)
    ax_inf_long.set_zorder(2)
    ax_inf_long.patch.set_visible(False)
    
    # 3. Total Diagnosed (Twin Axes)
    sns.histplot(stats_df['diag_std'], ax=ax_diag, color=PIE_COLORS['standard'], alpha=0.5, bins=10, stat="probability")
    ax_diag.set_ylabel("Prob. Normal", color=PIE_COLORS['standard'], fontsize=8)
    ax_diag.tick_params(axis='y', labelcolor=PIE_COLORS['standard'])
    ax_diag.set_xlabel("Diagnosed Count", fontsize=8)
    
    ax_diag_long = ax_diag.twinx()
    sns.histplot(stats_df['diag_long'], ax=ax_diag_long, color=PIE_COLORS['long'], alpha=0.5, bins=10, stat="probability")
    ax_diag_long.set_ylabel("Prob. Long", color=PIE_COLORS['long'], fontsize=8)
    ax_diag_long.tick_params(axis='y', labelcolor=PIE_COLORS['long'])
    ax_diag.set_title("Total Diagnosed", fontsize=9)
    
    # Bring twinx to foreground
    ax_diag.set_zorder(1)
    ax_diag_long.set_zorder(2)
    ax_diag_long.patch.set_visible(False)
    
    # 4. Survival Pie Chart
    died = (stats_df['duration'] < 365).sum()
    survived = (stats_df['duration'] >= 365).sum()
    if died + survived > 0:
        ax_pie.pie([died, survived], labels=['Died (<365d)', 'Survived'], colors=['#cccccc', '#4daf4a'], autopct='%1.1f%%', textprops={'fontsize': 7})
    else:
        ax_pie.axis('off')
    ax_pie.set_title("Simulation Survival", fontsize=9)
    
    for ax in [ax_dur, ax_inf, ax_diag]:
        ax.tick_params(axis='x', labelsize=7)
        ax.tick_params(axis='y', labelsize=7)
    ax_inf_long.tick_params(axis='y', labelsize=7) 
    ax_diag_long.tick_params(axis='y', labelsize=7)