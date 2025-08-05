#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 10:55:24 2025

@author: pietro
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import argparse

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.settings_manager as sm
import simplicity.clustering as cl

# --- Config ---
CONFIG = {
    'TAKEOVER_THRESH': 0.5,
    'MIN_DAYS': 21,
    
    'RISE_LOW': 0.01,
    'RISE_HIGH': 0.5,
    'PIE_COLORS': {
        'normal': 'blue',
        'long_shedder': 'orange'
    },
    'SOD_COLORS' : ['blue', 'orange', 'green']

}


# --- Helpers ---
def count_takeovers(freq_series, threshold=0.1, min_days=21):
    freq_series = np.asarray(freq_series)
    over = freq_series >= threshold
    groups = np.split(over, np.where(~over)[0])
    return sum(group.sum() >= min_days for group in groups if group.size > 0)

def time_to_rise(time_series, freq_series, low=0.01, high=0.1):
    time_series, freq_series = pd.Series(time_series), pd.Series(freq_series)
    try:
        t_start = time_series[freq_series <= low].min()
        t_end = time_series[(freq_series >= high) & (time_series > t_start)].min()
        return t_end - t_start if pd.notnull(t_end) else None
    except Exception:
        return None

def analyze_ssod(ssod, use_clusters=True, cluster_threshold=5):
    # Read phylogenetic and frequency data
    phylo_df = om.read_phylogenetic_data(ssod)
    lineage_df = om.read_lineage_frequency(ssod)

    # Round time
    lineage_df['Time_sampling'] = lineage_df['Time_sampling'].round(2)

    # Build host type map
    host_type_series_raw = phylo_df.set_index("Lineage_name")["Host_type"]

    if use_clusters:
        # Pivot lineage frequencies
        pivot = lineage_df.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t")

        # Build lineage -> mutation dictionary
        lin2mut = cl.build_lineage_to_mutation_dict(phylo_df)

        # Cluster and aggregate frequencies
        freq_df_clusters, parents = cl.build_clustered_freqs(lin2mut, pivot, cluster_threshold)

        # Convert back to long format
        cluster_df = freq_df_clusters.reset_index().melt(id_vars="Time_sampling",
                                                         var_name="Lineage_name",
                                                         value_name="Frequency_at_t")
        lineage_df = cluster_df  # override with clustered frequencies

        # Assign host type from parent lineages
        host_type_series = cl.assign_cluster_hosttypes_from_parents(parents, freq_df_clusters.columns, host_type_series_raw)
    else:
        host_type_series = host_type_series_raw

    grouped = lineage_df.groupby('Lineage_name')

    takeover_counts = 0
    growth_meta = []
    freq_meta = []
    takeover_hosttypes = []

    for lineage, group in grouped:
        group = group.sort_values("Time_sampling")
        freqs = group["Frequency_at_t"].values
        times = group["Time_sampling"].values

        max_freq = freqs.max()
        freq_meta.append((lineage, max_freq))

        count = count_takeovers(freqs, threshold=CONFIG['TAKEOVER_THRESH'], min_days=CONFIG['MIN_DAYS'])
        if count > 0:
            takeover_counts += count
            host_type = host_type_series.get(lineage)
            if host_type is not None:
                takeover_hosttypes.append(host_type)

        ttr = time_to_rise(times, freqs, CONFIG['RISE_LOW'], CONFIG['RISE_HIGH'])
        if ttr is not None:
            host_type = host_type_series.get(lineage)
            if host_type is not None:
                growth_meta.append((lineage, ttr, host_type))

    top_growth = sorted(growth_meta, key=lambda x: x[1])[:3]
    enriched_freq_meta = []
    for lineage, maxf in sorted(freq_meta, key=lambda x: x[1], reverse=True)[1:4]:
        host_type = host_type_series.get(lineage)
        if host_type is not None:
            enriched_freq_meta.append((lineage, maxf, host_type))

    return {
        "takeover_count": takeover_counts,
        "growth": top_growth,
        "freq": enriched_freq_meta,
        "takeover_hosttypes": takeover_hosttypes
    }


def process_sod(sod, sod_name, min_final_time, use_clusters):
    ssods = dm.get_seeded_simulation_output_dirs(sod)
    long_ratio = sm.get_parameter_value_from_simulation_output_dir(sod, 'long_shedders_ratio')

    takeover_counts = []
    growth_bins = [[], [], []]
    freq_bins = [[], [], []]
    growth_pies = [Counter() for _ in range(3)]
    freq_pies = [Counter() for _ in range(3)]
    takeover_pie = Counter()

    for ssod in ssods:
        try:
            final_time = om.read_final_time(ssod)
        except Exception:
            continue  # skip if no final_time file

        if final_time < min_final_time:
            continue  # exclude short simulations

        result = analyze_ssod(ssod, use_clusters) 
        takeover_counts.append(result["takeover_count"])

        for i, (lin, val, ht) in enumerate(result["growth"]):
            growth_bins[i].append(val)
            growth_pies[i][ht] += 1

        for i, (lin, val, ht) in enumerate(result["freq"]):
            freq_bins[i].append(val)
            freq_pies[i][ht] += 1

        for ht in result["takeover_hosttypes"]:
            takeover_pie[ht] += 1

    return {
        "sod_name": sod_name,
        "takeovers": takeover_counts,
        "growth_bins": growth_bins,
        "freq_bins": freq_bins,
        "growth_pies": growth_pies,
        "freq_pies": freq_pies,
        "takeover_pie": takeover_pie,
        "long_shedder_ratio": long_ratio
    }


# --- Plotting Helpers ---
def _plot_box(ax, data, title, ylabel, labels=None):
    ax.boxplot(data, patch_artist=True, labels=labels)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)

def _plot_pie(ax, pie_data, title):
    keys = list(pie_data.keys())
    values = list(pie_data.values())
    colors = [CONFIG['PIE_COLORS'].get(k, None) for k in keys]

    ax.pie(values, labels=keys, autopct='%1.0f%%',
           textprops={'fontsize': 8}, colors=colors)
    ax.set_title(title, fontsize=9)

def plot_comparative_sod_boxplots(sod_datas, experiment_name):
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Comparative Boxplots Across SODs", fontsize=12)

    sod_labels = [d['sod_name'] for d in sod_datas]
    num_sods = len(sod_datas)

    from matplotlib.patches import Patch
    legend_labels = [
        f"long_shedder_ratio = {d['long_shedder_ratio']}" for d in sod_datas
    ]
    legend_patches = [
        Patch(facecolor=CONFIG['SOD_COLORS'][i], edgecolor='black', label=legend_labels[i])
        for i in range(len(sod_datas))
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=len(sod_datas), fontsize=9)

    # Panel 1
    box = axs[0].boxplot([d["takeovers"] for d in sod_datas], patch_artist=True)
    for patch, color in zip(box['boxes'], CONFIG['SOD_COLORS'][:num_sods]):
        patch.set_facecolor(color)
    axs[0].set_title("Takeover Events", fontsize=10)
    axs[0].set_ylabel("Count", fontsize=9)
    axs[0].set_xticklabels(sod_labels)
    axs[0].tick_params(labelsize=8)

    # Panel 2
    labels = ["Fastest", "2nd Fastest", "3rd Fastest"]
    growth_data = [sod["growth_bins"][rank] for rank in range(3) for sod in sod_datas]
    box = axs[1].boxplot(growth_data, patch_artist=True)
    for i, patch in enumerate(box['boxes']):
        patch.set_facecolor(CONFIG['SOD_COLORS'][i % num_sods])
    axs[1].set_title("Top Growth Speeds", fontsize=10)
    axs[1].set_ylabel("Days", fontsize=9)
    axs[1].set_xticks([1.5 + i * num_sods for i in range(3)])
    axs[1].set_xticklabels(labels)
    axs[1].tick_params(labelsize=8)

    # Panel 3
    freq_data = [sod["freq_bins"][rank] for rank in range(3) for sod in sod_datas]
    box = axs[2].boxplot(freq_data, patch_artist=True)
    for i, patch in enumerate(box['boxes']):
        patch.set_facecolor(CONFIG['SOD_COLORS'][i % num_sods])
    axs[2].set_title("Top Max Frequencies", fontsize=10)
    axs[2].set_ylabel("Relative Freq", fontsize=9)
    axs[2].set_xticks([1.5 + i * num_sods for i in range(3)])
    axs[2].set_xticklabels(labels)
    axs[2].tick_params(labelsize=8)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plot_dir = dm.get_experiment_plots_dir(experiment_name)
    output_path = os.path.join(plot_dir, f"comparative_boxplots_{experiment_name}.tiff")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    print(f"Saved comparative boxplots to: {output_path}")
    plt.close(fig)


def plot_sod_pies(sod_data, experiment_name):
    import os
    sod_name = sod_data["sod_name"]
    fig, axs = plt.subplots(3, 3, figsize=(18, 12))
    fig.suptitle(f"Host Type Composition for {sod_name}", fontsize=12)

    _plot_pie(axs[0, 0], sod_data["takeover_pie"], "Host Type (Takeovers)")
    for j in [1, 2]:
        axs[0, j].axis("off")

    labels = ["Fastest", "2nd Fastest", "3rd Fastest"]
    for i in range(3):
        _plot_pie(axs[1, i], sod_data["growth_pies"][i], f"Host Type ({labels[i]} Growth)")
    for i in range(3):
        _plot_pie(axs[2, i], sod_data["freq_pies"][i], f"Host Type ({labels[i]} Freq)")

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plot_dir = dm.get_experiment_plots_dir(experiment_name)
    output_path = os.path.join(plot_dir, f"{sod_name}_pies_{experiment_name}.tiff")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    print(f"Saved pie charts for {sod_name} to: {output_path}")
    plt.close(fig)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('min_final_time', type=int, help="min_final_time")
    parser.add_argument('--use_clusters', action='store_true', help="Enable lineage clustering")
    parser.add_argument('--no_clusters', dest='use_clusters', action='store_false', help="Disable lineage clustering")
    parser.set_defaults(use_clusters=True)


    args = parser.parse_args()
    sod_dirs = dm.get_simulation_output_dirs(args.experiment_name)

    if len(sod_dirs) != 2:
        raise ValueError(f"Expected exactly 2 SODs for comparative plots, got {len(sod_dirs)}.")

    sod_datas = [
        process_sod(sod_dirs[0], "SOD 1", args.min_final_time, args.use_clusters),
        process_sod(sod_dirs[1], "SOD 2", args.min_final_time, args.use_clusters)
    ]

    plot_comparative_sod_boxplots(sod_datas, args.experiment_name)
    plot_sod_pies(sod_datas[0],args.experiment_name)
    plot_sod_pies(sod_datas[1],args.experiment_name)

