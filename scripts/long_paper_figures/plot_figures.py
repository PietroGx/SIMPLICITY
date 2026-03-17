#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import random
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec

import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm

from long_shedders_preprocess import *
from long_shedders_plots import *

pm.apply_plos_rcparams()

def build_figure_1(exp_num, M, R, ratio):
    print(f"[Figure 1] Building with slice M={M}, R={R}, long_ratio={ratio}")
    df = read_master_log(exp_num)
    subset = df[(np.isclose(df['M'], M)) & (np.isclose(df['R_long'], R)) & (np.isclose(df['long_shedders_ratio'], ratio))]
    taus = sorted(subset['tau_3_long'].unique())
    
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))
    
    color_map = plot_fig1_intra_host(axA, taus)
    axA.set_title("A. Theoretical Intra-Host Dynamics", loc="left")

    labels = ["Baseline"] + [f"{float(t):g}d" for t in taus]
    data_list = []
    
    base_sod = get_baseline_sod(exp_num)
    base_ssod = dm.get_seeded_simulation_output_dirs(base_sod)[0]
    data_list.append(extract_infection_durations(base_ssod, "standard"))
    
    for t in taus:
        sod = get_grid_sod(M, ratio, t, R, exp_num)
        ssod = dm.get_seeded_simulation_output_dirs(sod)[0]
        data_list.append(extract_infection_durations(ssod, "long_shedder"))
        
    plot_fig1_violins(axB, data_list, labels, color_map)
    axB.set_title("B. Empirical Infection Durations", loc="left")
    
    outpath = os.path.join(dm.get_data_dir(), f"Figure_1_M{M}_R{R}_lsr{ratio}.png")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    print(f"[Saved] {outpath}")
    plt.close(fig)

def build_figure_2(exp_num, M, R, ratio, tau, seed):
    print(f"[Figure 2] Extracting timelines for M={M}, R={R}, long_ratio={ratio}, tau={tau}, seed={seed}")
    sod = get_grid_sod(M, ratio, tau, R, exp_num)
    ssod = dm.get_ssod(sod, seed)
    
    polished_df, t_final = analyze_long_shedders_trajectories(ssod)
    colormap_df = pm.make_lineages_colormap(ssod)
    
    fig, axA = plt.subplots(1, 1, figsize=(16, 8))
    plot_fig2_timelines(axA, polished_df, colormap_df, t_final, only_long=False, label="A")
    
    axB = axA.inset_axes([0.05, 0.5, 0.5, 0.5])
    axB.set_facecolor("none")
    plot_fig2_timelines(axB, polished_df, colormap_df, t_final, only_long=True, label="B")
    
    for s in ("top", "right"): axA.spines[s].set_visible(False)
    
    outpath = os.path.join(dm.get_data_dir(), f"Figure_2_M{M}_R{R}_lsr{ratio}_tau{tau}_seed{seed}.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"[Saved] {outpath}")
    plt.close(fig)

def build_figure_3(exp_num, M, R, ratio, tau, seed, cluster_threshold):
    print(f"[Figure 3] Extracting lineage frequencies and Hamming distances for seed={seed}...")
    
    base_sod = get_baseline_sod(exp_num)
    grid_sod = get_grid_sod(M, ratio, tau, R, exp_num)
    
    base_ssod = dm.get_ssod(base_sod, seed)
    grid_ssod = dm.get_ssod(grid_sod, seed)
    
    df = read_master_log(exp_num)
    subset = df[(np.isclose(df['M'], M)) & (np.isclose(df['R_long'], R)) & (np.isclose(df['long_shedders_ratio'], ratio))]
    taus = sorted(subset['tau_3_long'].unique())
    dfE = cache_and_aggregate_hamming(exp_num, M, R, ratio, taus)
    labelsE = ["Baseline (Standard)"] + [f"{float(t):g}d" for t in taus]
    
    fig = plt.figure(figsize=(16, 12))
    gs = plt.GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1.1], hspace=0.35, wspace=0.25)
    
    axA = fig.add_subplot(gs[0, 0])
    cmap1, tf1 = plot_fig3_lineage_freq(axA, base_ssod, "A. Baseline Lineages", 0.05)
    
    axB = fig.add_subplot(gs[1, 0])
    plot_fig3_clustered_freq(axB, base_ssod, cluster_threshold, cmap1, tf1, "B. Baseline Clades")
    
    axC = fig.add_subplot(gs[0, 1])
    cmap2, tf2 = plot_fig3_lineage_freq(axC, grid_ssod, f"C. Long Shedders (t={tau}) Lineages", 0.05)
    
    axD = fig.add_subplot(gs[1, 1])
    plot_fig3_clustered_freq(axD, grid_ssod, cluster_threshold, cmap2, tf2, "D. Long Shedders Clades")
    
    axE = fig.add_subplot(gs[2, :])
    plot_fig3_violins(axE, dfE, labelsE)
    
    outpath = os.path.join(dm.get_data_dir(), f"Figure_3_M{M}_R{R}_lsr{ratio}_tau{tau}_seed{seed}.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"[Saved] {outpath}")
    plt.close(fig)

def build_figure_4(exp_num, M, cluster_threshold, min_days):
    print(f"[Figure 4] Building master clade grid for M={M}")
    df = read_master_log(exp_num)
    df = df[np.isclose(df['M'], M)]
    
    taus = sorted(df["tau_3_long"].unique(), reverse=True)
    ratios = sorted(df["long_shedders_ratio"].unique(), reverse=True)
    rs = sorted(df["R_long"].unique(), reverse=True)
    
    fig = plt.figure(figsize=(5 * len(rs) * len(ratios), 5 * len(taus)))
    outer = fig.add_gridspec(len(taus), len(rs) * len(ratios), wspace=0.3, hspace=0.4)
    
    for i_tau, tau in enumerate(taus):
        col_idx = 0
        for ratio in ratios:
            for R in rs:
                cell_spec = outer[i_tau, col_idx]
                inner = GridSpecFromSubplotSpec(2, 2, subplot_spec=cell_spec, wspace=0.05, hspace=0.2)
                axes = [fig.add_subplot(inner[0, 0]), fig.add_subplot(inner[0, 1]), 
                        fig.add_subplot(inner[1, 0]), fig.add_subplot(inner[1, 1])]
                
                try:
                    sod = get_grid_sod(M, ratio, tau, R, exp_num)
                    results, valid_seeds = summarize_sod_pies(sod, cluster_threshold, min_days)
                    plot_fig4_pies(axes, results)
                    axes[0].set_title(f"t={tau} | R={R} | long_ratio={ratio}\n({valid_seeds} runs)", fontsize=10, loc="left", pad=15)
                except Exception as e:
                    for ax in axes: ax.axis("off")
                    axes[0].text(0.5, 0.5, f"Missing Data\n{e}", ha="center", va="center", fontsize=8)
                col_idx += 1

    add_global_pie_legend(fig)
    outpath = os.path.join(dm.get_data_dir(), f"Figure_4_Grid_M{M}.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"[Saved] {outpath}")
    plt.close(fig)

def build_figure_5(exp_num, M):
    print(f"[Figure 5] Building epidemiological stats grid for M={M}")
    df = read_master_log(exp_num)
    df = df[np.isclose(df['M'], M)]
    
    taus = sorted(df["tau_3_long"].unique(), reverse=True)
    ratios = sorted(df["long_shedders_ratio"].unique(), reverse=True)
    rs = sorted(df["R_long"].unique(), reverse=True)
    
    # Matching Figure 4 dimensions
    fig = plt.figure(figsize=(4.5 * len(rs) * len(ratios), 4.5 * len(taus)))
    outer = fig.add_gridspec(len(taus), len(rs) * len(ratios), wspace=0.35, hspace=0.45)
    
    for i_tau, tau in enumerate(taus):
        col_idx = 0
        for ratio in ratios:
            for R in rs:
                cell_spec = outer[i_tau, col_idx]
                # Changed to 2x2 to perfectly match Fig 4 layout visually
                inner = GridSpecFromSubplotSpec(2, 2, subplot_spec=cell_spec, wspace=0.45, hspace=0.45)
                axes = [fig.add_subplot(inner[0, 0]), fig.add_subplot(inner[0, 1]), 
                        fig.add_subplot(inner[1, 0]), fig.add_subplot(inner[1, 1])]
                
                try:
                    sod = get_grid_sod(M, ratio, tau, R, exp_num)
                    stats_df = extract_epi_stats(sod)
                    plot_fig5_stats_cell(axes, stats_df)
                    
                    axes[0].set_title(f"t={tau} | R={R} | long_ratio={ratio}", fontsize=10, loc="left", pad=15)
                except Exception as e:
                    for ax in axes: ax.axis("off")
                    axes[0].text(0.5, 0.5, f"Missing Data", ha="center", va="center", fontsize=8)
                col_idx += 1

    outpath = os.path.join(dm.get_data_dir(), f"Figure_5_StatsGrid_M{M}.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"[Saved] {outpath}")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Master Plotting Script for Long Shedders Grid")
    parser.add_argument("--figure", type=str, choices=["1", "2", "3", "4", "5", "all"], default="all")
    parser.add_argument("--exp-name", type=str, help="Experiment name")
    parser.add_argument("--exp-num", type=int, default=1, help="Experiment run iteration (default 1)")
    
    parser.add_argument("--M", type=float, default=2.0)
    parser.add_argument("--R", type=float, default=1.0)
    parser.add_argument("--ratio", type=float, default=0.01)
    
    parser.add_argument("--tau", type=float, default=90.0, help="Specific tau for Fig 2 & 3")
    parser.add_argument("--seed", type=int, default=1, help="Specific seed to plot (if not using random)")
    parser.add_argument("--random-seeds", type=int, default=0, help="Number of random seeds to plot for Figs 2 & 3")
    
    parser.add_argument("--cluster-threshold", type=int, default=5, help="Clustering threshold for Fig 3 & 4")
    parser.add_argument("--min-days", type=int, default=100, help="Min days for Fig 4 evaluation")

    args = parser.parse_args()
    
    if args.figure in ["1", "all"]: build_figure_1(args.exp_num, args.M, args.R, args.ratio)
    if args.figure in ["4", "all"]: build_figure_4(args.exp_num, args.M, args.cluster_threshold, args.min_days)
    if args.figure in ["5", "all"]: build_figure_5(args.exp_num, args.M)
    
    seeds_to_plot = [args.seed]
    experiment_numbered_name = f'{args.exp_name}_#{args.exp_num}'
    
    if args.random_seeds > 0 and args.figure in ["2", "3", "all"]:
        exp_name = f"long_shedders_exp_M{args.M}_lsr{args.ratio}_tau{args.tau}_R{args.R}_#{args.exp_num}"
        try:
            total_seeds = int(sm.get_n_seeds_from_experiment_settings(experiment_numbered_name))
        except Exception as e:
            print(f"[Warning] Could not read 'seeds' from settings: {e}. Falling back to default (100).")
            total_seeds = 100
        
        num_to_pick = min(args.random_seeds, total_seeds)
        seeds_to_plot = random.sample(range(1, total_seeds + 1), num_to_pick)
        
        print(f"\n[Randomizer] Read total seeds = {total_seeds} from settings.")
        print(f"[Randomizer] Selected {num_to_pick} random seeds: {seeds_to_plot}\n")

    for s in seeds_to_plot:
        if args.figure in ["2", "all"]: build_figure_2(args.exp_num, args.M, args.R, args.ratio, args.tau, s)
        if args.figure in ["3", "all"]: build_figure_3(args.exp_num, args.M, args.R, args.ratio, args.tau, s, args.cluster_threshold)

if __name__ == "__main__":
    main()