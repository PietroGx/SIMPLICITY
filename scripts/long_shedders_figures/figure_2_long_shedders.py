#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:48:46 2025

@author: pietro
"""
import scripts.long_shedders_figures.preprocess_data as preprocess
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import argparse
import os
import numpy as np
import simplicity.intra_host_model as ih

import matplotlib.pyplot as plt
import seaborn as sns


def plot_segmented_infection_timeline_ax(ax, polished_data, colormap_df, t_final, only_long, label):
    """
    Plot IH lineage segments colored by lineage per individual.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Subplot axis to draw on.
    polished_data : pd.DataFrame
        Must contain 'IH_lineages_trajectory', 'type', and be preprocessed.
    colormap_df : pd.DataFrame
        Output from make_lineages_colormap(ssod)
    t_final : float
        Final time for vertical line marker.
    only_long : bool
        If True, plot only long shedders. If False, plot all, with normal faded.
    label : str
        Subplot label (e.g., "B")
    """
    jitter_step = 0.15

    # Pre-filter dataset for plotting
    if only_long:
        data = polished_data[polished_data['type'] == 'long_shedder']
    else:
        data = polished_data

    # Sort by infection time
    data = data.sort_values(by='t_infection').reset_index(drop=True)

    for idx, row in data.iterrows():
        is_long = row['type'] == 'long_shedder'
        traj = row['IH_lineages_trajectory']
        base_y = idx
        
        segments = [(lineage, times['ih_birth'], times['ih_death']) for lineage, times in traj.items()]
        segments.sort(key=lambda x: x[1])
        overlaps = any(segments[i][1] < segments[i - 1][2] for i in range(1, len(segments)))

        for j, (lineage, ih_start, ih_end) in enumerate(segments):
            y_jittered = base_y + (j * jitter_step - jitter_step / 2) if overlaps else base_y
            color = pm.get_lineage_color(lineage, colormap_df)
            alpha = 1.0 if is_long else 0.3
            zorder = 5 if is_long else 1

            ax.hlines(y_jittered, ih_start, ih_end, color=color, linewidth=2, alpha=alpha, zorder=zorder)
            if j == 0 and is_long:
                ax.plot(ih_start, y_jittered, marker='o', color=color, markersize=4, zorder=zorder + 1)

    ax.axvline(t_final, color='gray', linestyle='--', linewidth=1.5)
    ax.set_ylabel("Individuals")
    ax.set_xlabel("Time")
    ax.set_xlim(0)
    ax.set_ylim(0)

    ax.grid(False)

    ax.text(-0.05, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')


def plot_figure_2(ssod):
    polished_df, colormap_df, t_final = preprocess.prepare_figure_2_data(ssod)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    label = 'A'
    only_long=False
    plot_segmented_infection_timeline_ax(axes[0], polished_df, colormap_df, t_final, only_long, label)
    label = 'B'
    only_long=True
    plot_segmented_infection_timeline_ax(axes[1], polished_df, colormap_df, t_final, only_long, label)


    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()

    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    sim_out_name = dm.get_simulation_output_foldername_from_SSOD(ssod)
    seed = os.path.basename(ssod)
    output_path = os.path.join(experiment_plots_dir, f"Figure_2_{sim_out_name}_{seed}.tiff")

    print(f"Saving figure to: {output_path}")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # # Set up the argument parser
    # parser = argparse.ArgumentParser(description="Plot")
    # parser.add_argument('experiment_name', type=str, help="experiment name")
    # args = parser.parse_args()
    # # Run the script with the provided parameter
    
    experiment_name = 'test_long_shedders_r1_kv_#1'
    seed_number =  1
    
    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)
    
    for sim_out_dir in sim_out_dirs:
        
        ssod = dm.get_ssod(sim_out_dir, seed_number)
    
        plot_figure_2(ssod)
    
    
if __name__ == "__main__":
    metrics = main()
    
    
    