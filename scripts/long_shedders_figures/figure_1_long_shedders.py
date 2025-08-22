#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 18:04:59 2025

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



# === Subplot A: IH model ===

def plot_comparison_intra_host_models_ax(ax, label="A"):
    """
    Plot intra-host dynamics comparison onto provided subplot axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target subplot axis.
    label : str
        Subplot label (e.g. 'A').
    """
    # --- Baseline ('normal') individual ---
    intra_host_model = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8)
    t = np.arange(0, 300, 1)
    y = intra_host_model.data_plot_ih_solution(0, 300, 1)[0]
    ax.plot(t, y, label='Normal', linewidth=2)

    # --- Long-shedder scenarios (11/16/32/72 days) ---
    long_shed_durations = [11, 16, 32, 72]
    for d in long_shed_durations:
        ih_long = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=d, tau_4=8)
        y_long = ih_long.data_plot_ih_solution(0, 300, 1)[0]
        ax.plot(t, y_long, linestyle='--', label=f'Long shedder – {d}d')

    ax.set_xlabel("Time (d)")
    ax.set_ylabel("P(infectious after t)")
    ax.set_xlim(0, 300)
    ax.set_ylim(0)  # keep your original y lower bound; remove or set to (0,1) if preferred
    ax.legend(title="Intra-host profiles", frameon=False, ncol=2)

    pm.apply_standard_axis_style(ax)

    # Subplot label
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

    
# === Subplot B: Avg Duration ===

def plot_avg_duration_ax(ax, data, label="B"):
    sns.violinplot(
        data=data,
        x='type',
        y='infection_duration',
        hue='type',
        palette='muted',
        cut=0,
        inner='quartile',
        linewidth=1,
        ax=ax
    )
    ax.set_ylabel("Average Duration")
    ax.set_xlabel("")

    ax.text(-0.05, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


# === Main Figure ===

def plot_figure_1(ssod):
    duration_data = preprocess.prepare_figure_1_data(ssod)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [4, 2]})

    plot_comparison_intra_host_models_ax(axes[0])
    # plot_avg_duration_ax(axes[1], duration_data)

    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()

    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    sim_out_name = dm.get_simulation_output_foldername_from_SSOD(ssod)
    seed = os.path.basename(ssod)
    output_path = os.path.join(experiment_plots_dir, f"Figure_1_{sim_out_name}_{seed}_.tiff")

    print(f"Saving figure to: {output_path}")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close()

def plot(experiment_name, seed_number):
    
    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)

    ssod = dm.get_ssod(sim_out_dirs[1], seed_number)
    plot_data = preprocess.analyze_long_shedders(ssod)
    metrics = plot_data[0]

    plot_figure_1(ssod)

    return metrics

def main():
    # # Set up the argument parser
    # parser = argparse.ArgumentParser(description="Plot")
    # parser.add_argument('experiment_name', type=str, help="experiment name")
    # args = parser.parse_args()
    # # Run the script with the provided parameter
    
    experiment_name = 'test_long_shedders_r1_kv_#1'
    seed_number = 1
    metrics = plot(experiment_name, seed_number)
    return metrics
    
if __name__ == "__main__":
    metrics = main()
