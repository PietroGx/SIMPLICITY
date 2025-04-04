#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 09:18:31 2025

@author: pietro
"""

import simplicity.output_manager as om
import pandas as pd
import numpy as np
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def R_eff_trajectory(ssod, time_window):
    df = om.read_individuals_data(ssod)
    final_time = om.read_final_time(ssod)

    trajectory = []

    for t in range(time_window, final_time + 1):
        window_start = t - time_window
        window_end = t  # exclusive

        # Select individuals whose entire infectious period is within the window
        in_window = df[
            (df['t_infectious'] >= window_start) &
            (df['t_not_infectious'] <= window_end)
        ]

        if not in_window.empty:
            avg_secondary = in_window['new_infections'].apply(
                lambda x: len(x) if isinstance(x, (list, tuple)) else 0
            ).mean()
        else:
            avg_secondary = np.nan  # or 0, depending on how you want to handle gaps

        trajectory.append((t, avg_secondary))

    # Convert to DataFrame
    return pd.DataFrame(trajectory, columns=["t", "R_eff"])

def get_avg_simulation_R_eff(simulation_output_dir):
    
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    all_infection_counts = []
    
    for ssod in ssod_list:
        df = om.read_individuals_data(ssod)
        
        if df.empty or "new_infections" not in df.columns:
            print(f"Skipping {ssod}: missing or empty 'new_infections'")
            continue
        
        infection_counts = df["new_infections"].apply(lambda x: len(x) if isinstance(x, (list, tuple)) else 0)
        all_infection_counts.extend(infection_counts.tolist())  # Flatten into one list
        
    if not all_infection_counts:
        return None
    R_eff_avg = np.mean(all_infection_counts)
    R_eff_std = np.std(all_infection_counts)
    return R_eff_avg, R_eff_std

def get_avg_infectious_duration(simulation_output_dir, state_filter):
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    durations = []

    for ssod in ssod_list:
        df = om.read_individuals_data(ssod)

        if df.empty or "state" not in df.columns:
            continue

        filtered_df = df[df["state"] == state_filter]
        if not filtered_df.empty:
            duration = (filtered_df["t_not_infectious"] - filtered_df["t_infectious"]).dropna()
            durations.extend(duration.tolist())

    if not durations:
        return None, None

    return np.mean(durations), np.std(durations)

def get_avg_delta_t(simulation_output_dir):
    """
    Computes the average and std of delta_t (time step differences) across all runs.
    """
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    all_dts = []

    for ssod in ssod_list:
        df = om.read_simulation_trajectory(ssod)

        if df.empty or "time" not in df.columns:
            continue

        dts = df["time"].diff().dropna()
        all_dts.extend(dts.tolist())

    if not all_dts:
        return None, None

    return np.mean(all_dts), np.std(all_dts)

def extract_simulation_summary(experiment_name):
    rows = []
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)

     for sim_out_dir in tqdm(simulation_output_dirs, desc="Extracting summaries"):
        phenotype_model = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'phenotype_model')
        R = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'R')
        diagnosis_rate = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'diagnosis_rate')

        R_eff_avg, R_eff_std = get_avg_simulation_R_eff(sim_out_dir)
        recovered_avg, recovered_std = get_avg_infectious_duration(sim_out_dir, "recovered")
        diagnosed_avg, diagnosed_std = get_avg_infectious_duration(sim_out_dir, "diagnosed")
        delta_t_avg, delta_t_std = get_avg_delta_t(sim_out_dir)

        rows.append({
            "phenotype_model": phenotype_model,
            "R": R,
            "diagnosis_rate": diagnosis_rate,
            "R_eff_avg": R_eff_avg,
            "R_eff_std": R_eff_std,
            "recovered_avg_duration": recovered_avg,
            "recovered_std_duration": recovered_std,
            "diagnosed_avg_duration": diagnosed_avg,
            "diagnosed_std_duration": diagnosed_std,
            "delta_t_avg": delta_t_avg,
            "delta_t_std": delta_t_std,
        })

    return pd.DataFrame(rows)


def plot_R_eff_vs_R(df):
    sns.set(style="whitegrid")

    # Sort and map diagnosis rates to fixed jitter offsets
    diag_rates_sorted = sorted(df["diagnosis_rate"].unique())
    n_rates = len(diag_rates_sorted)
    jitter_spacing = 0.02
    jitter_offsets = {
        rate: (i - n_rates // 2) * jitter_spacing
        for i, rate in enumerate(diag_rates_sorted)
    }

    phenotype_models = sorted(df["phenotype_model"].unique())
    n_phenotypes = len(phenotype_models)

    # Create subplots (1 row, n_phenotypes columns)
    fig, axes = plt.subplots(1, n_phenotypes, figsize=(7 * n_phenotypes, 6), sharey=True)

    if n_phenotypes == 1:
        axes = [axes]  # ensure axes is always iterable

    for ax, phenotype in zip(axes, phenotype_models):
        subset = df[df["phenotype_model"] == phenotype]

        for diag_rate in diag_rates_sorted:
            
            # add jitter to better separate d_rates 
            diag_subset = subset[subset["diagnosis_rate"] == diag_rate]
            jitter = jitter_offsets[diag_rate]
            R_jittered = diag_subset["R"] + jitter

            ax.errorbar(
                R_jittered,
                diag_subset["R_eff_avg"],
                yerr=diag_subset["R_eff_std"],
                fmt='o',
                capsize=3,
                linestyle='None',
                label=f"diagnosis_rate = {diag_rate}"
            )

        ax.set_title(f"Phenotype: {phenotype}")
        ax.set_xlabel("Input R")

    axes[0].set_ylabel("Observed R_eff")
    axes[0].legend(title="Diagnosis Rate")
    plt.tight_layout()
    plt.show()

def plot_infectious_durations(df):
    sns.set(style="whitegrid")

    phenotype_models = sorted(df["phenotype_model"].unique())
    diagnosis_rates = sorted(df["diagnosis_rate"].unique())
    colors = sns.color_palette("tab10", n_colors=len(diagnosis_rates))
    color_map = {d: c for d, c in zip(diagnosis_rates, colors)}

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 10), sharex=True)

    for i, phenotype in enumerate(phenotype_models):
        pheno_df = df[df["phenotype_model"] == phenotype]

        for j, state in enumerate(["recovered", "diagnosed"]):
            ax = axes[i][j]
            y_col = f"{state}_avg_duration"
            yerr_col = f"{state}_std_duration"

            for d_rate in diagnosis_rates:
                d_sub = pheno_df[pheno_df["diagnosis_rate"] == d_rate]

                ax.errorbar(
                    d_sub["R"],
                    d_sub[y_col],
                    yerr=d_sub[yerr_col],
                    fmt='o',
                    label=f"d_rate={d_rate}",
                    capsize=3,
                    linestyle='None',
                    color=color_map[d_rate]
                )

            ax.set_title(f"{phenotype} — {state}")
            ax.set_ylabel("Avg infectious duration")
            ax.set_xlabel("R")

            if i == 0 and j == 1:
                ax.legend(title="Diagnosis rate")

    plt.tight_layout()
    plt.show()

def plot_avg_deltat_vs_R(df):
    sns.set(style="whitegrid")

    phenotype_models = sorted(df["phenotype_model"].unique())
    diagnosis_rates = sorted(df["diagnosis_rate"].unique())
    colors = sns.color_palette("tab10", n_colors=len(diagnosis_rates))
    color_map = {d: c for d, c in zip(diagnosis_rates, colors)}

    fig, axes = plt.subplots(1, len(phenotype_models), figsize=(7 * len(phenotype_models), 5), sharey=True)

    if len(phenotype_models) == 1:
        axes = [axes]

    for ax, phenotype in zip(axes, phenotype_models):
        pheno_df = df[df["phenotype_model"] == phenotype]

        for d_rate in diagnosis_rates:
            subset = pheno_df[pheno_df["diagnosis_rate"] == d_rate]

            ax.errorbar(
                subset["R"],
                subset["delta_t_avg"],
                yerr=subset["delta_t_std"],
                fmt='o',
                label=f"d_rate={d_rate}",
                capsize=3,
                linestyle='None',
                color=color_map[d_rate]
            )

        ax.set_title(f"Phenotype: {phenotype}")
        ax.set_xlabel("R")

    axes[0].set_ylabel("Average Δt (time step)")
    axes[0].legend(title="Diagnosis rate")
    plt.tight_layout()
    plt.show()


def plot_combined_summary(df,experiment_name):
    sns.set(style="whitegrid")

    # Get sorted values
    phenotype_models = sorted(df["phenotype_model"].unique())
    diagnosis_rates = sorted(df["diagnosis_rate"].unique())

    # Fixed jitter offsets per diagnosis_rate
    n_rates = len(diagnosis_rates)
    jitter_spacing = 0.02
    jitter_offsets = {
        rate: (i - n_rates // 2) * jitter_spacing
        for i, rate in enumerate(diagnosis_rates)
    }

    # Colors for diagnosis rates
    colors = sns.color_palette("tab10", n_colors=n_rates)
    color_map = {d: c for d, c in zip(diagnosis_rates, colors)}

    # Create subplot grid: 3 rows (R_eff, recovered, diagnosed) x 2 phenotypes
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(14, 16), sharex='col')

    for col, phenotype in enumerate(phenotype_models):
        pheno_df = df[df["phenotype_model"] == phenotype]

        for d_rate in diagnosis_rates:
            d_sub = pheno_df[pheno_df["diagnosis_rate"] == d_rate]
            jitter = jitter_offsets[d_rate]
            R_jittered = d_sub["R"] + jitter

            # --- Row 0: R_eff ---
            ax = axes[0][col]
            ax.errorbar(
                R_jittered,
                d_sub["R_eff_avg"],
                yerr=d_sub["R_eff_std"],
                fmt='o',
                capsize=3,
                linestyle='None',
                label=f"d_rate={d_rate}",
                color=color_map[d_rate]
            )
            ax.set_ylabel("Observed R_eff")

            # --- Row 1: Recovered durations ---
            ax = axes[1][col]
            ax.errorbar(
                R_jittered,
                d_sub["recovered_avg_duration"],
                yerr=d_sub["recovered_std_duration"],
                fmt='o',
                capsize=3,
                linestyle='None',
                color=color_map[d_rate]
            )
            ax.set_ylabel("Avg duration (recovered)")

            # --- Row 2: Diagnosed durations ---
            ax = axes[2][col]
            ax.errorbar(
                R_jittered,
                d_sub["diagnosed_avg_duration"],
                yerr=d_sub["diagnosed_std_duration"],
                fmt='o',
                capsize=3,
                linestyle='None',
                color=color_map[d_rate]
            )
            ax.set_ylabel("Avg duration (diagnosed)")
            ax.set_xlabel("Input R")
            
            # -- Row 3: Δt --
            ax = axes[3][col]
            ax.errorbar(
               R_jittered, d_sub["delta_t_avg"], yerr=d_sub["delta_t_std"],
               fmt='o', capsize=3, linestyle='None',
               color=color_map[d_rate]
            )
            ax.set_ylabel("Avg Δt (time step)")
            ax.set_xlabel("Input R")

        # Titles at top of each column
        axes[0][col].set_title(f"Phenotype: {phenotype}")

    # Legend in top-right subplot
    axes[0][1].legend(title="Diagnosis rate")

    plt.tight_layout()
    plt.savefig(f'R_eff_time_inf_{experiment_name}_check.png')
    plt.close()
    
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot R_eff and infectious durations from simulation outputs.")
    parser.add_argument("experiment_name", type=str, help="Name of the experiment")
    args = parser.parse_args()

    experiment_name = args.experiment_name

    print(f"Extracting simulation summary for: {experiment_name}")
    df_summary = extract_simulation_summary(experiment_name)

    print("Generating plots...")
    plot_combined_summary(df_summary,experiment_name)

if __name__ == "__main__":
    main()