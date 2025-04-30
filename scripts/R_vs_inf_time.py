#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 09:18:31 2025

@author: pietro

This script contains a set of functions used to analyse and plot the relationship
between the infection rate in the model and the observed average infection time.

"""

import simplicity.output_manager as om
import pandas as pd
import numpy as np
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def get_avg_infectious_duration(simulation_output_dir, state_filter, min_final_time=None):
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    durations = []

    for ssod in ssod_list:
        if min_final_time is not None:
            final_time = om.read_final_time(ssod)
            if final_time is None or final_time < min_final_time:
                continue

        df = om.read_individuals_data(ssod)
        if df.empty or "state" not in df.columns:
            continue

        if state_filter == 'recovered':
            filtered_df = df[(df["state"] == state_filter)]
        else:
            filtered_df = df[(df["state"] == state_filter)]

        if not filtered_df.empty:
            duration = (filtered_df["t_not_infectious"] - filtered_df["t_infectious"]).dropna()
            durations.extend(duration.tolist())

    if not durations:
        return None, None

    return np.mean(durations), np.std(durations)

def get_avg_infection_duration(simulation_output_dir, state_filter, min_final_time=None):
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    durations = []

    for ssod in ssod_list:
        if min_final_time is not None:
            final_time = om.read_final_time(ssod)
            if final_time is None or final_time < min_final_time:
                continue

        df = om.read_individuals_data(ssod)
        if df.empty or "state" not in df.columns:
            continue

        filtered_df = df[(df["state"] == state_filter)]
        if not filtered_df.empty:
            duration = (filtered_df["t_not_infected"] - filtered_df["t_infection"]).dropna()
            durations.extend(duration.tolist())

    if not durations:
        return None, None

    return np.mean(durations), np.std(durations)

def get_avg_delta_t(simulation_output_dir, min_final_time=None):
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    all_dts = []

    for ssod in ssod_list:
        if min_final_time is not None:
            final_time = om.read_final_time(ssod)
            if final_time is None or final_time < min_final_time:
                continue

        df = om.read_simulation_trajectory(ssod)
        if df.empty or "time" not in df.columns:
            continue

        dts = df["time"].diff().dropna()
        all_dts.extend(dts.tolist())

    if not all_dts:
        return None, None

    return np.mean(all_dts), np.std(all_dts)

def get_avg_debug_delta_t(simulation_output_dir, min_final_time=None):
    ssod_list = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    all_dts = []

    for ssod in ssod_list:
        if min_final_time is not None:
            final_time = om.read_final_time(ssod)
            if final_time is None or final_time < min_final_time:
                continue

        df = om.read_DEBUG_update_ih(ssod)
        if df.empty or "delta_t" not in df.columns:
            continue

        dts = df["delta_t"].dropna()
        all_dts.extend(dts.tolist())

    if not all_dts:
        return None, None

    return np.mean(all_dts), np.std(all_dts)

def extract_simulation_summary(experiment_name, min_final_time=None):
    rows = []
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)

    for sim_out_dir in tqdm(simulation_output_dirs, desc="Extracting summaries"):
        phenotype_model = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'phenotype_model')
        R = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'R')
        diagnosis_rate = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'diagnosis_rate')

        recovered_avg_infectious, recovered_std_infectious = get_avg_infectious_duration(sim_out_dir, "recovered", min_final_time)
        recovered_avg_infection, recovered_std_infection = get_avg_infection_duration(sim_out_dir, "recovered", min_final_time)
        diagnosed_avg_infectious, diagnosed_std_infectious = get_avg_infectious_duration(sim_out_dir, "diagnosed", min_final_time)
        delta_t_avg, delta_t_std = get_avg_delta_t(sim_out_dir, min_final_time)
        # delta_t_update_avg, delta_t_update_std = get_avg_debug_delta_t(sim_out_dir, min_final_time)

        rows.append({
            "phenotype_model": phenotype_model,
            "R": R,
            "diagnosis_rate": diagnosis_rate,

            "recovered_avg_infectious": recovered_avg_infectious,
            "recovered_std_infectious": recovered_std_infectious,
            "recovered_avg_infection": recovered_avg_infection,
            "recovered_std_infection": recovered_std_infection,

            "diagnosed_avg_infectious": diagnosed_avg_infectious,
            "diagnosed_std_infectious": diagnosed_std_infectious,

            "delta_t_avg": delta_t_avg,
            "delta_t_std": delta_t_std,

            # "delta_t_update_avg": delta_t_update_avg,
            # "delta_t_update_std": delta_t_update_std
        })

    return pd.DataFrame(rows)


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


def plot_combined_summary(df, experiment_name):
    sns.set(style="whitegrid")

    phenotype_models = sorted(df["phenotype_model"].unique())
    diagnosis_rates = sorted(df["diagnosis_rate"].unique())

    n_rates = len(diagnosis_rates)
    jitter_spacing = 0.02
    jitter_offsets = {
        rate: (i - n_rates // 2) * jitter_spacing
        for i, rate in enumerate(diagnosis_rates)
    }

    colors = sns.color_palette("tab10", n_colors=n_rates)
    color_map = {d: c for d, c in zip(diagnosis_rates, colors)}

    # 6 rows: R_eff, recovered infection, recovered infectious, diagnosed infectious, Δt
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(14, 20), sharex='col')

    for col, phenotype in enumerate(phenotype_models):
        pheno_df = df[df["phenotype_model"] == phenotype]

        for d_rate in diagnosis_rates:
            d_sub = pheno_df[pheno_df["diagnosis_rate"] == d_rate].copy()
            
            # Ensure all numeric and fill missing
            for col_name in [
                "R_eff_avg", "R_eff_std",
                "recovered_avg_infection", "recovered_std_infection",
                "recovered_avg_infectious", "recovered_std_infectious",
                "diagnosed_avg_infectious", "diagnosed_std_infectious",
                "delta_t_avg", "delta_t_std"
            ]:
                d_sub[col_name] = pd.to_numeric(d_sub[col_name], errors="coerce").fillna(0)

            jitter = jitter_offsets[d_rate]
            R_jittered = d_sub["R"] + jitter

            # Row 1: Recovered infection duration
            ax = axes[0][col]
            label = f"d_rate={d_rate}" if col == 0 else None  # Only assign label once per rate
            ax.errorbar(R_jittered, d_sub["recovered_avg_infection"],
                        yerr=d_sub["recovered_std_infection"], fmt='o', capsize=3,
                        linestyle='None', color=color_map[d_rate])
            ax.set_ylabel("Infection duration (recovered)")

            # Row 2: Recovered infectious duration
            ax = axes[1][col]
            ax.errorbar(R_jittered, d_sub["recovered_avg_infectious"],
                        yerr=d_sub["recovered_std_infectious"], fmt='o', capsize=3,
                        linestyle='None', color=color_map[d_rate], label=None)
            ax.set_ylabel("Infectious duration (recovered)")

            # Row 3: Diagnosed infectious duration
            ax = axes[2][col]
            ax.errorbar(R_jittered, d_sub["diagnosed_avg_infectious"],
                        yerr=d_sub["diagnosed_std_infectious"], fmt='o', capsize=3,
                        linestyle='None', color=color_map[d_rate], label=None)
            ax.set_ylabel("Infectious duration (diagnosed)")
            ax.set_xlabel("Input R")

            # Row 4: Δt (step size)
            ax = axes[3][col]
            ax.errorbar(R_jittered, d_sub["delta_t_avg"], yerr=d_sub["delta_t_std"],
                        fmt='o', capsize=3, linestyle='None', color=color_map[d_rate], label=None)
            ax.set_ylabel("Avg Δt (step size)")
            ax.set_xlabel("Input R")
            ax.set_yscale('log')
            y_vals = d_sub["delta_t_avg"] - d_sub["delta_t_std"]
            y_vals = y_vals[(y_vals > 0) & (~y_vals.isna())]
            
            if not y_vals.empty:
                min_y = max(y_vals.min() * 0.5, 1e-6)
            else:
                min_y = 1e-3  # fallback lower limit
                
            ax.set_ylim(bottom=min_y)

        axes[0][col].set_title(f"Phenotype: {phenotype}")
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        axes[0][1].legend(handles, labels, title="Diagnosis rate")
        
    plt.tight_layout()
    plt.savefig(f'R_eff_time_inf_{experiment_name}_check.png')
    plt.close()


def plot_avg_deltat_update_vs_R(df_summary,experiment_name):
    """
    Plot average delta_t from DEBUG_update_ih files vs R, grouped by phenotype_model.
    """
    sns.set(style="whitegrid")
    df = df_summary
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
                subset["delta_t_update_avg"],
                yerr=subset["delta_t_update_std"],
                fmt='o',
                label=f"d_rate={d_rate}",
                capsize=3,
                linestyle='None',
                color=color_map[d_rate]
            )

        ax.set_title(f"Phenotype: {phenotype}")
        ax.set_xlabel("R")

    axes[0].set_ylabel("Average Δt (DEBUG_update_ih)")
    plt.tight_layout()
    plt.savefig(f"avg_deltat_update_vs_R_{experiment_name}.png")
    plt.close()

def get_all_debug_update_ih(experiment_name):
    """
    Aggregates DEBUG_update_ih data from all seeded simulation output dirs
    in all experiments under the given experiment_name.
    """
    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)
    all_dfs = []

    for sim_out_dir in sim_out_dirs:
        ssods = dm.get_seeded_simulation_output_dirs(sim_out_dir)
        for ssod in ssods:
            df = om.read_DEBUG_update_ih(ssod)
            if not df.empty:
                df["ssod"] = ssod
                all_dfs.append(df)

    if not all_dfs:
        return pd.DataFrame()

    return pd.concat(all_dfs, ignore_index=True)

import os 

def plot_debug_deltat_timeseries(experiment_name):
    """
    Extracts and plots delta_t vs sim_time from DEBUG_update_ih across all simulations.
    - One subplot per (phenotype, R) in a grid: rows = R, cols = phenotype
    """
    print("Loading DEBUG_update_ih data...")
    df_debug = get_all_debug_update_ih(experiment_name)

    if df_debug.empty:
        print("No DEBUG data found. Skipping timeseries plot.")
        return
    
    df_debug["sim_out_dir"] = df_debug["ssod"].apply(os.path.dirname)
    # Add metadata
    df_debug["phenotype_model"] = df_debug["sim_out_dir"].apply(
        lambda s: sm.get_parameter_value_from_simulation_output_dir(s, "phenotype_model")
    )
    df_debug["R"] = df_debug["sim_out_dir"].apply(
        lambda s: sm.get_parameter_value_from_simulation_output_dir(s, "R")
    )

    phenotype_models = sorted(df_debug["phenotype_model"].unique())
    R_values = sorted(df_debug["R"].unique())
    n_rows = len(R_values)
    n_cols = len(phenotype_models)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 3.5 * n_rows), sharex=True)

    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for i_r, R in enumerate(R_values):
        for j_p, phenotype in enumerate(phenotype_models):
            ax = axes[i_r, j_p]
            subset = df_debug[(df_debug["phenotype_model"] == phenotype) & (df_debug["R"] == R)]
            top_10_ssods = sorted(subset["ssod"].unique())[:10]

            for ssod in top_10_ssods:
                sim_df = subset[subset["ssod"] == ssod]
                ax.plot(sim_df["sim_time"], sim_df["delta_t"], alpha=0.8, lw=1)

            if i_r == n_rows - 1:
                ax.set_xlabel("Simulation Time")
            if j_p == 0:
                ax.set_ylabel(f"Δt (R={R})")

            ax.set_title(f"Phenotype: {phenotype}" if i_r == 0 else "")

    fig.suptitle("Δt vs Simulation Time — First 10 Runs per R/Phenotype", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"delta_t_update_timeseries_all_{experiment_name}.png")
    plt.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot infectious durations from simulation outputs.")
    parser.add_argument("experiment_name", type=str, help="Name of the experiment")
    args = parser.parse_args()

    experiment_name = args.experiment_name

    print(f"Extracting simulation summary for: {experiment_name}")
    df_summary = extract_simulation_summary(experiment_name, min_final_time=50)

    print("Generating plots...")
    plot_combined_summary(df_summary,experiment_name)
    # plot_avg_deltat_update_vs_R(df_summary,experiment_name)
    # plot_debug_deltat_timeseries(experiment_name)

    print("All plots complete.")

if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    