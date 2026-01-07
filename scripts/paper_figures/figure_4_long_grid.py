#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build comparison grids for long-shedders experiments.

Each figure: rows = tau_3_long values, cols = (R_long x long_shedders_ratio) combos.
Each cell: 2x2 pies (Peak, Burden, Survival, Growth).
Global legend added once from the first pies block.

Output filenames:
  Data/fig4_compare_long_evo_rate_f=<VAL>_CL_<CLUSTER>_MD_<MINDAYS>.png
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

import simplicity.dir_manager as dm
import figure_4_long_shedders as f4


def read_index_csv(path):
    """Read CSV robustly, handling Windows encoding if needed."""
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp1252")


def plot_one_figure(df, evo_val, cluster_threshold, min_days, outdir, exp_number_suffix="#1"):
    taus = sorted(df["tau_3_long"].unique())
    r_vals = sorted(df["R_long"].unique())
    ratios = sorted(df["long_shedders_ratio"].unique())

    n_tau = len(taus)
    n_r = len(r_vals)
    n_ratio = len(ratios)

    # Outer grid: rows = taus, cols = R*ratio combos
    fig = plt.figure(figsize=(4 * n_r * n_ratio, 4 * n_tau))
    outer = GridSpec(n_tau, n_r * n_ratio, figure=fig, wspace=0.25, hspace=0.35)
    # fig.suptitle("Figure 4 - long_evo_rate_f = {}".format(evo_val), fontsize=16)

    legend_added = False

    for i_tau, tau in enumerate(taus):
        for i_r, R in enumerate(r_vals):
            for i_ratio, ratio in enumerate(ratios):
                cell_col = i_r * n_ratio + i_ratio
                cell_spec = outer[i_tau, cell_col]

                # Subgrid for the 4 pies
                inner = GridSpecFromSubplotSpec(2, 2, subplot_spec=cell_spec, wspace=0.05, hspace=0.05)
                inner_axes = [
                    fig.add_subplot(inner[0, 0]),
                    fig.add_subplot(inner[0, 1]),
                    fig.add_subplot(inner[1, 0]),
                    fig.add_subplot(inner[1, 1]),
                ]

                # Filter subset
                sub = df[(df["tau_3_long"] == tau) &
                         (df["R_long"] == R) &
                         (df["long_shedders_ratio"] == ratio)]
                if sub.empty:
                    for ax in inner_axes:
                        ax.axis("off")
                    inner_axes[0].text(0.5, 0.5, "Missing", ha="center", va="center", fontsize=9)
                    continue

                exp_name = str(sub.iloc[0]["generated_experiment_name"]) + "_" + exp_number_suffix
                sods = dm.get_simulation_output_dirs(exp_name)
                if not sods:
                    for ax in inner_axes:
                        ax.axis("off")
                    inner_axes[0].text(0.5, 0.5, "No SOD", ha="center", va="center", fontsize=9)
                    continue

                sod = sods[0]
                try:
                    f4.summarize_sod_to_pies(
                        experiment_name=exp_name,
                        sod_path=sod,
                        cluster_threshold=cluster_threshold,
                        min_days=min_days,
                        axes=inner_axes,
                        savefig=False,
                        showlegend=False
                    )
                    # Label block with parameters
                    inner_axes[0].set_title("tau={} R={} ratio={}".format(tau, R, ratio),
                                            fontsize=9, loc="left")
                except Exception as e:
                    for ax in inner_axes:
                        ax.axis("off")
                    inner_axes[0].text(0.5, 0.5, "Error", ha="center", va="center", fontsize=9)
                    print("[WARN] Failed pies for {} (R={}, ratio={}): {}".format(exp_name, R, ratio, e))
                    continue

                
    try:
        f4.add_clade_and_panel_legend(fig, fontsize=10, pad=0.012)
    except Exception as e:
        print("[DEBUG] Legend creation failed:", e)
        
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(
        outdir,
        "Figure_4_long_grid={}_CL_{}_MD_{}.png".format(evo_val, cluster_threshold, min_days)
    )
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[SAVED] {}".format(outfile))


def main():
    ap = argparse.ArgumentParser(description="Build Figure 4 comparison grids")
    ap.add_argument("--index-csv",
                    default="scripts/long_shedders_experiments/INDEX_exp_scripts.csv",
                    help="Path to INDEX_exp_scripts.csv")
    ap.add_argument("--cluster-threshold", type=int, default=5,
                    help="Clustering threshold (default: 5)")
    ap.add_argument("--min-days", type=int, default=0,
                    help="Only include SSODs with final_time >= this many days (default: 0)")
    ap.add_argument("--exp-number", type=int, default=1,
                    help="Experiment number suffix (default: 1)")
    ap.add_argument("--outdir", default="Data",
                    help="Output directory (default: Data)")
    args = ap.parse_args()

    df = read_index_csv(args.index_csv)
    for evo_val in sorted(df["long_evo_rate_f"].unique()):
        sub = df[df["long_evo_rate_f"] == evo_val]
        if sub.empty:
            continue
        plot_one_figure(sub, evo_val, args.cluster_threshold, args.min_days,
                        args.outdir, "#{}".format(args.exp_number))


if __name__ == "__main__":
    main()


