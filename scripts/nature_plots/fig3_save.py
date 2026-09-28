import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig3_preprocess_data as preproc
import fig3_plots as plots
from _scenarios import (BOUND_EXP_NAME, figure_path, resolve_scenarios,
                        scenario_meta)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_plots import plot_fig4_pies, add_global_pie_legend

SCENARIO_PALETTE = {
    "control": "#333333", "SOT": "#56B4E9", "HIV_low": "#D55E00",
    "HIV_high": "#E69F00", "edge_case": "#CC79A7",
}


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Figure 3 for the long-shedders paper")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number to plot. Required on purpose: the "
                            "old default of 4 pointed at one of the runs that "
                            "produced invalid science.")
    parser.add_argument('--exp-name', type=str, default=BOUND_EXP_NAME,
                        help=f"Pipeline arm (default: {BOUND_EXP_NAME}).")
    parser.add_argument('--cluster-threshold', type=int, default=5)
    parser.add_argument('--min-days', type=int, default=100)
    parser.add_argument('--format', type=str, choices=['pdf', 'png'], default='png')
    return parser.parse_args()


def set_nature_rcparams():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


def add_panel_label(ax, label, x_offset=-0.15, y_offset=1.15):
    ax.text(x_offset, y_offset, label, transform=ax.transAxes,
           fontsize=8, fontweight='bold', va='top', ha='right')


def build_figure_3(exp_num, exp_name, cluster_threshold, min_days, fmt):
    set_nature_rcparams()

    scenarios, _ = resolve_scenarios(exp_name, exp_num)
    if not scenarios:
        raise SystemExit(f"No scenario output found for {exp_name} #{exp_num}.")

    fig = plt.figure(figsize=(180 / 25.4, 150 / 25.4))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.42, hspace=0.62)

    # --- A: conversion efficiency, the non-circular measurement -------------
    print("Panel A: conversion efficiency...")
    eff_raw = preproc.get_efficiency_data(exp_num, scenarios, exp_name=exp_name,
                                          min_days=min_days)
    eff = preproc.efficiency_ratio(eff_raw, threshold=cluster_threshold)
    ax_a = fig.add_subplot(gs[0, 0])
    plots.plot_efficiency_bars(ax_a, eff, scenario_order=scenarios)
    add_panel_label(ax_a, "A")

    # --- B: what drives it --------------------------------------------------
    print("Panel B: efficiency vs duration...")
    ax_b = fig.add_subplot(gs[0, 1])
    plots.plot_efficiency_vs_duration(ax_b, eff, scenario_meta(exp_name))
    add_panel_label(ax_b, "B")

    # --- C: outcome ---------------------------------------------------------
    print("Panel C: clade metrics...")
    gs_c = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1, 0],
                                            wspace=0.45, hspace=0.75)
    c_axes = [fig.add_subplot(gs_c[0, 0]), fig.add_subplot(gs_c[0, 1]),
              fig.add_subplot(gs_c[1, 0]), fig.add_subplot(gs_c[1, 1])]
    metrics_df = preproc.get_panel_b_data(exp_num, scenarios, cluster_threshold,
                                          min_days, exp_name=exp_name)
    plots.plot_fig3_metrics(c_axes, metrics_df, palette=SCENARIO_PALETTE,
                            scenario_order=scenarios)
    add_panel_label(c_axes[0], "C")

    # --- D: robustness ------------------------------------------------------
    print("Panel D: threshold x criterion robustness...")
    ax_d = fig.add_subplot(gs[1, 1])
    plots.plot_efficiency_robustness(ax_d, preproc.efficiency_robustness(eff_raw),
                                     scenario_order=scenarios)
    add_panel_label(ax_d, "D")

    fig.suptitle(
        f"Do long shedders drive variant emergence?   "
        f"(clades at {cluster_threshold} substitutions; all {len(scenarios)} scenarios, "
        f"all seeds)", fontsize=7.5, y=0.995)

    output_filename = figure_path(3, exp_name, exp_num, fmt)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_3(args.exp_num, args.exp_name, args.cluster_threshold,
                  args.min_days, args.format)
