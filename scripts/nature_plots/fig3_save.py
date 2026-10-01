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
                        help=f"Pipeline to plot (default: {BOUND_EXP_NAME}).")
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

    fig = plt.figure(figsize=(180 / 25.4, 165 / 25.4))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.40, hspace=0.55,
                           top=0.90, bottom=0.07)

    eff_raw = preproc.get_efficiency_data(exp_num, scenarios, exp_name=exp_name,
                                          min_days=min_days)
    eff = preproc.efficiency_ratio(eff_raw, threshold=cluster_threshold)

    # --- A: which host type produced the standout clade --------------------
    print("Panel A: standout clade by origin...")
    wr = preproc.winrate_summary(
        preproc.get_winrate_data(exp_num, scenarios, exp_name=exp_name,
                                 cluster_threshold=cluster_threshold,
                                 min_days=min_days))
    gs_a = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0, 0],
                                            hspace=0.50)
    a_axes = [fig.add_subplot(gs_a[i, 0]) for i in range(3)]
    plots.plot_winrate(a_axes, wr, scenario_order=scenarios)
    # Header and legend positioned from the block's own bbox, so they sit
    # inside panel A rather than drifting into the suptitle or the panel below.
    box = gs[0, 0].get_position(fig)
    fig.text(box.x0 - 0.030, box.y1 + 0.052, "A",
             fontsize=8, fontweight="bold", ha="right", va="bottom")
    fig.text(box.x0 - 0.042, box.y0 + box.height / 2, plots.WIN_YLABEL,
             rotation=90, va="center", ha="center", fontsize=5.8)
    plots.winrate_legend(fig, loc="upper center",
                         bbox_to_anchor=(box.x0 + box.width / 2, box.y1 + 0.056),
                         ncol=2)

    # --- B: conversion rate, pooled ----------------------------------------
    print("Panel B: conversion rate, pooled...")
    ax_b = fig.add_subplot(gs[0, 1])
    plots.plot_efficiency_bars(ax_b, eff, scenario_order=scenarios)
    add_panel_label(ax_b, "B", y_offset=1.22)

    # --- C: conversion rate, per simulation --------------------------------
    print("Panel C: conversion rate, per simulation...")
    ax_c = fig.add_subplot(gs[1, 0])
    plots.plot_efficiency_vs_duration(ax_c, eff_raw, scenario_meta(exp_name),
                                      threshold=cluster_threshold,
                                      scenario_order=scenarios)
    add_panel_label(ax_c, "C", y_offset=1.30)

    # --- D: robustness to both analysis choices ----------------------------
    print("Panel D: robustness...")
    ax_d = fig.add_subplot(gs[1, 1])
    plots.plot_efficiency_robustness(ax_d, preproc.efficiency_robustness(eff_raw),
                                     scenario_order=scenarios,
                                     palette=SCENARIO_PALETTE)
    add_panel_label(ax_d, "D", y_offset=1.30)

    fig.suptitle(
        "Do long shedders drive variant emergence?"
        f"   (clades at {cluster_threshold} substitutions; "
        f"simulations \u2265 {preproc.MIN_FINAL_TIME} d, "
        f"measured after day {preproc.BURNIN_CUTOFF_DAYS})",
        fontsize=7.5, y=1.025)

    output_filename = figure_path(3, exp_name, exp_num, fmt)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_3(args.exp_num, args.exp_name, args.cluster_threshold,
                  args.min_days, args.format)
