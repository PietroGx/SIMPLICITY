import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig3_preprocess_data as preproc
import fig3_plots as plots
from _scenarios import BOUND_EXP_NAME, figure_path, resolve_scenarios

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
    parser.add_argument('--group', type=str, default="HIV_high",
                        help="Scenario for Panels A and C (default: HIV_high)")
    parser.add_argument('--seed', type=str, default="1", help="Representative seed for Panel C's scatter")
    parser.add_argument('--cluster-threshold', type=int, default=5)
    parser.add_argument('--min-days', type=int, default=100)
    parser.add_argument('--max-long-per-individual', type=int, default=None,
                        help="Cap sequenced lineages per long shedder (default: no subsampling)")
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


def build_figure_3(exp_num, exp_name, group, seed, cluster_threshold, min_days,
                   max_long_per_individual, fmt):
    set_nature_rcparams()

    scenarios, _ = resolve_scenarios(exp_name, exp_num)
    if not scenarios:
        raise SystemExit(f"No scenario output found for {exp_name} #{exp_num}.")
    # --group used to be validated by argparse choices=SCENARIOS, evaluated at
    # import against a hardcoded list -- which rejected edge_case outright.
    # Validate against what this pipeline actually has instead.
    if group not in scenarios:
        raise SystemExit(
            f"--group {group!r} is not available for {exp_name} #{exp_num}. "
            f"Choose from: {scenarios}")

    width_in = 180 / 25.4
    height_in = 150 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.5)

    print(f"--- Building Figure 3 (group={group}, exp #{exp_num}) ---")

    # --- Panel A: pies for the selected group, nested 2x2 ---
    gs_pies = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 0], wspace=0.3, hspace=0.4)
    pie_axes = [fig.add_subplot(gs_pies[0, 0]), fig.add_subplot(gs_pies[0, 1]),
               fig.add_subplot(gs_pies[1, 0]), fig.add_subplot(gs_pies[1, 1])]

    print("Panel A: pies...")
    pie_data, valid_seeds = preproc.get_panel_a_data(exp_num, group, cluster_threshold, min_days,
                                                  exp_name=exp_name)
    if valid_seeds:
        plot_fig4_pies(pie_axes, pie_data)
        add_global_pie_legend(fig)
    else:
        for ax in pie_axes:
            ax.axis("off")
        pie_axes[0].text(0.5, 0.5, "No Data", ha="center", va="center", fontsize=8)
    add_panel_label(pie_axes[0], "A")

    # --- Panel B: the four clade metrics, one sub-panel each ---
    print("Panel B: clade metrics...")
    gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 1],
                                            wspace=0.45, hspace=0.75)
    b_axes = [fig.add_subplot(gs_b[0, 0]), fig.add_subplot(gs_b[0, 1]),
              fig.add_subplot(gs_b[1, 0]), fig.add_subplot(gs_b[1, 1])]
    metrics_df = preproc.get_panel_b_data(exp_num, scenarios, cluster_threshold, min_days,
                                           exp_name=exp_name)
    plots.plot_fig3_metrics(b_axes, metrics_df, palette=SCENARIO_PALETTE,
                            scenario_order=scenarios)
    add_panel_label(b_axes[0], "B")

    # --- Panel C: sequence-space PCA, single seed of the selected group ---
    ax_c = fig.add_subplot(gs[1, 0])
    print(f"Panel C: sequence PCA (seed {seed})...")
    snp_df, labels = preproc.get_panel_c_scatter_data(exp_num, group, seed, max_long_per_individual,
                                                    exp_name=exp_name)
    plots.plot_fig3_sequence_pca(ax_c, snp_df, labels)
    add_panel_label(ax_c, "C")

    # --- Panel D: per-seed consistency check for the selected group ---
    ax_d = fig.add_subplot(gs[1, 1])
    print("Panel D: per-seed consistency...")
    consistency_df = preproc.get_panel_c_consistency_data(exp_num, group, max_long_per_individual,
                                                       exp_name=exp_name)
    plots.plot_fig3_consistency(ax_d, consistency_df)
    add_panel_label(ax_d, "D")

    # Three of the four panels describe ONE scenario; only B spans all of them.
    # Say so on the figure rather than leaving it to the caption.
    fig.suptitle(f"Scenario: {group}   —   A, C, D show this scenario only; "
                 f"B spans all {len(scenarios)}",
                 fontsize=7.5, y=0.998)
    fig.subplots_adjust(bottom=0.08, top=0.88, left=0.08, right=0.95)

    # group MUST be in the name: figure 3 is per-scenario, so without it
    # every scenario would overwrite the same file.
    output_filename = figure_path(3, exp_name, exp_num, fmt, seed=seed,
                                  group=group)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_3(args.exp_num, args.exp_name, args.group, args.seed,
                  args.cluster_threshold, args.min_days,
                  args.max_long_per_individual, args.format)
