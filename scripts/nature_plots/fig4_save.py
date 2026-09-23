import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig4_preprocess_data as preproc
import fig4_plots as plots
from _scenarios import BOUND_EXP_NAME, figure_path, resolve_scenarios

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_plots import add_global_pie_legend




def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate Figure 4 for the long-shedders paper (tree comparison)")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number to plot. Required on purpose: the "
                            "old default of 4 pointed at one of the runs that "
                            "produced invalid science.")
    parser.add_argument('--exp-name', type=str, default=BOUND_EXP_NAME,
                        help=f"Pipeline arm (default: {BOUND_EXP_NAME}).")
    parser.add_argument('--seed', type=str, default="1",
                        help="Seed used for all 5 trees (same seed across scenarios)")
    parser.add_argument('--cluster-threshold', type=int, default=5)
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


def build_figure_4(exp_num, exp_name, seed, cluster_threshold, fmt):
    set_nature_rcparams()

    scenarios, _ = resolve_scenarios(exp_name, exp_num)
    if not scenarios:
        raise SystemExit(f"No scenario output found for {exp_name} #{exp_num}.")

    width_in = 180 / 25.4
    height_in = 70 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    gs = gridspec.GridSpec(1, len(scenarios), figure=fig, wspace=0.3)

    print(f"--- Building Figure 4 (seed={seed}, exp #{exp_num}) ---")

    for col, scenario in enumerate(scenarios):
        print(f"Tree {col+1}/{len(scenarios)}: {scenario}...")
        ax = fig.add_subplot(gs[0, col])
        bt_tree = preproc.get_tree_for_scenario(exp_num, scenario, seed,
                                                cluster_threshold, exp_name=exp_name)
        # get_tree_for_scenario returns None when the seed has no tree; passing
        # that straight to the plotter crashed. Figures 2 and 3 already tolerate
        # a missing scenario -- this one did not.
        if bt_tree is None:
            print(f"  [warn] no tree for '{scenario}' seed {seed}; empty panel.")
            ax.axis("off")
            ax.text(0.5, 0.5, f"{scenario}\nno tree", ha="center", va="center",
                    fontsize=6.5, color="0.45")
            continue
        plots.plot_fig4_tree(ax, bt_tree, scenario)

    add_global_pie_legend(fig)
    fig.subplots_adjust(bottom=0.15, top=0.8, left=0.05, right=0.98)

    output_filename = figure_path(4, exp_name, exp_num, fmt, seed=seed)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_4(args.exp_num, args.exp_name, args.seed, args.cluster_threshold,
                  args.format)
