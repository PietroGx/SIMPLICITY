import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig4_preprocess_data as preproc
import fig4_plots as plots

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_plots import add_global_pie_legend

SCENARIOS = ["control", "SOT", "HIV_low", "HIV_high"]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate Figure 4 for the long-shedders paper (tree comparison)")
    parser.add_argument('--exp-num', type=int, default=4, help="Target experiment number (default: 4)")
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


def build_figure_4(exp_num, seed, cluster_threshold, fmt):
    set_nature_rcparams()

    width_in = 180 / 25.4
    height_in = 70 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    gs = gridspec.GridSpec(1, len(SCENARIOS), figure=fig, wspace=0.3)

    print(f"--- Building Figure 4 (seed={seed}, exp #{exp_num}) ---")

    for col, scenario in enumerate(SCENARIOS):
        print(f"Tree {col+1}/{len(SCENARIOS)}: {scenario}...")
        ax = fig.add_subplot(gs[0, col])
        bt_tree = preproc.get_tree_for_scenario(exp_num, scenario, seed, cluster_threshold)
        plots.plot_fig4_tree(ax, bt_tree, scenario)

    add_global_pie_legend(fig)
    fig.subplots_adjust(bottom=0.15, top=0.8, left=0.05, right=0.98)

    output_filename = f'Figure_4.{fmt}'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_4(args.exp_num, args.seed, args.cluster_threshold, args.format)
