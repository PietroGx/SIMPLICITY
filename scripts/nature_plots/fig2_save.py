import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig2_preprocess_data as preproc
from _scenarios import BOUND_EXP_NAME, figure_path, resolve_scenarios
import fig2_plots as plots

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Figure 2 for SIMPLICITY Nature Paper")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number to plot. Required on purpose: the "
                            "old default of 4 pointed at one of the runs that "
                            "produced invalid science.")
    parser.add_argument('--exp-name', type=str, default=BOUND_EXP_NAME,
                        help=f"Pipeline arm (default: {BOUND_EXP_NAME}).")
    parser.add_argument('--seed', type=str, default="1", help="Target seed for lineage frequency plots (Rows 1 & 2)")
    parser.add_argument('--format', type=str, choices=['pdf', 'png'], default='png', help="Output format: pdf or png (default: png)")
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
        'ps.fonttype': 42
    })

def build_figure_2(exp_num, exp_name, target_seed, fmt):
    set_nature_rcparams()
    # from the pipeline's own config; edge_case appears on the unbound arm
    scenarios, _ = resolve_scenarios(exp_name, exp_num)
    if not scenarios:
        raise SystemExit(f"No scenario output found for {exp_name} #{exp_num}.")
    palette = {"control": "#333333", "SOT": "#56B4E9", "HIV_low": "#D55E00",
               "HIV_high": "#E69F00", "edge_case": "#CC79A7"}

    # --- Phase A: Pre-calculation for strict alignments ---
    shared_violin_seeds = preproc.get_shared_valid_seeds(exp_num, scenarios, exp_name=exp_name)
    
    print(f"Calculating global time bounds for Seed #{target_seed}...")
    global_max_t = 0
    for scenario in scenarios:
        _, _, t_final = preproc.get_fig2_freq_data(exp_num, scenario, target_seed, exp_name=exp_name)
        if t_final > global_max_t:
            global_max_t = t_final
            
    if global_max_t == 0:
        print(f"Warning: Seed {target_seed} not found or produced empty timeframes.")
        global_max_t = 365 
        
    print(f"Global X-axis standardized to {global_max_t} days.")

    # --- Phase B: Build Grid ---
    width_in = 180 / 25.4
    height_in = 150 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    
    # One column per scenario, widely spaced so the scenarios read as separate
    # blocks rather than one continuous grid.
    gs = gridspec.GridSpec(3, len(scenarios), figure=fig, wspace=0.62, hspace=0.45)
    
    print(f"--- Building Figure 2 Grid (Exp #{exp_num}) ---")

    # Row 3 shares one y-scale across scenarios; per-column autoscaling made
    # control's jumps (max 3) look comparable to edge_case's (max 15).
    global_div_max = 0
    all_dists = {}
    for scenario in scenarios:
        d = preproc.get_fig2_divergence_data(exp_num, scenario, shared_violin_seeds,
                                             exp_name=exp_name)
        all_dists[scenario] = d
        if d:
            global_div_max = max(global_div_max, max(d))
    global_div_max = global_div_max * 1.08 if global_div_max else None
    print(f"Divergence-jump y-axis standardized to {global_div_max}")

    div_axes = []
    for col, scenario in enumerate(scenarios):
        display_name = preproc.get_clinical_label(scenario)
        print(f"Processing column {col+1}/{len(scenarios)}: {display_name}...")
        
        ax_freq = fig.add_subplot(gs[0, col])
        ax_clust = fig.add_subplot(gs[1, col])
        ax_div = fig.add_subplot(gs[2, col])
        
        lf, cmap_df, _ = preproc.get_fig2_freq_data(exp_num, scenario, target_seed, exp_name=exp_name)
        _, clade_to_lineages, clade_meta_df, _ = preproc.get_fig2_clustered_data(exp_num, scenario, target_seed, exp_name=exp_name)
        dists = all_dists[scenario]
        
        plots.plot_fig2_lineage_freq(ax_freq, lf, cmap_df, global_max_t, title=display_name)
        plots.plot_fig2_clustered_freq(ax_clust, lf, clade_to_lineages, clade_meta_df, cmap_df, global_max_t)
        plots.plot_fig2_divergence_violin(ax_div, dists, color=palette[scenario],
                                          y_max=global_div_max)
        div_axes.append(ax_div)

        # One letter per scenario, above its column.
        plots.add_column_label(ax_freq, chr(ord('A') + col))

        if col == 0:
            ax_freq.set_ylabel("Lineage Freq", fontsize=7)
            ax_clust.set_ylabel("Clade Freq", fontsize=7)
            ax_div.set_ylabel("distinct gen.pos.count", fontsize=7)
        else:
            ax_freq.set_ylabel("")
            ax_clust.set_ylabel("")
            ax_div.set_ylabel("")
            
        ax_freq.set_xlabel("")
        ax_clust.set_xlabel("Time (days)")
        
    output_filename = figure_path(2, exp_name, exp_num, fmt, seed=target_seed)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")

if __name__ == "__main__":
    args = parse_arguments()
    build_figure_2(args.exp_num, args.exp_name, args.seed, args.format)
