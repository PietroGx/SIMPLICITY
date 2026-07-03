import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig2_preprocess_data as preproc
import fig2_plots as plots

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Figure 2 for SIMPLICITY Nature Paper")
    parser.add_argument('--exp-num', type=int, default=4, help="Target experiment number (default: 4)")
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

def build_figure_2(exp_num, target_seed, fmt):
    set_nature_rcparams()
    scenarios = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
    palette = {"control": "#333333", "SOT": "#56B4E9", "HIV_low": "#D55E00", "HIV_high": "#E69F00", "edge_case": "#CC79A7"}

    # --- Phase A: Pre-calculation for strict alignments ---
    shared_violin_seeds = preproc.get_shared_valid_seeds(exp_num, scenarios)
    
    print(f"Calculating global time bounds for Seed #{target_seed}...")
    global_max_t = 0
    for scenario in scenarios:
        _, _, t_final = preproc.get_fig2_freq_data(exp_num, scenario, target_seed)
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
    
    gs = gridspec.GridSpec(3, 5, figure=fig, wspace=0.35, hspace=0.45)
    
    print(f"--- Building Figure 2 Grid (Exp #{exp_num}) ---")
    
    for col, scenario in enumerate(scenarios):
        display_name = preproc.get_clinical_label(scenario)
        print(f"Processing column {col+1}/5: {display_name}...")
        
        ax_freq = fig.add_subplot(gs[0, col])
        ax_clust = fig.add_subplot(gs[1, col])
        ax_div = fig.add_subplot(gs[2, col])
        
        lf, cmap_df, _ = preproc.get_fig2_freq_data(exp_num, scenario, target_seed)
        _, clade_to_lineages, clade_meta_df, _ = preproc.get_fig2_clustered_data(exp_num, scenario, target_seed)
        dists = preproc.get_fig2_divergence_data(exp_num, scenario, shared_violin_seeds)
        
        plots.plot_fig2_lineage_freq(ax_freq, lf, cmap_df, global_max_t, title=display_name)
        plots.plot_fig2_clustered_freq(ax_clust, lf, clade_to_lineages, clade_meta_df, cmap_df, global_max_t)
        plots.plot_fig2_divergence_violin(ax_div, dists, color=palette[scenario])
        
        if col == 0:
            ax_freq.set_ylabel("Lineage Freq", fontsize=7)
            ax_clust.set_ylabel("Clade Freq", fontsize=7)
            ax_div.set_ylabel("Divergence Jump", fontsize=7)
        else:
            ax_freq.set_ylabel("")
            ax_clust.set_ylabel("")
            ax_div.set_ylabel("")
            
        ax_freq.set_xlabel("")
        ax_clust.set_xlabel("Time (days)")
        
    output_filename = f'Figure_2.{fmt}'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {output_filename}")

if __name__ == "__main__":
    args = parse_arguments()
    build_figure_2(args.exp_num, args.seed, args.format)
