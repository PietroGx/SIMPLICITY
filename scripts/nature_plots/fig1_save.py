import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import fig1_preprocess_data as preproc
import fig1_plots as plots

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Figure 1")
    parser.add_argument('--format', type=str, choices=['pdf', 'png'], default='png',
                        help='Output format for the figure (default: png)')
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

def add_panel_label(ax, label, x_offset=-0.1, y_offset=1.15):
    ax.text(x_offset, y_offset, label, transform=ax.transAxes, 
            fontsize=8, fontweight='bold', va='top', ha='right')

def build_figure_1(fmt):
    set_nature_rcparams()
    
    width_in = 180 / 25.4
    height_in = 120 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    
    # Adjusted GridSpec with weight_ratios to favor regression panels
    gs = gridspec.GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 1.2], wspace=0.4, hspace=0.7)
    
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])
    ax_D = fig.add_subplot(gs[1, 0])
    ax_E = fig.add_subplot(gs[1, 1])
    
    # Increased hspace to add vertical space between Panels F and G
    gs_FG = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 2], hspace=1.0)
    ax_F = fig.add_subplot(gs_FG[0, 0])
    ax_G = fig.add_subplot(gs_FG[1, 0])
    
    # Enforce square aspect ratio ONLY for D and E
    for ax in [ax_D, ax_E]:
        ax.set_box_aspect(1)
    
    # === LOAD DATA ===
    TARGET_EXP_NUM = 4
    
    df_ab_theoretical = preproc.get_panel_a_data(exp_num=TARGET_EXP_NUM)
    df_ab_realized = preproc.get_panel_b_data(exp_num=TARGET_EXP_NUM)
    df_c = preproc.get_panel_c_data()
    df_d, df_e = preproc.get_panel_de_data()
    df_f, df_g = preproc.get_panel_fg_data(exp_num=TARGET_EXP_NUM)
    
    # === PLOT DATA ===
    plots.plot_fig1_intra_host(ax_A, df_ab_theoretical)
    plots.plot_fig1_violins(ax_B, df_ab_realized)
    plots.plot_infection_duration(ax_C, df_c)
    
    # Panel D: Standard (Dates)
    plots.plot_tempest_regression(ax_D, df_d, "Real Data (Standard)", 
                                     x_col='sampling_date', y_col='hamming_per_overlap', 
                                     is_date=True, scatter_color='#3b528b',
                                     force_zero_intercept=False, force_intercept_at_min_x=True)
    
    # Panel E: Long Shedder (Days since infection)
    plots.plot_tempest_regression(ax_E, df_e, "Real Data (Long Shedder)", 
                                     x_col='time_point', y_col='hamming_distance_per_overlap', 
                                     is_date=False, scatter_color='#e41a1c', 
                                     force_zero_intercept=True, force_intercept_at_min_x=False,
                                     x_label="Days since infection")
    
    # Panel F: Model Global Clock (Simulation Time)
    plots.plot_tempest_regression(ax_F, df_f, "Model (Standard)", scatter_color='#3b528b', 
                                     x_label="Simulation Time (days)")
    
    # Panel G: Model Intra-Host Clock (Days since infection)
    plots.plot_tempest_regression(ax_G, df_g, "Model (Long Shedder)", scatter_color='#e41a1c', 
                                     x_label="Days since infection")
    
    # Reduce font size for text annotations in regression boxes
    for ax in [ax_F, ax_G]:
        for text in ax.texts:
            text.set_fontsize(6)
    
    # === ANNOTATIONS ===
    axes_dict = {'A': ax_A, 'B': ax_B, 'C': ax_C, 'D': ax_D, 'E': ax_E, 'F': ax_F, 'G': ax_G}
    for label, ax in axes_dict.items():
        if label == 'F':
            # Sit over the y-axis label, raised so D/E/F letters align
            add_panel_label(ax, label, x_offset=-0.05, y_offset=1.4)
        elif label == 'G':
            # Sit over the y-axis label, raised above its own panel top
            add_panel_label(ax, label, x_offset=-0.05, y_offset=1.4)
        else:
            add_panel_label(ax, label)
        
    # Layout Adjustments
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9)
    
    output_filename = f'Figure_1.{fmt}'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Successfully generated {output_filename}")

if __name__ == "__main__":
    args = parse_arguments()
    build_figure_1(args.format)
