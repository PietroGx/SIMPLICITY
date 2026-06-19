import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import nature_figures_preprocess_data as preproc
import nature_figures_plots as nf_plots

def set_nature_rcparams():
    """Enforces Nature Microbiology typography and style standards."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,        # Standard text size
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'pdf.fonttype': 42,    # Editable text in Illustrator/Inkscape
        'ps.fonttype': 42
    })

def add_panel_label(ax, label):
    """Adds bold 8pt panel letters."""
    ax.text(-0.1, 1.15, label, transform=ax.transAxes, 
            fontsize=8, fontweight='bold', va='top', ha='right')

def build_figure_1():
    set_nature_rcparams()
    
    # Nature dimensions: Max width 180mm. Let's use 180mm x 120mm
    width_in = 180 / 25.4
    height_in = 120 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    
    # Create main 2x3 grid
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.3, hspace=0.4)
    
    # Initialize main axes
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])
    ax_D = fig.add_subplot(gs[1, 0])
    ax_E = fig.add_subplot(gs[1, 1])
    
    # Create nested grid for F and G in the bottom right corner
    gs_FG = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 2], hspace=0.1)
    ax_F = fig.add_subplot(gs_FG[0, 0])
    ax_G = fig.add_subplot(gs_FG[1, 0])
    
    # Load Data (Currently empty stubs)
    df_ab = preproc.get_panel_ab_data()
    df_c = preproc.get_panel_c_data()
    df_d, df_e = preproc.get_panel_de_data()
    df_f, df_g = preproc.get_panel_fg_data()
    


    # Plot Data
    nf_plots.plot_baseline_dynamics(ax_A, df_ab, "Model Dynamics")
    nf_plots.plot_baseline_dynamics(ax_B, df_ab, "Model Evolution")
    nf_plots.plot_infection_duration(ax_C, df_c)
    
    # Panel D: Standard (Dates) - Anchor intercept to the first sample
    nf_plots.plot_tempest_regression(ax_D, df_d, "Real Data (Standard)", 
                                     x_col='sampling_date', y_col='hamming_per_overlap', 
                                     is_date=True, scatter_color='#3b528b',
                                     force_zero_intercept=False,
                                     force_intercept_at_min_x=True) # <--- NEW FLAG SET TO TRUE
    
    # Panel E: Long Shedder (Days since infection) - Anchor intercept to strictly 0
    nf_plots.plot_tempest_regression(ax_E, df_e, "Real Data (Long Shedder)", 
                                     x_col='time_point', y_col='hamming_distance_per_overlap', 
                                     is_date=False, scatter_color='#e41a1c', 
                                     force_zero_intercept=True,
                                     force_intercept_at_min_x=False)
    
    # Model Panels F & G (Placeholders for now, using defaults)
    nf_plots.plot_tempest_regression(ax_F, df_f, "Model (Standard)", scatter_color='#3b528b')
    nf_plots.plot_tempest_regression(ax_G, df_g, "Model (Long Shedder)", scatter_color='#e41a1c')



    # Add standardized panel labels (Uppercase for Nature defaults)
    axes_dict = {'A': ax_A, 'B': ax_B, 'C': ax_C, 'D': ax_D, 'E': ax_E, 'F': ax_F, 'G': ax_G}
    for label, ax in axes_dict.items():
        add_panel_label(ax, label)
        
    # Save the figure
    output_filename = 'Figure_1_Skeleton.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Successfully generated {output_filename}")

if __name__ == "__main__":
    build_figure_1()
