import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import fig1_preprocess_data as preproc
import fig1_plots as plots
from _scenarios import EXP_NAME, figure_path

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Figure 1")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number to plot. Required on purpose: the "
                            "old default of 4 pointed at one of the runs that "
                            "produced invalid science.")
    parser.add_argument('--exp-name', type=str, default=EXP_NAME,
                        help=f"Pipeline to plot (default: {EXP_NAME}).")
    parser.add_argument('--format', type=str, choices=['pdf', 'png'], default='png',
                        help='Output format for the figure (default: png)')
    parser.add_argument('--no-fit-stats', action='store_true',
                        help="Hide the fitted rate / R^2 annotations in panels D "
                            "and F. The submission version should use this.")
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

def build_figure_1(exp_num, exp_name, fmt):
    set_nature_rcparams()

    # Six panels: the duration story (A-C), then the global clock (D) with
    # its distribution beside it (E), and the intra-host clock (F).
    # D, E and F each carry BOTH the real cohort and the model, overlaid, so
    # the comparison is read off one axes instead of across a data/model pair.
    # D runs to the end of the cohort rather than 400 d: its fit always used
    # every sequence, and the 400 d window hid the 46% beyond it -- the ones
    # that stop accumulating, and that put the second mode in E.
    # 180 mm is the double-column width and is not negotiable, so the legend
    # card's strip on the right is taken out of the panels, not added to the
    # canvas. plots.LEGEND_X is where the card's left edge sits.
    fig = plt.figure(figsize=(180 / 25.4, 100 / 25.4))
    gs = gridspec.GridSpec(2, 6, figure=fig, wspace=1.5, hspace=0.8,
                           right=plots.LEGEND_X - 0.045)

    ax_A = fig.add_subplot(gs[0, 0:2])
    ax_B = fig.add_subplot(gs[0, 2:4])
    ax_C = fig.add_subplot(gs[0, 4:6])
    ax_D = fig.add_subplot(gs[1, 0:2])
    ax_E = fig.add_subplot(gs[1, 2:4])
    ax_F = fig.add_subplot(gs[1, 4:6])

    # === LOAD DATA ===
    df_ab_theoretical = preproc.get_panel_a_data(exp_num=exp_num, exp_name=exp_name)
    df_ab_realized = preproc.get_panel_b_data(exp_num=exp_num, exp_name=exp_name)
    df_c = preproc.get_panel_c_data()
    df_real_std, df_real_long = preproc.get_patient_clock_data()
    df_model_global = preproc.get_model_global_clock(exp_num=exp_num, exp_name=exp_name)
    df_model_ih = preproc.get_model_intrahost_clock(exp_num=exp_num, exp_name=exp_name)

    # === PLOT ===
    # Panel A: 800 d only when a scenario actually runs that long (edge_case).
    # The bound pipeline's longest is HIV at 109 d, so 800 left it mostly empty.
    has_edge_scen = 'Edge Case' in set(df_ab_theoretical.get('cohort', []))
    plots.plot_fig1_intra_host(ax_A, df_ab_theoretical,
                               x_max=800 if has_edge_scen else 400)
    plots.plot_fig1_violins(ax_B, df_ab_realized)
    plots.plot_infection_duration(ax_C, df_c)

    # D -- global clock: standard patients vs the model's control standards.
    # days_since_ref is already in the real dataset; the old panel D plotted
    # absolute dates, which cannot be overlaid on simulation time.
    plots.plot_clock_overlay(
        ax_D, df_real_std, df_model_global,
        real_x='days_since_ref', real_y='hamming_per_overlap',
        title="Global clock (standard)", data_label="Patient data",
        x_label="Days since first sample", x_max=750, y_max=0.004)

    # E -- D's two series with time collapsed. D shows the rates agree to
    # within 5%; E shows the distributions do not.
    plots.plot_divergence_density(
        ax_E, df_real_std, df_model_global,
        real_x='days_since_ref', real_y='hamming_per_overlap',
        title="Divergence distribution", data_label="Patient data",
        x_max=0.004)

    # F -- intra-host clock: long-shedder patients vs the model's long shedders.
    plots.plot_clock_overlay(
        ax_F, df_real_long, df_model_ih,
        real_x='time_point', real_y='hamming_distance_per_overlap',
        title="Intra-host clock (long shedders)", data_label="Patient data",
        x_label="Days since infection", x_max=400, y_max=0.0035)

    # === ANNOTATIONS ===
    for label, ax in {'A': ax_A, 'B': ax_B, 'C': ax_C, 'D': ax_D, 'E': ax_E,
                      'F': ax_F}.items():
        add_panel_label(ax, label)

    # The edge-case colour also marks panel C's >300 d patients, so the legend
    # entry belongs whenever EITHER appears -- on the bound pipeline there is no
    # edge_case scenario but panel C still shows those patients in pink.
    has_edge = ('Edge Case' in set(df_ab_realized.get('cohort', [])))
    if not df_c.empty:
        has_edge = has_edge or bool(
            (df_c['duration'] > plots.LONG_SHEDDER_CUTOFF_DAYS).any())
    plots.add_figure_legend(fig, has_edge_case=has_edge,
                            has_patient_data=not (df_real_std.empty and df_real_long.empty))

    output_filename = figure_path(1, exp_name, exp_num, fmt)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Successfully generated {output_filename}")


if __name__ == "__main__":
    args = parse_arguments()
    if args.no_fit_stats:
        plots.SHOW_FIT_STATS = False
    build_figure_1(args.exp_num, args.exp_name, args.format)
