# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.image import BboxImage
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from matplotlib import colormaps
import argparse

from simplicity.phenotype.distance import hamming_iw
import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.clustering as cl

pm.apply_plos_rcparams()

class RainbowLegendLine:
    """Dummy handle for a rainbow‐gradient line in the legend."""
    pass

class HandlerRainbowLine(HandlerBase):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        """
        Draw a horizontal line that goes from left to right, 
        colored by a ‘gist_rainbow’ gradient.
        """
        # We will break the line into many small segments:
        n_segments = 100
        xs = np.linspace(xdescent, xdescent + width, n_segments)
        y_center = ydescent + height / 2.0

        # Build list of segment endpoints: [ [(x0,y), (x1,y)], [(x1,y),(x2,y)], ... ]
        segments = []
        for i in range(n_segments - 1):
            x0, x1 = xs[i], xs[i + 1]
            segments.append([(x0, y_center), (x1, y_center)])

        # Create a LineCollection from these segments
        # “array” controls the color along the line (0→1 over all segments)
        cmap = colormaps["gist_rainbow"]
        colors = np.linspace(0, 1, n_segments - 1)

        lc = LineCollection(
            segments,
            cmap=cmap,
            norm=Normalize(0, 1),
            linewidths=1.5,
            transform=trans,  
            zorder=2
        )
        lc.set_array(colors)

        return [lc]

class RainbowLegendBox:
    """A dummy handle for the rainbow‐gradient legend entry."""
    pass

class HandlerRainbowBox(HandlerBase):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        # Build a small horizontal gradient array
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        gradient = np.vstack([gradient] * 10)  # give it some vertical “thickness”

        # Define a bounding box in the legend’s coordinate system
        bbox = Bbox.from_bounds(xdescent, ydescent, width, height)
        tbbox = TransformedBbox(bbox, trans)

        # Create a BboxImage that will fill that box with the “gist_rainbow” colormap
        im = BboxImage(
            tbbox,
            cmap="gist_rainbow",
            interpolation="bilinear",
            norm=Normalize(0, 1),
            zorder=2
        )
        im.set_array(gradient)

        # Create a thin black‐border rectangle on top
        border = Rectangle(
            (xdescent, ydescent),
            width,
            height,
            transform=trans,
            edgecolor="black",
            facecolor="none",
            linewidth=0.5,
            zorder=3
        )

        return [im, border]

def add_custom_legend_right(fig, axs1, axs2=None):
    """
    Places all legend sections (abbreviations + subplot entries) vertically
    on the right-hand side of the figure.
    """

    # Section 1: Abbreviations
    subtitle_abbrev = Line2D([], [], linestyle="None", label=r"$\bf{Abbreviations:}$")
    l_abbrev = [
        "  N.i.  = N individuals",
        "  R.a.t.f. = Relative avg. transmission fitness",
        "  L.r.f. = Lineage relative frequency",
        "  R.e.  = R effective",
    ]
    abbrev_handles = [subtitle_abbrev] + [Line2D([], [], linestyle="None") for _ in l_abbrev]
    abbrev_labels  = [subtitle_abbrev.get_label()] + l_abbrev

    # Section 2: Subplot 1 – Infected trajectory
    subtitle1 = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 1\ –\ Infected\ trajectory}$")
    h1, l1 = axs1[0].get_legend_handles_labels()
    sub1_handles = [subtitle1] + h1
    sub1_labels  = [subtitle1.get_label()] + l1

    # Section 3: Subplot 2 – Relative Fitness
    subtitle2 = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 2\ –\ Relative\ Fitness}$")
    h2, l2 = axs1[1].get_legend_handles_labels()
    sub2_handles = [subtitle2] + h2
    sub2_labels  = [subtitle2.get_label()] + l2

    # Section 4: Subplot 3 – Lineage Frequency
    subtitle3    = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 3\ –\ Lineage\ Frequency}$")
    rainbow_box  = RainbowLegendBox()
    sub3_handles = [subtitle3, rainbow_box]
    sub3_labels  = [subtitle3.get_label(), "Lineage frequency"]

    # Section 5: Subplot 4 – Clustered Frequencies
    subtitle4     = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 4\ –\ Clustered\ Frequencies}$")
    cluster_line  = RainbowLegendLine()
    sub4_handles  = [subtitle4, cluster_line]
    sub4_labels   = [subtitle4.get_label(), "Cluster frequency"]

    # Combine all
    all_handles = abbrev_handles + [Line2D([], [], linestyle="None")] + \
                  sub1_handles + [Line2D([], [], linestyle="None")] + \
                  sub2_handles + [Line2D([], [], linestyle="None")] + \
                  sub3_handles + [Line2D([], [], linestyle="None")] + \
                  sub4_handles

    all_labels = abbrev_labels + [" "] + \
                 sub1_labels + [" "] + \
                 sub2_labels + [" "] + \
                 sub3_labels + [" "] + \
                 sub4_labels

    # Final right-side legend
    fig.legend(
        all_handles,
        all_labels,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=1.0,
        handletextpad=0.6,
        labelspacing=0.4,
        handler_map={
            RainbowLegendBox: HandlerRainbowBox(),
            RainbowLegendLine: HandlerRainbowLine(),
        }
    )

    # Leave more room on the right
    fig.subplots_adjust(right=0.9   )



def compute_transmission_distances(ssod):
    """
    Compute the Hamming genomic distance between the inherited lineage and each
    transmitted lineage for all individuals.

    Returns:
        pd.DataFrame with columns:
            - individual_id
            - inherited_lineage
            - transmitted_lineage
            - hamming_distance
    """
    
    # Load data
    individuals_data = om.read_individuals_data(ssod)
    lineage_df = om.read_phylogenetic_data(ssod)

    # Build lineage → genome map
    genome_map = dict(zip(lineage_df['Lineage_name'], lineage_df['Genome']))

    # data to collect
    data = []

    for individual_id, row in individuals_data.iterrows():
        inherited = row['inherited_lineage']
        new_infs = row['new_infections']
        ind_type = row['type']

        # Skip if no transmissions or inherited lineage missing
        if not new_infs or inherited is None:
            continue

        genome_inherited = genome_map.get(inherited)
        if genome_inherited is None:
            continue  # skip if missing genome

        for inf in new_infs:
            transmitted = inf.get('transmitted_lineage')
            if transmitted is None:
                continue

            genome_transmitted = genome_map.get(transmitted)
            if genome_transmitted is None:
                continue

            dist = hamming_iw(genome_inherited, genome_transmitted)

            data.append({
                'individual_id': individual_id,
                'inherited_lineage': inherited,
                'transmitted_lineage': transmitted,
                'hamming_distance': dist,
                'type': ind_type
            })
            
    return pd.DataFrame(data)

def plot_hamming_distance_boxplot_ax(ax, transmission_distance_data, label="C"):
    """
    Plot a boxplot of Hamming distances from transmission events on a given Axes,
    grouped by:
        - 'normal'
        - 'long'
        - 'normal+long' (combined)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on (as a subplot).
    df : pd.DataFrame
        Output of compute_transmission_distances(...)
    label : str
        Subplot label to annotate (e.g., "D")
    """
    import seaborn as sns
    # Build group column
    df_plot = transmission_distance_data.copy()
    df_plot['group'] = df_plot['type'].map({
        'normal': 'normal',
        'long_shedder': 'long'
    })

    # Plot to provided axis
    sns.boxplot(data=df_plot, x='group', y='hamming_distance',
                hue = 'group', palette='pastel', ax=ax)
    sns.stripplot(data=df_plot, x='group', y='hamming_distance',
                  color='black', alpha=0.3, size=4, jitter=True, ax=ax)

    ax.set_xlabel("Transmitting Individual Type")
    ax.set_ylabel("Hamming Distance")
    ax.set_title("")
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

    # Optional: remove top and right spines for aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
def plot_hamming_distance_violinplot_ax(ax, transmission_distance_data, label="C"):
    """
    Plot a violin plot of Hamming distances from transmission events on a given Axes,
    grouped by:
        - 'normal'
        - 'long'
        - 'normal+long' (combined)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on (as a subplot).
    transmission_distance_data : pd.DataFrame
        Output of compute_transmission_distances(...)
    label : str
        Subplot label to annotate (e.g., "D")
    """
    import seaborn as sns

    # Build group column
    df_plot = transmission_distance_data.copy()
    df_plot['group'] = df_plot['type'].map({
        'normal': 'normal',
        'long_shedder': 'long'
    })

    # Plot violin
    sns.violinplot(data=df_plot, x='group', y='hamming_distance',
                   hue='group', palette='pastel', ax=ax, inner=None)

    # Add individual data points
    # sns.stripplot(data=df_plot, x='group', y='hamming_distance',
    #               color='black', alpha=0.3, size=4, jitter=True, ax=ax)

    ax.set_xlabel("Transmitting Individual Type")
    ax.set_ylabel("Tras. Dist.")
    ax.set_title("")
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

    # Aesthetic cleanup
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)



def extended_simulation_trajectory(axs, ssod, experiment_name, min_freq_threshold, cluster_threshold, label):
    axs[0].text(-0.1, 1.05, label, transform=axs[0].transAxes,
                fontsize=16, fontweight='bold', va='top', ha='left')

    # --- Subplot 1: System trajectory ----------------------------------------
    trajectory_df = om.read_simulation_trajectory(ssod)
    time = trajectory_df['time'].tolist()
    infected = trajectory_df['infected'].tolist()

    axs[0].plot(time, infected, label='Infected individuals at time t', color='blue')
    try:
        long_shedders = trajectory_df['long_shedders'].tolist()
        ax2 = axs[0].twinx()
        ax2.plot(time, long_shedders, label='Long shedders at time t', color='orange')
    except:
        pass
    axs[0].set_ylim(0, max(infected) * 1.5)
    axs[0].set_ylabel('N. i.')

    # --- Subplot 2: Fitness trajectory ---------------------------------------
    fitness_trajectory_file_path = os.path.join(ssod, 'fitness_trajectory.csv')
    try:
        fitness_trajectory_df = pd.read_csv(fitness_trajectory_file_path)
    except Exception as e:
        print(f"[ERROR] Reading fitness_trajectory.csv failed: {e}")
        return

    times = fitness_trajectory_df['Time'].tolist()
    means = fitness_trajectory_df['Mean'].tolist()
    stds = fitness_trajectory_df['Std'].tolist()
    entropy = fitness_trajectory_df['Entropy'].tolist()

    axs[1].plot(times, means, linestyle='-', color='#2ca02c', label='Relative avg. transmission fitness')
    axs[1].fill_between(times,
                        [m - s for m, s in zip(means, stds)],
                        [m + s for m, s in zip(means, stds)],
                        color='#2ca02c', alpha=0.3)
    ax3 = axs[1].twinx()
    ax3.plot(times, entropy, linestyle='-', color='purple', label='Transmission fitness entropy')
    
    axs[1].set_ylabel("R.a.t.f.")

    # --- Subplot 3: Lineages frequency ---------------------------------------
    try:
        lineage_frequency_df = om.read_lineage_frequency(ssod)
        colormap_df = pm.make_lineages_colormap(ssod, cmap_name='gist_rainbow')
        pm.plot_lineages_colors_tab(ssod)
    except Exception as e:
        print(f"[ERROR] Lineage frequency or colormap failed: {e}")
        return

    filtered_df, _ = om.filter_lineage_frequency_df(lineage_frequency_df, min_freq_threshold)
    colors = [pm.get_lineage_color(name, colormap_df) for name in filtered_df.columns]
    filtered_df.plot(kind='area', stacked=False, color=colors, alpha=0.5, ax=axs[2], legend=False)

    time_file_path = os.path.join(ssod, 'final_time.csv')
    try:
        time_final = pd.read_csv(time_file_path, header=None).iloc[0, 0]
    except Exception as e:
        print(f"[ERROR] final_time.csv could not be read: {e}")
        time_final = max(time)
    axs[2].set_xlim([0, time_final])
    axs[2].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{100 * x:.0f}%"))
    axs[2].set_ylabel("L.r.f.")

    # --- Subplot 4: Clustered frequency ---------------------------------------
    try:
    
        # Read phylogenetic info
        phylo_df = om.read_phylogenetic_data(ssod)
    
        # Parse genome into mutation sets
        lineage_to_mutations = cl.build_lineage_to_mutation_dict(phylo_df)
        
        # Pivot frequency dataframe
        full_pivot_df = lineage_frequency_df.pivot(
            index="Time_sampling",
            columns="Lineage_name",
            values="Frequency_at_t"
        )

    
        # Apply clustering
        clustered_df, cluster_parents = cl.build_clustered_freqs(
            lineage_to_mutations,
            full_pivot_df,
            cluster_threshold
        )

        # Plot
        colors = [pm.get_lineage_color(parent, colormap_df) for parent in cluster_parents]
    
        clustered_df.plot(kind='area', stacked=True,
                          color=colors, alpha=0.6, ax=axs[3], legend=False)
        axs[3].set_xlim([0, time_final])
        axs[3].set_ylim(0, 1.0)
        axs[3].set_ylabel("Cluster freq.")
        
    except:
        import traceback
        print("[ERROR] Clustered frequency plot failed:")
        traceback.print_exc()

    
    
    # Apply consistent styling
    for ax in axs:
        pm.apply_standard_axis_style(ax)

    # y-axis limits
    max_y1 = max(infected) * 1.5
    max_y2 = max([m + s for m, s in zip(means, stds)]) * 1.2
    max_y3 = 1.0
    max_y4 = 1.0  # same for clustered frequency

    return [max_y1, max_y2, max_y3, max_y4]

def plot_secondary_infections_boxplot_ax(ax, ssod, label="D"):
    """
    Plot a boxplot of the number of secondary infections (offspring)
    per individual, grouped by type (normal vs long_shedder).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on (as a subplot).
    ssod : str
        Path to seeded simulation output directory.
    label : str
        Subplot label (e.g., "D").
    """
    import seaborn as sns
    import simplicity.output_manager as om

    df = om.read_individuals_data(ssod)
    df['group'] = df['type'].map({
        'normal': 'normal',
        'long_shedder': 'long'
    })

    # Compute number of infections caused
    df['n_infections'] = df['new_infections'].apply(len)
    
    # Get top 3
    top_infectors = df.sort_values('n_infections', ascending=False).head(3)
    
    # Print result
    for idx, row in top_infectors.iterrows():
        print(f"Index: {idx}, Type: {row['type']}, Infections caused: {row['n_infections']}")
        
    # Plot
    sns.boxplot(data=df, x='group', y='n_infections',
                palette='pastel', hue = 'group', ax=ax)
    # sns.stripplot(data=df, x='group', y='n_infections',
    #               color='black', alpha=0.3, size=4, jitter=True, ax=ax)

    ax.set_xlabel("Transmitting Individual Type")
    ax.set_ylabel("Inf./ind.")
    ax.set_title("")
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# def plot(experiment_name, seed, min_freq_threshold, cluster_threshold):
#     sim_dirs = dm.get_simulation_output_dirs(experiment_name)

#     if len(sim_dirs) < 2:
#         raise ValueError("Need at least 2 simulation output directories.")

#     # Get long_shedders_ratio for each sim_dir and sort
#     sim_dirs_sorted = sorted(sim_dirs, key=lambda d: sm.get_parameter_value_from_simulation_output_dir(d, 'long_shedders_ratio'))

#     ssod1 = dm.get_ssod(sim_dirs_sorted[0], seed)  # lower value → left
#     ssod2 = dm.get_ssod(sim_dirs_sorted[1], seed)  # higher value → right

#     # ── Setup figure layout ───────────────────────────────────────────────
#     fig = plt.figure(figsize=(16, 12))
#     gs = GridSpec(6, 7, figure=fig)

#     # Top left: ssod1 simulation outputs
#     axs1 = [fig.add_subplot(gs[0, 0:3])]
#     for i in range(1, 4):
#         axs1.append(fig.add_subplot(gs[i, 0:3], sharex=axs1[0]))
#     ylims1 = extended_simulation_trajectory(axs1, ssod1, experiment_name, min_freq_threshold, cluster_threshold, label='A')

#     # Top right: ssod2 simulation outputs
#     axs2 = [fig.add_subplot(gs[0, 4:7])]
#     for i in range(1, 4):
#         axs2.append(fig.add_subplot(gs[i, 4:7], sharex=axs2[0]))
#     ylims2 = extended_simulation_trajectory(axs2, ssod2, experiment_name, min_freq_threshold, cluster_threshold, label='B')

#     # Match y-axis limits across the two sides
#     for i in range(3):
#         ymax = max(ylims1[i], ylims2[i])
#         axs1[i].set_ylim(top=ymax)
#         axs2[i].set_ylim(top=ymax)

#     # ── Bottom row: transmission distance boxplot ─────────────────────────
#     axs3 = fig.add_subplot(gs[4, :])
#     transmission_distance_data = compute_transmission_distances(ssod2)
#     plot_hamming_distance_violinplot_ax(axs3, transmission_distance_data)

#     # ── Adjust layout and finalize ────────────────────────────────────────
#     fig.subplots_adjust(hspace=0.5, top=0.94, bottom=0.08)

#     axs1[-1].tick_params(labelbottom=True)
#     axs1[2].text(0.5, -0.3, "Time (d)", transform=axs1[2].transAxes,
#                  ha='center', va='top')
#     axs2[2].text(0.5, -0.3, "Time (d)", transform=axs2[2].transAxes,
#                  ha='center', va='top')

#     # Add combined legend using top panel (axs1)
#     add_custom_legend_right(fig, axs1, axs2)

#     # Save figure
#     experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
#     output_path = os.path.join(experiment_plots_dir, f"Figure_3_long_{experiment_name}_seed{seed}.tiff")
#     plt.savefig(output_path, format='tiff',  dpi=300, bbox_inches='tight')
#     print(f"Figure saved to: {output_path}")
#     plt.close(fig)

def plot(experiment_name, seed, min_freq_threshold, cluster_threshold):
    sim_dirs = dm.get_simulation_output_dirs(experiment_name)

    if len(sim_dirs) < 2:
        raise ValueError("Need at least 2 simulation output directories.")

    # Sort sim_dirs by long_shedders_ratio
    sim_dirs_sorted = sorted(sim_dirs, key=lambda d: sm.get_parameter_value_from_simulation_output_dir(d, 'long_shedders_ratio'))
    ssod1 = dm.get_ssod(sim_dirs_sorted[0], seed)
    ssod2 = dm.get_ssod(sim_dirs_sorted[1], seed)

    # ── Setup figure layout ───────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(6, 7, figure=fig)

    # Top left: ssod1 simulation outputs
    axs1 = [fig.add_subplot(gs[0, 0:3])]
    for i in range(1, 4):
        axs1.append(fig.add_subplot(gs[i, 0:3], sharex=axs1[0]))
    ylims1 = extended_simulation_trajectory(axs1, ssod1, experiment_name, min_freq_threshold, cluster_threshold, label='A')

    # Top right: ssod2 simulation outputs
    axs2 = [fig.add_subplot(gs[0, 4:7])]
    for i in range(1, 4):
        axs2.append(fig.add_subplot(gs[i, 4:7], sharex=axs2[0]))
    ylims2 = extended_simulation_trajectory(axs2, ssod2, experiment_name, min_freq_threshold, cluster_threshold, label='B')

    # Match y-axis limits
    for i in range(3):
        ymax = max(ylims1[i], ylims2[i])
        axs1[i].set_ylim(top=ymax)
        axs2[i].set_ylim(top=ymax)

    # ── Bottom row: violin plot (left) and placeholder (right) ─────────────
    ax_violin = fig.add_subplot(gs[4:6, 0:3])  
    ax_infections    = fig.add_subplot(gs[4:6, 4:7])  

    transmission_distance_data = compute_transmission_distances(ssod2)
    plot_hamming_distance_violinplot_ax(ax_violin, transmission_distance_data)
    
    plot_secondary_infections_boxplot_ax(ax_infections, ssod2, label="D")
    
    # ── Adjust layout and finalize ────────────────────────────────────────
    fig.subplots_adjust(hspace=0.5, wspace=0.3, top=0.94, bottom=0.08)

    axs1[-1].tick_params(labelbottom=True)
    axs1[2].text(0.5, -0.3, "Time (d)", transform=axs1[2].transAxes, ha='center', va='top')
    axs2[2].text(0.5, -0.3, "Time (d)", transform=axs2[2].transAxes, ha='center', va='top')

    add_custom_legend_right(fig, axs1, axs2)

    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    output_path = os.path.join(experiment_plots_dir, f"Figure_3_long_{experiment_name}_seed{seed}.tiff")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    plt.close(fig)


def main():
    
    # Set up the argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('max_seed', type=int, help="max_seed")
    parser.add_argument('cluster_threshold', type=int, help="cluster_threshold")
    
    args = parser.parse_args()
    

    min_freq_threshold = 0.05
    for seed in range (0,args.max_seed):
        plot(args.experiment_name, seed, min_freq_threshold, args.cluster_threshold)

if __name__ == "__main__":
    main()
    