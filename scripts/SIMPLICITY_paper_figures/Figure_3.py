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
from PIL import Image
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
from matplotlib.cm import get_cmap

from simplicity.tree import build_tree
import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm

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
        cmap = get_cmap("gist_rainbow")
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


def add_custom_legend(fig, axs1, axs2=None):
    """
    Adds a multi‐section legend beneath the figure, arranged in three columns:
      • Column 1: Abbreviations
      • Column 2: Subplot 1 + Subplot 2 entries
      • Column 3: Subplot 3 + Subplot 4 entries
    """
    # ── Build all the individual legend handles and labels ─────────────────────
    # Section 1: Abbreviations
    subtitle_abbrev = Line2D([], [], linestyle="None", label=r"$\bf{Abbreviations:}$")
    l_abbrev = [
        "  N.i.  = N individuals",
        "  R.a.t.f. = Relative avg. transmission fitness",
        "  L.r.f. = Lineage relative frequency",
        "  R.e.  = R effective",
    ]
    abbrev_handles = [subtitle_abbrev, Line2D([], [], linestyle="None")]
    abbrev_labels = [subtitle_abbrev.get_label(), " "]  # blank line for spacing
    
    for text in l_abbrev:
        abbrev_handles.append(Line2D([], [], linestyle="None"))
        abbrev_labels.append(text)

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
    subtitle3     = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 3\ –\ Lineage\ Frequency}$")
    rainbow_box   = RainbowLegendBox()
    sub3_handles  = [subtitle3, rainbow_box]
    sub3_labels   = [subtitle3.get_label(), "Lineage frequency"]

    # Section 5: Subplot 4 – R effective
    subtitle4    = Line2D([], [], linestyle="None", label=r"$\bf{Subplot\ 4\ –\ R\ effective}$")
    h4, l4       = axs1[3].get_legend_handles_labels()
    rainbow_line = RainbowLegendLine()
    sub4_handles = [subtitle4] + h4 + [rainbow_line]
    sub4_labels  = [subtitle4.get_label()] + l4 + ["R effective lineage"]

    # ── First Column: Abbreviations ─────────────────────────────────────────────
    # Place the "Abbreviations" block at x ≈ 0.16 (in figure fraction coordinates)
    fig.legend(
        abbrev_handles,
        abbrev_labels,
        loc='lower left',
        bbox_to_anchor=(0.16, -0.02),
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=0,       
        handletextpad=0.5,
    )

    # ── Second Column: Subplot 1 + Subplot 2 ────────────────────────────────────
    # Combine the two sub‐blocks into one list, keeping subtitles at top of each
    sub12_handles = sub1_handles + [Line2D([], [], linestyle="None", label=" ")] + sub2_handles
    sub12_labels  = sub1_labels  + [" "] + sub2_labels

    # Place that combined block at x ≈ 0.50
    fig.legend(
        sub12_handles,
        sub12_labels,
        loc='lower center',
        bbox_to_anchor=(0.50, -0.02),
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=1.0,
        handletextpad=0.5,
        labelspacing=0.4,
        handler_map={
            RainbowLegendBox: HandlerRainbowBox(),
            RainbowLegendLine: HandlerRainbowLine(),
        }
    )

    # ── Third Column: Subplot 3 + Subplot 4 ────────────────────────────────────
    sub34_handles = sub3_handles + [Line2D([], [], linestyle="None", label=" ")] + sub4_handles
    sub34_labels  = sub3_labels  + [" "] + sub4_labels

    # Place that combined block at x ≈ 0.84
    fig.legend(
        sub34_handles,
        sub34_labels,
        loc='lower right',
        bbox_to_anchor=(0.84, -0.02),
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=1.0,
        handletextpad=0.5,
        labelspacing=0.4,
        handler_map={
            RainbowLegendBox: HandlerRainbowBox(),
            RainbowLegendLine: HandlerRainbowLine(),
        }
    )

    # ── Make room at the bottom for all three legend columns ───────────────────
    fig.subplots_adjust(bottom=0.20)


def extended_simulation_trajectory(axs, ssod, experiment_name, threshold, label):
    axs[0].text(-0.1, 1.05, label, transform=axs[0].transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')
    # --- Subplot 1: System trajectory ----------------------------------------
    trajectory_df = om.read_simulation_trajectory(ssod)
    time = trajectory_df['time'].tolist()
    infected = trajectory_df['infected'].tolist()
    # diagnosed = trajectory_df['diagnosed'].tolist()

    axs[0].plot(time, infected, label='Infected individuals at time t', color='red')
    # axs[0].plot(time, susceptibles, label='Susceptibles individuals at time t', color='green')
    # axs[0].plot(time, diagnosed, label='Diagnosed individuals (cumulative)', color='yellow')
    axs[0].set_ylim(0, max(infected) * 1.5)
    axs[0].set_ylabel('N. i.')
    # axs[0].legend(loc='upper left')
    
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

    axs[1].plot(times, means, linestyle='-', color='#2ca02c', label='Relative avg. transmission fitness')
    axs[1].fill_between(times,
                        [m - s for m, s in zip(means, stds)],
                        [m + s for m, s in zip(means, stds)],
                        color='#2ca02c', alpha=0.3, label='Std.')
    # axs[1].legend(loc='upper left')
    axs[1].set_ylabel("R.a.t.f.")

    # --- Subplot 3: Lineages frequency ---------------------------------------
    try:
        lineage_frequency_df = om.read_lineage_frequency(ssod)
        colormap_df = pm.make_lineages_colormap(ssod, cmap_name='gist_rainbow')
        pm.plot_lineages_colors_tab(ssod)  # Save color key figure
    except Exception as e:
        print(f"[ERROR] Lineage frequency or colormap failed: {e}")
        return

    filtered_df, _ = om.filter_lineage_frequency_df(lineage_frequency_df, threshold)
    colors = [pm.get_lineage_color(name, colormap_df) for name in filtered_df.columns]
    filtered_df.plot(kind='area', stacked=False, color=colors, alpha=0.5, ax=axs[2], legend=False)

    time_file_path = os.path.join(ssod, 'final_time.csv')
    try:
        time_final = pd.read_csv(time_file_path, header=None).iloc[0, 0]
    except Exception as e:
        print(f"[ERROR] final_time.csv could not be read: {e}")
        time_final = max(time)
    axs[2].set_xlim([0, time_final])

    def to_percent(x, _): return f"{100 * x:.0f}%"
    axs[2].yaxis.set_major_formatter(FuncFormatter(to_percent))
    axs[2].set_ylabel("L.r.f.")

    # --- Subplot 4: R Effective by Lineage ----------------------------------
    try:
        colormap_df = pm.make_lineages_colormap(ssod)
        r_eff_pop, r_eff_lineages = om.read_r_effective_trajs_csv(experiment_name, ssod, time_window=21, threshold=threshold)
        lineage_freq_df = om.read_lineage_frequency(ssod)
        _, filtered_lineages = om.filter_lineage_frequency_df(lineage_freq_df, threshold)
    except Exception as e:
        print(f"[ERROR] Failed to process R Effective data: {e}")
        return

    for lineage in filtered_lineages:
        if lineage not in r_eff_lineages:
            continue
        r_series = r_eff_lineages[lineage]
        color = pm.get_lineage_color(lineage, colormap_df)
        axs[3].plot(r_series.index, r_series.values, color=color)

    axs[3].plot(r_eff_pop.index, r_eff_pop.values, color='black', label='Average R_effective', zorder=100)
    axs[3].set_ylabel("Re")
    axs[3].set_xlabel("Time (d)")
    
    
    for ax in axs:
        # ax.tick_params(labelbottom=True)
        pm.apply_standard_axis_style(ax)
    
    # set y lims --------------------------------------------------------------
    max_y1 = max(infected) * 1.5
    axs[0].set_ylim(0, max_y1)
    
    max_y2 = max([m + s for m, s in zip(means, stds)])
    axs[1].set_ylim(0, max_y2 * 1.2)
 
    max_y3 = 1.0  
    axs[2].set_ylim(0, max_y3)
    
    # # Collect max values from lineage series
    # lineage_max_values = []
    # for r in r_eff_lineages.values():
    #     if isinstance(r, pd.Series):
    #         rmax = r.max()
    #         if pd.notna(rmax) and np.isfinite(rmax):
    #             lineage_max_values.append(rmax)
    
    # # Add population-wide R_eff if valid
    # if isinstance(r_eff_pop, pd.Series):
    #     r_eff_pop_max = r_eff_pop.max()
    #     if pd.notna(r_eff_pop_max) and np.isfinite(r_eff_pop_max):
    #         lineage_max_values.append(r_eff_pop_max)
    
    # # Final max value with fallback
    # if lineage_max_values:
    #     max_y4 = max(lineage_max_values) * 1.2
    # else:
    #     print(f"[WARNING] No valid R_effective values found in {ssod}")
    #     max_y4 = 1.0  # default
    max_y4 = 12
    axs[3].set_ylim(0, max_y4 * 1.2)
    
    return [max_y1, max_y2, max_y3, max_y4]



def load_tree_to_ax(ax, ssod, experiment_name, label):
    try:
        # Get the path to the circular phylogenetic tree image
        image_path = om.get_tree_plot_filepath(
            experiment_name,
            ssod,
            tree_type='phylogenetic',
            tree_subtype='circular'
        )

        # Open and display it on the axis
        img = Image.open(image_path)
        ax.text(-0.1, 1.05, label, transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='left')
        ax.imshow(np.array(img))
        ax.axis('off')

    except Exception as e:
        print(f"[ERROR] Could not load tree image from {image_path}: {e}")
        ax.text(0.5, 0.5, "Tree Image Missing", ha='center', va='center')
        ax.axis('off')

def figure_3(experiment_name, seed, threshold):
    sim_dirs = dm.get_simulation_output_dirs(experiment_name)

    if len(sim_dirs) < 2:
        raise ValueError("Need at least 2 simulation output directories.")

    ssod1 = dm.get_ssod(sim_dirs[0], seed)
    ssod2 = dm.get_ssod(sim_dirs[1], seed)

    # Create figure
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(9, 7, figure=fig)

    # Top half (SSOD1)
    axs1 = [fig.add_subplot(gs[0, 0:4])]
    for i in range(1, 4):
        axs1.append(fig.add_subplot(gs[i, 0:4], sharex=axs1[0]))
    ylims1 = extended_simulation_trajectory(axs1, ssod1, experiment_name, threshold,
                                            label= 'A')
   
    # ------------------------------------------------------------------------- 
    ax_tree1 = fig.add_subplot(gs[0:4, 4:])
    build_tree.draw_tree(ssod1, experiment_name, 'phylogenetic')
    load_tree_to_ax(ax_tree1, ssod1, experiment_name, label='B')

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # Spacer row
    fig.add_subplot(gs[4, :]).axis('off')
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # Create axs2 subplots with shared x-axis
    axs2 = [fig.add_subplot(gs[5, 0:4])]
    for i in range(6, 9):
        axs2.append(fig.add_subplot(gs[i, 0:4], sharex=axs2[0]))
    ylims2 = extended_simulation_trajectory(axs2, ssod2, experiment_name, threshold,
                                            label= 'C')

    # set max y axis across corresponding plots
    for i in range(4):
        ymax = max(ylims1[i], ylims2[i])
        axs1[i].set_ylim(top=ymax)
        axs2[i].set_ylim(top=ymax)
        
    # --------------------------------------------------------------------------
    
    ax_tree2 = fig.add_subplot(gs[5:9, 4:])
    build_tree.draw_tree(ssod2, experiment_name, 'phylogenetic')
    load_tree_to_ax(ax_tree2, ssod2, experiment_name, label='D')
    
    # --------------------------------------------------------------------------
    fig.subplots_adjust(hspace=0.5, top=0.94, bottom=0.05)
    
    axs1[-1].tick_params(labelbottom=True)  
    axs1[3].text(0.5, -0.3, "Time (d)", transform=axs1[3].transAxes,
             ha='center', va='top' )
    
    add_custom_legend(fig, axs1, axs2)
    
    output_path = os.path.join("Data", f"Figure_3_{experiment_name}_seed{seed}.tiff")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    output_path = os.path.join("Data", "Figure_3.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    plt.close(fig)
    
def main():
    experiment_name = 'SIMPLICITY_exp_output'
    seed = 7
    threshold = 0.05
    figure_3(experiment_name, seed, threshold)

if __name__ == "__main__":
    main()
    