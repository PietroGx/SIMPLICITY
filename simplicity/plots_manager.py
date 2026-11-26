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
import matplotlib
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
# matplotlib.use('Agg')
import anytree
from anytree.exporter import DotExporter
from ete3 import NodeStyle, TreeStyle, faces, AttrFace
import matplotlib.pyplot as plt
import seaborn           as sns
import pandas            as pd
import numpy             as np
import simplicity.dir_manager      as dm
import simplicity.settings_manager as sm
import simplicity.output_manager   as om

def apply_plos_rcparams():
    text_size = 16
    plt.rcParams.update({
        'savefig.dpi': 300,
        'font.size': text_size,
        'axes.labelsize': text_size,
        'axes.titlesize': text_size,
        'xtick.labelsize': text_size,
        'ytick.labelsize': text_size,
        'legend.fontsize': text_size,
#        'font.family': 'sans-serif',
#        'font.sans-serif': ['Arial'],  
        'pdf.fonttype': 42,            
        'ps.fonttype': 42,
        'figure.dpi': 300
    })

apply_plos_rcparams()

def apply_standard_axis_style(ax, has_secondary_y=False):
    """
    Apply PLOS-compliant style to an axis:
    - 8pt font (inherited via rcParams)
    - No top/right spines unless `has_secondary_y` is True
    """
    ax.spines['top'].set_visible(False)
    if not has_secondary_y:
        ax.spines['right'].set_visible(False)
        
    text_size = 16
    # Ensure ticks and labels are 8pt (in case rcParams missed it)
    ax.tick_params(labelsize=text_size)
    ax.title.set_fontsize(text_size)
    ax.xaxis.label.set_fontsize(text_size)
    ax.yaxis.label.set_fontsize(text_size)


# -----------------------------------------------------------------------------
#                        Simplicity simulation plots
# -----------------------------------------------------------------------------


    
def plot_trajectory(seeded_simulation_output_dir):
    """
    Plot the population trajectory of a SARS-CoV-2 simulation.

    This figure shows the evolution of three compartments over time:
    - Susceptibles
    - Infected individuals
    - Diagnosed individuals (cumulative)

    Parameters
    ----------
    seeded_simulation_output_dir : str
        Path to the seeded simulation output directory.
    """

    # Load simulation trajectory data
    simulation_trajectory = om.read_simulation_trajectory(seeded_simulation_output_dir)
    time = simulation_trajectory['time'].tolist()
    infected = simulation_trajectory['infected']
    diagnosed = simulation_trajectory['diagnosed']
    susceptibles = simulation_trajectory['susceptibles']

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot population compartments
    ax.plot(time, susceptibles, label='Susceptibles', color='black')
    ax.plot(time, infected, label='Infected individuals', color='red')
    ax.plot(time, diagnosed, label='Diagnosed individuals (cumulative)', color='green')

    # Set labels, limits, and legend
    ax.set_title('Simulation of SARS-CoV-2 outbreak')
    ax.set_xlabel('Time - days')
    ax.set_ylabel('Number of individuals in compartment')
    ax.set_xlim(0, time[-1])
    ax.set_ylim(0)
    ax.legend()
    plt.tight_layout()

    # Apply axis styling
    apply_standard_axis_style(ax)

    # Determine save path
    experiment_name = dm.get_experiment_foldername_from_SSOD(seeded_simulation_output_dir)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    seed = os.path.basename(seeded_simulation_output_dir)
    figure_output_path = os.path.join(experiment_plots_dir, f"{experiment_name}_{seed}_trajectory.tiff")

    # Save the plot 
    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_simulation(seeded_simulation_output_dir, threshold):
    """
    Plot a combined simulation summary figure.

    This multi-panel figure includes:
    1. System trajectory: Number of infected individuals over time.
    2. Lineage frequencies over time (filtered by threshold).
    3. Average fitness score with standard deviation.

    Parameters
    ----------
    seeded_simulation_output_dir : str
        Path to the seeded simulation output directory.
    threshold : float
        Relative frequency threshold for filtering lineages (e.g., 0.05 = 5%).
    """
    # Set up figure and layout
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(3, 1, height_ratios=[1, 1, 1])  
    ax2 = fig.add_subplot(gs[0, 0])  
    ax1 = fig.add_subplot(gs[2, 0], sharex=ax2)
    ax3 = fig.add_subplot(gs[1, 0], sharex=ax2)

    fig.subplots_adjust(
        top=0.95, bottom=0.085,
        left=0.07, right=0.81,
        hspace=0.1, wspace=0.085
    )

    # --- Subplot 1: System trajectory ----------------------------------------
    trajectory_df = om.read_simulation_trajectory(seeded_simulation_output_dir)
    time = trajectory_df['time'].tolist()
    infected = trajectory_df['infected'].tolist()

    ax2.plot(time, infected, label='Infected individuals at time t', color='red')
    ax2.set_ylim(0, max(infected) * 1.5)
    ax2.set_ylabel('N individuals')
    ax2.legend(loc='upper left')
    apply_standard_axis_style(ax2)

    # --- Subplot 2: Lineages frequency ---------------------------------------
    lineage_frequency_df = om.read_lineage_frequency(seeded_simulation_output_dir)
    colormap_df = make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')
    plot_lineages_colors_tab(seeded_simulation_output_dir)  # Save color key figure

    filtered_df, _ = om.filter_lineage_frequency_df(lineage_frequency_df, threshold)
    colors = [get_lineage_color(name, colormap_df) for name in filtered_df.columns]
    filtered_df.plot(kind='area', stacked=True, color=colors, alpha=0.5, ax=ax1)

    time_file_path = os.path.join(seeded_simulation_output_dir, 'final_time.csv')
    time_final = pd.read_csv(time_file_path, header=None).iloc[0, 0]
    ax1.set_xlim([0, time_final])

    def to_percent(x, _): return f"{100 * x:.0f}%"
    ax1.yaxis.set_major_formatter(FuncFormatter(to_percent))
    ax1.set_ylabel("Relative frequency of lineage")
    ax1.legend().remove()
    apply_standard_axis_style(ax1)

    # --- Subplot 3: Fitness trajectory ---------------------------------------
    fitness_trajectory_file_path = os.path.join(seeded_simulation_output_dir, 'fitness_trajectory.csv')
    fitness_trajectory_df = pd.read_csv(fitness_trajectory_file_path)

    times = fitness_trajectory_df['Time'].tolist()
    means = fitness_trajectory_df['Mean'].tolist()
    stds = fitness_trajectory_df['Std'].tolist()

    ax3.plot(times, means, linestyle='-', color='#2ca02c', label='Average fitness score')
    ax3.fill_between(times,
                     [m - s for m, s in zip(means, stds)],
                     [m + s for m, s in zip(means, stds)],
                     color='#2ca02c', alpha=0.3, label='Fitness score std')

    ax3.set_xlabel('Time - days')
    ax3.set_ylabel('Fitness score')
    ax3.set_xticks(np.arange(0, time_final, 50))
    ax3.set_xticklabels(np.arange(0, time_final, 50).astype(int))
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    ax3.legend(loc='upper left')
    apply_standard_axis_style(ax3)

    # --- Final formatting and saving -----------------------------------------
    fig.align_ylabels([ax2, ax3])
    plt.tight_layout()

    experiment_name = dm.get_experiment_foldername_from_SSOD(seeded_simulation_output_dir)
    so_foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    experiment_simulations_plots_dir = dm.get_experiment_simulations_plots_dir(experiment_name)
    figure_output_path = os.path.join(
        experiment_simulations_plots_dir,
        f"{so_foldername}_{seed}_simulation_trajectory.tiff"
    )

    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close(fig)


# -----------------------------------------------------------------------------
#                       Tempest regression and OSR plots
# -----------------------------------------------------------------------------


def plot_tempest_regression(sequencing_data_df,
                            fitted_tempest_regression,
                            ax):
    """
    Plot the tempest regression for observed substitution rates (OSR).

    This plot shows the linear regression of root-to-tip distance over time,
    with the observed substitution rate (OSR) shown as the slope. The OSR value
    is also included in the legend.

    Parameters
    ----------
    sequencing_data_df : pandas.DataFrame
        DataFrame with 'Sequencing_time' and 'Distance_from_root' columns.

    fitted_tempest_regression : sklearn.linear_model.LinearRegression
        Fitted regression model for OSR estimation.

    ax : matplotlib.axes.Axes
        The matplotlib axis object to draw the plot on.
    """

    # Get linear regression prediction
    x = sequencing_data_df['Sequencing_time'].values.reshape(-1, 1)
    y_pred = fitted_tempest_regression.predict(x)
    osr = fitted_tempest_regression.coef_[0]

    # Plot observed points
    sequencing_data_df.plot(
        kind='scatter',
        x='Sequencing_time',
        y='Distance_from_root',
        color='blue',
        ax=ax
    )

    # Plot regression line
    ax.plot(x, y_pred, color='orange', linewidth=2, label='y = OSR · x')

    # Axis labels and limits
    ax.set_xlabel('Simulation time (y)')
    ax.set_ylabel('Genetic distance from root (ssy)') #(substitutions/site, normed)
    ax.set_xlim(left=0, right=3)
    ax.set_ylim(0)

    # Add OSR value to legend
    extra_text = f'OSR = {osr:.5f}'
    handles, labels = ax.get_legend_handles_labels()
    handles.append(plt.Line2D([0], [0], color='white', label=extra_text))
    ax.legend(handles=handles)

    # Apply standard styling
    apply_standard_axis_style(ax)
    


def plot_OSR_fit(experiment_name, 
                 fit_result, 
                 model_type,
                 min_seq_number,
                 min_sim_lenght, group=None):
    """
    Plot nucleotide substitution rate vs. observed substitution rate (OSR)
    including three panel views (linear, semilog, and log-log).

    Each panel shows:
    - OSR estimates (single simulations)
    - Fitted curve
    - Combined regression estimates
    - Mean OSR values

    Parameters
    ----------
    experiment_name : str
        Name of the experiment, used to load and save data.

    fit_result : lmfit.model.ModelResult
        Fitted model result containing `.best_fit`.

    model_type : str
        Model used (used in output file name).

    min_seq_number : int
        Minimum number of sequences required for inclusion.

    min_sim_lenght : int
        Minimum simulation duration required for inclusion.
    """
    parameter = 'nucleotide_substitution_rate'

    line_color = 'black'
    scatter_color = '#DE8F05'            # Orange
    combined_OSR_marker_color = '#E64B9D'  # Pink
    scatter_color_2 = '#0173B2'          # Blue

    fig, axes = plt.subplots(3, 1, figsize=(8, 12))
    
    # import combined regression data
    combined_data = om.read_combined_OSR_vs_parameter_csv(
        experiment_name, parameter, min_seq_number, min_sim_lenght, group=group)
    # import single simulations regression data
    data = om.read_OSR_vs_parameter_csv(
        experiment_name, parameter, min_seq_number, min_sim_lenght, group=group)
    # Group by nucleotide_substitution_rate and compute mean and standard deviation for OSR
    data_mean_df = om.get_mean_std_OSR(
        experiment_name, parameter, min_seq_number, min_sim_lenght, group=group)

    x_data = data['nucleotide_substitution_rate']

    for ax in axes:
        sns.scatterplot(
            x=parameter, y='observed_substitution_rate',
            label='Estimated OSR - single simulation', data=data,
            color=scatter_color, alpha=0.5, ax=ax, zorder=0
        )
        sns.lineplot(
            x=x_data, y=fit_result.best_fit,
            label=f'Fitted {model_type} curve',
            color=line_color, linewidth=1, ax=ax, zorder=1
        )
        sns.scatterplot(
            x='nucleotide_substitution_rate', y='observed_substitution_rate', marker='X',
            label='Combined tempest regression estimate of OSR', data=combined_data,
            color=combined_OSR_marker_color, alpha=1, ax=ax, zorder=2
        )
        sns.scatterplot(
            x=parameter, y='mean', marker='X',
            label='Mean of estimated OSR (single simulations)', data=data_mean_df,
            color=scatter_color_2, alpha=1, ax=ax, zorder=3
        )
        ax.set_ylim(0.000009, 0.02)
        apply_standard_axis_style(ax)

    # Linear scale
    axes[0].set_xlabel("Nucleotide substitution rate (site/year)")
    axes[0].set_ylabel("Observed Substitution Rate (site/year)")

    # Semi-log x-axis
    axes[1].set_xlabel("Nucleotide substitution rate (site/year)")
    axes[1].set_ylabel("Observed Substitution Rate (site/year)")
    axes[1].set_xscale("log")
    axes[1].legend_.remove()

    # Log-log
    axes[2].set_xlabel("Nucleotide substitution rate (site/year)")
    axes[2].set_ylabel("Observed Substitution Rate (site/year)")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].legend_.remove()

    plt.tight_layout()

    figure_output_path = os.path.join(
        dm.get_experiment_plots_dir(experiment_name),
        f"{experiment_name}_observed_substitution_rates_{model_type}_fit.tiff"
    )
    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close(fig)

# -----------------------------------------------------------------------------
#                          Phenotype model plots
# -----------------------------------------------------------------------------
    

# -----------------------------------------------------------------------------
#                             Intra host model plots
# -----------------------------------------------------------------------------

    

    
# -----------------------------------------------------------------------------
#                           Runtime analysis plot
# -----------------------------------------------------------------------------

    
# -----------------------------------------------------------------------------
#                          Diagnosis rate plots
# -----------------------------------------------------------------------------
    


# -----------------------------------------------------------------------------
#                      Intra host lineage distribution plot
# -----------------------------------------------------------------------------



def plot_IH_lineage_distribution_grouped_by_simulation(experiment_name, group=None):
    """
    Plot intra-host lineage distributions grouped by tau_3 parameter across simulations.

    Each bar represents the proportion of individuals with a given number of distinct
    intra-host lineages, grouped and color-coded by `tau_3` values.

    Parameters
    ----------
    experiment_name : str
        Name of the experiment used to retrieve simulation data and save the plot.
    """
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name, group)
    plot_output_dir = dm.get_experiment_plots_dir(experiment_name)

    all_data = []

    for sim_dir in simulation_output_dirs:
        df = om.get_IH_lineages_data_simulation(sim_dir)[['IH_unique_lineages_number', 'ih_lineage_count']]
        tau_3 = sm.get_parameter_value_from_simulation_output_dir(sim_dir, 'tau_3')
        df['tau_3'] = tau_3

        total = df['ih_lineage_count'].sum()
        df['ih_lineage_count'] = df['ih_lineage_count'] / total if total > 0 else 0

        all_data.append(df)

    full_df = pd.concat(all_data, axis=0)
    full_df['IH_unique_lineages_number'] = pd.Categorical(
        full_df['IH_unique_lineages_number'],
        categories=[1, 2, 3, 4, 5],
        ordered=True
    )

    tau_3_order = sorted(full_df['tau_3'].unique())

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        x='IH_unique_lineages_number',
        y='ih_lineage_count',
        hue='tau_3',
        data=full_df,
        palette=sns.color_palette("coolwarm", n_colors=len(tau_3_order)),
        hue_order=tau_3_order,
        dodge=True,
        alpha=0.8,
        ax=ax
    )

    ax.set_xlabel("Intra host distinct lineages number")
    ax.set_ylabel("Proportion of individuals with x distinct lineages")
    ax.legend(title='tau_3', bbox_to_anchor=(1.05, 1), loc='upper left')

    apply_standard_axis_style(ax)
    plt.tight_layout()

    output_path = os.path.join(plot_output_dir, f'{experiment_name}_IH_unique_lineages_by_tau_3.tiff')
    print(f"Saving figure to: {output_path}")
    plt.savefig(output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close(fig)

# -----------------------------------------------------------------------------
#                       Simulations final times histogram
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#                              Lineages colors 
# -----------------------------------------------------------------------------

def make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow'):
    """
      Generate a colormap for lineages based on their order of emergence.
    
      This function reads phylogenetic data for a seeded simulation and assigns
      each lineage a unique color from a specified matplotlib colormap, ordered
      by time of emergence. The output is a DataFrame that maps lineage names
      to hexadecimal color codes.
    
      Parameters
      ----------
      seeded_simulation_output_dir : str
          Path to the directory containing simulation output and phylogenetic data.
      cmap_name : str, optional
          Name of the matplotlib colormap to use for color assignment (default is 'gist_rainbow').
    
      Returns
      -------
      pd.DataFrame
          DataFrame with columns:
          - 'Lineage_name': str, lineage identifiers
          - 'Color': str, hexadecimal color codes
    """
    def lineage_to_colormap(index, lineages, cmap_name):
        """Map a lineage to a unique color in a gradient colormap."""
        cmap = plt.get_cmap(cmap_name)
        rgb =  cmap(index / lineages)[:3]  # Ensure it returns an RGB tuple without alpha
        return mcolors.rgb2hex(rgb) #convert to HEX format
    
    df = om.read_phylogenetic_data(seeded_simulation_output_dir)
    # Sort lineages by time of emergence
    sorted_lineages = df.sort_values(by='Time_emergence')['Lineage_name'].tolist()
    # Generate colors
    colors = [lineage_to_colormap(i, len(sorted_lineages), cmap_name) for i in range(len(sorted_lineages))]

    # Create a new DataFrame mapping lineages to colors
    lineage_color_mapping_df = pd.DataFrame({'Lineage_name': sorted_lineages, 'Color': colors})

    return lineage_color_mapping_df

def get_lineage_color(lineage_name, colormap_df, cmap_name='gist_rainbow'):
    """
    Retrieve the color assigned to a specific lineage.

    Given a lineage name and a lineage-to-color mapping DataFrame (as generated
    by `make_lineages_colormap`), this function returns the corresponding color
    as a hexadecimal string.

    Parameters
    ----------
    lineage_name : str
        Name of the lineage to retrieve the color for.
    colormap_df : pd.DataFrame
        DataFrame containing 'Lineage_name' and 'Color' columns, as returned by `make_lineages_colormap`.
    cmap_name : str, optional
        Unused in this function (retained for compatibility), default is 'gist_rainbow'.

    Returns
    -------
    str
        Hexadecimal color code (e.g., '#aabbcc') corresponding to the given lineage.

    Raises
    ------
    ValueError
        If `lineage_name` is None or not found in the colormap.
    """
    df = colormap_df
    if lineage_name is None:
        raise ValueError('Lineage cannot be NoneType!')
    elif lineage_name not in df['Lineage_name'].values:
        raise ValueError(f'Lineage {lineage_name} not in lineages list! Cannot compute color.')
        
    lineage_color = df.loc[df['Lineage_name'] == lineage_name, 'Color'].values[0]
    return lineage_color  

def plot_lineages_colors_tab(seeded_simulation_output_dir):
    """
    Plot a labeled color reference for all lineages in the simulation.

    Each lineage is represented by a colored circle and its label, displayed in a
    vertical table format. Colors are assigned based on the order of emergence using
    a specified colormap.

    The figure is saved to the experiment's simulation plots directory using the
    seeded simulation folder and seed as part of the filename.

    Parameters
    ----------
    seeded_simulation_output_dir : str
        Path to the seeded simulation output directory. Used to read lineage data
        and determine the output save location.
    """
    colormap_df = make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')

    fig, ax = plt.subplots(num=2, figsize=(5, len(colormap_df) * 0.5))

    for i, lineage_name in enumerate(reversed(colormap_df['Lineage_name'].values)):
        ax.scatter(0, i, color=get_lineage_color(lineage_name, colormap_df), s=200, marker='o')
        ax.text(0.2, i, lineage_name, va='center', fontsize=8)

    ax.set_xlim(-0.1, 1)
    ax.set_ylim(-1, len(colormap_df))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    apply_standard_axis_style(ax)

    experiment_name = dm.get_experiment_foldername_from_SSOD(seeded_simulation_output_dir)
    so_foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    experiment_simulations_plots_dir = dm.get_experiment_simulations_plots_dir(experiment_name)

    figure_output_path = os.path.join(
        experiment_simulations_plots_dir,
        f'{so_foldername}_{seed}_lineages_colors_tab.tiff'
    )

    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close()

# -----------------------------------------------------------------------------
#                              Trees plots
# -----------------------------------------------------------------------------

def get_fitness_color(fitness_score, nodes_data):
    cmap = matplotlib.pyplot.get_cmap('cool')
    
    # Normalize the value
    norm = matplotlib.colors.Normalize(vmin=nodes_data.fitness_score.min(), vmax=nodes_data.fitness_score.max())
    normalized_value = norm(fitness_score)
    hexcolor = matplotlib.colors.rgb2hex(cmap(normalized_value))
    
    return hexcolor

def get_state_color(state):
    if state == 'infected':
        return 'red'
    if state == 'recovered':
        return 'green'
    if state == 'diagnosed':
        return 'orange'
    if state == 'deceased':
        return 'black'
    
def get_node_color(node, 
                   coloring, 
                   tree_data, 
                   colormap_df): 
    if coloring == 'state':
        color = get_state_color(node.state)
        return color
    if coloring == 'fitness':
        color = get_fitness_color(node.fitness_score, tree_data)
        return color
    if coloring == 'lineage':
        color = get_lineage_color(node.lineage, colormap_df)
        return color

def tree_fitness_legend(tree_data, tree_type, tree_plot_filepath): 
    '''
    Create legend for fitness color in plotted trees 
    
    tree_data: 
        phylogenetic data OR individuals data
    tree: str 
        'infection' or 'phylogenetic' 
    '''
    # Create a dummy invisible image.
    d = np.linspace(0, tree_data.fitness_score.max(),
                    int(tree_data.fitness_score.max())).reshape(1, -1)
    d = np.vstack((d, d))

    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 2))
    fig.subplots_adjust(bottom=0.5, top=0.99, left=0.01, right=0.8)

    # Set the extent to cover the range of data
    extent = [0, tree_data.fitness_score.max(), 0, 1]

    # The imshow plots the dummy data.
    ax.imshow(d, aspect='auto',
                    cmap=matplotlib.pyplot.get_cmap('cool'),
                    extent=extent)

    # Set the ticks at the beginning and end of the color bar
    ax.set_xticks([0, tree_data.fitness_score.max()])
    # ax.set_xticks(np.arange(0,data.fitness.max(),10))
    # Set the labels "low" and "high" for the ticks
    ax.set_xticklabels(["low", "high"])
    
    # Remove y-ticks and their labels
    ax.set_yticklabels([])
    ax.set_yticks([])

    # Set the title for the plot
    ax.set_title("Relative Fitness")
    
    # Extract filename without extension
    tree_filename = os.path.splitext(os.path.basename(tree_plot_filepath))[0]
    
    # Construct the new legend file path
    legend_filename = tree_filename + "_legend.png"
    legend_path = os.path.join(os.path.dirname(tree_plot_filepath), legend_filename)
    
    matplotlib.pyplot.savefig(legend_path,
                  dpi=600, bbox_inches='tight')
    matplotlib.pyplot.close()


def plot_infection_tree(root,
                           infection_tree_data,
                           tree_subtype,
                           coloring,
                           colormap_df,
                           tree_plot_filepath):
    tree_data = infection_tree_data
    
    def nodeattrfunc(node):
        
        if node.infection_type == 'standard':
            return 'color="{}", label="{}",'.format(get_node_color(node, 
                                                                  coloring, 
                                                                  tree_data, 
                                                                  colormap_df), 
                                                   node.name)
        else:
            return 'color="{}", label="{}",shape=diamond, style=filled'.format(get_node_color(node, 
                                                                  coloring, 
                                                                  tree_data, 
                                                                  colormap_df), 
                                                   node.label)
            
    if tree_subtype == 'binary':
        
        DotExporter(root,
                    nodeattrfunc=nodeattrfunc,
                    ).to_picture(tree_plot_filepath)
     
    elif tree_subtype == 'compact':
        
        def edgeattrfunc(node, child):
            if node.name == child.name:
                return 'color=transparent'
            
        DotExporter(root,
                    nodeattrfunc=nodeattrfunc,
                    edgeattrfunc=edgeattrfunc,
                    ).to_picture(tree_plot_filepath)
        
def plot_phylogenetic_tree(root,
                           phylogenetic_data,
                           tree_subtype,
                           coloring,
                           colormap_df,
                           tree_plot_filepath):
    tree_data = phylogenetic_data
    
    # define nodeattrfunc for DotExporter to format the graph picture
    def nodeattrfunc(node):
        if coloring == 'state':
            return 'color="{}", label="{}"'.format('black', 
                                                   node.name)
        else:    
            return 'color="{}", label="{}"'.format(get_node_color(node, 
                                                                  coloring, 
                                                                  tree_data, 
                                                                  colormap_df), 
                                                   node.name)
    
    # format tree for binary tree
    if tree_subtype == 'binary':
    
        # save tree picture
        anytree.exporter.DotExporter(root,
                    nodeattrfunc=nodeattrfunc,
                    ).to_picture(tree_plot_filepath)
        
    # format tree for compact tree
    elif tree_subtype == 'compact':
        
        # define nodeattrfunc for DotExporter for it to not display edges 
        # that connect collapsed nodes
        def edgeattrfunc(node, child):
            if node.name == child.name:
                return 'color=transparent'
            
        # save tree picture
        anytree.exporter.DotExporter(root,
                    nodeattrfunc=nodeattrfunc,
                    edgeattrfunc=edgeattrfunc,
                    ).to_picture(tree_plot_filepath)

def plot_circular_tree(ete_root,
                       tree_type,
                       colormap_df,
                       individuals_lineages,
                       file_path):
    
    # function to set layout of branches
    def color_branches_black(tree):
        """
        Make all branch lines black and set style.
        """
        for node in tree.traverse():
            ns = NodeStyle()
            ns["hz_line_color"] = "black"
            ns["vt_line_color"] = "black"
            ns["hz_line_width"] = 1
            ns["hz_line_type"] = 0
            ns["vt_line_width"] = 1
            # Keep node circles invisible
            ns["size"] = 0
            node.set_style(ns)
            
    # functions to color tree nodes background
    def blend_with_white(hex_color, factor=0.5):
        """
        Blend the input hex color with white by the given factor.
        
        Parameters:
            hex_color (str): Color in the format "#RRGGBB" (e.g., "#ff0000").
            factor (float): Blend ratio between 0 and 1. 
                            0 returns the original color, 1 returns white.
        
        Returns:
            str: The blended color as a hex string in "#RRGGBB" format.
        """
        # Remove '#' if it exists.
        hex_color = hex_color.lstrip('#')
        
        # Convert hex components to integers.
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
        
        # Blend each component with white (255) using the factor.
        r_new = int(r + (255 - r) * factor)
        g_new = int(g + (255 - g) * factor)
        b_new = int(b + (255 - b) * factor)
        
        # Return the new hex color string.
        return f"#{r_new:02x}{g_new:02x}{b_new:02x}"
    
    def get_clade_color(lineage, colormap_df, factor=0.5):
        lineage_color = get_lineage_color(lineage, colormap_df)
        clade_color = blend_with_white(lineage_color, factor)
        return clade_color
    
    def color_clades_by_lineage(tree, colormap_df, individuals_lineages):
        """
        For each node, color the background wedge in circular layout. 
        If plotting phylo tree, color white lineages not present in infection tree.
        """
       
        for node in tree.traverse():
            if node.lineage in individuals_lineages:
                if node.lineage is None:
                    raise ValueError(f'Lineage is NoneType! Check {node}.')
                node.img_style["bgcolor"] = get_clade_color(node.lineage, colormap_df, factor=0.55)
            else:
                node.img_style["bgcolor"] = 'white'
    # tree layout (faces = labels of branches outside tree)
    def layout(node):
        if node.is_leaf():
            N = AttrFace("name", fsize=14)
            faces.add_face_to_node(N, node, 0, position="aligned")
    def layout_inf(node):
        pass
    # color tree
    color_branches_black(ete_root)
    color_clades_by_lineage(ete_root,colormap_df, individuals_lineages)
            
    # setup circular tree
    ts = TreeStyle()
    if tree_type == 'infection':
        ts.layout_fn = layout_inf
        ts.scale = 20
    elif tree_type == 'phylogenetic':
        ts.layout_fn = layout
        ts.scale = 20
        
    ts.mode = "c"  # circular layout    
    ts.show_leaf_name = False
    ts.show_branch_length = False
    # ts.complete_branch_lines_when_necessary = True
    ts.draw_guiding_lines = True
    ts.root_opening_factor = 1.0  
    ts.show_scale = False

    
    # Render tree
    ete_root.render(file_path, w=800, h=800, tree_style=ts)
    print(f"Saved plot: {file_path}.")

# -----------------------------------------------------------------------------
#                              R effective plots
# -----------------------------------------------------------------------------



def plot_OSR_fit_figure(experiment_name,
                parameter,
                fit_result,
                OSR_single_sim_data,
                OSR_combined_sim_data,
                OSR_mean_data,
                model_type,
                min_seq_number,
                min_sim_lenght):
    ''' plot fit of nucleotide substitution rate / observed substitution rates curve
    '''
    
    line_color = 'black' #'#DE8F05' # orange
    scatter_color = '#DE8F05'# orange
    combined_OSR_marker_color = '#E64B9D' # pink
    scatter_color_2 = '#0173B2' # blue '#029E73' # green
    
    # Create figure and axes
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    
    x_data = OSR_single_sim_data[parameter]

    # scatterplot Estimated OSR - single simulation
    sns.scatterplot(x=parameter, y='observed_substitution_rate', 
                    label='Estimated OSR - single simulation', data=OSR_single_sim_data,
                    color=scatter_color, alpha=0.5, ax=ax,
                    zorder=0)
    # plot fitted curve
    sns.lineplot(x=x_data, y=fit_result.best_fit, 
                 label=f'Fitted {model_type} curve', 
                 color=line_color, linewidth=1, ax=ax,
                 zorder=1)
    
    # scatterplot combined regression points (as comparison)
    sns.scatterplot(x=parameter, y='observed_substitution_rate', marker='X',
        label='Combined tempest regression estimate of OSR', data=OSR_combined_sim_data,
        color=combined_OSR_marker_color,alpha=1, ax=ax,
        zorder=2)
    
    # plot mean of observed_substitution_rate
    sns.scatterplot(x=parameter, y='mean', marker = 'X',
                    label='Mean of estimated OSR (single simulations)', data=OSR_mean_data,
                    color=scatter_color_2, alpha=1, ax=ax,
                    zorder=3)

    # Set axis (log log scale) -----------------------------------------------
    ax.set_xlabel(f'{parameter}')
    ax.set_ylabel('observed substitution rate')
    ax.set_xscale('log')
    ax.set_yscale('log') 
    ax.set_xlim(x_data.min(),x_data.max()*1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(dm.get_experiment_plots_dir(experiment_name), 
        f"OSR_{model_type}_fit.png"))

# ---------------------------------------------------------------- clade pies
# Moved here from scripts/long_paper_figures/long_shedders_plots.py in v2.5.6,
# when that directory was archived. figures 2, 2S and 4 all draw these, and a
# plotting helper three figure scripts share belongs with the other plotting
# helpers, not in a figure directory of its own.

PIE_COLORS = {"standard": "#4daf4a", "long": "#e41a1c", "mixed": "#377eb8",
              "founder": "#8172B2"}


def plot_fig4_pies(axes_2x2, pie_data_dict):
    """One pie per metric: which origin class won it, over the seeds of an SOD."""
    metrics = ["Peak", "Burden", "Survival", "Growth"]
    label_order = ["standard", "long", "mixed", "founder"]

    for idx, metric in enumerate(metrics):
        ax = axes_2x2[idx]
        counts = pie_data_dict[metric]

        if sum(counts.values()) == 0:
            ax.text(0.5, 0.5, "No Data", ha="center", va="center", fontsize=8)
            ax.axis("off")
            continue

        values = [counts.get(k, 0) for k in label_order]
        colors = [PIE_COLORS[k] for k in label_order]

        ax.pie(values, colors=colors,
               wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        ax.set_title(metric, fontsize=9)


def add_global_pie_legend(fig):
    """One legend for every pie on the figure, so each pie carries none."""
    handles = [Patch(facecolor=PIE_COLORS[k], edgecolor="white", label=k)
               for k in ["standard", "long", "mixed", "founder"]]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98),
               frameon=False, ncol=4)
