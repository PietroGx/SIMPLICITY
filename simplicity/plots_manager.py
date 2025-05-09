#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:00:16 2025

@author: pietro
"""
import os
import math
import matplotlib
import matplotlib.colors as mcolors
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
import simplicity.tuning.diagnosis_rate    as dr
import simplicity.tuning.evolutionary_rate as er
import simplicity.phenotype.weight         as pheno_weight

# -----------------------------------------------------------------------------
#                        Simplicity simulation plots
# -----------------------------------------------------------------------------

def plot_fitness(simulation_output):
    """
    Plots the fitness trajectory with mean and standard deviation.
    
    """
    # Extract time, mean, and standard deviation data
    times = [coord[0] for coord in simulation_output.fitness_trajectory]
    means = [coord[1][0] for coord in simulation_output.fitness_trajectory]
    stds  = [coord[1][1] for coord in simulation_output.fitness_trajectory]

    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(times, means, linestyle='-', color='g', label='Mean fitness')

    # Add shaded area for standard deviation
    plt.fill_between(times, 
                     [m - s for m, s in zip(means, stds)], 
                     [m + s for m, s in zip(means, stds)], 
                     color='green', alpha=0.3, label='Standard Deviation')

    # Add labels and title
    plt.xlabel('Time')
    plt.ylabel('Average fitness')
    plt.title('Average fitness during simulation')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    # Add legend
    plt.legend()
    # Save the plot
    # plt.savefig(figpath)
    plt.close()
    
def plot_trajectory(seeded_simulation_output_dir):
    '''
    Plots the simulation trajectory and the evolution of infectivity
    during the simulation. (In paper Figure 2B)

    '''
    simulation_trajectory = om.read_simulation_trajectory(seeded_simulation_output_dir)
    time       = simulation_trajectory['time'].tolist()
    infected   = simulation_trajectory['infected']
    diagnosed  = simulation_trajectory['diagnosed']
    susceptibles  = simulation_trajectory['susceptibles'] 

    # fig = plt.figure(0,figsize=(16, 9))
    plt.title('Simulation of SARS-CoV-2 outbreak')
    
    
    plt.plot(time,susceptibles,
        label='Susceptibles', 
        color= 'black')
    plt.plot(time,infected,
            label='Infected individuals', 
            color= 'red')
    plt.plot(time,diagnosed,
            label='Diagnosed individuals (cumulative)', 
            color= 'green')
    
    plt.xlabel('Time - days')
    plt.ylabel('Number of individuals in compartment')
    plt.xlim(0,time[-1])
    plt.ylim(0)
    plt.tight_layout()
    plt.legend()

def plot_simulation(seeded_simulation_output_dir, threshold):
    '''
    joins plots of population trajectory and lineages frequency in a single 
    plot

    '''
    
    from matplotlib.ticker import FuncFormatter
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(num=1,figsize=(12, 8))

    # Create a gridspec object with 3 rows. 
    # The height_ratios argument determines the relative heights of the plots.
    # Here, 1:2 ratio for the plots, and an empty space with '0' ratio.
    gs = GridSpec(3, 1, height_ratios=[1,1,1])  
    
    # Create the subplots using the gridspec indices
    ax2 = fig.add_subplot(gs[0, 0])  
    ax1 = fig.add_subplot(gs[2, 0],sharex=ax2) 
    ax3 = fig.add_subplot(gs[1, 0],sharex=ax2) 
    
    # Adjust subplot parameters
    fig.subplots_adjust(
        top=0.95,
        bottom=0.085,
        left=0.07,
        right=0.81,
        hspace=0.1,
        wspace=0.085
    )
    
    # -------------------- system trajectory subplot  -------------------------
    trajectory_df = om.read_simulation_trajectory(seeded_simulation_output_dir)
    
    time = trajectory_df['time'].tolist()
    infected   = trajectory_df['infected'].tolist()
            
    ax2.plot(time,infected,
            label='Infected individuals at time t', 
            color= 'red')
    
    ax2.set_ylim(0,max(infected)*1.5)
    ax2.tick_params(axis='both', labelsize=15)
    ax2.set_ylabel('N individuals', fontsize=15)
    ax2.legend(loc='upper left', fontsize=15)
    
    # ----------------- lineages frequency subplot ----------------------------
    lineage_frequency_df = om.read_lineage_frequency(seeded_simulation_output_dir)
    
    # get lineages colormap
    colormap_df = make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')
    # save lineages colors tab
    plot_lineages_colors_tab(seeded_simulation_output_dir)
   
    # filter lineages_df by threshold
    filtered_df,_ = om.filter_lineage_frequency_df(lineage_frequency_df, threshold)
    # get colors
    colors = [get_lineage_color(lineage_name, colormap_df) for lineage_name in filtered_df.columns]
    # plot
    filtered_df.plot(kind='area', stacked=False, color=colors,
                         alpha=0.5, ax=ax1)
    time_file_path = os.path.join(seeded_simulation_output_dir, 'final_time.csv')
    time_final = pd.read_csv(time_file_path, header=None).iloc[0, 0]
    ax1.set_xlim([0,time_final])
    # Define custom formatter function
    def to_percent(x, _):
        return f"{100 * x:.0f}%"
    ax1.yaxis.set_major_formatter(FuncFormatter(to_percent))
    # Set tick label sizes for both axes
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_ylabel("Relative frequency of lineage", fontsize=15)
    # ax1.set_xlabel("Time (days)", fontsize=20)
    # ax2.legend(loc='upper left', bbox_to_anchor=(1, 2.95))
    ax1.legend().remove()
    
    # ------------------------ fitness subplot --------------------------------
    fitness_trajectory_file_path = os.path.join(seeded_simulation_output_dir, 'fitness_trajectory.csv')
    fitness_trajectory_df = pd.read_csv(fitness_trajectory_file_path)
    
    # Extract time, mean, and standard deviation data
    times = fitness_trajectory_df['Time'].tolist()
    means = fitness_trajectory_df['Mean'].tolist()
    stds = fitness_trajectory_df['Std'].tolist()
    
    ax3.plot(times, means, linestyle='-', color='#2ca02c', label='Average fitness score')

    # Add shaded area for standard deviation
    ax3.fill_between(times, 
                     [m - s for m, s in zip(means, stds)], 
                     [m + s for m, s in zip(means, stds)], 
                     color='#2ca02c', alpha=0.3, label='Fitness score std')

    # Add labels and title
    # ax3.xlabel('Time')
    ax3.tick_params(axis='both', labelsize=15)
    ax3.set_xlabel('Time - days', fontsize=15)
    ax3.set_xticks(np.arange(0,time_final,50))
    ax3.set_xticklabels(np.arange(0,time_final,50).astype(int))
    ax3.set_ylabel('Fitness score',fontsize=15)
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    # ax3.yaxis.set_label_coords(-0.1, 0)
    # Add legend
    ax3.legend(loc='upper left', fontsize=15)
    
    # Align the y-axis labels
    fig.align_ylabels([ax2, ax3])
    # save plot
    plt.tight_layout()
    experiment_name = dm.get_experiment_foldername_from_SSOD(seeded_simulation_output_dir)
    so_foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    experiment_simulations_plots_dir = dm.get_experiment_simulations_plots_dir(experiment_name)
    figure_output_path = os.path.join(experiment_simulations_plots_dir,f"{so_foldername}_{seed}_simulation_trajectory.png")
    plt.savefig(figure_output_path)
    plt.close()

# -----------------------------------------------------------------------------
#                       Tempest regression and OSR plots
# -----------------------------------------------------------------------------

def ideal_subplot_grid(num_plots):
    # Calculate the number of columns (the ceiling of the square root of the number of plots)
    num_cols = math.ceil(math.sqrt(num_plots))
    
    # Calculate the number of rows needed
    num_rows = math.ceil(num_plots / num_cols)
    
    return num_rows, num_cols        

def plot_tempest_regression(sequencing_data_df,
                            fitted_tempest_regression,
                            ax):
    # get linear regression points
    x = sequencing_data_df['Sequencing_time'].values.reshape(-1, 1)
    y_pred = fitted_tempest_regression.predict(x)
    # get observed subswtitution rate (OSR) value
    observed_substitution_rate = fitted_tempest_regression.coef_[0] # substitution rate per site per year
    # plot data points
    sequencing_data_df.plot(kind='scatter', x='Sequencing_time', y='Distance_from_root', color='blue', ax=ax)
    # plot linear regression
    ax.plot(x, y_pred, color='red', linewidth=2, 
            label='y = OSR · x')
    ax.set_xlabel('Simulation time in years')
    ax.set_ylabel('Distance from root (normed substitions/site)')
    ax.set_xlim(left=0,right=3)
    ax.set_ylim(0)
    ax.grid(True)
    # Adding OSR value to the legend
    extra_text = f'OSR = {observed_substitution_rate:.5f}'
    handles, labels = ax.get_legend_handles_labels()
    handles.append(plt.Line2D([0], [0], color='white', label=extra_text))  # Invisible line for text
    # Show legend
    ax.legend(handles=handles)
   
def plot_figure_tempest_regression():
    experiment_name = 'generate_data_OSR_fit_#1'
    parameter = 'nucleotide_substitution_rate'
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    simulation_output_dir = simulation_output_dirs[7]
    param = sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir, parameter)
    sequencing_data_df = om.create_combined_sequencing_df(simulation_output_dir, min_sim_lenght=0)
    
    fig, ax = plt.subplots(1, 1, figsize=(8,8))
    
    if sequencing_data_df is None: 
        print('no data to plot')
        print(simulation_output_dir)
        pass
    else:
        fitted_tempest_regression = er.tempest_regression(sequencing_data_df)
        
        plot_tempest_regression(sequencing_data_df,
                                   fitted_tempest_regression,
                                   ax)
        ax.set_title(f'Regression - {parameter}: {param}')
        # ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(os.path.join(experiment_plots_dir, f"{experiment_name}_figure_combined_tempest_regression.png"))
    
def plot_combined_tempest_regressions(experiment_name, parameter, min_sim_lenght=0, y_axis_max=0.1):
    # Get sorted simulation output directories for experiment
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
  
    sorted_simulation_output_dirs = sorted(simulation_output_dirs, 
                     key=lambda dir: sm.get_parameter_value_from_simulation_output_dir(dir, parameter))
    
    # Determine the number of rows and columns for subplots
    num_rows, num_cols = ideal_subplot_grid(len(sorted_simulation_output_dirs))

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    
    # If there's only one subplot, axs is not an array, so we convert it to a 2D array for consistency
    if num_rows * num_cols == 1:
        axs = [[axs]]
    elif num_rows == 1 or num_cols == 1:
        axs = axs.reshape(1, -1)
    else:
        axs = axs
    
    for i, simulation_output_dir in enumerate(sorted_simulation_output_dirs):
        param = sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir, parameter)
        sequencing_data_df = om.create_combined_sequencing_df(simulation_output_dir, min_sim_lenght)
        if sequencing_data_df is None: 
            pass
        else:
            fitted_tempest_regression = er.tempest_regression(sequencing_data_df)
            ax = axs[i // num_cols][i % num_cols]
            plot_tempest_regression(sequencing_data_df,
                                       fitted_tempest_regression,
                                       ax)
            ax.set_title(f'Regression - {parameter}: {param}')
            ax.set_ylim(0, y_axis_max)
    
    plt.tight_layout()
    plt.savefig(os.path.join(experiment_plots_dir, f"{experiment_name}_combined_regression.png"))

def plot_combined_OSR_vs_parameter(experiment_name, 
                                    parameter,  
                                    min_seq_number=0,
                                    min_sim_lenght=0):
    ''' Plot observed substituion rate (tempest regression) against desired parameter values
    '''
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)

    df = om.read_OSR_vs_parameter_csv(experiment_name, 
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)

    # Boxplot for each X position
    sns.boxplot(x=df[parameter], y=df['observed_substitution_rate'], hue=df[parameter], 
                palette="coolwarm", width=0.6, showfliers=True)
    
    plt.xlabel(parameter)
    plt.ylabel('Observed substitution rate')
    plt.title(f'{parameter} vs observed substitution rate')
    plt.grid(True)
    plt.tight_layout()
    plt.legend().remove()
    plt.savefig(os.path.join(experiment_plots_dir, 
                             f"{experiment_name}_{parameter}_vs_observed_substitution_rate.png"))
        
def plot_combined_OSR_fit(experiment_name, 
                                                 fit_result, 
                                                 model_type, 
                                                 min_seq_number,
                                                 min_sim_lenght):
    ''' plot fit of nucleotide substitution rate / observed substitution rate curve
    '''

    data = om.read_combined_OSR_vs_parameter_csv(experiment_name, 
                                                 'nucleotide_substitution_rate', 
                                                 min_seq_number,
                                                 min_sim_lenght)
    x_data = data['nucleotide_substitution_rate'] 
    y_data = data['observed_substitution_rate']  
    
    # Create figure and axes
    fig, ax = plt.subplots(3,1, figsize=(8, 10))
    # scatterplot data (u)
    ax[0].scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
    ax[1].scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
    ax[2].scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
    
    # linear scale
    ax[0].plot(x_data, fit_result.best_fit, label=f'Fitted {model_type} curve', 
               color='red', linewidth=2)
    ax[0].set_ylim(0)
    ax[0].set_xlim(0)
    # ax[0].set_xlabel('Substitution Rate (site/year) (e)')
    # ax[0].set_ylabel('observed substitution rate (site/year) (u)')
    ax[0].legend()
    # semilog scale
    ax[1].plot(x_data, fit_result.best_fit, label=f'Fitted {model_type} curve',
               color='red', linewidth=2)
    ax[1].set_xscale("log")
    ax[1].set_ylim(0)
    # ax[1].set_xlabel('Substitution Rate (site/year) (e)')
    ax[1].set_ylabel('observed substitution rate (site/year) (u)')
    ax[1].legend()
    # loglog scale
    ax[2].plot(x_data, fit_result.best_fit, label=f'Fitted {model_type} curve',
               color='red', linewidth=2)
    ax[2].set_yscale("log")
    ax[2].set_xscale("log")
    ax[2].set_xlabel('Substitution Rate (site/year) (e)')
    # ax[2].set_ylabel('observed substitution rate (site/year) (u)')
    ax[2].legend()
    
    # fig.legend(loc="upper center", ncol=2, fontsize=10)
    
    plt.title('')
    # plt.subplots_adjust(top=0.6)
    plt.tight_layout()
    file_path = os.path.join(dm.get_experiment_plots_dir(experiment_name),
     f'{experiment_name}_combined_observed_substitution_rate_{model_type}_fit.png')
    plt.savefig(file_path)

def plot_OSR_fit(experiment_name, 
                fit_result, 
                model_type,
                min_seq_number,
                min_sim_lenght):
    ''' plot fit of nucleotide substitution rate / observed substitution rates curve
    '''
    parameter = 'nucleotide_substitution_rate'
    
    line_color = 'black' #'#DE8F05' # orange
    scatter_color = '#DE8F05'# orange
    combined_OSR_marker_color = '#E64B9D' # pink
    scatter_color_2 = '#0173B2' # blue '#029E73' # green
    
    # Create figure and axes
    fig, ax = plt.subplots(3,1, figsize=(8, 12))
    
    # import combined regression data
    combined_data = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                                          parameter,
                                                          min_seq_number,
                                                          min_sim_lenght)
    
    # import single simulations regression data
    data = om.read_OSR_vs_parameter_csv(experiment_name, 
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    # Group by nucleotide_substitution_rate and compute mean and standard deviation for OSR
    data_mean_df = om.get_mean_std_OSR(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    # get lower and upper confidence interval for fit results
    x_data = data['nucleotide_substitution_rate']
    # x, lower_curve, upper_curve = confidence_interval_fit(model_type, fit_result, x_data)

    
    # plot on each axis of subplots in a loop
    for a in ax:
        # # Fill between the upper and lower curves for the confidence interval region
        # a.fill_between(x,lower_curve, upper_curve, 
        #                color=line_color, alpha=0.3, label='95% Confidence Interval',
        #                zorder=-1)
        # scatterplot Estimated OSR - single simulation
        sns.scatterplot(x=parameter, y='observed_substitution_rate', 
                        label='Estimated OSR - single simulation', data=data,
                        color=scatter_color, alpha=0.5, ax=a,
                        zorder=0)
        # plot fitted curve
        sns.lineplot(x=x_data, y=fit_result.best_fit, 
                     label=f'Fitted {model_type} curve', 
                     color=line_color, linewidth=1, ax=a,
                     zorder=1)
        # scatterplot combined regression points (as comparison)
        sns.scatterplot(x='nucleotide_substitution_rate', y='observed_substitution_rate', marker='X',
            label='Combined tempest regression estimate of OSR', data=combined_data,
            color=combined_OSR_marker_color,alpha=1, ax=a,
            zorder=2)
        # plot mean of observed_substitution_rate from data_mean_std
        sns.scatterplot(x=parameter, y='mean', marker = 'X',
                        label='Mean of estimated OSR (single simulations)', data=data_mean_df,
                        color=scatter_color_2, alpha=1, ax=a,
                        zorder=3)
        # Add horizontal lines 
        # a.hlines(y=[1e-5, 1e-2], xmin=0, xmax=x_data.max(), colors=['r', 'r'], linestyles='--')
        # Set y-axis limits
        a.set_ylim(0.000009, 0.02)
        
        # minsimlenghts = [0,100,200,300]
        # palette = sns.color_palette("tab10", len(minsimlenghts))
        # for i, minsimlenght in enumerate(minsimlenghts):
        #     df = om.get_mean_std_observed_substitution_rates(experiment_name,
                                                                # parameter,
                                                                # min_seq_number,
                                                                # minsimlenght)
        #     # plot mean of observed_substitution_rate from data_mean_std
        #     sns.scatterplot(x=parameter, y='mean', marker = 'X',
        #                     label=f'Mean of estimated OSR - min {minsimlenght} d', data=df,
        #                     color=palette[i], alpha=0.8, ax=a,
        #                     zorder=4)

    # First plot (linear scale) -----------------------------------------------
    ax[0].set_xlabel(f'{parameter}')
    ax[0].set_ylabel('observed substitution rate')
    
    
    # Second plot (semilog scale) -----------------------------------------------
    
    ax[1].set_xlabel(f'{parameter}')
    ax[1].set_ylabel('observed substitution rate')
    ax[1].set_xscale('log')
    ax[1].legend_.remove()
    
    # Third plot (log log scale) -----------------------------------------------
    ax[2].set_xlabel(f'{parameter}')
    ax[2].set_ylabel('observed substitution rate')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].legend_.remove()
    
    plt.tight_layout()
    plt.savefig(os.path.join(dm.get_experiment_plots_dir(experiment_name), 
        f"{experiment_name}_observed_substitution_rates_{model_type}_fit.png"))
    
def plot_OSR_fit_figure(experiment_name, 
                fit_result, 
                model_type,
                min_seq_number,
                min_sim_lenght):
    ''' plot fit of nucleotide substitution rate / observed substitution rates curve
    '''
    parameter = 'nucleotide_substitution_rate'
    
    line_color = 'black' #'#DE8F05' # orange
    scatter_color = '#DE8F05'# orange
    combined_OSR_marker_color = '#E64B9D' # pink
    scatter_color_2 = '#0173B2' # blue '#029E73' # green
    
    # Create figure and axes
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    
    # import combined regression data
    combined_data = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                                          parameter,
                                                          min_seq_number,
                                                          min_sim_lenght)
    
    # import single simulations regression data
    data = om.read_OSR_vs_parameter_csv(experiment_name, 
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    # Group by nucleotide_substitution_rate and compute mean and standard deviation for OSR
    data_mean_df = om.get_mean_std_OSR(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    # get lower and upper confidence interval for fit results
    x_data = data['nucleotide_substitution_rate']
    # x, lower_curve, upper_curve = confidence_interval_fit(model_type, fit_result, x_data)
    
    
    # # Fill between the upper and lower curves for the confidence interval region
    # ax.fill_between(x,lower_curve, upper_curve, 
    #                color=line_color, alpha=0.3, label='95% Confidence Interval',
    #                zorder=-1)
    # scatterplot Estimated OSR - single simulation
    sns.scatterplot(x=parameter, y='observed_substitution_rate', 
                    label='Estimated OSR - single simulation', data=data,
                    color=scatter_color, alpha=0.5, ax=ax,
                    zorder=0)
    # plot fitted curve
    sns.lineplot(x=x_data, y=fit_result.best_fit, 
                 label=f'Fitted {model_type} curve', 
                 color=line_color, linewidth=1, ax=ax,
                 zorder=1)
    # scatterplot combined regression points (as comparison)
    sns.scatterplot(x='nucleotide_substitution_rate', y='observed_substitution_rate', marker='X',
        label='Combined tempest regression estimate of OSR', data=combined_data,
        color=combined_OSR_marker_color,alpha=1, ax=ax,
        zorder=2)
    # plot mean of observed_substitution_rate from data_mean_std
    sns.scatterplot(x=parameter, y='mean', marker = 'X',
                    label='Mean of estimated OSR (single simulations)', data=data_mean_df,
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
        f"Figure4_OSR_{model_type}_fit.png"))
    
# def confidence_interval_fit(model_type, fit_result, x):
#     params_lower = {}
#     params_upper = {}
    
    
#     # print(f'{model_type}')
#     for param in fit_result.params:
#         param_value = fit_result.params[param].value
#         param_stderr = fit_result.params[param].stderr if fit_result.params[param].stderr else 0
        
#         # print(f'CI std : {param} - {param_value} - {param_stderr}')
#         ci_lower = param_value - 1.96 * param_stderr
#         ci_upper = param_value + 1.96 * param_stderr
        
#         params_lower[param] = ci_lower
#         params_upper[param] = ci_upper
    
#     def remove_duplicates_array(arr):
#         _, idx = np.unique(arr, return_index=True)
#         return arr[np.sort(idx)]

#     x = remove_duplicates_array(x)

#     # Create the upper and lower bound curves for confidence intervals
#     upper_curve = er.evaluate_model(model_type, params_upper, x)
#     lower_curve = er.evaluate_model(model_type, params_lower, x)
#     return x, lower_curve, upper_curve


# -----------------------------------------------------------------------------
#                          Phenotype model plots
# -----------------------------------------------------------------------------

def plot_w_t(t_max, t_eval): # plot weights of phenotype model
    '''plot the phenotype model weight function for calculating fitness score 
    (distance from weighted consensus sequence)
    '''
    
    params = pheno_weight.w_t_params()
    k_e = params[0]
    k_a = params[1]
    
    # Time points 
    t = np.linspace(0,300,1000)

    # Calculate concentration for each time point
    w_t = pheno_weight.weights(t, t_eval, k_e, k_a, t_max)

    # Plotting the concentration-time profile
    plt.plot(t, w_t, label=f'w(t) at simulation time {t_eval}',color='black')
    plt.xlabel('Time (days)',fontsize=20)
    plt.ylabel('Normalized Antibody Concentration (weight)',fontsize=20)
    plt.ylim(0,1.1)
    plt.xlim(left=0)
    plt.tick_params(axis='both', labelsize=20)
    # plt.title('Concentration-Time Profile')
    plt.legend(loc ='upper left',fontsize=20)
    plt.show()
    
# -----------------------------------------------------------------------------
#                             Intra host model plots
# -----------------------------------------------------------------------------
def plot_intra_host(intra_host_model,time,step):
    '''
    Plot intra-host model solutions.

    Parameters
    ----------
    time : float
        Time for the intra-host model solution

    Returns
    -------
    Plot p_inf/p_dia/p_rec.

    '''
    import matplotlib.colors as mcolors
    from matplotlib import cm
    matplotlib.rcParams.update({'font.size': 22})

    t = np.arange(0,time,step)
    states = np.arange(0,21,1) #[0] # 
    # setup the normalization and the colormap
    normalize = mcolors.Normalize(vmin=min(states), vmax=max(states))
    colormap = cm.brg
    
    # plot curve for each starting state
    for state in states:
        plt.plot(t,intra_host_model._data_plot_model(state,time,step), 
                 color = colormap(normalize(state)))
        
    # setup the colorbar
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(states)
    cb = plt.colorbar(scalarmappaple)
    cb.set_ticks(np.arange(0,21,1))
    plt.text(300.5,1.1,'Initial state')
   
    # show the figure
    plt.xlabel('Time (d)')
    plt.ylabel('Probability of being inf after t days')
    plt.xlim(0,time)
    plt.ylim(0)
    plt.show()    
    
def plot_comparison_intra_host_models(intra_host_model):
    '''
    Plot intra-host model probability of being infectious after t days for
    normal and long shedders individuals

    Parameters
    ----------
    time : float
        Time for the intra-host model solution

    Returns
    -------
    Plot p_inf

    '''
    import scipy
    
    matplotlib.rcParams.update({'font.size': 22})

    t = np.arange(0,50,0.1)

    plt.plot(t,intra_host_model._data_plot_model(0,50,.1), color = 'orange', 
             label = 'normal')
    
    # setup intra-host model matrix and calculate matrix exponential
    intra_host_model.A = intra_host_model._matrix(133.5)
    intra_host_model.A_ex = scipy.linalg.expm(intra_host_model.A)
    
    t = np.arange(0,300,1)
    plt.plot(t,intra_host_model._data_plot_model(0,300,1), color = 'blue', 
             label = 'immunocompromised')
    # show the figure
    plt.xlabel('Time (d)')
    plt.ylabel('Probability of being inf. after t days')
    plt.xlim(0,300)
    plt.ylim(0)
    plt.legend()
    plt.show()    
    
# -----------------------------------------------------------------------------
#                           Runtime analysis plot
# -----------------------------------------------------------------------------

def plot_extrande_pop_runtime(extrande_pop_runtime_csv):
    ''' Plots infected number vs simulation RUNTIME.
    extrande_pop_runtime_csv is generated in extrande for memory profiling (commented out)
    '''
    import csv
    # Read the CSV file
    x_values = []
    y_values = []
    
    with open(extrande_pop_runtime_csv, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            x, y = map(float, row)  
            x_values.append(x)
            y_values.append(y)
            
    # Sort based on x_values 
    x_values, y_values = zip(*sorted(zip(x_values, y_values)))
    # Plotting the data
    plt.plot(x_values, y_values) 
    plt.xlabel('Runtime (s)')
    plt.ylabel('Infected number')
    plt.title('Plot infected number over simulation runtime')
    plt.grid(True)
    plt.savefig('extrande_pop_runtime.png')
    plt.show()
    
# -----------------------------------------------------------------------------
#                          Diagnosis rate plots
# -----------------------------------------------------------------------------

def plot_effective_theoretical_diagnosis_rate(experiment_name):
    ''' plot scatter and regression line of effective vs theoretical diagnosis rate
    '''
    import simplicity.tuning.diagnosis_rate as dr
    from sklearn.linear_model import LinearRegression
    # get theoretical and effective diagnosis rates 
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name) 
    diagnosis_rates_coord = []
    std_effective_rates = []
    for simulation_output_dir in simulation_output_dirs:
        diagnosis_rates, std_effective_rate = dr.get_diagnosis_rates(simulation_output_dir)
        diagnosis_rates_coord.append(diagnosis_rates) 
        std_effective_rates.append(std_effective_rate)
        
    x, y = zip(*diagnosis_rates_coord) # x is theoretical rate, y is effective
    
    # Convert x and y to numpy arrays for compatibility with sklearn
    x = np.array(x).reshape(-1, 1)  # Reshape for sklearn
    y = np.array(y)
    
    # Fit linear regression
    model = LinearRegression()
    model.fit(x, y)
    
    y_pred = model.predict(x)
    
    # Get the regression slope and intercept
    slope = model.coef_[0]
    intercept = model.intercept_

    # Creating the scatter plot
    plt.scatter(x, y)
    # Plot the regression line
    plt.plot(x, y_pred, color='red', 
             label=f'Regression Line: y = {slope:.2f}x + {intercept:.2f}')
    # Plot the standard deviation as vertical lines
    for i in range(len(x)):
        plt.errorbar(x[i], y[i], yerr=std_effective_rate, fmt='o', color='black', capsize=5)
   
    # Adding labels
    plt.title('Linear regression: theoretical vs effective diagnosis rate')
    plt.xlabel('Theoretical diagnosis rate')
    plt.xlim(0,max(x)*1.1)
    plt.ylabel('Effective diagnosis rate')
    plt.ylim(0,max(x)*1.1)
    plt.legend()
    # Display the plot
    plt.show()

def plot_heatmap_R_diagnosis_rate(experiment_name):
    ''' Plot heatnmap to see relationship of R and diagnosis rate.  
    y axis is R 
    x axis is theoretical diagnosis rate (user input)
    score for heatmap is effective/theoretical diagnosis rate
    '''
    # get theoretical and effective diagnosis rates 
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name) 
    # create df
    heatmap_data = pd.DataFrame()
    # get data
    theoretical_diagnosis_rates = []
    R_values = []
    diagnosis_rate_scores = [] 
    for simulation_output_dir in simulation_output_dirs:
        diagnosis_rates, _ = dr.get_diagnosis_rates(simulation_output_dir)
        
        theoretical_diagnosis_rates.append(sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir,
                                                                             'diagnosis_rate'))
        R_values.append(sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir,
                                                                             'R'))
        diagnosis_rate_scores.append(diagnosis_rates[1]/diagnosis_rates[0])
    # fill the df iteratively
    for diagnosis_rate, R, score in zip(theoretical_diagnosis_rates, 
                                        R_values, 
                                        diagnosis_rate_scores):
        heatmap_data.loc[diagnosis_rate, R] = score
    # sort the df by labels values
    heatmap_data = heatmap_data.sort_index(axis=1)
    heatmap_data = heatmap_data.sort_index(axis=0)
    # pivot df 
    heatmap_data = heatmap_data.T
    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap='viridis')
    plt.title('Heatmap of diagnosis rates vs R')
    plt.ylabel('R')
    plt.xlabel('Theoretical (user input) diagnosis rates')
    
# -----------------------------------------------------------------------------
#                      Intra host lineage distribution plot
# -----------------------------------------------------------------------------

def plot_IH_lineage_distribution(experiment_name):
    fig, axes  = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
    df = om.get_IH_lineages_data_experiment(experiment_name)
    
    ax0 = sns.barplot(x='IH_lineages_number', y='ih_virus_count', data=df, hue='IH_virus_emergence_rate', 
                ax=axes[0], dodge=True, palette="Set2", alpha=0.7)
    ax1 = sns.barplot(x='IH_unique_lineages_number', y='ih_lineage_count', data=df, hue='IH_virus_emergence_rate',
                ax=axes[1], dodge=True, palette="Set2",alpha=0.7)
    
    # Loop through each bar and set the edgecolor to match the bar color
    for patch in ax0.patches:
        # Get the color of the bar
        color = patch.get_facecolor()
        # Set the edgecolor to the same as the facecolor
        patch.set_edgecolor(color)
        patch.set_linewidth(1)
        
    for patch in ax1.patches:
        # Get the color of the bar
        color = patch.get_facecolor()
        # Set the edgecolor to the same as the facecolor
        patch.set_edgecolor(color)
        patch.set_linewidth(1)
    
    axes[0].set_title('Count of IH_lineages_number')
    axes[0].set_xlabel('IH_lineages_number')
    axes[0].set_ylabel('Count')

    axes[1].set_title('Count of unique lineages_number')
    axes[1].set_xlabel('Unique lineages_number')
    axes[1].set_ylabel('Count')
    
    
    plt.suptitle('Histogram of IH virus and lineage distribution')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    fig_path = os.path.join(experiment_output_dir,'IH variability plot.png')
    plt.savefig(fig_path)

def plot_IH_lineage_distribution_simulation(experiment_name):
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    num_sims = len(simulation_output_dirs)
    ncols = 2
    nrows = math.ceil(num_sims / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(7 * ncols, 5 * nrows))
    axes = axes.flatten()  # Flatten in case of a single row/col

    for i, simulation_output_dir in enumerate(simulation_output_dirs):
        df = om.get_IH_lineages_data_simulation(simulation_output_dir)
        ax = axes[i]

        sns.barplot(
            x='IH_unique_lineages_number', y='ih_lineage_count', 
            data=df,
            ax=ax, dodge=True, palette="Set2", alpha=0.7
        )

        # Beautify bars
        for patch in ax.patches:
            color = patch.get_facecolor()
            patch.set_edgecolor(color)
            patch.set_linewidth(1)

        sim_id = os.path.basename(simulation_output_dir)
        ax.set_title(f'Simulation: {sim_id}')
        ax.set_xlabel('IH_lineages_number')
        ax.set_ylabel('Proportion')

        # Optional: shrink legend
        ax.legend(title='Intra-host lineages', fontsize=9, title_fontsize=10)

    # Hide any unused axes
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle('IH Lineage Distribution by Simulation', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    fig_path = os.path.join(experiment_output_dir, 'IH_variability_by_simulation.png')
    plt.savefig(fig_path)
    plt.close(fig)

def plot_IH_lineage_distribution_grouped_by_simulation(experiment_name):
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    plot_output_dir = dm.get_experiment_plots_dir(experiment_name)

    all_data = []

    for sim_dir in simulation_output_dirs:
        df = om.get_IH_lineages_data_simulation(sim_dir)[['IH_unique_lineages_number', 'ih_lineage_count']]
        tau_3 = sm.get_parameter_value_from_simulation_output_dir(sim_dir, 'tau_3')
        df['tau_3'] = tau_3

        # Normalize ih_lineage_count
        total = df['ih_lineage_count'].sum()
        df['ih_lineage_count'] = df['ih_lineage_count'] / total if total > 0 else 0

        all_data.append(df)

    full_df = pd.concat(all_data, axis=0)

    # Ensure IH_unique_lineages_number is categorical with fixed order 1–5
    full_df['IH_unique_lineages_number'] = pd.Categorical(
        full_df['IH_unique_lineages_number'],
        categories=[1, 2, 3, 4, 5],
        ordered=True
    )

    # Sort hue order
    tau_3_order = sorted(full_df['tau_3'].unique())

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        x='IH_unique_lineages_number',
        y='ih_lineage_count',
        hue='tau_3',
        data=full_df,
        palette=sns.color_palette("coolwarm", n_colors=len(tau_3_order)),
        hue_order=tau_3_order,
        dodge=True,
        alpha=0.8
    )

    ax.set_title('Unique IH Lineage Distribution by Simulation (Normalized)', fontsize=14)
    ax.set_xlabel('IH_unique_lineages_number')
    ax.set_ylabel('Proportion')
    ax.legend(title='tau_3', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    output_path = os.path.join(plot_output_dir, f'{experiment_name}_IH_unique_lineages_by_tau_3.png')
    plt.savefig(output_path)
    plt.close()



# -----------------------------------------------------------------------------
#                       Simulations final times histogram
# -----------------------------------------------------------------------------

def plot_histograms(experiment_name, final_times_data_frames, r_order=None):
    ''' 
    Plot histograms of simulation final times.
    Plots are ordered by R if r_order is provided.
    All subplots share the same x-axis limits.
    '''
    num_folders = len(final_times_data_frames.columns)

    # Determine global min and max for x-axis
    all_values = []
    for col in final_times_data_frames.columns:
        all_values.extend(final_times_data_frames[col].dropna().values)
    if not all_values:
        print("No data to plot.")
        return
    global_min, global_max = min(all_values), max(all_values)

    # Apply R-ordering if provided
    if r_order:
        ordered_columns = list(final_times_data_frames.columns)
        col_to_r = {col: r for col, r in zip(ordered_columns, r_order)}
        ordered_columns = sorted(ordered_columns, key=lambda col: col_to_r[col])
        final_times_data_frames = final_times_data_frames[ordered_columns]

    fig, axes = plt.subplots(num_folders, 1, figsize=(10, 4.5 * num_folders), squeeze=False)

    for ax, (folder_name, data) in zip(axes.flatten(), final_times_data_frames.items()):
        ax.hist(data.dropna(), bins=30, edgecolor='black')
        ax.set_xlim(global_min, global_max)
        ax.set_title(f'Histogram for {folder_name}')
        ax.set_xlabel('Final Time Value')
        ax.set_ylabel('Frequency')

    plt.tight_layout()
    out_path = os.path.join(dm.get_experiment_dir(experiment_name), 'simulations_length_histogram.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved histogram plot to {out_path}")


# -----------------------------------------------------------------------------
#                              Lineages colors 
# -----------------------------------------------------------------------------

def make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow'):
    """
    Generate a DataFrame mapping lineages to unique colors based on their time order.

    Returns:
    pd.DataFrame: DataFrame with 'Lineage_name' and 'Color' columns.
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
    """Return the corresponding color for a given lineage name.
    lineage: str 
        string with lineage name
    color mapping: colormap
        output of make_lineages_colormap
    """
    df = colormap_df
    if lineage_name is None:
        raise ValueError('Lineage cannot be NoneType!')
    elif lineage_name not in df['Lineage_name'].values:
        raise ValueError(f'Lineage {lineage_name} not in lineages list! Cannot compute color.')
        
    lineage_color = df.loc[df['Lineage_name'] == lineage_name, 'Color'].values[0]
    return lineage_color  

def plot_lineages_colors_tab(seeded_simulation_output_dir):
    
    colormap_df = make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')
    # Plot lineages colors in a table format with circles
    fig, ax = plt.subplots(num=2,figsize=(5, len(colormap_df) * 0.5))
    for i, lineage_name in enumerate(reversed(colormap_df['Lineage_name'].values)):
        ax.scatter(0, i, color=get_lineage_color(lineage_name, colormap_df), s=200, marker='o')  # Draw a circle
        ax.text(0.2, i, lineage_name, va='center', fontsize=12)  # Add text next to it
    
    ax.set_xlim(-0.1, 1)
    ax.set_ylim(-1, len(colormap_df))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    experiment_name = dm.get_experiment_foldername_from_SSOD(seeded_simulation_output_dir)
    so_foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    experiment_simulations_plots_dir = dm.get_experiment_simulations_plots_dir(experiment_name)
    plt.savefig(os.path.join(experiment_simulations_plots_dir,
                             f'{so_foldername}_{seed}_lineages_colors_tab.png'))
    plt.close()

# -----------------------------------------------------------------------------
#                              Trees plots
# -----------------------------------------------------------------------------

def get_fitness_color(fitness, nodes_data):
    cmap = matplotlib.pyplot.get_cmap('cool')
    
    # Normalize the value
    norm = matplotlib.colors.Normalize(vmin=nodes_data.fitness.min(), vmax=nodes_data.fitness.max())
    normalized_value = norm(fitness)
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
        color = get_fitness_color(node.fitness, tree_data)
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
    d = np.linspace(0, tree_data.fitness.max(),
                    int(tree_data.fitness.max())).reshape(1, -1)
    d = np.vstack((d, d))

    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 2))
    fig.subplots_adjust(bottom=0.5, top=0.99, left=0.01, right=0.8)

    # Set the extent to cover the range of data
    extent = [0, tree_data.fitness.max(), 0, 1]

    # The imshow plots the dummy data.
    ax.imshow(d, aspect='auto',
                    cmap=matplotlib.pyplot.get_cmap('cool'),
                    extent=extent)

    # Set the ticks at the beginning and end of the color bar
    ax.set_xticks([0, tree_data.fitness.max()])
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
        
        if node.infection_type == 'normal':
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
    
    # Render tree
    ete_root.render(file_path, w=800, h=800, tree_style=ts)
    print(f"Saved plot: {file_path}.")

# -----------------------------------------------------------------------------
#                              R effective plots
# -----------------------------------------------------------------------------

def plot_infections_hist(individuals_df, ax, colormap_df, bin_size):
    """
    Plot stacked histogram of infections per lineage over time using `inherited_lineage`.
    """
    df = individuals_df[individuals_df['t_infection'].notnull()].copy()
    df = df[df['inherited_lineage'].notnull()]

    max_time = df['t_infection'].max()
    bin_edges = np.arange(0, max_time + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = (bin_edges[1] - bin_edges[0]) * 0.97

    bottom = np.zeros(len(bin_centers))

    for lineage in sorted(df['inherited_lineage'].unique()):
        lineage_df = df[df['inherited_lineage'] == lineage]
        counts, _ = np.histogram(lineage_df['t_infection'], bins=bin_edges)
        color = get_lineage_color(lineage, colormap_df)
        ax.bar(bin_centers, counts, width=bin_width, bottom=bottom, color=color, alpha=0.7, label=lineage)
        bottom += counts

    ax.set_ylim(0, bottom.max() * 1.2)
    ax.set_ylabel("Infections in time window")

def plot_combined_r_effective(experiment_name, seeded_simulation_output_dir, window_size, threshold):
    """
    Create a two-panel plot:
    - Top: population-level Rₑ and infection histogram
    - Bottom: lineage-specific Rₑ (filtered by frequency threshold)

    Data is loaded or computed from CSV using output_manager functions.
    """
    # Load individuals
    individuals_df = om.read_individuals_data(seeded_simulation_output_dir)

    # Load or compute R_effective trajectories
    try:
        r_effective_population_traj, r_effective_lineages_traj = om.read_r_effective_trajs_csv(
            experiment_name, seeded_simulation_output_dir, window_size, threshold
        )
        print("Rₑ data loaded from CSV.")
    except FileNotFoundError:
        print("Rₑ CSVs not found. Computing and saving...")
        om.write_r_effective_trajs_csv(
            experiment_name,
            seeded_simulation_output_dir,
            window_size,
            threshold
        )
        r_effective_population_traj, r_effective_lineages_traj = om.read_r_effective_trajs_csv(
            experiment_name, seeded_simulation_output_dir, window_size, threshold
        )
    # -------------------------------------------------------------------------
    # Prepare colormap
    colormap_df = make_lineages_colormap(seeded_simulation_output_dir)

    # Set up figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # -------------------------------------------------------------------------
    # Population R_eff + histogram
    ax1b = ax1.twinx()
    ax1.set_zorder(ax1b.get_zorder() + 1)
    ax1.patch.set_visible(False)
    
    bin_size = 7 # days
    plot_infections_hist(individuals_df, ax1b, colormap_df, bin_size)

    label_avg = f"Rₑ (window={window_size}d)"
    ax1.plot(r_effective_population_traj.index, r_effective_population_traj.values,
             color='black', label=label_avg)

    R_avg = r_effective_population_traj.mean()
    ax1.axhline(R_avg, color='red', linestyle='--', label=f"Overall avg Rₑ = {R_avg:.2f}")

    ax1.set_ylabel("Population Rₑ")
    ax1.legend(loc='upper left')
    ax1.set_title("Average Rₑ and Infections Histogram")
    # -------------------------------------------------------------------------
    # Lineage-specific R_eff
    lineage_freq_df = om.read_lineage_frequency(seeded_simulation_output_dir)
    _, filtered_lineages = om.filter_lineage_frequency_df(lineage_freq_df, threshold)

    for lineage in filtered_lineages:
        if lineage not in r_effective_lineages_traj:
            continue
        r_series = r_effective_lineages_traj[lineage]
        color = get_lineage_color(lineage, colormap_df)
        ax2.plot(r_series.index, r_series.values, label=lineage, color=color)

    # Overplot average
    ax2.plot(r_effective_population_traj.index, r_effective_population_traj.values,
             color='black', label=label_avg, zorder=100)

    ax2.set_title(f"Lineage Rₑ (filtered, > {threshold * 100:.0f}% prevalence)")
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Lineage Rₑ")
    ax2.legend(loc='upper left', fontsize='small', ncol=2)

    # Formatting
    ax1.set_xlim(0, r_effective_population_traj.index.max())
    ax2.set_ylim(bottom=0)
    ax1.spines['top'].set_visible(False)
    ax1b.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()

    # Save to file
    output_dir = dm.get_experiment_plots_dir(experiment_name)
    foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    filename = f"{experiment_name}_{foldername}_{seed}_combined_R_effective.png"
    filepath = os.path.join(output_dir, filename)

    plt.savefig(filepath)
    plt.close()

    print(f"Combined Rₑ plot saved to:\n- {filepath}")




















