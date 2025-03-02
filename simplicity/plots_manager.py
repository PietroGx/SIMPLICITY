#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:00:16 2025

@author: pietro
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import math
from simplicity.tuning.evolutionary_rate import create_joint_sequencing_df
from simplicity.tuning.evolutionary_rate import tempest_regression
import glob
import simplicity.tuning.diagnosis_rate as dr
import simplicity.tuning.evolutionary_rate as er
import simplicity.settings_manager as sm
import seaborn as sns

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

    fig = plt.figure(0,figsize=(16, 9))
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
    

def plot_lineage_frequency(simulation_output,threshold):
    '''
    Plot relative frequency of lineages during the simulation

    Parameters
    ----------
    threshold : float
        Do not plot variants under this frequency threshold
    t_final : int
        X axis upper limit.

    Returns
    -------
    None.

    '''
    df = pd.DataFrame(simulation_output.variant_inf_count)
    df.fillna(0,inplace=True)
    df_normalized = df.div(df.sum(axis=0), axis=1)
    # # remove last value (for correct plot visualization)
    # df_normalized = df_normalized.drop(df_normalized.columns[-1], axis=1)
    # set last column values for better visualization
    df_normalized.iloc[:, -1] = df_normalized.iloc[:, -2]
    
    df_cut = df_normalized[df_normalized.gt(threshold).any(axis=1)]
    df_cut = df_cut.transpose()
    # ax = df_cut.plot(kind='area', stacked=True, colormap='gist_rainbow')
    
    # Reverse the column order
    df_cut_reversed = df_cut[df_cut.columns[::-1]]
    
    # Plot with transparency
    ax = df_cut_reversed.plot(kind='area', stacked=True, colormap='gist_rainbow', alpha=0.5)
            
    ax.set_xlim([0,simulation_output.time])
    # Set tick label sizes for both axes
    ax.tick_params(axis='both', labelsize=20)
    plt.show()
    
def plot_simulation(seeded_simulation_output_dir, threshold):
    '''
    joins plots of population trajectory and lineages frequency in a single 
    plot

    '''
    threshold = 0.05 # dont show lineages under this frequency
    from matplotlib.ticker import FuncFormatter
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(12, 8))

    # Create a gridspec object with 3 rows. 
    # The height_ratios argument determines the relative heights of the plots.
    # Here, 1:2 ratio for the plots, and an empty space with '0' ratio.
    gs = GridSpec(3, 1, height_ratios=[1,1,1])  
    
    # Create the subplots using the gridspec indices
    ax2 = fig.add_subplot(gs[0, 0])  
    ax1 = fig.add_subplot(gs[1, 0],sharex=ax2) 
    ax3 = fig.add_subplot(gs[2, 0],sharex=ax2) 
    
    # Adjust subplot parameters
    fig.subplots_adjust(
        top=0.95,
        bottom=0.085,
        left=0.07,
        right=0.81,
        hspace=0.1,
        wspace=0.085
    )
    
    # system trajectory subplot
    ######################################################################
    trajectory_df = om.read_simulation_trajectory(seeded_simulation_output_dir)
    
    time = trajectory_df['time'].tolist()
    infected   = trajectory_df['infected'].tolist()
            
    ax2.plot(time,infected,
            label='Infected individuals', 
            color= 'red')
    
    ax2.set_ylim(0,max(infected)*1.5)
    ax2.tick_params(axis='both', labelsize=15)
    ax2.set_ylabel('N individuals', fontsize=15)
    ax2.legend(loc='upper left', fontsize=15)
    
    ######################################################################
    # variants frequency subplot
    lineage_frequency_file_path = os.path.join(seeded_simulation_output_dir, 'lineage_frequency.csv')
    lineage_frequency_df = pd.read_csv(lineage_frequency_file_path)
    
    lineage_frequency_df.fillna(0,inplace=True)
    lineage_frequency_df_normalized = lineage_frequency_df.div(lineage_frequency_df.sum(axis=0), axis=1)
    # set last column values for better visualization
    lineage_frequency_df_normalized.iloc[:, -1] = lineage_frequency_df_normalized.iloc[:, -2]
    
    df_cut = lineage_frequency_df_normalized[lineage_frequency_df_normalized.gt(threshold).any(axis=1)]
    df_cut = df_cut.transpose()   
    # Reverse the column order
    df_cut_reversed = df_cut[df_cut.columns[::-1]]
    # Plot with transparency
    df_cut_reversed.plot(kind='area', stacked=True, colormap='gist_rainbow',
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
    
    #######################################################################
    # fitness subplot
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
                     color='#2ca02c', alpha=0.3, label='Std')

    # Add labels and title
    # ax3.xlabel('Time')
    ax3.tick_params(axis='both', labelsize=15)
    ax3.set_xlabel('Time - days', fontsize=15)
    ax3.set_xticks(np.arange(0,time_final,10))
    ax3.set_xticklabels(np.arange(0,time_final,10).astype(int))
    ax3.set_ylabel('Fitness score',fontsize=15)
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    # ax3.yaxis.set_label_coords(-0.1, 0)
    # Add legend
    ax3.legend(loc='upper left', fontsize=15)
    
    #######################################################################
    # Align the y-axis labels
    fig.align_ylabels([ax2, ax3])
    plt.tight_layout()
    figure_output_path = os.path.join(seeded_simulation_output_dir,"simulation_trajectory.png")
    plt.savefig(figure_output_path)
    plt.close()

###############################################################################
###############################################################################
######################## PLOT TEMPEST REGRESSION ##############################
###############################################################################
###############################################################################

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
    # get observed_evolutionary_rate value
    observed_evolutionary_rate = fitted_tempest_regression.coef_[0] # substitution rate per site per year
    # plot data points
    sequencing_data_df.plot(kind='scatter', x='Sequencing_time', y='Distance_from_root', color='blue', ax=ax)
    # plot linear regression
    ax.plot(x, y_pred, color='red', linewidth=2, 
            label=f'observed_evolutionary_rate = {observed_evolutionary_rate:.5f}')
    ax.set_xlabel('Time [y]')
    ax.set_ylabel('Distance from root [#S/site]')
    ax.set_xlim(left=0)
    ax.grid(True)
    ax.legend()
 
def plot_combined_tempest_regressions(experiment_name, parameter, min_sim_lenght=0, y_axis_max=0.1):
    # Get sorted simulation output directories for experiment
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
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
        sequencing_data_df = create_joint_sequencing_df(simulation_output_dir, min_sim_lenght)
        if sequencing_data_df is None: 
            pass
        else:
            fitted_tempest_regression = tempest_regression(sequencing_data_df)
            ax = axs[i // num_cols][i % num_cols]
            plot_tempest_regression(sequencing_data_df,
                                       fitted_tempest_regression,
                                       ax)
            ax.set_title(f'Regression - {parameter}: {param}')
            ax.set_ylim(0, y_axis_max)
    
    plt.tight_layout()
    plt.savefig(os.path.join(experiment_output_dir, f"{experiment_name}_combined_regression.png"))

def plot_combined_observed_evolutionary_rate_vs_parameter(experiment_name, parameter, min_sim_lenght=0):
    ''' Plot observed evolutionary rate (tempest regression) against desired parameter values
    '''
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    # build combined dataframe, filtered by min simulation lenght
    om.build_combined_observed_evolutionary_rate_vs_parameter_df(experiment_name, 
                                                       parameter, 
                                                       min_sim_lenght)
    df = om.read_combined_observed_evolutionary_rate_csv(experiment_name, parameter,min_sim_lenght)
    # Plot target parameter vs u as a line plot with points
    plt.figure(figsize=(10, 6))
    plt.plot(df[parameter],  df['observed_evolutionary_rate'], 
             marker='o', 
             color='black', 
             linestyle='-', 
             label=f'{experiment_name}_{parameter} vs observed_evolutionary_rate')
    plt.xlabel(parameter)
    plt.ylabel('u')
    plt.title(f'{parameter} vs observed_evolutionary_rate')
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(experiment_output_dir, 
                             f"{experiment_name}_{parameter}_vs_observed_evolutionary_rate.png"))

def plot_observed_evolutionary_rates_vs_parameter_scatter(experiment_name, parameter, min_sim_lenght=0):
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    # build combined dataframe, filtered by min simulation lenght
    om.build_observed_evolutionary_rates_vs_parameter_df(experiment_name, 
                                                       parameter, 
                                                       min_sim_lenght)
    observed_evolutionary_rate_vs_parameter_df = om.read_observed_evolutionary_rates_csv(experiment_name, 
                                                                                         parameter,
                                                                                         min_sim_lenght)
    # Create figure and axes
    fig, ax = plt.subplots(3,1, figsize=(8, 10))
    # First scatter plot
    sns.scatterplot(x=parameter, y='observed_evolutionary_rate', label='Data', 
                    color='blue', alpha=0.5, ax=ax[0],
                    data=observed_evolutionary_rate_vs_parameter_df)
    ax[0].set_xlabel(f'{parameter}')
    ax[0].set_ylabel('Observed Evolutionary Rate')
    
    # Second scatter plot
    sns.scatterplot(x=parameter, y='observed_evolutionary_rate', label='Data', 
                    color='blue', alpha=0.5, ax=ax[1],
                    data=observed_evolutionary_rate_vs_parameter_df)
    ax[1].set_xlabel(f'{parameter}')
    ax[1].set_ylabel('Observed Evolutionary Rate')
    ax[1].set_xscale('log')
    
    # Third scatter plot 
    sns.scatterplot(x=parameter, y='observed_evolutionary_rate', label='Data', 
                    color='blue', alpha=0.5, ax=ax[2],
                    data=observed_evolutionary_rate_vs_parameter_df)
    ax[2].set_xlabel(f'{parameter}')
    ax[2].set_ylabel('Observed Evolutionary Rate')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(experiment_output_dir, 
            f"{experiment_name}_{parameter}_vs_observed_evolutionary_rate_scatter.png"))
    
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(experiment_output_dir,))
    
def export_tempest_regression_plots(experiment_name): 
    ''' move tempest regression plots from experiment folder to plots folder
    '''
    # get experiment output dir
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    # get plots directories
    plots = glob.glob(os.path.join(experiment_output_dir, '*.png'))
    # create target directories to move the plots to
    plots_folder_dir = os.path.join(dm.get_data_dir(), '00_Tempest_regression_plots')
    os.makedirs(plots_folder_dir,exist_ok=True)
    # get plots filenames
    plot_filenames = [os.path.basename(plot) for plot in plots]
    # move the plots
    for plot,plot_filename in zip(plots,plot_filenames):
        os.replace(plot, os.path.join(plots_folder_dir,plot_filename))
        

def plot_combined_observed_evolutionary_rate_fit(experiment_name, fit_result, model_type, min_sim_lenght):
    ''' plot fit of evolutionary rate / observed evolutionary rate curve
    '''

    data = om.read_combined_observed_evolutionary_rate_csv(experiment_name, 'evolutionary_rate', min_sim_lenght)
    x_data = data['evolutionary_rate'] 
    y_data = data['observed_evolutionary_rate']  
    
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
    # ax[0].set_ylabel('Observed evolutionary rate (site/year) (u)')
    ax[0].legend()
    # semilog scale
    ax[1].plot(x_data, fit_result.best_fit, label=f'Fitted {model_type} curve',
               color='red', linewidth=2)
    ax[1].set_xscale("log")
    ax[1].set_ylim(0)
    # ax[1].set_xlabel('Substitution Rate (site/year) (e)')
    ax[1].set_ylabel('Observed evolutionary rate (site/year) (u)')
    ax[1].legend()
    # loglog scale
    ax[2].plot(x_data, fit_result.best_fit, label=f'Fitted {model_type} curve',
               color='red', linewidth=2)
    ax[2].set_yscale("log")
    ax[2].set_xscale("log")
    ax[2].set_xlabel('Substitution Rate (site/year) (e)')
    # ax[2].set_ylabel('Observed evolutionary rate (site/year) (u)')
    ax[2].legend()
    
    # fig.legend(loc="upper center", ncol=2, fontsize=10)
    
    plt.title('')
    # plt.subplots_adjust(top=0.6)
    plt.tight_layout()
    file_path = os.path.join(dm.get_experiment_dir(experiment_name),
     f'{experiment_name}_combined_observed_evolutionary_rate_{model_type}_fit.png')
    plt.savefig(file_path)

def plot_observed_evolutionary_rates_fit(experiment_name, fit_result, model_type,min_sim_lenght):
    ''' plot fit of evolutionary rate / observed evolutionary rates curve
    '''
    parameter = 'evolutionary_rate'
    
    line_color = 'black' #'#DE8F05' # orange
    scatter_color = '#DE8F05'# orange
    combined_OER_marker_color = '#E64B9D' # pink
    scatter_color_2 = '#0173B2' # blue '#029E73' # green
    
    # Create figure and axes
    fig, ax = plt.subplots(3,1, figsize=(8, 12))
    
    # import combined regression data
    combined_data = om.read_combined_observed_evolutionary_rate_csv(experiment_name,parameter,min_sim_lenght)
    
    # import single simulations regression data
    data = om.read_observed_evolutionary_rates_csv(experiment_name, parameter,min_sim_lenght)
    # Group by evolutionary_rate and compute mean and standard deviation for OER
    data_mean_df = om.get_mean_std_observed_evolutionary_rates(experiment_name,parameter,min_sim_lenght)
    # get lower and upper confidence interval for fit results
    x_data = data['evolutionary_rate']
    x, lower_curve, upper_curve = confidence_interval_fit(model_type, fit_result, x_data.to_numpy())
    
    # plot on each axis of subplots in a loop
    for a in ax:
        # Fill between the upper and lower curves for the confidence interval region
        a.fill_between(x,lower_curve, upper_curve, 
                       color=line_color, alpha=0.3, label='95% Confidence Interval',
                       zorder=-1)
        # scatterplot Estimated OER - single simulation
        sns.scatterplot(x=parameter, y='observed_evolutionary_rate', 
                        label='Estimated OER - single simulation', data=data,
                        color=scatter_color, alpha=0.5, ax=a,
                        zorder=0)
        # plot fitted curve
        sns.lineplot(x=x_data, y=fit_result.best_fit, 
                     label=f'Fitted {model_type} curve', 
                     color=line_color, linewidth=1, ax=a,
                     zorder=1)
        # scatterplot combined regression points (as comparison)
        sns.scatterplot(x='evolutionary_rate', y='observed_evolutionary_rate', marker='X',
            label='Combined tempest regression estimate of OER', data=combined_data,
            color=combined_OER_marker_color,alpha=1, ax=a,
            zorder=2)
        # plot mean of observed_evolutionary_rate from data_mean_std
        sns.scatterplot(x=parameter, y='mean', marker = 'X',
                        label='Mean of estimated OER (single simulations)', data=data_mean_df,
                        color=scatter_color_2, alpha=1, ax=a,
                        zorder=3)
        
        minsimlenghts = [0]
        palette = sns.color_palette("tab10", len(minsimlenghts))
        for i, minsimlenght in enumerate(minsimlenghts):
            df = om.get_mean_std_observed_evolutionary_rates(experiment_name,parameter,minsimlenght)
            # plot mean of observed_evolutionary_rate from data_mean_std
            sns.scatterplot(x=parameter, y='mean', marker = 'X',
                            label=f'Mean of estimated OER - min {minsimlenght} d', data=df,
                            color=palette[i], alpha=0.8, ax=a,
                            zorder=4)

    # First plot (linear scale) -----------------------------------------------
    ax[0].set_xlabel(f'{parameter}')
    ax[0].set_ylabel('Observed Evolutionary Rate')
    
    # Second plot (semilog scale) -----------------------------------------------
    
    ax[1].set_xlabel(f'{parameter}')
    ax[1].set_ylabel('Observed Evolutionary Rate')
    ax[1].set_xscale('log')
    ax[1].legend_.remove()
    
    # Third plot (log log scale) -----------------------------------------------
    ax[2].set_xlabel(f'{parameter}')
    ax[2].set_ylabel('Observed Evolutionary Rate')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].legend_.remove()
    
    plt.tight_layout()
    plt.savefig(os.path.join(dm.get_experiment_dir(experiment_name), 
        f"{experiment_name}_observed_evolutionary_rates_{model_type}_fit.png"))
    
def confidence_interval_fit(model_type, fit_result, x):
    params_lower = {}
    params_upper = {}
    
    for param in fit_result.params:
        param_value = fit_result.params[param].value
        param_stderr = fit_result.params[param].stderr if fit_result.params[param].stderr else 0
        ci_lower = param_value - 1.96 * param_stderr
        ci_upper = param_value + 1.96 * param_stderr
        
        params_lower[param] = ci_lower
        params_upper[param] = ci_upper
        
    def remove_duplicates_array(arr):
        _, idx = np.unique(arr, return_index=True)
        return arr[np.sort(idx)]

    x = remove_duplicates_array(x)

    # Create the upper and lower bound curves for confidence intervals
    upper_curve = er.evaluate_model(model_type, params_upper, x)
    lower_curve = er.evaluate_model(model_type, params_lower, x)
    return x, lower_curve, upper_curve
 
###############################################################################
###############################################################################   
def plot_w_t(weights, t_max, params, t_eval): # plot weights of phenotype model
    '''plot the weight function
    from weight.py file in phenotype:
    # params = w_t_params() 
    # w_t = weights(t, t_eval, k_e, k_a, t_max)
    '''
    k_e = params[0]
    k_a = params[1]
    
    # Time points 
    t = np.linspace(0,300,1000)

    # Calculate concentration for each time point
    w_t = weights(t, t_eval, k_e, k_a, t_max)

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
    
###############################################################################
###############################################################################
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
    
###############################################################################
###############################################################################

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
    
###############################################################################
###############################################################################

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
    
###############################################################################
###############################################################################

def plot_IH_lineage_distribution(experiment_name):
    fig, axes  = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
    df = om.get_IH_lineages_data_experiment(experiment_name)
    
    ax0 = sns.barplot(x='IH_virus_number', y='ih_virus_count', data=df, hue='IH_virus_emergence_rate', 
                ax=axes[0], dodge=True, palette="Set2", alpha=0.7)
    ax1 = sns.barplot(x='lineages_number', y='ih_lineage_count', data=df, hue='IH_virus_emergence_rate',
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
    
    axes[0].set_title('Count of IH_virus_number')
    axes[0].set_xlabel('IH_virus_number')
    axes[0].set_ylabel('Count')

    axes[1].set_title('Count of lineages_number')
    axes[1].set_xlabel('lineages_number')
    axes[1].set_ylabel('Count')
    
    
    plt.suptitle('Histogram of IH virus and lineage distribution')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    fig_path = os.path.join(experiment_output_dir,'IH variability plot.png')
    plt.savefig(fig_path)

###############################################################################
###############################################################################
def plot_histograms(experiment_name, final_times_data_frames):
    ''' plot histograms of simulation final times.
    Called from statistics_simulation_lenght.py
    '''
    num_folders = len(final_times_data_frames.columns)
    fig, axes = plt.subplots(num_folders, 1, figsize=(10, 5 * num_folders), squeeze=False)
    
    for ax, (folder_name, data) in zip(axes.flatten(), final_times_data_frames.items()):
        ax.hist(data, bins=30, edgecolor='black')
        ax.set_title(f'Histogram for {folder_name}')
        ax.set_xlabel('Last Time Value')
        ax.set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(os.path.join(dm.get_experiment_dir(experiment_name),
                             'simulations_lenght_histogram.png'))
   





























