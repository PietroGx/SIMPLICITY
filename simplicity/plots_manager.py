#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:00:16 2025

@author: pietro
"""

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import json
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import math
from simplicity.tuning.evolutionary_rate import create_joint_sequencing_df
from simplicity.tuning.evolutionary_rate import tempest_regression
import glob


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
    
def plot_trajectory(simulation_output):
    '''
    Plots the simulation trajectory and the evolution of infectivity
    during the simulation.

    '''
    time       = [i[0] for i in simulation_output.trajectory]
    infected   = [i[1] for i in simulation_output.trajectory]
    diagnosed  = [i[2] for i in simulation_output.trajectory]
    recovered  = [i[3] for i in simulation_output.trajectory]
    
    fig = plt.figure(0,figsize=(16, 9))
    plt.title('Simulation of SARS-CoV-2 outbreak')
    
    ax = fig.add_subplot(1,1,1)
    
    
    ax.plot(time,recovered,
        label='Recovered (cumulative)', 
        color= 'black')
    ax.plot(time,infected,
            label='Infected individuals', 
            color= 'red')
    ax.plot(time,diagnosed,
            label='Diagnosed individuals (cumulative)', 
            color= 'green')
    
    ax.set_xlabel('Time - days')
    ax.set_ylabel('Number of individuals in compartment')
    ax.set_xlim(0,time[-1])
    ax.set_ylim(0)
    plt.tight_layout()
    plt.legend()
    
    plt.close()

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
    
def plot_simulation(output_directory, threshold):
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
    trajectory_file_path = os.path.join(output_directory, 'simulation_trajectory.csv')
    trajectory_df = pd.read_csv(trajectory_file_path)
    
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
    lineage_frequency_file_path = os.path.join(output_directory, 'lineage_frequency.csv')
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
    time_file_path = os.path.join(output_directory, 'final_time.csv')
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
    fitness_trajectory_file_path = os.path.join(output_directory, 'fitness_trajectory.csv')
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
    figure_output_path = os.path.join(output_directory,"simulation_trajectory.png")
    plt.savefig(figure_output_path)
    plt.close()
    # fig.show()

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
 
def plot_combined_regressions(experiment_name, parameter, min_sim_lenght=0, y_axis_max=0.1):
    # Get experiment output dir
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    # Get seeded simulations output subfolders
    seeeded_simulations_output_directories = [os.path.join(experiment_output_dir, 
                                    f.name) for f in os.scandir(experiment_output_dir
                                                                ) if f.is_dir()]
   
    def extract_parameter(seeeded_simulations_output_directory,parameter):
        json_file = seeeded_simulations_output_directory.replace(
                    '04_Output','02_Simulation_parameters') +'.json'
        try:
            with open(json_file, 'r') as file:
                data = json.load(file)
                return data.get(parameter)
        except Exception as e:
            print(f"Error reading JSON file: {e}")
            return None
        
    # Sort the subdirectories based on the evolutionary rate
    sorted_dirs = sorted(seeeded_simulations_output_directories, 
                     key=lambda dir: extract_parameter(dir, parameter))
    # Determine the number of rows and columns for subplots
    num_rows, num_cols = ideal_subplot_grid(len(seeeded_simulations_output_directories))

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    
    # If there's only one subplot, axs is not an array, so we convert it to a 2D array for consistency
    if num_rows * num_cols == 1:
        axs = [[axs]]
    elif num_rows == 1 or num_cols == 1:
        axs = axs.reshape(1, -1)
    else:
        axs = axs
    
    for i, subdir in enumerate(sorted_dirs):
        param = extract_parameter(subdir, parameter)
        combined_df = create_joint_sequencing_df(subdir, min_sim_lenght)
        if combined_df is None: 
            pass
        else:
            u, model = tempest_regression(combined_df)
            
            x = combined_df['Sequencing_time'].values.reshape(-1, 1)
            y_pred = model.predict(x)
            
            ax = axs[i // num_cols][i % num_cols]
            combined_df.plot(kind='scatter', x='Sequencing_time', y='Distance_from_root', color='blue', ax=ax)
            ax.plot(x, y_pred, color='red', linewidth=2, label=f'u = {u:.5f}')
            ax.set_xlabel('Time [y]')
            ax.set_ylabel('Distance from root [#S/site]')
            ax.set_xlim(left=0)
            ax.set_ylim(0, y_axis_max)
            ax.set_title(f'Regression - {parameter}: {param}')
            ax.grid(True)
            ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(experiment_output_dir, f"{experiment_name}_combined_regression.png"))

def plot_u_vs_parameter(experiment_name, parameter, min_sim_lenght=0):
    ''' Plot observed evolutionary rate (tempest regression) against desired parameter values
    '''
    # get simulation parameter files for the selected experiment
    simulation_parameters_dir = dm.get_simulation_parameters_dir(experiment_name)
    simulation_parameters_files = glob.glob(os.path.join(simulation_parameters_dir, '*.json'))
    # get output directory
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    # Get seeded simulations output subfolders
    seeeded_simulations_output_directories = [os.path.join(experiment_output_dir, 
                                    f.name) for f in os.scandir(experiment_output_dir
                                                                ) if f.is_dir()]
    # for each simulation set in the experiment perfom the tempest regression
    results = []
    for simulation_parameters_file, seeeded_simulations_output_directory in zip(
            simulation_parameters_files,seeeded_simulations_output_directories):
        # Read the parameter from the settings file
        with open(simulation_parameters_file, 'r') as file:
            data = json.load(file)
        parameter_value = data.get(parameter)
        
        # Perform regression for the settings output folder
        combined_df = create_joint_sequencing_df(seeeded_simulations_output_directory, min_sim_lenght)
        if combined_df is None: 
            pass
        else:
            u, _ = tempest_regression(combined_df)
            results.append({str(parameter): parameter_value, 'u': u})
    # add results to df
    results_df = pd.DataFrame(results)
    
    # Sort the results by the 'parameter' values
    results_df = results_df.sort_values(by=str(parameter))
    
    # Save the results to a CSV file
    csv_file_path = os.path.join(experiment_output_dir, 
                                 f'{experiment_name}_u_vs_{parameter}_values.csv')
    results_df.to_csv(csv_file_path, index=False)
    
    # Plot target parameter vs u as a line plot with points
    plt.figure(figsize=(10, 6))
    plt.plot(results_df[str(parameter)], results_df['u'], 
             marker='o', 
             color='black', 
             linestyle='-', 
             label=f'{experiment_name}_{parameter} vs u')
    plt.xlabel(parameter)
    plt.ylabel('u')
    plt.title(f'{parameter} vs u')
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(experiment_output_dir, 
                             f"{experiment_name}_{parameter}_vs_u.png"))
    
def export_u_regression_plots(experiment_name): 
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
        
###############################################################################
###############################################################################
def plot_histograms(experiment_name, final_times_data_frames):
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

def plot_u_fit(experiment_name,fit_result,scale:str):
    if scale == 'loglog' or 'semilog' or 'lin':
        pass
    else:
        raise ValueError('invalid scale settigs.')
        
    data = om.read_u_e_values(experiment_name)
    x_data = data['evolutionary_rate'] 
    y_data = data['u']  
    
    # Create figure and axes
    fig, ax = plt.subplots()
    ax.scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
    # print(x_data)
    # print(y_data)
    # print(fit_result.best_fit)
    # Set log scales
    
    if scale == 'loglog':
        ax.plot(x_data, fit_result.best_fit, label='Fitted curve', color='red', linewidth=2)
        ax.set_yscale("log")
        ax.set_xscale("log")
    elif scale == 'semilog':
        ax.plot(x_data, fit_result.best_fit, label='Fitted curve', color='red', linewidth=2)
        ax.set_xscale("log")
        ax.set_ylim(0)
    else:
        ax.plot(x_data, fit_result.best_fit, label='Fitted curve', color='red', linewidth=2)
        ax.set_ylim(0)
        ax.set_xlim(0)
    ax.set_xlabel('Evolutionary Rate')
    ax.set_ylabel('u')
    ax.legend()
    plt.title('Logarithmic Fit to Data')
    file_path = os.path.join(dm.get_experiment_dir(experiment_name),f'ue_fitting_{scale}.png')
    plt.savefig(file_path)
 
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







































