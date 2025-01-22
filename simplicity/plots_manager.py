#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:00:16 2025

@author: pietro
"""
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

def plot_fitness(simulation_output):
    """
    Plots the fitness trajectory with mean and standard deviation.
    
    """
    # Extract time, mean, and standard deviation data
    times = [coord[0] for coord in simulation_output.fitness_trajectory]
    means = [coord[1][0] for coord in simulation_output.fitness_trajectory]
    stds = [coord[1][1] for coord in simulation_output.fitness_trajectory]

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
    
def plot(simulation_output):
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

def plot_var_frequency(simulation_output,threshold):
    '''
    Plot relative frequency of variants during the simulation

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