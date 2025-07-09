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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import argparse
import os
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


# Plot 1: Timeline of infection for each long shedder
def plot_long_shedders_timeline(ssod):
    
    individuals_data = om.read_individuals_data(ssod)
    # Cleaned data: remove entries without valid infection period
    data = individuals_data.dropna(subset=['t_infection', 't_not_infected'])
    
    long_df = data[data['type'] == 'long_shedder'].copy()
    long_df = long_df.sort_values(by='t_infection')
    plt.figure(figsize=(10, 6))
    
    for idx, row in long_df.iterrows():
        plt.hlines(
            y=idx, xmin=row['t_infection'], xmax=row['t_not_infected'],
            color='red', linewidth=2
        )
        plt.plot(row['t_infection'], idx, 'o', color='black')  # start point
        plt.plot(row['t_not_infected'], idx, 'x', color='gray')  # end point
    
    plt.xlabel("Time")
    plt.ylabel("Long shedder individuals")
    plt.title("Infection duration for long shedders")
    plt.grid(True)
    plt.tight_layout()

    # Determine save path
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    seed = os.path.basename(ssod)
    figure_output_path = os.path.join(experiment_plots_dir, f"{experiment_name}_{seed}_long_shedders_trajectory.tiff")

    # Save the plot 
    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close()

def plot_avg_infection_duration_by_type(ssod):
    
    individuals_data = om.read_individuals_data(ssod)
    # Cleaned data: remove entries without valid infection period
    data = individuals_data.dropna(subset=['t_infection', 't_not_infected'])
    
    df_durations = data.copy()
    df_durations['infection_duration'] = df_durations['t_not_infected'] - df_durations['t_infection']
    df_durations = df_durations.dropna(subset=['infection_duration', 'type'])

    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=df_durations, 
        x='type', 
        y='infection_duration', 
        estimator=np.mean,  # Corrected from 'mean' to np.mean
        ci='sd', 
        palette='muted'
    )
    plt.ylabel("Average Infection Duration")
    plt.title("Average Infection Time: Normal vs Long Shedders")
    plt.tight_layout()
    
    # Determine save path
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    seed = os.path.basename(ssod)
    figure_output_path = os.path.join(experiment_plots_dir, f"{experiment_name}_{seed}_avg_inf_time.tiff")

    # Save the plot 
    print(f"Saving figure to: {figure_output_path}")
    plt.savefig(figure_output_path, format='tiff', dpi=300, bbox_inches='tight')
    plt.close()

def plot(experiment_name, seed_number):
    
    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)
    for sim_out_dir in sim_out_dirs:
        ssod = dm.get_ssod(sim_out_dir, seed_number)
        pm.plot_simulation(ssod, threshold=0)
        plot_long_shedders_timeline(ssod)
        plot_avg_infection_duration_by_type(ssod)

def main():
    # # Set up the argument parser
    # parser = argparse.ArgumentParser(description="Plot")
    # parser.add_argument('experiment_name', type=str, help="experiment name")
    # args = parser.parse_args()
    # # Run the script with the provided parameter
    
    experiment_name = 'test_long_shedders_#2'
    seed_number = 4
    plot(experiment_name, seed_number)
    
if __name__ == "__main__":
    main()
    