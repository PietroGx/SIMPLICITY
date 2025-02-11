#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 13:28:55 2025

@author: pietro
"""
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import simplicity.config as c
import argparse

def load_data(experiment_name):
    
    experiment_output_dir = c.get_experiment_output_dir(experiment_name)
    dic = {}
    
    seeded_simulations_folders = glob.glob(os.path.join(experiment_output_dir, '*/'))
    # iterate over each seeded simulation folder and extract final times 
    for seeded_simulations_folder_path in seeded_simulations_folders:
        final_times = []
        for seeded_simulation_path in glob.glob(os.path.join(seeded_simulations_folder_path, '*/')):
            final_time_file = os.path.join(seeded_simulation_path, 'final_time.csv')
            # Check if the final_time.csv exists in the subfolder
            if os.path.exists(final_time_file):
                # Load the value from final_time.csv 
                try:
                    with open(final_time_file, 'r') as f:
                       final_time_value = f.readline().strip()
                       final_times.append(float(final_time_value))  
                       
                except Exception as e:
                    print(f"Error reading {final_time_file}: {e}")
        
        dic[os.path.basename(os.path.normpath(seeded_simulations_folder_path))] = final_times
    print(dic)
 
    df = pd.DataFrame(dic)
    
    return df

def plot_histograms(experiment_name):
    final_times_data_frames = load_data(experiment_name)
    num_folders = len(final_times_data_frames)
    fig, axes = plt.subplots(num_folders, 1, figsize=(10, 5 * num_folders), squeeze=False)
    
    for ax, (folder_name, data) in zip(axes.flatten(), final_times_data_frames.items()):
        ax.hist(data, bins=30, edgecolor='black')
        ax.set_title(f'Histogram for {folder_name}')
        ax.set_xlabel('Last Time Value')
        ax.set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(os.path.join(c.get_experiment_dir(experiment_name),
                             'simulations_lenght_histogram.png'))

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to plot simulation lenght histogram")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    # Run the script with the provided parameter
    plot_histograms(args.experiment_name)
    
if __name__ == "__main__":
    main()