#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 17:15:49 2025

@author: pietro
"""
import simplicity.dir_manager as dm
import os
import shutil
import argparse

def move_seed_plots(seed_number, search_folder, destination_folder):
    # Convert the number to a four-digit string
    seed_str = f"seed_{seed_number:04d}"
    
    # Walk through all subdirectories and move plots that contain seed_str
    for root, _, files in os.walk(search_folder):
       for file in files:
           if file.endswith(".png") and seed_str in file:
               source_path = os.path.join(root, file)
               dest_path = os.path.join(destination_folder, file)
               shutil.move(source_path, dest_path)
               print(f"Moved: {source_path} -> {dest_path}")

def search_plots_and_trees(experiment_name, seed_number):
    
    # Create destination folder path
    experiment_dir = dm.get_experiment_dir(experiment_name)
    ordered_plots_dir = os.path.join(experiment_dir, '07_Ordered_seeds_plots')
    os.makedirs(ordered_plots_dir, exist_ok=True)
    destination_folder = os.path.join(ordered_plots_dir, f"plots_seed_{seed_number}")
    os.makedirs(destination_folder, exist_ok=True)
    
    # Define search folders
    plots_search_folder = dm.get_experiment_plots_dir(experiment_name)
    tree_search_folder = dm.get_experiment_tree_dir(experiment_name)
    
    # Move matching files
    move_seed_plots(seed_number, plots_search_folder, destination_folder)
    move_seed_plots(seed_number, tree_search_folder, destination_folder)

def main():
    parser = argparse.ArgumentParser(description="Move PNG files containing specific seed numbers in their filename to 07_Ordered_seeds_plots.")
    parser.add_argument("experiment_name", type=str, help="Experiment name")
    parser.add_argument("seeds_list", type=int, nargs='+', help="List of seed numbers")
    args = parser.parse_args()
    
    seeds_list = args.seeds_list
    
    for seed in seeds_list:
        search_plots_and_trees(args.experiment_name, seed)

if __name__ == "__main__":
    main()
    

