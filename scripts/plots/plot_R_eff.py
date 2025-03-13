#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import argparse
import os

def select_and_plot_R_effective(experiment_name, selected_seeds):
    SSODs = []
    for simulation_output_dir in dm.get_simulation_output_dirs(experiment_name):
        SSODs += (dm.get_seeded_simulation_output_dirs(simulation_output_dir))
        
        for ssod in SSODs:
            seed = os.path.basename(ssod)
            for selected_seed in selected_seeds:
                if seed == f"seed_{selected_seed:04d}":
                    seeded_simulation_output_dir = ssod
                    # set window size and lineage occurence threshold
                    window_size = 14
                    threshold = 0.1
                    # process R effective data and write csv files needed for plot
                    om.write_R_effective_dfs_csv(experiment_name, seeded_simulation_output_dir, window_size, threshold)
                    pm.plot_R_effective(experiment_name, seeded_simulation_output_dir, window_size, threshold)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('selected_seeds', nargs='+', type=int, help="List of selected seeds (plot relative simulation)")
    args = parser.parse_args()
    # Run the script 
    select_and_plot_R_effective(args.experiment_name, args.selected_seeds)
    
if __name__ == "__main__":
    main()
    