#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import argparse
import os

def select_and_plot_R_effective(experiment_name, window_size, lin_threshold, selected_seeds):
    SSODs = []
    # get all the seeded simulation output directories
    for simulation_output_dir in dm.get_simulation_output_dirs(experiment_name):
        SSODs += (dm.get_seeded_simulation_output_dirs(simulation_output_dir))
        # loop over the ssods and only plot the selected seeds folders
        for ssod in SSODs:
            seed = os.path.basename(ssod)
            for selected_seed in selected_seeds:
                if seed == f"seed_{selected_seed:04d}":
                    seeded_simulation_output_dir = ssod
                    print('')
                    print(f'Processing {ssod} R_effective data...')
                    # process R effective data and write csv files needed for plot
                    om.write_R_effective_dfs_csv(experiment_name, 
                                                 seeded_simulation_output_dir, 
                                                 window_size, lin_threshold)
                    pm.plot_R_effective(experiment_name, 
                                        seeded_simulation_output_dir, 
                                        window_size, lin_threshold)
                    print('Plot R effective: completed.')
                    print('')

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('window_size', type=int, help="time window_size")
    parser.add_argument('lin_threshold', type=float, help="lineage occurence threshold")
    parser.add_argument('selected_seeds', nargs='+', type=int, help="List of selected seeds (plot relative simulation)")
    args = parser.parse_args()
    
    # Run the script 
    select_and_plot_R_effective(args.experiment_name, 
                                args.window_size,
                                args.lin_threshold,
                                args.selected_seeds)
    
if __name__ == "__main__":
    main()
    