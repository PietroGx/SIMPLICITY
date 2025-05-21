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
                    try:
                        pm.plot_R_effective(experiment_name, 
                                            seeded_simulation_output_dir, 
                                            window_size, lin_threshold)
                        print('Plot R effective: success.')
                        print('')
                    except Exception as e:
                        raise ValueError(f'Plotting failed, check dataframes! Error: {e}')
                        break
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
    