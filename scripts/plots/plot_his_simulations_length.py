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


import os
import glob
import pandas as pd
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.plots_manager as pm
import argparse

def load_data(experiment_name):
    """
    Loads final_time.csv values for each seeded simulation.
    Sorts the columns by R value and returns a DataFrame.
    """
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    data = []

    # Get all sim_out_dirs (e.g., phenotype_model_R folders)
    sim_out_dirs = glob.glob(os.path.join(experiment_output_dir, '*/'))

    for sim_out_dir in sim_out_dirs:
        final_times = []
        # Get all ssods within this sim_out_dir
        for ssod in glob.glob(os.path.join(sim_out_dir, '*/')):
            final_time_file = os.path.join(ssod, 'final_time.csv')
            if os.path.exists(final_time_file):
                try:
                    with open(final_time_file, 'r') as f:
                        final_time_value = f.readline().strip()
                        final_times.append(float(final_time_value))  
                except Exception as e:
                    print(f"Error reading {final_time_file}: {e}")
        
        # Try to get R value from sim_out_dir
        try:
            R_val = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, "R")
        except:
            R_val = None
        
        sim_dir_name = os.path.basename(os.path.normpath(sim_out_dir))
        data.append({
            "sim_dir": sim_dir_name,
            "R": R_val,
            "final_times": final_times
        })

    # Sort the data by R
    sorted_data = sorted(data, key=lambda x: x["R"] if x["R"] is not None else float('inf'))

    # Build the DataFrame (columns are sim_dir names)
    df = pd.DataFrame({
        entry["sim_dir"]: entry["final_times"]
        for entry in sorted_data
    })

    # Store R values for ordering in plotting
    df.attrs["R_order"] = [entry["R"] for entry in sorted_data]

    return df

def plot_histograms(experiment_name):
    """
    Loads data and passes to plots_manager to create histograms ordered by R.
    """
    final_times_df = load_data(experiment_name)
    r_order = final_times_df.attrs.get("R_order", None)
    pm.plot_histograms(experiment_name, final_times_df, r_order=r_order)

def main():
    parser = argparse.ArgumentParser(description="Plot histogram of simulation lengths.")
    parser.add_argument('experiment_name', type=str, help="Name of the experiment")
    args = parser.parse_args()

    plot_histograms(args.experiment_name)

if __name__ == "__main__":
    main()
