# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om

def check_output_file(directory, filename):
    return os.path.isfile(os.path.join(directory, filename))

def check_seeded_simulation_output(ssod):
    required_files = [
        'final_time.csv', 'fitness_trajectory.csv', 'individuals_data.csv',
        'lineage_frequency.csv', 'phylogenetic_data.csv', 'sequencing_data_regression.csv',
        'sequencing_data.fasta', 'simulation_trajectory.csv',
    ]
    return all(check_output_file(ssod, f) for f in required_files)

def categorize_simulations_by_end_time(ssod_list, final_time_limit, pop_size):
    early, mid, late = 0, 0, 0
    early_died, early_sat = 0, 0
    mid_died, mid_sat = 0, 0
    late_died, late_sat = 0, 0
    early_seeds, first_third_seeds, second_third_seeds, last_third_seeds = [], [], [], []

    for ssod in ssod_list:
        try:
            end_time = om.read_final_time(ssod)
        except Exception:
            continue

        seed = dm.get_seed_from_SSOD(ssod)
        
        is_died_out = False
        is_saturated = False

        # Autopsy ONLY for early terminations
        if end_time < final_time_limit:
            try:
                df = om.read_simulation_trajectory(ssod)
                final_S = df.iloc[-1]['susceptibles']
                threshold = 0.05 * pop_size
                
                if (pop_size - final_S) <= threshold:
                    is_died_out = True
                elif final_S <= threshold:
                    is_saturated = True
            except Exception:
                pass

        if end_time <= 1.0:
            early_seeds.append(seed)

        # Binning into thirds and logging the specific autopsy results
        if end_time < final_time_limit / 3:
            early += 1
            if is_died_out: early_died += 1
            if is_saturated: early_sat += 1
            first_third_seeds.append(seed)
            
        elif end_time < (2 * final_time_limit / 3):
            mid += 1
            if is_died_out: mid_died += 1
            if is_saturated: mid_sat += 1
            second_third_seeds.append(seed)
            
        else:
            late += 1
            # Late bin includes fully completed runs, so we only track autopsies if it ended early
            if end_time < final_time_limit:
                if is_died_out: late_died += 1
                if is_saturated: late_sat += 1
            last_third_seeds.append(seed)

    return (early, mid, late, 
            early_died, early_sat, 
            mid_died, mid_sat, 
            late_died, late_sat, 
            early_seeds, first_third_seeds, second_third_seeds, last_third_seeds)

def format_seed_list(seeds, max_display=5):
    return ', '.join(seeds[:max_display]) if seeds else 'None'

def count_and_analyze_simulations(experiment_name):
    print(f"\n{'#' * 80}")
    print(f"Experiment: {experiment_name}\n")

    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)
    n_seeds = sm.read_n_seeds_file(experiment_name)['n_seeds']

    for sim_out_dir in sim_out_dirs:
        ssod_list = dm.get_seeded_simulation_output_dirs(sim_out_dir)
        valid_ssods = [s for s in ssod_list if check_seeded_simulation_output(s)]
        
        try:
            theoretical_final_time = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'final_time')
            pop_size = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'population_size')
        except Exception as e:
            folder_name = os.path.basename(sim_out_dir)
            print(f"  [!] Could not retrieve parameters for {folder_name}: {e}")
            continue

        (early, mid, late, 
         early_died, early_sat, 
         mid_died, mid_sat, 
         late_died, late_sat, 
         early_seeds, first_third_seeds, second_third_seeds, last_third_seeds) = (
            categorize_simulations_by_end_time(valid_ssods, theoretical_final_time, pop_size)
        )

        folder_name = os.path.basename(sim_out_dir)
        total_valid = len(valid_ssods)

        print(f"  ┌─ {folder_name}")
        print(f"  │  Valid simulations: {total_valid}/{n_seeds}")
        if total_valid > 0:
            # First Third
            print(f"  │  Ended before 1/3  of theoretical end: {early:3d} ({early/total_valid*100:5.1f}%)")
            if early > 0:
                print(f"  │     ↳ Died Out (Extinction):          {early_died:3d}")
                print(f"  │     ↳ Infected Everyone (Saturation): {early_sat:3d}")
                
            # Second Third
            print(f"  │  Ended before 2/3  of theoretical end: {mid:3d} ({mid/total_valid*100:5.1f}%)")
            if mid > 0:
                print(f"  │     ↳ Died Out (Extinction):          {mid_died:3d}")
                print(f"  │     ↳ Infected Everyone (Saturation): {mid_sat:3d}")
                
            # Last Third
            print(f"  │  Ended ≥ 2/3       of theoretical end: {late:3d} ({late/total_valid*100:5.1f}%)")
            # Only show autopsy stats here if there were early terminations in the final stretch
            if (late_died + late_sat) > 0:
                print(f"  │     ↳ (Ended early) Died Out:         {late_died:3d}")
                print(f"  │     ↳ (Ended early) Saturated:        {late_sat:3d}")
                
        print(f"  │  Seeds ending before first day:       {format_seed_list(early_seeds)}")
        print(f"  │  Seeds ending in first third:         {format_seed_list(first_third_seeds)}")
        print(f"  │  Seeds ending in second third:        {format_seed_list(second_third_seeds)}")
        print(f"  │  Seeds ending in last third:          {format_seed_list(last_third_seeds)}")
        print(f"  └{'─' * (len(folder_name) + 2)}")

    print(f"\n{'#' * 80}")

def main():
    parser = argparse.ArgumentParser(description="Check simulations in experiment and get stats")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    count_and_analyze_simulations(args.experiment_name)

if __name__ == "__main__":
    main()
