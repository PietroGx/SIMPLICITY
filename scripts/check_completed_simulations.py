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

# -*- coding: utf-8 -*-

import os
import sys
import argparse
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om

def check_output_file(directory, filename):
    return os.path.isfile(os.path.join(directory, filename))

def check_seeded_simulation_output(ssod):
    required_files = [
        'final_time.csv',
        'fitness_trajectory.csv',
        'individuals_data.csv',
        'lineage_frequency.csv',
        'phylogenetic_data.csv',
        'sequencing_data_regression.csv',
        'sequencing_data.fasta',
        'simulation_trajectory.csv',
    ]
    return all(check_output_file(ssod, f) for f in required_files)

def categorize_simulations_by_end_time(ssod_list, final_time_limit):
    early, mid, late = 0, 0, 0
    early_seeds = []
    first_third_seeds, second_third_seeds, last_third_seeds = [], [], []

    for ssod in ssod_list:
        if not check_seeded_simulation_output(ssod):
            continue

        try:
            end_time = om.read_final_time(ssod)
        except Exception:
            continue

        seed = os.path.basename(ssod)

        if end_time <= 1.0:
            early_seeds.append(seed)

        if end_time < final_time_limit / 3:
            early += 1
            first_third_seeds.append(seed)
        elif end_time < (2 * final_time_limit / 3):
            mid += 1
            second_third_seeds.append(seed)
        else:
            late += 1
            last_third_seeds.append(seed)

    return early, mid, late, early_seeds, first_third_seeds, second_third_seeds, last_third_seeds

def format_seed_list(seeds, max_display=5):
    """Extract numeric suffix from seed names and return a comma-separated string."""
    import re
    numbers = [re.sub(r'\D+', '', s) for s in seeds[:max_display]]
    return ', '.join(numbers) if numbers else 'None'

def count_and_analyze_simulations(experiment_name):
    print(f"\n{'#' * 80}")
    print(f"Experiment: {experiment_name}\n")

    n_seeds = sm.read_n_seeds_file(experiment_name)['n_seeds']
    sim_out_dirs = dm.get_simulation_output_dirs(experiment_name)

    for sim_out_dir in sim_out_dirs:
        ssod_list = dm.get_seeded_simulation_output_dirs(sim_out_dir)
        valid_ssods = [s for s in ssod_list if check_seeded_simulation_output(s)]

        try:
            theoretical_final_time = sm.get_parameter_value_from_simulation_output_dir(sim_out_dir, 'final_time')
        except Exception as e:
            print(f"  [!] Could not retrieve final_time for {os.path.basename(sim_out_dir)}: {e}")
            continue

        early, mid, late, early_seeds, first_third_seeds, second_third_seeds, last_third_seeds = (
            categorize_simulations_by_end_time(valid_ssods, theoretical_final_time)
        )

        folder_name = os.path.basename(sim_out_dir)
        total_valid = len(valid_ssods)

        print(f"  ┌─ {folder_name}")
        print(f"  │  Valid simulations: {total_valid}/{n_seeds}")
        print(f"  │  Ended before 1/3  of theoretical end: {early:3d} ({early/total_valid*100:5.1f}%)")
        print(f"  │  Ended before 2/3  of theoretical end: {mid:3d} ({mid/total_valid*100:5.1f}%)")
        print(f"  │  Ended ≥ 2/3       of theoretical end: {late:3d} ({late/total_valid*100:5.1f}%)")
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
    experiment_name = args.experiment_name
    
    try:
        experiment_dir = dm.get_experiment_dir(experiment_name)
        if not os.path.isdir(experiment_dir):
            raise ValueError("Provided experiment name does not map to a valid directory.")
        count_and_analyze_simulations(experiment_name)
    except Exception as e:
        print(f"\n[!] Error processing experiment '{experiment_name}': {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
