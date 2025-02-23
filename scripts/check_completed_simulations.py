# -*- coding: utf-8 -*-

import os 
import simplicity.dir_manager as dm 
import simplicity.settings_manager as sm

def list_experiment_folders():
    # List all folders and exclude those starting with '0'
    folders = [f for f in os.listdir(dm.get_data_dir())
               if os.path.isdir(os.path.join(dm.get_data_dir(), f)) and not f.startswith('0')]
    return sorted(folders)

def check_output_file(directory, filename): 
    file_path = os.path.join(directory, filename)
    if os.path.isfile(file_path):
        return True
    else:
        return False
    
def check_seeded_simulation_output(seeded_simulation_output_directory):
    
    a = check_output_file(seeded_simulation_output_directory, 'final_time.csv')
    b = check_output_file(seeded_simulation_output_directory, 'fitness_trajectory.csv')
    c = check_output_file(seeded_simulation_output_directory, 'individuals_data.csv')
    d = check_output_file(seeded_simulation_output_directory, 'lineage_frequency.csv')
    e = check_output_file(seeded_simulation_output_directory, 'phylogenetic_data.csv')
    f = check_output_file(seeded_simulation_output_directory, 'sequencing_data_regression.csv')
    g = check_output_file(seeded_simulation_output_directory, 'sequencing_data.fasta')
    h = check_output_file(seeded_simulation_output_directory, 'simulation_trajectory.csv')
    
    if a and b and c and d and e and f and g and h:
        return True 
    else:
        return False 

def count_completed_simulations(experiment_name):
    print('')
    print('##########################################')
    print(f'Simulation summary of {experiment_name}:')
    print('')
    n_seeds = sm.read_n_seeds_file(experiment_name)['n_seeds']
    
    # loop over each simulation output folder
    for simulation_output_dir in dm.get_simulation_output_dirs(experiment_name):
        counter = 0
        seeded_simulation_output_dirs = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
        # loop over each seeded simulation output folder
        for seeded_simulation_output_dir in seeded_simulation_output_dirs:
            if check_seeded_simulation_output(seeded_simulation_output_dir):
                counter += 1
                
        folder_name = os.path.basename(simulation_output_dir)  
        print(f'In {folder_name}:    {counter}/{n_seeds} simulations were run successfully.')
    print('##########################################')

def main():
    for folder in list_experiment_folders():
        try:
            count_completed_simulations(folder)
        except Exception as e:
            print('')
            print('')
            print('##########################################')
            print(f'In {folder}:    something went wrong.')
            print(e)
            print('##########################################')
            print('')
    
if __name__ == "__main__":
    main()