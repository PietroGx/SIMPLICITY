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
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    
    # loop over each simulation output folder
    for folder in os.listdir(experiment_output_dir):
        counter = 0
        simulation_directory = os.path.join(experiment_output_dir,folder)
        # loop over each seeded simulation output folder
        for subfolder in os.listdir(simulation_directory):
            seeded_simulation_output_directory = os.path.join(simulation_directory,subfolder)
            if check_seeded_simulation_output(seeded_simulation_output_directory):
                counter += 1
        
        print(f'In {folder}:    {counter}/{n_seeds} simulations were run successfully.')
    print('##########################################')

def main():
    for folder in list_experiment_folders():
        try:
            count_completed_simulations(folder)
        except:
            print('')
            print('')
            print('##########################################')
            print(f'In {folder}:    The experiment FAILED.')
            print('##########################################')
            print('')
    
if __name__ == "__main__":
    main()