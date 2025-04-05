#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:07:15 2025

@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 

 
If you want to change any, you can specify them in the parameters dictionary below. 
For each parameter, specify a list of values that you would like to use for the 
simulation. If you want to change more than one parameter at the time, consider 
that you need to enter the same number of values for each parameter, e.g. :
    par 1 = [value1, value2]
    par 2 = [value3, value4]
This will run a simulation with par 1 = value1 and par 2 = value 3, and a simulation
with par 1 = value2 and par 2 = value4. 

Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""
from simplicity.runme import run_experiment
from simplicity.dir_manager import get_experiment_output_dir
import simplicity.settings_manager as sm
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import os
import pandas as pd
import json

def fixture_experiment_settings():
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    # varying_params = {
    #     'phenotype_model': ['distance from wt', 'immune waning']
        
    # }
    # # parameters to keep fixed (but different from standard_value) across combinations
    # fixed_params = {
    #     'infected_individuals_at_start': 10,
    #     'final_time': 365
    # }
    varying_params = {
        'phenotype_model': ['linear', 'immune_waning']
    }

    fixed_params = {
        'infected_individuals_at_start': 10,
        # 'R' : 1.3,
        'final_time': 365
    }
    
    # ---------- OR import them from file -------------------------------------
    
    # # leave empty if you only want to import parameters values from file
    # varying_params = {}
    # # import fixed parameters from user geenerated file. You can either create 
    # # it manually or use the provided script: generate_user_set_parameters_file.py
    # filename = 'standard_values.json'
    # fixed_params = sm.read_user_set_parameters_file(filename)
   
    # -------------------------------------------------------------------------
    n_seeds = 3
    
    return (varying_params,fixed_params,n_seeds)

def check_output_file(directory, filename): 
    file_path = os.path.join(directory, filename)
    if os.path.isfile(file_path):
        pass
    else:
        raise ValueError(f"'{filename}' does not exist in the directory.")
        
##### <actual test>
def test_run_experiment_local(runner:str, test_number:int):
    if runner == 'serial':
        runner_module = simplicity.runners.serial
    elif runner == 'multiprocessing':
        runner_module = simplicity.runners.multiprocessing
    else:
        raise ValueError('Runner must be either "serial" or "multiprocessing" ')
    print('')
    print('----------------------------------------------------')
    print(f'testing local {runner} runner')
    print('----------------------------------------------------')
    experiment_name = f'test_local_experiment_{runner}_#{test_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = True,
                       archive_experiment = False)
    except Exception as e:
        raise RuntimeError(f'The code did not pass the running test. Error: {e}')
    print('')
    print(f'TEST LOCAL {runner} RUNNER {test_number} -- SUCCESS.')
    print('############################################')
    return experiment_name

def test_experiment_output(experiment_name):
    # get experiment output directory
    directory = get_experiment_output_dir(experiment_name)
    # loop over each simulation output folder
    for folder in os.listdir(directory):
        simulation_directory = os.path.join(directory,folder)
        # loop over each seeded simulation output folder
        for subfolder in os.listdir(simulation_directory):
            seed_directory = os.path.join(simulation_directory,subfolder)
            check_output_file(seed_directory, 'final_time.csv')
            check_output_file(seed_directory, 'fitness_trajectory.csv')
            check_output_file(seed_directory, 'individuals_data.csv')
            check_output_file(seed_directory, 'lineage_frequency.csv')
            check_output_file(seed_directory, 'phylogenetic_data.csv')
            check_output_file(seed_directory, 'sequencing_data_regression.csv')
            check_output_file(seed_directory, 'sequencing_data.fasta')
            check_output_file(seed_directory, 'simulation_trajectory.csv')
    print('')
    print(f'TEST {experiment_name} OUTPUT -- SUCCESS.')
    print('----------------------------------------------------')
    
def safe_parse_list_of_dicts(val):
    """Parse stringified list of dicts and round time_infection."""
    if pd.isna(val) or str(val).strip() == "":
        return []
    try:
        items = json.loads(val.replace("'", '"'))  
        for item in items:
            if "time_infection" in item:
                item["time_infection"] = round(item["time_infection"], 10)
        return items
    except Exception:
        return val  # fallback to raw string if it breaks

def compare_csv_files(file1, file2, float_tol=1e-6, label=None):
    try:
        if os.path.getsize(file1) == 0 and os.path.getsize(file2) == 0:
            return True
        elif os.path.getsize(file1) == 0 or os.path.getsize(file2) == 0:
            tag = label or f'{file1} vs {file2}'
            print(f'[EMPTY FILE MISMATCH] {tag}')
            return False

        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
    except Exception as e:
        tag = label or f'{file1} vs {file2}'
        print(f'[ERROR] Failed to read files: {tag}: {e}')
        return False

    if df1.shape != df2.shape:
        tag = label or f'{file1} vs {file2}'
        print(f'[SHAPE MISMATCH] {tag}: {df1.shape} != {df2.shape}')
        return False

    # Float rounding for numeric columns
    for col in df1.columns:
        if pd.api.types.is_float_dtype(df1[col]):
            df1[col] = df1[col].round(10)
            df2[col] = df2[col].round(10)

    # Special handling for 'new_infections' column (rounding time)
    if "new_infections" in df1.columns:
        df1["new_infections"] = df1["new_infections"].apply(safe_parse_list_of_dicts)
        df2["new_infections"] = df2["new_infections"].apply(safe_parse_list_of_dicts)

    try:
        pd.testing.assert_frame_equal(
            df1, df2,
            rtol=float_tol, atol=float_tol,
            check_dtype=False
        )
        return True
    except AssertionError as e:
        tag = label or f'{file1} vs {file2}'
        print(f'[DATA MISMATCH] {tag}')
        print(str(e)[:500])  # Trim output
        return False

    
def compare_experiment_outputs(experiment_name_1, experiment_name_2):
    print('')
    print('----------------------------------------------------')
    print(f'Comparing {experiment_name_1} vs {experiment_name_2}')
    print('----------------------------------------------------')

    dir1 = get_experiment_output_dir(experiment_name_1)
    dir2 = get_experiment_output_dir(experiment_name_2)

    files_to_compare = [
        'final_time.csv',
        'fitness_trajectory.csv',
        'individuals_data.csv',
        'lineage_frequency.csv',
        'phylogenetic_data.csv',
        'sequencing_data_regression.csv',
        'sequencing_data.fasta',  
        'simulation_trajectory.csv'
    ]

    differences_found = False

    for folder in os.listdir(dir1):
        sim_dir1 = os.path.join(dir1, folder)
        sim_dir2 = os.path.join(dir2, folder)

        for subfolder in os.listdir(sim_dir1):
            seed_dir1 = os.path.join(sim_dir1, subfolder)
            seed_dir2 = os.path.join(sim_dir2, subfolder)

            for filename in files_to_compare:
                file1 = os.path.join(seed_dir1, filename)
                file2 = os.path.join(seed_dir2, filename)

                # Build relative path from root of experiment dir
                relative_path = os.path.join(folder, subfolder, filename)
                
                if not os.path.exists(file2):
                    print(f'[MISSING FILE] {relative_path}')
                    differences_found = True
                    continue

                if filename.endswith('.fasta'):
                    with open(file1) as f1, open(file2) as f2:
                        if f1.read() != f2.read():
                            print(f'[FASTA MISMATCH] {relative_path}')
                            differences_found = True
                else:
                    if not compare_csv_files(file1, file2, label=relative_path):
                        print(f'[DATA MISMATCH] {relative_path}')
                        differences_found = True

    if not differences_found:
        print('No differences found. Outputs match perfectly!')
    else:
        print('Differences detected in experiment outputs.')
    print('----------------------------------------------------')
    
def main(runner:str, test_number:int, compare_to: int = None):
    test_experiment_output(test_run_experiment_local(runner, test_number))
    if compare_to:
        current_experiment_name = f'test_local_experiment_{runner}_#{test_number}'
        benchmark_experiment_name = f'test_local_experiment_{runner}_#{compare_to}'
        compare_experiment_outputs(current_experiment_name, benchmark_experiment_name)

##### </actual test>

if __name__ == "__main__":
    import time
    start = time.time()
    main('serial',30)#,compare_to=4)
    elapsed = time.time() - start
    mins, secs = divmod(elapsed, 60)
    print(f"Test completed in {int(mins)} min {secs:.2f} sec")
    