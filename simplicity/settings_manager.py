#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 19:38:43 2024

@author: pietro
"""

import os
import json
import simplicity.dir_manager as dm
import pandas as pd
import itertools
import copy

_data_dir = dm.get_data_dir()

def get_standard_parameters_values_file_path():
    standard_parameters_values_file_path = os.path.join(dm.get_reference_parameters_dir(), "standard_values.json")
    return standard_parameters_values_file_path

def get_parameter_specs_file_path():
    parameter_specs_file_path = os.path.join(dm.get_reference_parameters_dir(), "parameter_specs.json")
    return parameter_specs_file_path

def write_standard_parameters_values():
    filename= get_standard_parameters_values_file_path()
    standard_values = {
        "population_size": 1000,
        "infected_individuals_at_start": 10,
        "final_time": 365,
        "tau_3": 7.5,
        "R": 1.1,
        "diagnosis_rate": 0.1, # in percentage, will be converted to kd in model 
        "IH_virus_emergence_rate": 0,      # k_v in theoretical model equations
        "nucleotide_substitution_rate":  4.0112730277780845e-05,  # e in theoretical model equations
        "phenotype_model": 'immune_waning',  # or 'linear'
        "sequencing_rate": 0.05,
        "max_runtime": 86000, 
        "seed": None
    }
    with open(filename, "w") as file:
        json.dump(standard_values, file, indent=4)
    print(f"Standard values written to {filename}")

def write_parameter_specs():
    filename= get_parameter_specs_file_path()
    parameter_specs = {
        "population_size":               {"type": "int", "min": 0, "max": 10000},
        "tau_3":                         {"type": "float", "min": 0, "max": 1000},
        "infected_individuals_at_start": {"type": "int", "min": 0},
        "R":                             {"type": "float", "min": 0, "max": 20},
        "diagnosis_rate":                {"type": "float", "min": 0, "max": 1},
        "IH_virus_emergence_rate":       {"type": "float", "min": 0},
        "nucleotide_substitution_rate":   {"type": "float", "min": 0, "max": 1},
        "final_time":                    {"type": "int", "min": 0},    
        "max_runtime":                   {"type": "int", "min": 0},
        "phenotype_model":               {"type": "str"},
        "sequencing_rate":               {"type": "float", "min": 0, "max": 1}
        }

    with open(filename, "w") as file:
        json.dump(parameter_specs, file, indent=4)
    print(f"Parameter specifications written to {filename}")

def read_standard_parameters_values():
    filename= get_standard_parameters_values_file_path()
    try:
        with open(filename, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default standard values.")
        write_standard_parameters_values()
        return read_standard_parameters_values()

def read_parameter_specs():
    filename= get_parameter_specs_file_path()
    try:
        with open(filename, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default parameter specifications.")
        write_parameter_specs()
        return read_parameter_specs()
    
def write_user_set_parameters_file(user_set_parameters, filename):
    file_path = os.path.join(dm.get_reference_parameters_dir(),filename)
    with open(file_path, "w") as file:
        json.dump(user_set_parameters, file, indent=4)
    print(f"user_set_parameters saved to {file_path}")

def read_user_set_parameters_file(filename):
    file_path = os.path.join(dm.get_reference_parameters_dir(),filename)
    try:
        with open(file_path, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default standard values.")
        write_standard_parameters_values()
        return read_standard_parameters_values()
    
def get_experiment_settings_file_path(experiment_name):
    return os.path.join(_data_dir,
                        f'{experiment_name}',
                        '01_Experiments_settings', 
                        f'{experiment_name}_settings.json')

def get_n_seeds_file_path(experiment_name):
    return os.path.join(_data_dir,
                        f'{experiment_name}',
                        '01_Experiments_settings', 
                        f'{experiment_name}_n_seeds.json')

def check_parameters_names(parameters_dic):
    STANDARD_VALUES = read_standard_parameters_values()
    for key in parameters_dic.keys():
        if key not in STANDARD_VALUES.keys():
            raise ValueError(f'Parameter {key} is not a valid parameter')

def read_experiment_settings(experiment_name):
    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    with open(experiment_settings_file_path, 'r') as json_file:
        experiment_settings = json.load(json_file)
    return experiment_settings

def read_n_seeds_file(experiment_name):
    n_seeds_file_path = get_n_seeds_file_path(experiment_name)
    with open(n_seeds_file_path, 'r') as json_file:
        n_seeds_dic = json.load(json_file)
    return n_seeds_dic


def generate_experiment_settings(varying_params: dict, fixed_params: dict = None):
    """
    Generates a list of parameter combinations from varying and fixed parameters.
    
    Args:
        varying (dict): Parameters for which all combinations should be generated.
        fixed (dict): Parameters that should have the same value across all combinations.
    
    Returns:
        List[dict]: A list of dictionaries with combined parameter sets.
    """
    fixed_params = fixed_params or {}

    keys, values = zip(*varying_params.items()) if varying_params else ([], [])
    combinations = list(itertools.product(*values)) if values else [()]

    experiment_settings = []
    for combo in combinations:
        setting = dict(zip(keys, combo))
        setting.update(copy.deepcopy(fixed_params))  # Avoid mutation
        experiment_settings.append(setting)
    
    return experiment_settings

def write_experiment_settings(experiment_name: str, experiment_settings: list, n_seeds: int):
    """
    Writes experiment settings (a list of parameter dictionaries) to a JSON file.

    Args:
        experiment_name (str): Name of the experiment (used for output folder).
        experiment_settings (list): List of parameter dictionaries.
        n_seeds (int): Number of random seeds to be stored separately.
    """
    # check parameter names validity
    for param_set in experiment_settings:
        check_parameters_names(param_set)

    # Write settings to JSON
    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    n_seeds_file_path = get_n_seeds_file_path(experiment_name)

    with open(experiment_settings_file_path, 'w') as settings_file:
        json.dump(experiment_settings, settings_file, indent=4)

    with open(n_seeds_file_path, 'w') as n_seeds_file:
        json.dump({'n_seeds': n_seeds}, n_seeds_file, indent=4)

    print(f"Experiment settings file written to {experiment_settings_file_path}")

def write_simulation_parameters(file_path,
                                population_size, 
                                tau_3,
                                infected_individuals_at_start, 
                                R,
                                diagnosis_rate,
                                IH_virus_emergence_rate,
                                nucleotide_substitution_rate,
                                final_time, 
                                max_runtime, 
                                phenotype_model,
                                sequencing_rate,
                                seed,
                                F
                                ):
    settings = {
        "population_size": population_size,
        "tau_3": tau_3,
        "infected_individuals_at_start": infected_individuals_at_start,
        "R": R,
        "diagnosis_rate": diagnosis_rate, 
        "IH_virus_emergence_rate" : IH_virus_emergence_rate,
        "nucleotide_substitution_rate": nucleotide_substitution_rate,
        "t_0": 0,
        "final_time": final_time,
        "max_runtime": max_runtime,
        "phenotype_model": phenotype_model,
        "sequencing_rate": sequencing_rate,
        "seed": seed
    }
    
    # Serialize the dictionary to a JSON-formatted string and write it to a file
    with open(file_path, "w") as json_file:
        json.dump(settings, json_file, indent=4)
        
def generate_filename_from_params(params: dict):
    
    abbreviations = {
    "population_size": "N",
    "tau_3": "tau3",
    "infected_individuals_at_start": "init",
    "R": "R",
    "diagnosis_rate": "kd",
    "IH_virus_emergence_rate": "kv",
    "nucleotide_substitution_rate": "NSR",
    "final_time": "T",
    "phenotype_model": "pheno"
    # excluded: max_runtime, sequencing_rate, seed, F
}
    exclude  = {"max_runtime", "sequencing_rate", "seed"}
    parts = []
    for key, value in params.items():
        if key in exclude:
            continue
        abbrev = abbreviations.get(key, key)
        if isinstance(value, float):
            value_str = f"{value:.2g}".replace('.', 'p')  # e.g., 0.01 → 1p0
        elif isinstance(value, int):
            value_str = str(value)
        elif isinstance(value, str):
            value_str = value.replace(' ', '')
        else:
            value_str = str(value)
        parts.append(f"{abbrev}_{value_str}")
    
    file_name = "_".join(parts) + ".json"
    return file_name


def read_settings_and_write_simulation_parameters(experiment_name):
    """
    Reads an experiment settings file in JSON format and generates individual simulation parameter files
    based on the combinations of parameters in the settings. The function will create a separate JSON file
    for each parameter combination within a directory named after the experiment.

    Parameters:
    -----------
    experiment_name : str
        The name of the experiment. This is used to locate the settings file and to create the corresponding
        simulation parameters directory.

    """
    
    # Define the path to the experiment settings file
    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    STANDARD_VALUES = read_standard_parameters_values()
    # Read the experiment settings file
    with open(experiment_settings_file_path, 'r') as settings_file:
        all_experiment_settings = json.load(settings_file)
  
    # Loop over each set of experiment settings and create a simulation parameters file
    for i, experiment_settings in enumerate(all_experiment_settings):
        # Only include parameters that differ from standard
        modified_params = {
            key: value for key, value in experiment_settings.items()
            if key in STANDARD_VALUES and value != STANDARD_VALUES[key]
        }
        # Fallback: if all params match standard, name it by index
        if not modified_params:
            file_name = 'standard_values.json'
        else:
            file_name = generate_filename_from_params(modified_params)
            
        simulation_parameters_file_path = os.path.join(
                       dm.get_simulation_parameters_dir(experiment_name), 
                                                 file_name)
        
        # Merge the standard values with the current experiment settings
        settings = {**STANDARD_VALUES, **experiment_settings}
        
        # Write the simulation parameters to a JSON file
        write_simulation_parameters(simulation_parameters_file_path, 
                                    settings["population_size"],
                                    settings["tau_3"],
                                    settings["infected_individuals_at_start"],
                                    settings["R"],
                                    settings["diagnosis_rate"],
                                    settings["IH_virus_emergence_rate"],
                                    settings["nucleotide_substitution_rate"],
                                    settings["final_time"],
                                    settings["max_runtime"],
                                    settings["phenotype_model"],
                                    settings["sequencing_rate"],
                                    settings["seed"]
                                    )

    print(f"Simulation parameters written to directory: {simulation_parameters_file_path}")

def write_seeded_simulation_parameters(experiment_name: str):
    """
    Generates multiple JSON files with different seeds for each simulation parameter file 
    within a specified experiment. The function reads the original simulation parameter 
    files, adds a 'seed' field, and writes the modified files to subdirectories named 
    after the original files.

    Parameters:
    -----------
    experiment_name : str
        The name of the experiment. 
    
    n_seeds : int
        The number of seeded JSON files to generate for each simulation parameter file. 

    Directory Structure:
    --------------------
    
    Data/
    └── experiment_name/
        ├── 02_Simulation_parameters/
        │   ├── param_file_1.json
        │   ├── param_file_2.json
        │   └── ...
        └── 03_Seeded_simulation_parameters/
            ├── param_file_1/
            │   ├── seed_0.json
            │   ├── seed_1.json
            │   └── ...
            ├── param_file_2/
            │   ├── seed_0.json
            │   ├── seed_1.json
            │   └── ...
            └── ...
    """
    simulation_parameters_dir = dm.get_simulation_parameters_dir(experiment_name)
    seeded_simulation_parameters_dir = dm.get_seeded_simulation_parameters_dir(experiment_name)

    # Iterate over all files in the simulation parameters directory
    for filename in os.listdir(simulation_parameters_dir):
        filepath = os.path.join(simulation_parameters_dir, filename)
        
        # Read the original JSON file
        with open(filepath, 'r') as file:
            simulation_parameters = json.load(file)
        
        # Create a subdirectory for the seeded files
        subdir_name = filename.replace(".json", "")
        subdir_path = os.path.join(seeded_simulation_parameters_dir, subdir_name)
        os.makedirs(subdir_path, exist_ok=True)
        
        # Generate multiple files with different seeds
        n_seeds = read_n_seeds_file(experiment_name)['n_seeds']
        for i in range(n_seeds):
            # seed = random.randint(0, 1000000)
            simulation_parameters['seed'] = i

            seeded_file_path = os.path.join(subdir_path, f"seed_{i:04}.json")
            
            # Write the new JSON file with the added seed
            with open(seeded_file_path, 'w') as seeded_file:
                json.dump(simulation_parameters, seeded_file, indent=4)


def get_seeded_simulation_parameters_paths(experiment_name):
    """
    Retrieves all seeded simulation parameter file paths for a given experiment.

    This function searches through the directory structure of the provided
    experiment name, looking for all JSON seeded simulation parameters files. 
    It returns a list of full paths to these files.

    Args:
        experiment_name (str): The name of the experiment for which seeded simulation 
                               parameter paths are to be retrieved.

    Returns:
        list of str: A list of file paths, each pointing to a JSON file containing
                     seeded simulation parameters.

    Example:
        If `experiment_name` is 'Experiment_1', and the directory structure contains
        multiple JSON files under:
            '/data_dir/Experiment_1/03_Seeded_simulation_parameters',
        the function will return a list like:
            [
                '/data_dir/Experiment_1/03_Seeded_simulation_parameters/subdir/file1.json',
                '/data_dir/Experiment_1/03_Seeded_simulation_parameters/subdir/file2.json'
            ]
    """
    # Define the base path for the simulation parameters
    base_dir = os.path.join(_data_dir, experiment_name, "03_Seeded_simulation_parameters")
    
    # List to store all seed file paths
    seeded_simulation_parameters_paths = []
    
    # Walk through the directory structure to find all .json files
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".json"):
                # Construct the full path to the file and add it to the list
                file_path = os.path.join(root, file)
                seeded_simulation_parameters_paths.append(file_path)
    
    return seeded_simulation_parameters_paths

def read_seeded_simulation_parameters(experiment_name, seeded_simulation_parameters_path):
    with open(seeded_simulation_parameters_path, 'r') as seeded_file:
        seeded_simulation_parameters = json.load(seeded_file)
    return seeded_simulation_parameters

def get_simulation_parameters_filepath_of_simulation_output_dir(simulation_output_dir):
    # get filepath of simulation parameters file from which the simulation_output_dir was generated
    from pathlib import Path
    simulation_output_dir_path = Path(simulation_output_dir)
    parts = simulation_output_dir_path.parts
    experiment_name = parts[-3]
    simulation_output_folder_name = parts[-1]
    # get simulation parameters dir for that experiment
    simulation_parameters_dir = dm.get_simulation_parameters_dir(experiment_name)
    simulation_parameters_file_path = os.path.join(simulation_parameters_dir,
                                              simulation_output_folder_name +'.json')

    return simulation_parameters_file_path

def get_parameter_value_from_simulation_output_dir(simulation_output_dir, parameter):
    # read and return desired parameter value for specific simulation output directory
    simulation_parameters_file_path = get_simulation_parameters_filepath_of_simulation_output_dir(
                       simulation_output_dir)
    with open(simulation_parameters_file_path, 'r') as file:
        parameters_dict = json.load(file)
    
    return parameters_dict[parameter]

def read_OSR_NSR_regressor_parameters():
    file_path = os.path.join(dm.get_reference_parameters_dir(),
              'OSR_NSR_regressor_parameters_for_standard_parameter_values_exp.csv')
    df = pd.read_csv(file_path,index_col=0)
    best_fit_df = pd.to_numeric(df['Best Fit'], errors='coerce')
    return best_fit_df






















