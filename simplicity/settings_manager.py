#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 19:38:43 2024

@author: pietro
"""

import os
import json
import simplicity.config as config

_data_dir = config.get_data_dir()

# Define the standard values
STANDARD_VALUES = {
    "population_size": 1000,
    "infected_individuals_at_start": 100,
    "R": 1.5,
    "k_d": 0.0055,
    "k_v": 0.0085,
    "e": 0.0017,
    "final_time": 365*3 ,
    "max_runtime": 300, 
    "phenotype_model": 'immune waning',  # or 'distance from wt'
    "sequencing_rate": 0.05,
    "seed": None,
    "F": 1.25
}

def write_settings(parameters: dict, 
                   n_seeds: int,
                   target_value=None) -> None:
    """
    Writes settings to a JSON file.

    Args:
        parameters (dict): A dictionary where keys are parameter names and values are lists 
                           of possible values for those parameters. All lists must be of 
                           the same length and contain at least one value.
        n_seeds (int): Number of random seeds for the experiment.
        target_value (float): The target value if running a parameter search.
    """
    
    # Check if target_value is not None and parameters contain more than one entry
    if target_value is not None:
        if len(parameters) != 1 or len(parameters[list(parameters.keys())[0]]) != 1:
            raise ValueError("If target_value is specified, parameters dict must contain exactly one entry with one value'.")

    settings = {
        "parameters": parameters,
        "n_seeds": n_seeds,
        "target_value": target_value,
        "learning_rate": 0.0003,
        "error": None,
        "experiment_number": 0
    }

    with open(os.path.join(_data_dir,'settings.json'), 'w') as json_file:
        json.dump(settings, json_file, indent=4)

    print("Settings written to Data/settings.json")
    
def read_experiment_settings(experiment_name):
    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    with open(experiment_settings_file_path, 'r') as json_file:
        experiment_settings = json.load(json_file)
    return experiment_settings
    
def write_experiment_settings(experiment_name):
    """
    Generates and writes experiment settings to a JSON file.

    This function generates combinations of parameters. The generated
    configurations are written to a JSON file, organized under a directory specific
    to the given experiment name.

    Args:
        experiment_name (str): The name of the experiment, which is used to define the 
                               directory and filename for the output JSON file.

    Raises:
        ValueError: If any parameter list is empty or if the parameter lists do not all 
                    have the same length.

    Example:
        If `parameters` is:
            {'param1': [1, 2, 3], 'param2': [4, 5, 6]}
        and `experiment_name` is 'Experiment_1', the function will create a JSON file
        with the following content:
            [
                {"param1": 1, "param2": 4},
                {"param1": 2, "param2": 5},
                {"param1": 3, "param2": 6}
            ]
        The file will be saved under:
            '/data_dir/Experiment_1/01_Experiments_settings/Experiment_1_settings.json'
    
    """
    # read settings and get the parameters dict
    with open(os.path.join(_data_dir,'settings.json'), 'r') as json_file:
        settings = json.load(json_file)
    parameters = settings['parameters']
    
    # Validate that all lists have the same length
    list_lengths = [len(v) for v in parameters.values()]
    if any(length == 0 for length in list_lengths):
        raise ValueError("All parameter lists must contain at least one value.")
    if len(set(list_lengths)) > 1:
        raise ValueError("All parameter lists must have the same length.")
    
    # Create a list to hold all paired parameter sets
    all_experiment_settings = []

    # Zip the values together to form pairs
    keys = list(parameters.keys())
    for combination in zip(*parameters.values()):
        experiment_settings = dict(zip(keys, combination))
        all_experiment_settings.append(experiment_settings)
    
    # Write the generated configurations to a JSON file
    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    with open(experiment_settings_file_path, 'w') as config_file:
        json.dump(all_experiment_settings, config_file, indent=4)
    
    print(f"Experiment settings file written to {experiment_settings_file_path}")

def get_experiment_settings_file_path(experiment_name):
    # Define the path to the output file
    return os.path.join(_data_dir,
                        f'{experiment_name}',
                        '01_Experiments_settings', 
                        f'{experiment_name}_settings.json')
    
def write_simulation_parameters(file_path,
                                population_size, 
                                infected_individuals_at_start, 
                                R,
                                k_d,
                                k_v,
                                e,
                                final_time, 
                                max_runtime, 
                                phenotype_model,
                                sequencing_rate,
                                seed,
                                F
                                ):
    settings = {
        "population_size": population_size,
        "infected_individuals_at_start": infected_individuals_at_start,
        "R": R,
        "diagnosis rate": k_d,
        "IH virus emergence rate" : k_v,
        "evolutionary rate": e,
        "t_0": 0,
        "final_time": final_time,
        "max_runtime": max_runtime,
        "phenotype_model": phenotype_model,
        "sequencing_rate": sequencing_rate,
        "seed": seed,
        "F": F
    }
    
    # Serialize the dictionary to a JSON-formatted string and write it to a file
    with open(file_path, "w") as json_file:
        json.dump(settings, json_file, indent=4)

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
    experiment_settings_file_path = os.path.join(_data_dir,
                                                 f'{experiment_name}',
                                                 '01_Experiments_settings', 
                                                 f'{experiment_name}_settings.json')
    
    # Read the experiment settings file
    with open(experiment_settings_file_path, 'r') as config_file:
        all_experiment_settings = json.load(config_file)
    # handle case of no changes from standard values
    if all_experiment_settings == []:
        # Create a filename
        file_name = 'standard_values.json'
        
        # Ensure the directory for simulation parameters exists
        simulation_parameters_file_path = os.path.join(_data_dir,
                                                 f'{experiment_name}',
                                                 '02_Simulation_parameters', 
                                                 file_name)
        
        # Write the simulation parameters to a JSON file
        write_simulation_parameters(simulation_parameters_file_path, 
                                    STANDARD_VALUES["population_size"],
                                    STANDARD_VALUES["infected_individuals_at_start"],
                                    STANDARD_VALUES["R"],
                                    STANDARD_VALUES["k_d"],
                                    STANDARD_VALUES["k_v"],
                                    STANDARD_VALUES["e"],
                                    STANDARD_VALUES["final_time"],
                                    STANDARD_VALUES["max_runtime"],
                                    STANDARD_VALUES["phenotype_model"],
                                    STANDARD_VALUES["sequencing_rate"],
                                    STANDARD_VALUES["seed"],
                                    STANDARD_VALUES["F"]
                                    )
    else: 
        # Loop over each set of experiment settings and create a simulation parameters file
        for i, experiment_settings in enumerate(all_experiment_settings):
            modified_params = []
            modified_values = []
            
            for key, value in experiment_settings.items():
                modified_params.append(key)
                modified_values.append(value)
            
            # Create a filename based on the parameters and their values
            file_name_parts = [f"{param}_{value}" for param, value in zip(modified_params, modified_values)]
            file_name = "_".join(file_name_parts) + '.json'
            
            # Ensure the directory for simulation parameters exists
            simulation_parameters_file_path = os.path.join(_data_dir,
                                                     f'{experiment_name}',
                                                     '02_Simulation_parameters', 
                                                     file_name)
            
            # Merge the standard values with the current experiment settings
            settings = {**STANDARD_VALUES, **experiment_settings}
            
            # Write the simulation parameters to a JSON file
            write_simulation_parameters(simulation_parameters_file_path, 
                                        settings["population_size"],
                                        settings["infected_individuals_at_start"],
                                        settings["R"],
                                        settings["k_d"],
                                        settings["k_v"],
                                        settings["e"],
                                        settings["final_time"],
                                        settings["max_runtime"],
                                        settings["phenotype_model"],
                                        settings["sequencing_rate"],
                                        settings["seed"],
                                        settings["F"]
                                        )

    print(f"Simulation parameters written to directory: {simulation_parameters_file_path}")

def write_seeded_simulation_parameters(experiment_name):
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
    # Define paths
    base_dir = os.path.join(_data_dir, f"{experiment_name}")
    simulation_parameters_dir = os.path.join(base_dir, "02_Simulation_parameters")
    seeded_simulation_parameters_dir = os.path.join(base_dir, "03_Seeded_simulation_parameters")

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
        
        # read settings file and get n_seeds
        with open(os.path.join(_data_dir,'settings.json'), 'r') as json_file:
            settings = json.load(json_file)
        n_seeds = settings['n_seeds']
        
        # Generate multiple files with different seeds
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
    base_dir = os.path.join(_data_dir, experiment_name, "03_Seeded_simulation_parameters")
    seeded_simulation_parameters_full_path = os.path.join(base_dir, seeded_simulation_parameters_path)
    with open(seeded_simulation_parameters_full_path, 'r') as seeded_file:
        seeded_simulation_parameters = json.load(seeded_file)
    return seeded_simulation_parameters
    

