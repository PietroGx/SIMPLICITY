#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 19:36:24 2024

@author: pietro
"""

import os

# Default data directory 
DEFAULT_DATA_DIR = os.path.join(os.getcwd(), "Data")
os.makedirs(DEFAULT_DATA_DIR,exist_ok=True)

# Global variable to store the data directory path
_data_dir = DEFAULT_DATA_DIR

# set env variables 
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_MULTIPROCESS"] = str(15)
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"] = str(150)

def set_data_dir(path):
    """Set the data directory path."""
    global _data_dir
    _data_dir = path

def get_data_dir():
    """Get the current data directory path."""
    return _data_dir

def create_directories(experiment_name):
    """Create necessary subdirectories within the data directory."""
    
    experiment_dir = os.path.join(_data_dir,experiment_name)
    subdirectories = ["01_Experiments_settings",
                      "02_Simulation_parameters", 
                      "03_Seeded_simulation_parameters",
                      "04_Output", 
                      ]
    
    for subdir in subdirectories:
        path = os.path.join(experiment_dir, subdir)
        os.makedirs(path, exist_ok=True)

def get_experiment_dir(experiment_name):
    """Get the experiment_name directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return experiment_dir
    else:
        raise ValueError('No experiment with that name!')

def get_simulation_parameters_dir(experiment_name):
    """Get the experiment_name simulation parameters directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, '02_Simulation_parameters')
    else:
         raise ValueError('No experiment with that name!')
         
def get_seeded_simulation_parameters_dir(experiment_name):
    """Get the experiment_name seeded simulation parameters directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, '03_Seeded_simulation_parameters')
    else:
         raise ValueError('No experiment with that name!')
         
def get_experiment_output_dir(experiment_name):
    """Get the experiment_name output directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, '04_Output')
    else:
         raise ValueError('No experiment with that name!')

def get_slurm_logs_dir(experiment_name):
    slurm_logs_dir = os.path.join(get_experiment_dir(experiment_name), 'slurm','slurm_logs')
    os.makedirs(slurm_logs_dir, exist_ok=True)
    return slurm_logs_dir

def get_slurm_id_map_dir(experiment_name):
    slurm_id_map_dir  = os.path.join(get_experiment_dir(experiment_name), 'slurm','job_id_mapping')
    os.makedirs(slurm_id_map_dir , exist_ok=True)
    return slurm_id_map_dir

def get_simulation_output_dirs(experiment_name):
    simulation_output_dirs = []
    experiment_output_dir = get_experiment_output_dir(experiment_name)
    for folder in os.listdir(experiment_output_dir):
        simulation_output_dir = os.path.join(experiment_output_dir,folder)
        if os.path.isdir(simulation_output_dir):
            simulation_output_dirs.append(simulation_output_dir)
    return simulation_output_dirs

def get_seeded_simulation_output_dirs(simulation_output_dir):
    seeded_simulation_output_dirs = []
    for subfolder in os.listdir(simulation_output_dir):
        seeded_simulation_output_dir = os.path.join(simulation_output_dir,subfolder)
        if os.path.isdir(seeded_simulation_output_dir):
            seeded_simulation_output_dirs.append(seeded_simulation_output_dir)
    return seeded_simulation_output_dirs









