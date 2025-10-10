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
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_MULTIPROCESS"] = str(10)
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"] = str(500)
os.environ["QT_QPA_PLATFORM"] = "offscreen"

def set_data_dir(path):
    """Set the data directory path."""
    global _data_dir
    _data_dir = path

def get_data_dir():
    """Get the current data directory path."""
    return _data_dir

def get_reference_parameters_dir():
    reference_parameters_dir = os.path.join(get_data_dir(), '00_Reference_parameters')
    os.makedirs(reference_parameters_dir, exist_ok=True)
    return reference_parameters_dir

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
         
def get_experiment_plots_dir(experiment_name):
    experiment_dir = os.path.join(_data_dir,experiment_name)
    experiment_plots_dir = os.path.join(experiment_dir, '05_Plots')
    os.makedirs(experiment_plots_dir, exist_ok=True) 
    return experiment_plots_dir

def get_experiment_simulations_plots_dir(experiment_name):
    experiment_plots_dir = get_experiment_plots_dir(experiment_name)
    experiment_simulations_plots_dir = os.path.join(experiment_plots_dir, 'Simulations')
    os.makedirs(experiment_simulations_plots_dir, exist_ok=True) 
    return experiment_simulations_plots_dir

def get_experiment_tree_dir(experiment_name):
    experiment_dir = get_experiment_dir(experiment_name)
    experiment_tree_dir = os.path.join(experiment_dir, '06_Trees')
    os.makedirs(experiment_tree_dir, exist_ok=True) 
    return experiment_tree_dir

def get_experiment_fit_result_dir(experiment_name):
    fit_result_dir = os.path.join(get_experiment_dir(experiment_name), 'Fit_results')
    os.makedirs(fit_result_dir, exist_ok=True)
    return fit_result_dir

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

# SSOD = seeded_simulation_output_dir
def get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir):
    ''' Get the simulation_output folder_name of Data/experiment/04_Output/simulation_output/seed_nr'''
    simulation_output_foldername = os.path.basename(os.path.dirname(seeded_simulation_output_dir))
    return simulation_output_foldername

def get_experiment_foldername_from_SSOD(seeded_simulation_output_dir):
    # Split the path into parts
    path_parts = os.path.normpath(seeded_simulation_output_dir).split(os.sep)
    # Extract the experiment folder name
    experiment_foldername = path_parts[-4]
    return experiment_foldername

def get_experiment_tree_simulation_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_dir = get_experiment_tree_dir(experiment_name)
    simulation_output_foldername = get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    experiment_tree_simulation_dir = os.path.join(experiment_tree_dir,simulation_output_foldername)
    os.makedirs(experiment_tree_simulation_dir,exist_ok=True)
    return experiment_tree_simulation_dir

def get_experiment_tree_simulation_files_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_simulation_dir = get_experiment_tree_simulation_dir(experiment_name,
                                           seeded_simulation_output_dir)
   
    experiment_tree_simulation_files_dir = os.path.join(experiment_tree_simulation_dir,
                                                       'files')
    os.makedirs(experiment_tree_simulation_files_dir,exist_ok=True)
    return experiment_tree_simulation_files_dir

def get_experiment_tree_simulation_plots_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_simulation_dir = get_experiment_tree_simulation_dir(experiment_name,
                                           seeded_simulation_output_dir)
    experiment_tree_simulation_plots_dir = os.path.join(experiment_tree_simulation_dir,
                                                       'plots')
    os.makedirs(experiment_tree_simulation_plots_dir,exist_ok=True)
    return experiment_tree_simulation_plots_dir

def get_ssod(sim_out_dir, seed_number):
    """
    Returns the seeded simulation output directory (SSOD) for the given seed number.

    Args:
        sim_out_dir (str): Path to a simulation output directory (e.g., one from dm.get_simulation_output_dirs()).
        seed_number (int): The seed number (e.g., 7).

    Returns:
        str: Full path to the matching seed_XXXX directory.

    Raises:
        ValueError: If the seed directory cannot be found.
    """
    all_seed_dirs = get_seeded_simulation_output_dirs(sim_out_dir)
    target = f"seed_{seed_number:04d}"
    
    for path in all_seed_dirs:
        if os.path.basename(path) == target:
            return path

    raise ValueError(f"Seed folder '{target}' not found in {sim_out_dir}")  

def get_figure_source_data_dir(experiment_name: str, figure_name: str):
    """
    Create (if missing) and return the path to the figure-specific source data folder.

    The base path is the experiment folder, e.g.:
        Data/<experiment_name>/source_data/<figure_name>/

    Parameters
    ----------
    experiment_name : str
        Name of the experiment folder.
    figure_name : str
        Name of the figure

    Returns
    -------
    str
        Absolute path to the figure source data directory.
    """
    exp_dir = get_experiment_dir(experiment_name)
    source_data_dir = os.path.join(exp_dir, "source_data", figure_name)
    os.makedirs(source_data_dir, exist_ok=True)
    print(f"[dir_manager] Source data directory ready: {source_data_dir}")
    return source_data_dir





