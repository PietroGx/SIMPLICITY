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
Created on Thu Aug 29 13:56:20 2024

@author: pietro
"""

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import os
import shutil
import tarfile
import csv
from simplicity.evolution.decoder import decode_genome
import simplicity.phenotype.distance  as dis
import pandas as pd
import simplicity.tuning.evolutionary_rate as er
import glob
import ast
import numpy as np
from tqdm import tqdm
import warnings
from pathlib import Path

def setup_output_directory(experiment_name, seeded_simulation_parameters_path):
    """
    Sets up the output directory structure for a simulation.

    This function constructs the output directory path based on the provided
    seeded simulation parameters file path and the experiment name.

    Args:
        seeded_simulation_parameters_path (str): The file path to the seeded simulation
                                                 parameters.
        experiment_name (str): The name of the experiment for which the output
                               directory is being set up.

    Returns:
        str: The path to the main output directory created for this experiment.

    Example:
        If the `seeded_simulation_parameters_path` is 
        '/path/to/seed/files/seed_file.txt' and the `experiment_name` is 
        'Experiment_1', the function will create and return the path:
        '/data_dir/Experiment_1/04_Output/files/seed_file'
    """
    # Split the path into components
    path_components = os.path.split(seeded_simulation_parameters_path)
    
    # Extract the necessary components
    root_folder_of_seed_file = os.path.split(path_components[0])[-1]
    
    # Get the seed file name without the extension
    seed_file_name = os.path.splitext(path_components[-1])[0]  
    # Construct the output directory path
    output_dir = os.path.join(dm.get_experiment_output_dir(experiment_name), 
                              root_folder_of_seed_file, 
                              seed_file_name)
    try:
        os.makedirs(output_dir)
    except:
        raise RuntimeError('You already run an experiment with the same name!')
    # Return the output directory path
    return output_dir

def archive_experiment(experiment_name):
    """
    Archives the specified experiment folder as a .tar.gz file and deletes the original folder.

    This script:
    1. Checks if the 'Archive' directory exists within the data directory, and creates it if not.
    2. Compresses the experiment folder into a .tar.gz archive that is compatible with Windows.
    3. Deletes the original experiment folder after archiving.

    Args:
        experiment_name (str): The name of the experiment to be archived.

    Raises:
        FileNotFoundError: If the experiment folder does not exist.
    """
    # Get the data directory from the dm
    data_dir = dm.get_data_dir()

    # Define the path to the experiment folder
    experiment_folder_path = os.path.join(data_dir, experiment_name)

    # Check if the experiment folder exists
    if not os.path.exists(experiment_folder_path) or experiment_name=='': 
        raise FileNotFoundError(f"The experiment folder '{experiment_folder_path}' does not exist.")

    # Define the path to the Archive directory
    archive_dir = os.path.join(data_dir, "Archive")

    # Create the Archive directory if it does not exist
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
        print(f"Created archive directory: {archive_dir}")

    # Define the path for the tar.gz archive file
    archive_file_path = os.path.join(archive_dir, f"{experiment_name}.tar.gz")

    # Create the tar.gz archive
    with tarfile.open(archive_file_path, "w:gz") as tar:
        tar.add(experiment_folder_path, arcname=os.path.basename(experiment_folder_path))
        print(f"Archived '{experiment_folder_path}' to '{archive_file_path}'")

    # Delete the original experiment folder
    shutil.rmtree(experiment_folder_path)
    print(f"Deleted the original experiment folder: {experiment_folder_path}")

# -----------------------------------------------------------------------------
#                            read and write simulation output
# -----------------------------------------------------------------------------


def save_sequencing_dataset(simulation_output, output_path, sequence_long_shedders=False):
    """
    Writes the simulated sequencing data to FASTA files and a regression CSV.
    
    Parameters:
    - simulation_output: The object containing simulation results.
    - output_path: Directory to save files.
    - sequence_long_shedders (bool): If True, extracts and saves all lineages 
      from long shedders (at end of infection/sim) to a separate FASTA 
      and includes them in the regression CSV.
    
    Outputs:
    1. sequencing_data.fasta: Random surveillance sequences 
    2. sequencing_data_long.fasta: All lineages from long-shedders 
    3. sequencing_data_regression.csv: Combined metrics for both groups.
    4. sequencing_data.csv: Raw surveillance metadata 
    """
    
    # Ensure output directory exists
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    try:
        # 1. Setup Paths
        fasta_file_path = os.path.join(output_path, 'sequencing_data.fasta')
        fasta_long_path = os.path.join(output_path, 'sequencing_data_long.fasta')
        data_regression_file_path = os.path.join(output_path, 'sequencing_data_regression.csv')
        data_file_path = os.path.join(output_path, 'sequencing_data.csv')

        # List to collect data for the regression CSV (combine surveillance + long shedders)
        sequencing_data_regression_dic = []

        # ==============================================================================
        # Process Random Surveillance Data 
        # ==============================================================================
        sequencing_data = simulation_output.sequencing_data
        
        if sequencing_data:
            # Open standard FASTA
            with open(fasta_file_path, 'w') as fasta_file:
                for individual_sequencing_data in sequencing_data:
                    # Extract details
                    ind_index    = individual_sequencing_data['individual_index']
                    ih_lin_index = individual_sequencing_data['intra-host_lineage_index']
                    seq_time     = individual_sequencing_data['sequencing_time']
                    lin_name     = individual_sequencing_data['lineage_name']
                    sequence     = individual_sequencing_data['sequence']
                    ind_type     = individual_sequencing_data['individual_type']
                    
                    # Write to FASTA
                    sequences_id = f'>{ind_index}_{ih_lin_index}_time_{seq_time:.2f}_lin_{lin_name}'
                    fasta_file.write(f"{sequences_id}\n")
                    
                    decoded_genome = decode_genome(sequence)
                    genome_str = decoded_genome.replace("'", "")
                    fasta_file.write(f"{genome_str}\n")      
                    
                    # Add to Regression Dict
                    sequencing_data_regression_dic.append({
                        'Sequencing_time': round(seq_time / 365.25, 4), # Years
                        'Distance_from_root': dis.hamming(sequence) / len(dis.reference),
                        'individual_type': ind_type
                    })
            
            # Save Sequencing Metadata CSV
            with open(data_file_path, mode='w', newline='') as file:
                fieldnames = sequencing_data[0].keys()  
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(sequencing_data)

        # ==============================================================================
        # Process All Long Shedders (Optional)
        # ==============================================================================
        long_shedder_entries = [] # Buffer
        
        if sequence_long_shedders:
            
            # Iterate through all individuals to find long shedders
            for i, ind_data in simulation_output.individuals.items():
                
                # Filter: Only Long Shedders who were infected and completed their infection course before the end of a simulation
                if ind_data['type'] == 'long_shedder' and ind_data['state'] == 'recovered':
                    
                    # 1. Determine Sequencing Time
                    if ind_data['t_not_infectious'] is not None:
                        seq_time = ind_data['t_not_infectious']
                    else:
                        seq_time = simulation_output.time 
                    
                    # 2. Iterate through ALL lineages
                    if not ind_data['IH_lineages']:
                        continue 
                        
                    for lineage_name in ind_data['IH_lineages']:
                        # Retrieve the genome sequence 
                        sequence_list = simulation_output.get_lineage_genome(lineage_name)
                        
                        # Prepare FASTA entry string
                        header = f">{i}_long_shedder_time_{seq_time:.2f}_lin_{lineage_name}"
                        decoded_genome = decode_genome(sequence_list)
                        genome_str = decoded_genome.replace("'", "")
                        
                        # Store in buffer
                        long_shedder_entries.append(f"{header}\n{genome_str}\n")
                        
                        # Add to Regression Dict 
                        sequencing_data_regression_dic.append({
                            'Sequencing_time': round(seq_time / 365.25, 4),
                            'Distance_from_root': dis.hamming(sequence_list) / len(dis.reference),
                            'individual_type': 'long_shedder'
                        })

            # Only write the Long Shedder FASTA file if we actually found data
            if long_shedder_entries:
                with open(fasta_long_path, 'w') as fasta_long:
                    fasta_long.writelines(long_shedder_entries)

        # ==============================================================================
        # PART C: Save Regression CSV (Combined)
        # ==============================================================================
        if sequencing_data_regression_dic:
            with open(data_regression_file_path, mode='w', newline='') as file:
                fieldnames = sequencing_data_regression_dic[0].keys()  
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(sequencing_data_regression_dic)

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Warning: An error occurred while saving sequencing datasets: {e}")

def read_sequencing_data_regression(seeded_simulation_output_dir):
    sequencing_data_file_path = os.path.join(seeded_simulation_output_dir,
                                        "sequencing_data_regression.csv")
    df = pd.read_csv(sequencing_data_file_path)
    return df
        
def save_simulation_trajectory(simulation_output, seeded_simulation_output_dir):
    df = pd.DataFrame(simulation_output.trajectory, columns= 
                      ['time','infected','diagnosed','recovered',
                       'infectious','detectables','susceptibles','long_shedders'])
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "simulation_trajectory.csv")
    df.to_csv(trajectory_file_path, index=False)

def read_simulation_trajectory(seeded_simulation_output_dir):
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "simulation_trajectory.csv")
    df = pd.read_csv(trajectory_file_path)
    return df
    
def save_lineage_frequency(simulation_output, seeded_simulation_output_dir):
    df = simulation_output.lineage_frequency_to_df()
    lineage_frequency_file_path = os.path.join(seeded_simulation_output_dir,
                                               "lineage_frequency.csv")
    df.to_csv(lineage_frequency_file_path, index=False)

def read_lineage_frequency(seeded_simulation_output_dir):
    lineage_frequency_file_path = os.path.join(seeded_simulation_output_dir,
                                        "lineage_frequency.csv")
    df = pd.read_csv(lineage_frequency_file_path)
    return df

def save_individuals_data(simulation_output, seeded_simulation_output_dir):
    individuals_data = simulation_output.individuals_data_to_df()
    individuals_data_file_path = os.path.join(seeded_simulation_output_dir,
                                               "individuals_data.csv")
    individuals_data.to_csv(individuals_data_file_path)
    
def read_individuals_data(seeded_simulation_output_dir):
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "individuals_data.csv")
    df = pd.read_csv(trajectory_file_path, index_col=0,low_memory=False)
    df['type'] = df['type'].astype(str)
    df['IH_lineages'] = df['IH_lineages'].apply(ast.literal_eval)
    df['IH_lineages_fitness_score'] = df['IH_lineages_fitness_score'].apply(ast.literal_eval) 
    df['new_infections'] = df['new_infections'].apply(ast.literal_eval)
    df['IH_lineages_trajectory'] = df['IH_lineages_trajectory'].apply(ast.literal_eval)
    return df

def save_phylogenetic_data(simulation_output, seeded_simulation_output_dir):
    phylogenetic_data = simulation_output.phylogenetic_data_to_df()
    phylogenetic_data_file_path = os.path.join(seeded_simulation_output_dir,
                                               "phylogenetic_data.csv")
    phylogenetic_data.to_csv(phylogenetic_data_file_path)
    
def read_phylogenetic_data(seeded_simulation_output_dir):
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "phylogenetic_data.csv")
    
    df = pd.read_csv(trajectory_file_path, index_col=0)
    df['Genome'] = df['Genome'].apply(ast.literal_eval)
    return df
    
def save_fitness_trajectory(simulation_output, seeded_simulation_output_dir):
    fitness_trajectory = simulation_output.fitness_trajectory_to_df()
    fitness_trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                               "fitness_trajectory.csv")
    fitness_trajectory.to_csv(fitness_trajectory_file_path)
    
def save_final_time(simulation_output, seeded_simulation_output_dir):
    final_time = round(simulation_output.time,6)
    final_time_file_path = os.path.join(seeded_simulation_output_dir,
                                               "final_time.csv")
    with open(final_time_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([final_time]) 
        
def read_final_time(seeded_simulation_output_dir):
    final_time_file_path = os.path.join(seeded_simulation_output_dir,
                                               "final_time.csv")
    with open(final_time_file_path, 'r') as f:
        final_time = float(f.read().strip())
    return final_time

def save_R_effective_trajectory(simulation_output, seeded_simulation_output_dir):
    R_effective_trajectory = simulation_output.R_effective_trajectory_to_df()
    R_effective_trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                               "R_effective_trajectory.csv")
    R_effective_trajectory.to_csv(R_effective_trajectory_file_path)

def read_R_effective_trajectory(seeded_simulation_output_dir):
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "R_effective_trajectory.csv")
    df = pd.read_csv(trajectory_file_path,index_col=0)
    return df

# -----------------------------------------------------------------------------
#                              Intra host lineages
# -----------------------------------------------------------------------------

def get_all_individuals_data_for_simulation_output_dir(simulation_output_dir):
    seeded_simulation_output_dirs = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
    all_individuals_data = pd.DataFrame()
    for seeded_simulation_output_dir in seeded_simulation_output_dirs:
        individuals_data = os.path.join(seeded_simulation_output_dir,'individuals_data.csv')
        df = pd.read_csv(individuals_data, index_col=0)
        all_individuals_data = pd.concat([all_individuals_data, df], axis=0)
    
    return all_individuals_data

def get_IH_lineages_data_simulation(simulation_output_dir):
    # get individuals dataframe 
    data = get_all_individuals_data_for_simulation_output_dir(simulation_output_dir)
    # keep only needed columns
    data = data[['IH_lineages_number', 'IH_unique_lineages_number']]

    # DataFrames with all possible IH_lineages_number and IH_unique_lineages_numbervalues (1 to 5)
    IH_lineages_number_values = pd.DataFrame({'IH_lineages_number': [1, 2, 3, 4, 5]})
    IH_unique_lineages_number_values = pd.DataFrame({'IH_unique_lineages_number': [1, 2, 3, 4, 5]})

    # Calculate the counts of each value for each type
    ih_virus_count = data.groupby(['IH_lineages_number']).size().reset_index(name='ih_virus_count')
    ih_lineage_count = data.groupby(['IH_unique_lineages_number']).size().reset_index(name='ih_lineage_count')

    # Merge with the IH_lineages_number_values DataFrame, filling NaNs with 0
    ih_virus_count = pd.merge(IH_lineages_number_values, ih_virus_count, on='IH_lineages_number', how='left').fillna({'ih_virus_count': 0})
    ih_lineage_count = pd.merge(IH_unique_lineages_number_values, ih_lineage_count, on='IH_unique_lineages_number', how='left').fillna({'ih_lineage_count': 0})

    # Normalize counts by the total count
    ih_virus_count['ih_virus_count'] = ih_virus_count['ih_virus_count'] / ih_virus_count['ih_virus_count'].sum()
    ih_lineage_count['ih_lineage_count'] = ih_lineage_count['ih_lineage_count'] / ih_lineage_count['ih_lineage_count'].sum()
    
    df = pd.concat([ih_virus_count,ih_lineage_count],axis=1)
    # add hue column to df - hue is IH_lineage_emergence_rate
    df['IH_virus_emergence_rate'] = sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir, 'IH_virus_emergence_rate')

    return df
 
def get_IH_lineages_data_experiment(experiment_name):
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    df = pd.DataFrame()
    for simulation_output_dir in simulation_output_dirs:
        new_df = get_IH_lineages_data_simulation(simulation_output_dir)
        df = pd.concat([df,new_df],axis=0)
    return df
    
# -----------------------------------------------------------------------------
#                              OSR fitting
# ------------------------------------------------------------------------------

def filter_sequencing_files_by_simulation_lenght(files, min_sim_lenght):
    """
    Filters sequencing files by keeping only the ones from simulation that 
    lasted at least min_sim_lenght.
    """
    filtered_files = []
    
    for file in files:
        seeded_simulation_output_dir = os.path.dirname(file)
        try: 
            final_time = read_final_time(seeded_simulation_output_dir)
            if final_time >= min_sim_lenght:
                filtered_files.append(file)
        except Exception as e:
            print(f"Error reading final_time: {e}")
    print(f'Keeping files with simulation lenght >= {min_sim_lenght}')
    # print(filtered_files)
    return filtered_files

def create_combined_sequencing_df(source, 
                                  min_seq_number=0,
                                  min_sim_lenght=0,
                                  individual_type=None):
    '''
    Join sequencing_data_regression.csv files.
    
    Args:
        source: Either a string (path to simulation directory) or a list of SSOD paths.
                - If str: Globs all seeds in that directory.
                - If list: Uses only the provided seed paths.
    '''
    
    # 1. Determine the list of files to process
    if isinstance(source, str):
        # Case A: Source is a directory 
        print(f"processing simulation batch: {os.path.basename(source)}")
        csv_files = glob.glob(os.path.join(source, '**', 'sequencing_data_regression.csv'), recursive=True)
    elif isinstance(source, list):
        # Case B: Source is a list of specific SSOD paths (filtered behavior)
        # print(f"processing filtered list of {len(source)} seeds")
        csv_files = [os.path.join(ssod, 'sequencing_data_regression.csv') for ssod in source]
    else:
        return None

    if individual_type:
        print(f"Filtering for individual_type: {individual_type}")

    # 2. Filter files by simulation length
    filtered_csv_files = filter_sequencing_files_by_simulation_lenght(csv_files, min_sim_lenght)
    print(f'Keeping files with at least {min_seq_number} sequences')
    print('')

    # 3. Read and concatenate DataFrames
    data_frames = []
    for csv_file in filtered_csv_files:
       try:
           df = pd.read_csv(csv_file)
           
           # filter by host type
           if individual_type:
               if 'individual_type' in df.columns:
                   df = df[df['individual_type'] == individual_type]
           
           # Check sequence count AFTER filtering 
           if len(df) >= min_seq_number:
               data_frames.append(df)
       except: pass 
       
    print('Total kept files: ', len(data_frames))
    try:
        combined_df = pd.concat(data_frames, ignore_index=True)
        return combined_df
    except:
        print('No sequencing data available! Check filter settings!')
        print('')
        return None

def detect_sod_outliers(sod_df, threshold=1.5, hard_floor=1e-9):
    """
    Identifies outliers using a Hard Floor for failures AND IQR for statistical deviants.
    """
    sod_df['is_outlier'] = 0
    
    # -----------------------------------------------------------
    # HARD FLOOR
    # Any rate below 1e-9 is physically impossible/noise -> Outlier
    # -----------------------------------------------------------
    sod_df.loc[sod_df['observed_substitution_rate'] < hard_floor, 'is_outlier'] = 1

    # IQR CHECK
    valid_df = sod_df[sod_df['is_outlier'] == 0]
    
    if len(valid_df) < 3:
        return sod_df 

    # Log-transform only the valid points
    log_rates = np.log10(valid_df['observed_substitution_rate'] + 1e-15)
    
    q1 = log_rates.quantile(0.25)
    q3 = log_rates.quantile(0.75)
    iqr = q3 - q1
    
    lower_bound = q1 - (threshold * iqr)
    upper_bound = q3 + (threshold * iqr)
    
    # Find indices of statistical outliers within the valid set
    stat_outliers = valid_df[
        (log_rates < lower_bound) | (log_rates > upper_bound)
    ].index
    
    # Mark them in the main DF
    sod_df.loc[stat_outliers, 'is_outlier'] = 1
    
    return sod_df

def get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
                                                parameter,
                                                min_seq_number,
                                                min_sim_lenght,
                                                individual_type=None):
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    type_suffix = f"_{individual_type}" if individual_type else ""
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_combined_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}{type_suffix}.csv')
    return csv_file_path

def get_OSR_vs_parameter_csv_file_path(experiment_name,
                                       parameter,
                                       min_seq_number,
                                       min_sim_lenght,
                                       individual_type=None):
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    type_suffix = f"_{individual_type}" if individual_type else ""
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}{type_suffix}.csv')
    return csv_file_path

def write_OSR_vs_parameter_csv(experiment_name, parameter, min_seq_number=0, min_sim_lenght=0, individual_type=None):
    """
    Calculates OSR for every individual seed and detects outliers per parameter set.
    """
    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name, parameter, min_seq_number, min_sim_lenght, individual_type)
    
    if os.path.exists(csv_file_path):
        return

    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    all_sod_dfs = []

    for sod in simulation_output_dirs:
        parameter_value = sm.get_parameter_value_from_simulation_output_dir(sod, parameter)
        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        
        sod_results = []
        for ssod in seeded_dirs:
            try:
                final_time = read_final_time(ssod)
                sequencing_data = read_sequencing_data_regression(ssod)
                
                if individual_type and 'Individual_type' in sequencing_data.columns:
                    sequencing_data = sequencing_data[sequencing_data['Individual_type'] == individual_type]
                
                seq_number = len(sequencing_data)
                
                if final_time >= min_sim_lenght and seq_number >= min_seq_number:
                    observed_rate = er.tempest_regression(sequencing_data).coef_[0]
                    sod_results.append({
                        parameter: parameter_value,
                        'observed_substitution_rate': observed_rate,
                        'ssod_path': ssod,
                        'simulation_id': os.path.basename(sod),
                        'is_outlier': 0
                    })
            except Exception:
                continue
        
        # Detect outliers for this SOD
        if sod_results:
            sod_df = pd.DataFrame(sod_results)
            sod_df = detect_sod_outliers(sod_df)
            all_sod_dfs.append(sod_df)

    if all_sod_dfs:
        final_df = pd.concat(all_sod_dfs, ignore_index=True)
        final_df.sort_values(by=str(parameter)).to_csv(csv_file_path, index=False)
    else:
        print("No results to write for individual OSR.")

def read_OSR_vs_parameter_csv(experiment_name, 
                              parameter,
                              min_seq_number,
                              min_sim_lenght,
                              individual_type=None,
                              include_outliers=False): 
    
    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name,
                                                       parameter,
                                                       min_seq_number,
                                                       min_sim_lenght,
                                                       individual_type)
    
    if not os.path.exists(csv_file_path):
        print(f"[Warning] File not found: {csv_file_path}")
        return pd.DataFrame()
        
    df = pd.read_csv(csv_file_path)

    # Filter out outliers if include_outliers is False
    if not include_outliers and 'is_outlier' in df.columns:
        initial_count = len(df)
        df = df[df['is_outlier'] == 0]
        # removed_count = initial_count - len(df)
        # print(f"Excluded {removed_count} outliers.")
        
    return df

def write_combined_OSR_vs_parameter_csv(experiment_name, 
                                        parameter, 
                                        min_seq_number=0,
                                        min_sim_lenght=0,
                                        individual_type=None,
                                        include_outliers=False):
    ''' 
    Create df of observed substitution rate (tempest regression on joint data).
    Uses individual OSR results to filter out outliers before combining, unless include_outliers=True.
    '''
    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name, 
                                                                parameter, 
                                                                min_seq_number,
                                                                min_sim_lenght,
                                                                individual_type)
    
    if os.path.exists(csv_file_path):
        return

    # Load individual OSR results 
    individual_osr_df = read_OSR_vs_parameter_csv(
        experiment_name, parameter, min_seq_number, min_sim_lenght, individual_type,
        include_outliers=include_outliers
    )

    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    results = []

    for sod in simulation_output_dirs:
        parameter_value = sm.get_parameter_value_from_simulation_output_dir(sod, parameter)
        sod_id = os.path.basename(sod)
        combined_df = None
        
        if not individual_osr_df.empty:
            # Get list of valid paths from the individual DF
            clean_paths = individual_osr_df[
                individual_osr_df['simulation_id'] == sod_id
            ]['ssod_path'].tolist()
            
            if clean_paths:
                 # Pass LIST of paths to the function
                 combined_df = create_combined_sequencing_df(clean_paths, 
                                                             min_seq_number, 
                                                             min_sim_lenght, 
                                                             individual_type)
        
        # Fallback Strategy 
        if combined_df is None or combined_df.empty:
             # Pass directory path to the function
             combined_df = create_combined_sequencing_df(sod, 
                                                         min_seq_number, 
                                                         min_sim_lenght, 
                                                         individual_type)

        # Perform Regression
        if combined_df is not None and not combined_df.empty: 
            try:
                observed_substitution_rate = er.tempest_regression(combined_df).coef_[0]
                results.append({str(parameter): parameter_value, 
                                'observed_substitution_rate': observed_substitution_rate})
            except Exception as e:
                print(f"Regression failed for {sod}: {e}")

    results_df = pd.DataFrame(results)
    if not results_df.empty:
        df = results_df.sort_values(by=str(parameter))
        df.to_csv(csv_file_path, index=False)
    else:
        print("No results to write for combined OSR.")

def read_combined_OSR_vs_parameter_csv(experiment_name, 
                                       parameter,
                                       min_seq_number,
                                       min_sim_lenght,
                                       individual_type=None):
    # Reads the file generated by write_combined_OSR...
    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
                                                                parameter,
                                                                min_seq_number,
                                                                min_sim_lenght,
                                                                individual_type)
    return pd.read_csv(csv_file_path)

def get_mean_std_OSR(experiment_name,
                     parameter,
                     min_seq_number,
                     min_sim_lenght,
                     individual_type=None,
                     include_outliers=False): 
    
    # Ensure the CSV exists 
    write_OSR_vs_parameter_csv(experiment_name, 
                               parameter, 
                               min_seq_number, 
                               min_sim_lenght,
                               individual_type)
    
    # Read the data
    df = read_OSR_vs_parameter_csv(experiment_name,
                                   parameter,
                                   min_seq_number,
                                   min_sim_lenght,
                                   individual_type,
                                   include_outliers=include_outliers)
    
    if df.empty:
        return pd.DataFrame()

    # Calculate statistics 
    df_mean_std = df.groupby(parameter)['observed_substitution_rate'].agg(['mean','std']).reset_index()
    return df_mean_std

def get_fit_results_filepath(experiment_name, model_type, individual_type=None):
    experiment_fit_result_dir = dm.get_experiment_fit_result_dir(experiment_name)
    
    type_suffix = f"_{individual_type}" if individual_type else ""
    filename = f'{experiment_name}_{model_type}_fit_results{type_suffix}.csv'
    
    return os.path.join(experiment_fit_result_dir, filename)

def write_fit_results_csv(experiment_name, model_type, fit_result, individual_type=None):
    # Pass individual_type to get unique filepath
    fit_results_filepath = get_fit_results_filepath(experiment_name, model_type, individual_type)
    
    # Save best-fit parameters to CSV
    param_dict = {name: param.value for name, param in fit_result.params.items()}
    df = pd.DataFrame.from_dict(param_dict, orient='index', columns=['Best Fit'])
    df.to_csv(fit_results_filepath, index=True, header=True)
    print(f'Saved fit results to: {fit_results_filepath}')

def read_fit_results_csv(experiment_name, model_type, individual_type=None):
    fit_results = get_fit_results_filepath(experiment_name, model_type, individual_type)
    
    if not os.path.exists(fit_results):
        print(f"[Warning] Fit results file not found: {fit_results}")
        return None

    df = pd.read_csv(fit_results, index_col=0)
    best_fit_df = pd.to_numeric(df['Best Fit'], errors='coerce')
    return best_fit_df

## -----------------------------------------------------------------------------
##                              OSR fitting
## ------------------------------------------------------------------------------
#
#def get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
#                                                parameter,
#                                                min_seq_number,
#                                                min_sim_lenght,
#                                                individual_type=None):
#    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
#    
#    # Add a suffix if a specific host type filter is used
#    type_suffix = f"_{individual_type}" if individual_type else ""
#    
#    csv_file_path = os.path.join(experiment_output_dir, 
#      f'{experiment_name}_combined_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}{type_suffix}.csv')
#    return csv_file_path
#
#def get_OSR_vs_parameter_csv_file_path(experiment_name,
#                                       parameter,
#                                       min_seq_number,
#                                       min_sim_lenght,
#                                       individual_type=None):
#    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
#    
#    type_suffix = f"_{individual_type}" if individual_type else ""
#    
#    csv_file_path = os.path.join(experiment_output_dir, 
#      f'{experiment_name}_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}{type_suffix}.csv')
#    return csv_file_path
#
#def filter_sequencing_files_by_simulation_lenght(files, min_sim_lenght):
#    """
#    Filters sequencing files by keeping only the ones from simulation that 
#    lasted at least min_sim_lenght.
#    """
#    filtered_files = []
#    
#    for file in files:
#        seeded_simulation_output_dir = os.path.dirname(file)
#        try: 
#            final_time = read_final_time(seeded_simulation_output_dir)
#            if final_time >= min_sim_lenght:
#                filtered_files.append(file)
#        except Exception as e:
#            print(f"Error reading final_time: {e}")
#    print(f'Keeping files with simulation lenght >= {min_sim_lenght}')
#    # print(filtered_files)
#    return filtered_files
#
#def detect_sod_outliers(sod_df, threshold=1.5):
#    """
#    Identifies outliers within a single Simulation Output Directory (SOD) 
#    using the Interquartile Range (IQR) method on log-transformed rates.
#    """
#    # Default everything to not being an outlier (0)
#    sod_df['is_outlier'] = 0
#    
#    # We need at least 3 points to define a distribution (Q1, Median, Q3)
#    if len(sod_df) < 3:
#        return sod_df
#
#    # Log-transform for the calculation (base 10)
#    # We add a tiny epsilon (1e-15) to prevent log(0) 
#    log_rates = np.log10(sod_df['observed_substitution_rate'] + 1e-15)
#    
#    # Calculate Quartiles
#    q1 = log_rates.quantile(0.25)
#    q3 = log_rates.quantile(0.75)
#    iqr = q3 - q1
#    
#    # Bounds
#    lower_bound = q1 - (threshold * iqr)
#    upper_bound = q3 + (threshold * iqr)
#    
#    # Mark rows that fall outside the bounds
#    sod_df['is_outlier'] = ((log_rates < lower_bound) | (log_rates > upper_bound)).astype(int)
#    
#    return sod_df
#
#def create_combined_sequencing_df(filtered_ssod, individual_type=None):
#    '''
#    Join sequencing files only from a specific list of filtered ssods.
#    '''
#    data_frames = []
#    for ssod in filtered_ssod:
#        csv_file = os.path.join(ssod, 'sequencing_data_regression.csv')
#        try:
#            df = pd.read_csv(csv_file)
#            # Apply host filter if specified
#            if individual_type and 'individual_type' in df.columns:
#                df = df[df['individual_type'] == individual_type]
#            data_frames.append(df)
#        except Exception:
#            continue
#            
#    return pd.concat(data_frames, ignore_index=True) if data_frames else None
#    
#def write_OSR_vs_parameter_csv(experiment_name, parameter, min_seq_number=0, min_sim_lenght=0, individual_type=None):
#    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name, parameter, min_seq_number, min_sim_lenght, individual_type)
#    
#    if os.path.exists(csv_file_path):
#        return
#
#    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
#    all_sod_dfs = []
#
#    for sod in simulation_output_dirs:
#        parameter_value = sm.get_parameter_value_from_simulation_output_dir(sod, parameter)
#        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
#        
#        sod_results = []
#        for ssod in seeded_dirs:
#            try:
#                final_time = read_final_time(ssod)
#                sequencing_data = read_sequencing_data_regression(ssod)
#                
#                if individual_type and 'individual_type' in sequencing_data.columns:
#                    sequencing_data = sequencing_data[sequencing_data['individual_type'] == individual_type]
#                
#                seq_number = len(sequencing_data)
#                
#                if final_time >= min_sim_lenght and seq_number >= min_seq_number:
#                    # Perform regression
#                    observed_rate = er.tempest_regression(sequencing_data).coef_[0]
#                    
#                    sod_results.append({
#                        parameter: parameter_value,
#                        'observed_substitution_rate': observed_rate,
#                        'ssod_path': ssod, # The "Breadcrumb"
#                        'simulation_id': os.path.basename(sod), # The Group ID
#                        'is_outlier': 0 # Default value
#                    })
#            except Exception:
#                continue
#        
#        # Create DF for this specific parameter set (SOD)
#        if sod_results:
#            sod_df = pd.DataFrame(sod_results)
#            # Detect outliers within this SOD logic
#            sod_df = detect_sod_outliers(sod_df)
#            all_sod_dfs.append(sod_df)
#
#    if all_sod_dfs:
#        final_df = pd.concat(all_sod_dfs, ignore_index=True)
#        final_df.sort_values(by=str(parameter)).to_csv(csv_file_path, index=False)
#    else:
#        print("No results to write for individual OSR.")
#
#def write_combined_OSR_vs_parameter_csv(experiment_name, 
#                                        parameter, 
#                                        min_seq_number=0,
#                                        min_sim_lenght=0,
#                                        individual_type=None):
#    ''' 
#    Create df of observed substitution rate (tempest regression on joint data).
#    Intelligently filters out outliers if individual OSR data is available.
#    '''
#    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name, 
#                                                                parameter, 
#                                                                min_seq_number,
#                                                                min_sim_lenght,
#                                                                individual_type)
#    
#    if os.path.exists(csv_file_path):
#        return
#
#    # 1. Try to load individual OSR results to identify outliers
#    individual_osr_df = read_OSR_vs_parameter_csv(
#        experiment_name, parameter, min_seq_number, min_sim_lenght, individual_type
#    )
#
#    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
#    results = []
#
#    for sod in simulation_output_dirs:
#        parameter_value = sm.get_parameter_value_from_simulation_output_dir(sod, parameter)
#        sod_id = os.path.basename(sod)
#        combined_df = None
#
#        # 2. Smart Selection Strategy
#        if not individual_osr_df.empty and 'is_outlier' in individual_osr_df.columns:
#            # Filter for clean SSOD paths for THIS specific parameter value (SOD)
#            clean_paths = individual_osr_df[
#                (individual_osr_df['simulation_id'] == sod_id) & 
#                (individual_osr_df['is_outlier'] == 0)
#            ]['ssod_path'].tolist()
#            
#            if clean_paths:
#                 combined_df = create_combined_sequencing_df_from_paths(clean_paths, individual_type)
#        
#        # 3. Fallback Strategy (if no OSR file or no clean paths found)
#        if combined_df is None or combined_df.empty:
#             # Use the original broad method
#             combined_df = create_combined_sequencing_df(sod, min_seq_number, min_sim_lenght, individual_type)
#
#        # 4. Perform Regression
#        if combined_df is not None and not combined_df.empty: 
#            try:
#                observed_substitution_rate = er.tempest_regression(combined_df).coef_[0]
#                results.append({str(parameter): parameter_value, 
#                                'observed_substitution_rate': observed_substitution_rate})
#            except Exception as e:
#                print(f"Regression failed for {sod}: {e}")
#
#    results_df = pd.DataFrame(results)
#    if not results_df.empty:
#        df = results_df.sort_values(by=str(parameter))
#        df.to_csv(csv_file_path, index=False)
#    else:
#        print("No results to write for combined OSR.")
#        
#def read_combined_OSR_vs_parameter_csv(experiment_name, 
#                                       parameter,
#                                       min_seq_number,
#                                       min_sim_lenght,
#                                       individual_type=None):
#    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
#                                                                parameter,
#                                                                min_seq_number,
#                                                                min_sim_lenght,
#                                                                individual_type)
#    combined_observed_substitution_rates_df = pd.read_csv(csv_file_path)
#    return combined_observed_substitution_rates_df
#
#def write_OSR_vs_parameter_csv(experiment_name, parameter, min_seq_number=0, min_sim_lenght=0, individual_type=None):
#    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name, parameter, min_seq_number, min_sim_lenght, individual_type)
#    
#    if os.path.exists(csv_file_path):
#        return
#
#    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
#    all_sod_dfs = []
#
#    for sod in simulation_output_dirs:
#        parameter_value = sm.get_parameter_value_from_simulation_output_dir(sod, parameter)
#        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
#        
#        sod_results = []
#        for ssod in seeded_dirs:
#            try:
#                final_time = read_final_time(ssod)
#                sequencing_data = read_sequencing_data_regression(ssod)
#                
#                if individual_type and 'individual_type' in sequencing_data.columns:
#                    sequencing_data = sequencing_data[sequencing_data['individual_type'] == individual_type]
#                
#                seq_number = len(sequencing_data)
#                
#                if final_time >= min_sim_lenght and seq_number >= min_seq_number:
#                    # Perform regression
#                    observed_rate = er.tempest_regression(sequencing_data).coef_[0]
#                    
#                    sod_results.append({
#                        parameter: parameter_value,
#                        'observed_substitution_rate': observed_rate,
#                        'ssod_path': ssod, # The "Breadcrumb"
#                        'simulation_id': os.path.basename(sod),
#                        'is_outlier': 0 # Default value
#                    })
#            except Exception:
#                continue
#        
#        # Create DF for this specific parameter set (SOD)
#        if sod_results:
#            sod_df = pd.DataFrame(sod_results)
#            # Detect outliers within this SOD
#            sod_df = detect_sod_outliers(sod_df)
#            all_sod_dfs.append(sod_df)
#
#    if all_sod_dfs:
#        final_df = pd.concat(all_sod_dfs, ignore_index=True)
#        final_df.sort_values(by=str(parameter)).to_csv(csv_file_path, index=False)
#    else:
#        print("No results to write.")
#    
#def read_OSR_vs_parameter_csv(experiment_name, 
#                              parameter,
#                              min_seq_number,
#                              min_sim_lenght,
#                              individual_type=None,
#                              exclude_outliers=False): 
#    
#    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name,
#                                                       parameter,
#                                                       min_seq_number,
#                                                       min_sim_lenght,
#                                                       individual_type)
#    
#    if not os.path.exists(csv_file_path):
#        print(f"[Warning] File not found: {csv_file_path}")
#        return pd.DataFrame()
#        
#    observed_substitution_rates_df = pd.read_csv(csv_file_path)
#
#    # Filter: only keep rows where is_outlier is 0
#    if exclude_outliers and 'is_outlier' in observed_substitution_rates_df.columns:
#        initial_count = len(observed_substitution_rates_df)
#        observed_substitution_rates_df = observed_substitution_rates_df[
#            observed_substitution_rates_df['is_outlier'] == 0
#        ]
#        removed_count = initial_count - len(observed_substitution_rates_df)
#        print(f"Excluded {removed_count} outliers from the data.")
#        
#    return observed_substitution_rates_df
#
#def get_mean_std_OSR(experiment_name,
#                     parameter,
#                     min_seq_number,
#                     min_sim_lenght,
#                     individual_type=None):
#    
#    write_OSR_vs_parameter_csv(experiment_name, 
#                               parameter, 
#                               min_seq_number, 
#                               min_sim_lenght,
#                               individual_type)
#    
#    df = read_OSR_vs_parameter_csv(experiment_name,
#                                   parameter,
#                                   min_seq_number,
#                                   min_sim_lenght,
#                                   individual_type)
#    
#    if df.empty:
#        return pd.DataFrame()
#
#    df_mean_std = df.groupby(parameter)['observed_substitution_rate'].agg(['mean','std']).reset_index()
#    return df_mean_std
#
#def get_fit_results_filepath(experiment_name, model_type, individual_type=None):
#    experiment_fit_result_dir = dm.get_experiment_fit_result_dir(experiment_name)
#    
#    type_suffix = f"_{individual_type}" if individual_type else ""
#    filename = f'{experiment_name}_{model_type}_fit_results{type_suffix}.csv'
#    
#    return os.path.join(experiment_fit_result_dir, filename)
#
#def write_fit_results_csv(experiment_name, model_type, fit_result, individual_type=None):
#    # Pass individual_type to get unique filepath
#    fit_results_filepath = get_fit_results_filepath(experiment_name, model_type, individual_type)
#    
#    # Save best-fit parameters to CSV
#    param_dict = {name: param.value for name, param in fit_result.params.items()}
#    df = pd.DataFrame.from_dict(param_dict, orient='index', columns=['Best Fit'])
#    df.to_csv(fit_results_filepath, index=True, header=True)
#    print(f'Saved fit results to: {fit_results_filepath}')
#
#def read_fit_results_csv(experiment_name, model_type, individual_type=None):
#    fit_results = get_fit_results_filepath(experiment_name, model_type, individual_type)
#    
#    if not os.path.exists(fit_results):
#        print(f"[Warning] Fit results file not found: {fit_results}")
#        return None
#
#    df = pd.read_csv(fit_results, index_col=0)
#    best_fit_df = pd.to_numeric(df['Best Fit'], errors='coerce')
#    return best_fit_df

# -----------------------------------------------------------------------------
#                              Tree exporter
# -----------------------------------------------------------------------------

def get_tree_filename(experiment_name,
                      seeded_simulation_output_dir,
                      tree_type,
                      tree_subtype,
                      file_type):
    seed_folder = os.path.basename(seeded_simulation_output_dir)
    ext = {'json'  :'.json',
           'newick': '_newick.txt',
           'nexus' : '_nexus.txt',
           'img'   : '.png',
           }
    if file_type not in ext.keys():
        raise ValueError('Invalid file_type selection for tree export.')
    foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    filename = f'{experiment_name}_{foldername}_{seed_folder}_{tree_type}_tree_{tree_subtype}{ext[file_type]}'
    return filename

def get_tree_file_filepath(experiment_name,
                      seeded_simulation_output_dir,
                      tree_type,
                      tree_subtype,
                      file_type):
    experiment_tree_simulation_files_dir = dm.get_experiment_tree_simulation_files_dir(
                                                                    experiment_name,
                                                                    seeded_simulation_output_dir)
    tree_filename = get_tree_filename(experiment_name,
                                seeded_simulation_output_dir,
                                tree_type,
                                tree_subtype,
                                file_type)
    return os.path.join(experiment_tree_simulation_files_dir,tree_filename)

def get_tree_plot_filepath(experiment_name,
                      seeded_simulation_output_dir,
                      tree_type,
                      tree_subtype,
                      file_type='img'):
    experiment_tree_simulation_plots_dir = dm.get_experiment_tree_simulation_plots_dir(
                                                                    experiment_name,
                                                                    seeded_simulation_output_dir)
    tree_filename = get_tree_filename(experiment_name,
                                      seeded_simulation_output_dir,
                                      tree_type,
                                      tree_subtype,
                                      file_type)
    return os.path.join(experiment_tree_simulation_plots_dir,tree_filename)

def export_tree(tree,
                experiment_name,
                seeded_simulation_output_dir,
                tree_type,
                tree_subtype,
                file_type):
    """Export an anytree tree to file"""
    import simplicity.tree.newick as nwk
    import simplicity.tree.nexus as nx
    import simplicity.tree.json_tree as jt
    root = tree[0]
    
    if file_type =='json':
        json_tree_filepath = get_tree_file_filepath(experiment_name,
                              seeded_simulation_output_dir,
                              tree_type,
                              tree_subtype,
                              file_type)
        print(f'saving json file: {json_tree_filepath}')
        jt.write_json_tree_file(root, json_tree_filepath)
        
    if file_type =='nexus':
        nexus_filepath = get_tree_file_filepath(experiment_name,
                              seeded_simulation_output_dir,
                              tree_type,
                              tree_subtype,
                              file_type)
        print(f'saving nexus file: {nexus_filepath}')
        nx.write_nexus_file(tree, nexus_filepath)
       
    if file_type == 'newick':
        newick_filepath = get_tree_file_filepath(experiment_name,
                              seeded_simulation_output_dir,
                              tree_type,
                              tree_subtype,
                              file_type)
        print(f'saving newick file: {newick_filepath}')
        nwk.write_newick_file(root, newick_filepath)
    print('')
    print('Tree file exported successfully.')

# -----------------------------------------------------------------------------
#                                R effective 
# -----------------------------------------------------------------------------
def filter_lineage_frequency_df(lineage_frequency_df, threshold):
    """
    Filters the lineage frequency DataFrame by threshold.

    Args:
        lineage_frequency_df: raw DataFrame from read_lineage_frequency()
        threshold: float, minimum frequency a lineage must reach at any time

    Returns:
        filtered_df: pivoted and filtered DataFrame (Time_sampling x Lineage)
    """
    pivot_df = lineage_frequency_df.pivot(
        index='Time_sampling',
        columns='Lineage_name',
        values='Frequency_at_t'
    )

    selected_lineages = pivot_df.columns[pivot_df.max() >= threshold].tolist()
    filtered_df = pivot_df[selected_lineages]

    return filtered_df, selected_lineages

def get_r_effective_population_csv_filepath(experiment_name, 
                                            seeded_simulation_output_dir,
                                            window_size, threshold):
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    filename = f"{experiment_name}_{foldername}_{seed}_r_effective_population_win_{window_size}_ts_{threshold}.csv"
    return os.path.join(experiment_output_dir, filename)

def get_r_effective_lineages_csv_filepath(experiment_name, 
                                          seeded_simulation_output_dir, 
                                          window_size, threshold):
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)
    foldername = dm.get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    seed = os.path.basename(seeded_simulation_output_dir)
    filename = f"{experiment_name}_{foldername}_{seed}_r_effective_lineages_win_{window_size}_ts_{threshold}.csv"
    return os.path.join(experiment_output_dir, filename)

def get_r_effective_population_traj(ssod, time_window):
    """
    Compute and return the population-level R_effective trajectory as a pandas Series.
    """
    df = read_individuals_data(ssod)
    # df = df[df['t_infection'] > 0]
    final_time = read_final_time(ssod)

    timeline = np.arange(0, final_time, 1)
    r_effective = []

    for t in tqdm(timeline, desc="Computing R_effective (population)"):
        window_start = t - time_window

        births = df[
            (df['t_infectious'] > window_start) & (df['t_infectious'] <= t)
        ].shape[0]

        deaths = df[
            (df['t_not_infectious'] > window_start) & (df['t_not_infectious'] <= t)
        ].shape[0]

        r = births / deaths if deaths > 0 else float('nan')
        r_effective.append((t, r))

    return pd.Series([r for _, r in r_effective], index=[t for t, _ in r_effective], name='R_effective')

def get_r_effective_lineages_traj(ssod, time_window, threshold):
    """
    Compute and return the lineage-level R_effective trajectories as a dictionary of Series.
    """
    
    def flatten_lineage_events(individuals_df):
        """
        Convert nested IH_lineages_trajectory dicts into a flat DataFrame.
        Returns: DataFrame with columns: individual_id, lineage, ih_birth, ih_death
        """
        rows = []

        for ind_id, row in individuals_df.iterrows():
            traj = row.get('IH_lineages_trajectory', {})
            for lineage, event in traj.items():
                rows.append({
                    'individual_id': ind_id,
                    'lineage': lineage,
                    'ih_birth': event.get('ih_birth'),
                    'ih_death': event.get('ih_death')
                })

        return pd.DataFrame(rows)
    
    individuals_df = read_individuals_data(ssod)
    # individuals_df = individuals_df[individuals_df['t_infection'] > 0]
    final_time = read_final_time(ssod)
    # Read and filter lineage frequency
    lineage_freq_df = read_lineage_frequency(ssod)
    _, filtered_lineages = filter_lineage_frequency_df(lineage_freq_df, threshold)
    
    flat_df = flatten_lineage_events(individuals_df)
    timeline = np.arange(0, final_time, 1)
    lineage_r_dict = {}

    for lineage in tqdm(filtered_lineages, desc="Computing R_effective (lineages)"):
        r_series = []
        sub = flat_df[flat_df['lineage'] == lineage]

        for t in timeline:
            window_start = t - time_window

            births = sub[
                (sub['ih_birth'] > window_start) & (sub['ih_birth'] <= t)
            ].shape[0]

            deaths = sub[
                (sub['ih_death'] > window_start) & (sub['ih_death'] <= t)
            ].shape[0]

            r = births / deaths if deaths > 0 else float('nan')
            r_series.append((t, r))

        lineage_r_dict[lineage] = pd.Series(
            [r for _, r in r_series],
            index=[t for t, _ in r_series],
            name=f'R_eff_{lineage}'
        )

    return lineage_r_dict

def write_r_effective_trajs_csv(experiment_name, seeded_simulation_output_dir, 
                                time_window, threshold):
    """
    Compute and write R_effective (population + lineage) trajectories to CSV.
    """
    # Compute R_effective
    r_effective_population_traj = get_r_effective_population_traj(seeded_simulation_output_dir, time_window)
    r_effective_lineages_traj = get_r_effective_lineages_traj(seeded_simulation_output_dir, time_window, threshold)

    # Get filepaths
    pop_fp = get_r_effective_population_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)
    line_fp = get_r_effective_lineages_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)

    # Save
    r_effective_population_traj.to_csv(pop_fp, header=True)
    pd.DataFrame(r_effective_lineages_traj).to_csv(line_fp)

    print(f"R_effective CSVs saved to:\n- {pop_fp}\n- {line_fp}")
    
def read_r_effective_trajs_csv(experiment_name, seeded_simulation_output_dir, 
                               time_window, threshold):
    """
    Read R_effective (population + lineage) trajectories from CSV.
    If files are missing, issues a warning and generates them.

    Returns:
        - r_effective_population_traj: pd.Series
        - r_effective_lineages_traj: dict of {lineage_name: pd.Series}
    """
    pop_fp = get_r_effective_population_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)
    line_fp = get_r_effective_lineages_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)

    # Check if files exist
    if not os.path.exists(pop_fp) or not os.path.exists(line_fp):
        warnings.warn("R_effective CSVs missing. Recomputing and writing them.")
        write_r_effective_trajs_csv(experiment_name, seeded_simulation_output_dir, time_window, threshold)

    # Read population R_eff
    r_eff_pop = pd.read_csv(pop_fp, index_col=0)
    if isinstance(r_eff_pop, pd.DataFrame):
        r_eff_pop = r_eff_pop.iloc[:, 0]

    # Read lineage R_eff
    df_line = pd.read_csv(line_fp, index_col=0)
    r_eff_lin_dict = {lineage: df_line[lineage] for lineage in df_line.columns}

    return r_eff_pop, r_eff_lin_dict

# -----------------------------------------------------------------------------

def get_procomputed_matrix_table_filepath(tau_1,tau_2,tau_3,tau_4):
    file_name = f"precomputed_A_exponentials_tau1={tau_1}_tau2={tau_2}_tau3={tau_3}_tau4={tau_4}.pkl"
    data_dir = dm.get_data_dir()
    return os.path.join(data_dir,file_name)

# -----------------------------------------------------------------------------

def get_clustering_table_filepath(experiment_name: str, ssod: str) -> str:
    """
    gets the filepath: 07_Clustering/{SOD}_seed_{SEED}_clustering.csv
    """
    cluster_dir = dm.get_experiment_cluster_dir(experiment_name)
    sod = dm.get_simulation_output_foldername_from_SSOD(ssod)  
    seed = dm.get_seed_from_SSOD(ssod)                        
    fname = f"{sod}_seed_{seed}_clustering.csv"
    return os.path.join(cluster_dir, fname)


def write_clustering_table(experiment_name: str, ssod: str, df: pd.DataFrame):
    """
    Writes the clustering CSV. 
    """
    path = get_clustering_table_filepath(experiment_name, ssod)
    df.to_csv(path, index=False)

def read_clustering_table(experiment_name: str, ssod: str) -> pd.DataFrame:
    """
    Reads the single clustering CSV and literal-evals nested columns.
    """
    path = get_clustering_table_filepath(experiment_name, ssod)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Clustering table not found: {path}")

    # peek header to decide which converters to apply
    import csv
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)

    converters = {}
    for col in ("lineages", "defining_mutations"):
        if col in header:
            converters[col] = ast.literal_eval

    return pd.read_csv(path, converters=converters) 

def get_nextstrain_dataset_paths(experiment_name: str, ssod: str):
    """
    Returns (dataset_json_path, metadata_tsv_path, base_name) for a given SSOD.
    """
    # Base dir 
    ns_root = Path(dm.get_nextstrain_dir(experiment_name))

    # Get the original name parts
    sod  = dm.get_simulation_output_foldername_from_SSOD(ssod)
    seed = dm.get_seed_from_SSOD(ssod)

    # Clean the long simulation name from underscores 
    clean_sod_name = sod.replace("_", "")
    
    # Create the new base filename
    base = f"seed{seed}_{clean_sod_name}"

    dataset_json_path = ns_root / f"{base}-nextstrain.json"
    metadata_tsv_path = ns_root / f"{base}-metadata.tsv"

    return dataset_json_path, metadata_tsv_path, base

    