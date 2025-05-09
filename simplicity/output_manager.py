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
def save_sequencing_dataset(simulation_output, output_path):
    """
    Writes the simulated sequencing data to a FASTA file.
    A genome from simulation_output.sequencing_data contains
    [
       patient id,
       time of sequencing,
       single encoded genome,
       len(genome) (number of substitutions),
       individuals[patient]['type'],
       infection lenght
       intra-host genome id (integer)
    ]
    """
    try:
        # dataframe for phylogenetic regression
        sequencing_data_dic = []
        # paths for files to save
        fasta_file_path = os.path.join(output_path,'sequencing_data.fasta')
        csv_file_path = os.path.join(output_path,'sequencing_data_regression.csv')
        # open fasta file to write on
        with open(fasta_file_path, 'w') as fasta_file:
             # loop over the list of sequencing data, getting the identifiers
             # and the sequences to write them in the fasta file
             for individual_sequencing_data in simulation_output.sequencing_data:
                # get seq ID
                sequences_id = f'>{individual_sequencing_data[0]}_{individual_sequencing_data[6]}_time_{individual_sequencing_data[1]:.2f}'
                # Write the sequence identifier
                fasta_file.write(f"{sequences_id}\n")
                # get full genome
                decoded_genome = decode_genome(individual_sequencing_data[2])
                # remove ' from genome string and write to fasta
                genome = decoded_genome.replace("'","")
                # Write the sequence data
                fasta_file.write(f"{genome}\n")     
                
                # add sequence data to df for regression 
                sequencing_data_dic.append({
                    'Sequencing_time': round(individual_sequencing_data[1]/365.25,4), #seq time in years
                    'Distance_from_root': dis.hamming_true(individual_sequencing_data[2])/len(dis.reference) #Subst/site - normed
                           })
        # save dictionary of sequencing data for regression
        with open(csv_file_path, mode='w', newline='') as file:
            # Create a DictWriter object and pass the fieldnames (keys of the dictionary)
            # Use the keys of the first dictionary to get the column names
            fieldnames = sequencing_data_dic[0].keys()  
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            
            # Write the header (column names)
            writer.writeheader()
            
            # Write the rows (the actual data)
            writer.writerows(sequencing_data_dic)
    except:
            if simulation_output.sequencing_data == []:
                print('No individuals sequenced during simulation, not saving FASTA file')
                print('')
            else: 
                raise ValueError('There is something wrong with the sequencing data')

def read_sequencing_data(seeded_simulation_output_dir):
    sequencing_data_file_path = os.path.join(seeded_simulation_output_dir,
                                        "sequencing_data_regression.csv")
    df = pd.read_csv(sequencing_data_file_path)
    return df
        
def save_simulation_trajectory(simulation_output, seeded_simulation_output_dir):
    df = pd.DataFrame(simulation_output.trajectory, columns= 
                      ['time','infected','diagnosed','recovered','deceased',
                       'infectious_normal','detectables','susceptibles'])
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

# def save_DEBUG_update_ih(simulation_output, seeded_simulation_output_dir):
#     df = simulation_output.DEBUG_update_ih_to_df()
#     DEBUG_update_ih_file_path = os.path.join(seeded_simulation_output_dir,
#                                                "DEBUG_update_ih.csv")
#     df.to_csv(DEBUG_update_ih_file_path, index=False)

# def read_DEBUG_update_ih(seeded_simulation_output_dir):
#     DEBUG_update_ih_file_path = os.path.join(seeded_simulation_output_dir,
#                                                "DEBUG_update_ih.csv")
#     df = pd.read_csv(DEBUG_update_ih_file_path)
#     return df

def save_individuals_data(simulation_output, seeded_simulation_output_dir):
    individuals_data = simulation_output.individuals_data_to_df()
    individuals_data.drop('model',axis=1,inplace=True)
    individuals_data_file_path = os.path.join(seeded_simulation_output_dir,
                                               "individuals_data.csv")
    individuals_data.to_csv(individuals_data_file_path)
    
def read_individuals_data(seeded_simulation_output_dir):
    trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                        "individuals_data.csv")
    df = pd.read_csv(trajectory_file_path, index_col=0)
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
    print('')
    print(filtered_files)
    return filtered_files

def filter_files_by_min_sequences(files, min_seq_number):
    """
    Filters CSV files by min_seq_number sequences.
    """
    filtered_files = []
    
    for file in files:
        try:
            df = pd.read_csv(file)
            if len(df) >= min_seq_number:
                filtered_files.append(file)
        except Exception as e:
            pass
            # print(f"Error reading CSV file {file}: {e}")
    
    print(f'Keeping files with at least {min_seq_number} sequences')
    print('')
    return filtered_files

def create_combined_sequencing_df(seeeded_simulations_output_directory, 
                                  min_seq_number=0,
                                  min_sim_lenght=0):
    '''
    seeeded_simulations_output_directory ==> path to subfolder of 
                                             experiment_name/04_Output/
    
    Join all sequencing_data_regression.csv files of different 
    seeded simulation runs (SAME PARAMETERS, different seeds) 
    in a single df and returns it (for tempest regression).
    '''
    print('##################################################################')
    print(f"processing simulation batch: {os.path.basename(seeeded_simulations_output_directory)}")
    csv_files = glob.glob(os.path.join(seeeded_simulations_output_directory,'**',
                                       'sequencing_data_regression.csv'),
                                        recursive=True)
    filtered_csv_files = filter_sequencing_files_by_simulation_lenght(csv_files, min_sim_lenght)
    filtered_csv_files = filter_files_by_min_sequences(filtered_csv_files, min_seq_number)
    
    data_frames = []
    # Concatenate all DataFrames into one
    try:
        combined_df = pd.concat(data_frames, ignore_index=True)
        return combined_df
    except:
        print('No sequencing data available to plot! Check filter settings!')
        print('')
        return None

def get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
                                                parameter,
                                                min_seq_number,
                                                min_sim_lenght):
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_combined_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}.csv')
    return csv_file_path

def write_combined_OSR_vs_parameter_csv(experiment_name, 
                                        parameter, 
                                        min_seq_number=0,
                                        min_sim_lenght=0):
    ''' Create df of observed substitution rate (tempest regression on joint data) and parameter values
    '''
    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name, 
                                                                parameter, 
                                                                min_seq_number,
                                                                min_sim_lenght)
    if os.path.exists(csv_file_path):
        # print(f'{csv_file_path} already exists, continuing...')
        pass
    else:
        # Get seeded simulations output subfolders
        simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
     
        results = []
        for simulation_output_dir in simulation_output_dirs:
            # Read the parameter from the settings file
            parameter_value = sm.get_parameter_value_from_simulation_output_dir(
                                                  simulation_output_dir, parameter)
            
            # Perform regression for the settings output folder
            combined_df = create_combined_sequencing_df(simulation_output_dir, 
                                                        min_seq_number,
                                                        min_sim_lenght)
            if combined_df is None: 
                pass
            else:
                observed_substitution_rate = er.tempest_regression(combined_df).coef_[0] # substitution rate per site per year
                results.append({str(parameter): parameter_value, 
                                'observed_substitution_rate': observed_substitution_rate})
        # add results to df
        results_df = pd.DataFrame(results)
        
        # Sort the results by the 'parameter' values
        df = results_df.sort_values(by=str(parameter))
    
        # Save the results to a CSV file
        df.to_csv(csv_file_path, index=False)
        
def read_combined_OSR_vs_parameter_csv(experiment_name, 
                                       parameter,
                                       min_seq_number,
                                       min_sim_lenght):
    csv_file_path = get_combined_OSR_vs_parameter_csv_file_path(experiment_name,
                                                                parameter,
                                                                min_seq_number,
                                                                min_sim_lenght)
    combined_observed_substitution_rates_df = pd.read_csv(csv_file_path)
    return combined_observed_substitution_rates_df

def get_OSR_vs_parameter_csv_file_path(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght):
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_OSR_vs_{parameter}_sim_lenght_{min_sim_lenght}_seq_n_{min_seq_number}.csv')
    return csv_file_path
    
def write_OSR_vs_parameter_csv(experiment_name, 
                                parameter, 
                                min_seq_number=0,
                                min_sim_lenght=0):
    ''' Create df of observed substitution rate (tempest regression) and parameter values
    '''
    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name,
                                                      parameter,
                                                      min_seq_number,
                                                      min_sim_lenght)
    # if the file exists already, does nothing
    if os.path.exists(csv_file_path):
        # print(f'{csv_file_path} already exists, continuing...')
        pass
    # create the df and saves it to csv 
    else:
        # Get seeded simulations output subfolders
        simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
        results = []
        # loop over each simulation_output_dir
        for simulation_output_dir in simulation_output_dirs:
            # Read the parameter from the settings file
            parameter_value = sm.get_parameter_value_from_simulation_output_dir(
                                                  simulation_output_dir, parameter)
            seeded_simulation_output_dirs = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
            # loop over each seeded_simulation_output_dir
            for seeded_simulation_output_dir in seeded_simulation_output_dirs:
                try:
                    # get simulation final time
                    final_time = read_final_time(seeded_simulation_output_dir)
                    # read sequencing_data
                    sequencing_data = read_sequencing_data(seeded_simulation_output_dir)
                    seq_number = len(sequencing_data)
                    # filter by simulation lenght and seq_number
                    if final_time >= min_sim_lenght and seq_number >= min_seq_number:
                        # Perform regression for each sequencing file
                        observed_substitution_rate = er.tempest_regression(sequencing_data).coef_[0] # substitution rate per site per year
                        results.append({
                            parameter: parameter_value, 
                            'observed_substitution_rate': observed_substitution_rate,
                            'simulation_final_time':read_final_time(seeded_simulation_output_dir),
                            'settings_final_time': sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir, 'final_time')
                                       })
                except: pass # only add files that are present
        
        # add fit results to df
        df = pd.DataFrame(results)
        
        # Sort the df by the 'parameter' values
        df = df.sort_values(by=str(parameter))
        
        # Save to a CSV file
        df.to_csv(csv_file_path, index=False)

def read_OSR_vs_parameter_csv(experiment_name, 
                              parameter,
                              min_seq_number,
                              min_sim_lenght):
    csv_file_path = get_OSR_vs_parameter_csv_file_path(experiment_name,
                                                      parameter,
                                                      min_seq_number,
                                                      min_sim_lenght)
    observed_substitution_rates_df = pd.read_csv(csv_file_path)
    return observed_substitution_rates_df

def get_mean_std_OSR(experiment_name,
                     parameter,
                     min_seq_number,
                     min_sim_lenght):
    write_OSR_vs_parameter_csv(experiment_name, 
                               parameter, 
                               min_seq_number, 
                               min_sim_lenght)
    df = read_OSR_vs_parameter_csv(experiment_name,
                                   parameter,
                                   min_seq_number,
                                   min_sim_lenght)
    df_mean_std = df.groupby(parameter)['observed_substitution_rate'].agg(['mean','std']).reset_index()
    return df_mean_std

def get_fit_results_filepath(experiment_name, model_type):
    experiment_fit_result_dir = dm.get_experiment_fit_result_dir(experiment_name)
    filename = f'{experiment_name}_{model_type}_fit_results.csv'
    return os.path.join(experiment_fit_result_dir,filename)

def write_fit_results_csv(experiment_name, model_type, fit_result):
    fit_results_filepath = get_fit_results_filepath(experiment_name, model_type)
    # Save best-fit parameters to CSV
    param_dict = {name: param.value for name, param in fit_result.params.items()}
    df = pd.DataFrame.from_dict(param_dict, orient='index', columns=['Best Fit'])
    df.to_csv(fit_results_filepath, index=True, header=True)

def read_fit_results_csv(experiment_name, model_type):
    fit_results = get_fit_results_filepath(experiment_name, model_type)
    df = pd.read_csv(fit_results,index_col=0)
    best_fit_df = pd.to_numeric(df['Best Fit'], errors='coerce')
    return best_fit_df
    
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
    
    Returns:
        - r_effective_population_traj: pd.Series
        - r_effective_lineages_traj: dict of {lineage_name: pd.Series}
    """
    pop_fp = get_r_effective_population_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)
    line_fp = get_r_effective_lineages_csv_filepath(experiment_name, seeded_simulation_output_dir, time_window, threshold)

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
    
    



    