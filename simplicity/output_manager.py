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
                print('')
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
    df = pd.DataFrame(simulation_output.lineage_frequency)
    lineage_frequency_file_path = os.path.join(seeded_simulation_output_dir,
                                               "lineage_frequency.csv")
    df.to_csv(lineage_frequency_file_path, index=False)
        
def save_individuals_data(simulation_output, seeded_simulation_output_dir):
    individuals_data = simulation_output.data()
    individuals_data.drop('model',axis=1,inplace=True)
    individuals_data_file_path = os.path.join(seeded_simulation_output_dir,
                                               "individuals_data.csv")
    individuals_data.to_csv(individuals_data_file_path)

def save_phylogenetic_data(simulation_output, seeded_simulation_output_dir):
    phylogenetic_data = simulation_output.phylogeny()
    phylogenetic_data_file_path = os.path.join(seeded_simulation_output_dir,
                                               "phylogenetic_data.csv")
    phylogenetic_data.to_csv(phylogenetic_data_file_path)
    
def save_fitness_trajectory(simulation_output, seeded_simulation_output_dir):
    fitness_trajectory = simulation_output.fitness_trajectory_to_df()
    fitness_trajectory_file_path = os.path.join(seeded_simulation_output_dir,
                                               "fitness_trajectory.csv")
    fitness_trajectory.to_csv(fitness_trajectory_file_path)
    
def save_final_time(simulation_output, seeded_simulation_output_dir):
    final_time = simulation_output.time
    final_time_file_path = os.path.join(seeded_simulation_output_dir,
                                               "final_time.csv")
    with open(final_time_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([final_time]) 
                          
def extract_u_e_values(experiment_name):
    ''' Perform tempest regression (get u) and saves the values vs evolutionary rate
    '''
    import glob
    import json
    from simplicity.tuning.evolutionary_rate import tempest_regression, create_joint_sequencing_df
    # get simulation parameter files for the selected experiment
    simulation_parameters_dir = dm.get_simulation_parameters_dir(experiment_name)
    simulation_parameters_files = glob.glob(os.path.join(simulation_parameters_dir, '*.json'))
    # get output directory
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    # Get seeded simulations output subfolders
    seeeded_simulations_output_directories = [os.path.join(experiment_output_dir, 
                                    f.name) for f in os.scandir(experiment_output_dir
                                                                ) if f.is_dir()]
    # for each simulation set in the experiment perfom the tempest regression
    results = []
    for simulation_parameters_file, seeeded_simulations_output_directory in zip(
            simulation_parameters_files,seeeded_simulations_output_directories):
        # Read the parameter from the settings file
        with open(simulation_parameters_file, 'r') as file:
            data = json.load(file)
        parameter_value = data.get('evolutionary_rate')
        
        # Perform regression for the settings output folder
        combined_df, _ = create_joint_sequencing_df(seeeded_simulations_output_directory)
        u, _ = tempest_regression(combined_df)
        
        results.append({str('evolutionary_rate'): parameter_value, 'u': u})
    # add results to df
    results_df = pd.DataFrame(results)
    
    # Sort the results by the 'parameter' values
    results_df = results_df.sort_values(by=str('evolutionary_rate'))
    
    # Save the results to a CSV file
    csv_file_path = os.path.join(experiment_output_dir, 
                                 f'{experiment_name}_u_vs_evolutionary_rate_values.csv')
    results_df.to_csv(csv_file_path, index=False)
    
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
    data = data[['IH_virus_number', 'lineages_number']]

    # DataFrames with all possible IH_virus_number and lineages_numbervalues (1 to 5)
    IH_virus_number_values = pd.DataFrame({'IH_virus_number': [1, 2, 3, 4, 5]})
    lineages_number_values = pd.DataFrame({'lineages_number': [1, 2, 3, 4, 5]})

    # Calculate the counts of each value for each type
    ih_virus_count = data.groupby(['IH_virus_number']).size().reset_index(name='ih_virus_count')
    ih_lineage_count = data.groupby(['lineages_number']).size().reset_index(name='ih_lineage_count')

    # Merge with the IH_virus_number_values DataFrame, filling NaNs with 0
    ih_virus_count = pd.merge(IH_virus_number_values, ih_virus_count, on='IH_virus_number', how='left').fillna({'ih_virus_count': 0})
    ih_lineage_count = pd.merge(lineages_number_values, ih_lineage_count, on='lineages_number', how='left').fillna({'ih_lineage_count': 0})

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

def get_combined_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter):
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_combined_observed_evolutionary_rate_vs_{parameter}_values.csv')
    return csv_file_path

def build_combined_observed_evolutionary_rate_vs_parameter_df(experiment_name, parameter, min_sim_lenght=0):
    ''' Create df of observed evolutionary rate (tempest regression on joint data) and parameter values
    '''
    
    # Get seeded simulations output subfolders
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
 
    results = []
    for simulation_output_dir in simulation_output_dirs:
        # Read the parameter from the settings file
        parameter_value = sm.get_parameter_value_from_simulation_output_dir(
                                              simulation_output_dir, parameter)
        
        # Perform regression for the settings output folder
        combined_df = er.create_joint_sequencing_df(simulation_output_dir, min_sim_lenght)
        if combined_df is None: 
            pass
        else:
            observed_evolutionary_rate = er.tempest_regression(combined_df).coef_[0] # substitution rate per site per year
            results.append({str(parameter): parameter_value, 
                            'observed_evolutionary_rate': observed_evolutionary_rate})
    # add results to df
    results_df = pd.DataFrame(results)
    
    # Sort the results by the 'parameter' values
    results_df = results_df.sort_values(by=str(parameter))
    
    # Save the results to a CSV file
    csv_file_path = get_combined_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name, parameter)
    results_df.to_csv(csv_file_path, index=False)

def get_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter):
    experiment_output_dir     = dm.get_experiment_output_dir(experiment_name)
    csv_file_path = os.path.join(experiment_output_dir, 
      f'{experiment_name}_observed_evolutionary_rate_vs_{parameter}_values.csv')
    return csv_file_path

def build_observed_evolutionary_rates_vs_parameter_df(experiment_name, parameter, min_sim_lenght=0):
    ''' Create df of observed evolutionary rate (tempest regression) and parameter values
    '''
    # Get seeded simulations output subfolders
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
 
    results = []
    for simulation_output_dir in simulation_output_dirs:
        # Read the parameter from the settings file
        parameter_value = sm.get_parameter_value_from_simulation_output_dir(
                                              simulation_output_dir, parameter)
        seeded_simulation_output_dirs = dm.get_seeded_simulation_output_dirs(simulation_output_dir)
        
        for seeded_simulation_output_dir in seeded_simulation_output_dirs:
            # Perform regression for each sequencing file
            try:
                sequencing_data = read_sequencing_data(seeded_simulation_output_dir)
                observed_evolutionary_rate = er.tempest_regression(sequencing_data).coef_[0] # substitution rate per site per year
                results.append({parameter: parameter_value, 
                            'observed_evolutionary_rate': observed_evolutionary_rate})
            except: pass # only add files that are present
    
    # add fit results to df
    df = pd.DataFrame(results)
    
    # Sort the df by the 'parameter' values
    df = df.sort_values(by=str(parameter))
    
    # Save to a CSV file
    csv_file_path = get_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter)
    df.to_csv(csv_file_path, index=False)

def read_observed_evolutionary_rates_csv(experiment_name, parameter):
    csv_file_path = get_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter)
    observed_evolutionary_rates_df = pd.read_csv(csv_file_path)
    return observed_evolutionary_rates_df

def read_combined_observed_evolutionary_rate_csv(experiment_name, parameter):
    csv_file_path = get_combined_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter)
    combined_observed_evolutionary_rates_df = pd.read_csv(csv_file_path)
    return combined_observed_evolutionary_rates_df







    