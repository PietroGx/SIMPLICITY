'''
In this file we perform the TempEst linear regression to estimate the observed
evolutionary rate u from the simulated data of a SIMPLICITY run. There are also
the functions to plot E (model evolutionary rate) vs u (observed evolutionary rate)
or vs any other simulation parameter.
'''

import pandas as pd
import math
import sklearn.linear_model 
import matplotlib.pyplot as plt
import os
import glob
import json
import simplicity.config as config

def create_joint_df(experiment_name):
    # Retrieve all Sequencing_data.csv files in the output folder
    experiment_output_dir = config.get_experiment_output_dir(experiment_name)
    csv_files = glob.glob(os.path.join(experiment_output_dir,'**','Sequencing_data.csv'), recursive=True)
    
    # List to store individual DataFrames
    data_frames = []
    
    for csv_file in csv_files:
        # Read each CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        data_frames.append(df)
    
    # Concatenate all DataFrames into one
    combined_df = pd.concat(data_frames, ignore_index=True)
    
    return combined_df, csv_files

def perform_regression(df):
    x = df['Sequencing_time'].values.reshape(-1, 1)
    y = df['Distance_from_root'].values
    
    model = sklearn.linear_model.LinearRegression(fit_intercept=False)
    model.fit(x, y)
    
    u = model.coef_[0]
    
    return u, model

def ideal_subplot_grid(num_plots):
    # Calculate the number of columns (the ceiling of the square root of the number of plots)
    num_cols = math.ceil(math.sqrt(num_plots))
    
    # Calculate the number of rows needed
    num_rows = math.ceil(num_plots / num_cols)
    
    return num_rows, num_cols        
 
def plot_combined_regressions(experiment_name):
    # Get list of subdirectories
    experiment_output_dir = config.get_experiment_output_dir(experiment_name)
    # subfolders = [os.path.relpath(f.path, OUTPUT_DIR) for f in os.scandir(experiment_output_folder) if f.is_dir()]
    simulations_subfolders = [os.path.join(experiment_output_dir, f.name) for f in os.scandir(experiment_output_dir) if f.is_dir()]
    
    # Determine the number of rows and columns for subplots
    num_rows, num_cols = ideal_subplot_grid(len(simulations_subfolders))

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    
    # If there's only one subplot, axs is not an array, so we convert it to a 2D array for consistency
    if num_rows * num_cols == 1:
        axs = [[axs]]
    elif num_rows == 1 or num_cols == 1:
        axs = axs.reshape(1, -1)
    else:
        axs = axs

    for i, subfolder in enumerate(simulations_subfolders):
        combined_df, _ = create_joint_df(experiment_name)
        u, model = perform_regression(combined_df)
        
        x = combined_df['Sequencing_time'].values.reshape(-1, 1)
        y_pred = model.predict(x)
        
        ax = axs[i // num_cols][i % num_cols]
        combined_df.plot(kind='scatter', x='Sequencing_time', y='Distance_from_root', color='blue', ax=ax)
        ax.plot(x, y_pred, color='red', linewidth=2, label=f'u = {u:.5f}')
        ax.set_xlabel('Time [y]')
        ax.set_ylabel('Distance from root [#S/site]')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_title(f'Regression - {os.path.basename(subfolder)}')
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(experiment_output_dir, "combined_regression.png"))

def plot_u_vs_parameter(experiment_name, parameter):
    results = []
    
    settings_folder = config.get_simulation_parameters_dir(experiment_name)
    experiment_output_dir   = config.get_output_dir(experiment_name)
    
    settings_files = glob.glob(os.path.join(settings_folder, '*.json'))
    for settings_file in settings_files:
        settings_name = os.path.splitext(os.path.basename(settings_file))[0]
        experiment_output_folder = os.path.join(experiment_output_dir, settings_name)
        
        # Read the parameter from the settings file
        with open(settings_file, 'r') as file:
            data = json.load(file)
        
        target_value = data.get(parameter)
        
        # Perform regression for the settings output folder
        try:
            combined_df, _ = create_joint_df(experiment_name)
        except ValueError:
            continue
        u, _ = perform_regression(combined_df)
        
        results.append({str(parameter): target_value, 'u': u})
    
    results_df = pd.DataFrame(results)
    
    # Sort the results by the 'target_parameter' values
    results_df = results_df.sort_values(by=str(parameter))
    
    # Save the results to a CSV file
    csv_file_path = os.path.join(experiment_output_folder, f'u_vs_{parameter}_values.csv')
    results_df.to_csv(csv_file_path, index=False)
    
    # Plot target parameter vs u as a line plot with points
    plt.figure(figsize=(10, 6))
    plt.plot(results_df[str(parameter)], results_df['u'], marker='o', color='black', linestyle='-', label=f'{target_parameter} vs u')
    plt.xlabel(target_parameter)
    plt.ylabel('u')
    plt.title(f'{parameter} vs u')
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(experiment_output_folder, f"{target_parameter}_vs_u.png"))

if __name__ == "__main__":    
    
    import sys
    experiment_name = sys.argv[1] 
    target_parameter = sys.argv[2] 
    # target_parameter = 'evolutionary rate'
    # experiment_name = 'test_refactoring'
    plot_combined_regressions(experiment_name)
    plot_u_vs_parameter(experiment_name,target_parameter)
