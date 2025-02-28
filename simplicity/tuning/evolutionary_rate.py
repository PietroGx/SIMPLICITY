'''
In this file we perform the TempEst linear regression to estimate the observed
evolutionary rate u from the simulated data of a SIMPLICITY run. There are also
the functions to plot E (model evolutionary rate) vs u (observed evolutionary rate)
or vs any other simulation parameter.
'''

import pandas as pd
import sklearn.linear_model 
import os
import glob
import numpy as np
from lmfit import Model
from lmfit.models import SplineModel
import simplicity.output_manager as om

def filter_sequencing_files_by_simulation_lenght(files, min_sim_lenght):
    """
    Filters sequencing files by keeping only the ones from simulation that 
    lasted at least min_sim_lenght.
    """
    filtered_files = []
    
    for file in files:
        directory = os.path.dirname(file)
        csv_path = os.path.join(directory, 'final_time.csv')
        try:
            with open(csv_path, 'r') as f:
                csv_value = float(f.read().strip())
                if csv_value >= min_sim_lenght:
                    filtered_files.append(file)
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")
    print(f'Keeping files with simulation lenght >= {min_sim_lenght}')
    print('')
    return filtered_files

def create_joint_sequencing_df(seeeded_simulations_output_directory, min_sim_lenght=0):
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
    # List to store individual DataFrames
    data_frames = []
    for csv_file in filtered_csv_files:
        # Read each CSV file into a DataFrame
        try:
            df = pd.read_csv(csv_file)
            data_frames.append(df)
        except:
            pass
    # Concatenate all DataFrames into one
    try:
        combined_df = pd.concat(data_frames, ignore_index=True)
        return combined_df
    except:
        print('No sequencing data available to plot! Check filter settings!')
        print('')
        return None
    
def tempest_regression(sequencing_data_df):
    '''
    perform TempEst regression on dataframe of sequencing data

    Parameters
    ----------
    df : pandas df 
        output of create_joint_sequencing_df.

    Returns
    -------
    u : TYPE
        observed evolutionary rate.
    model : func
        fitted model (sklearn linear regression(.

    '''
    x = sequencing_data_df['Sequencing_time'].values.reshape(-1, 1)
    y = sequencing_data_df['Distance_from_root'].values
    model = sklearn.linear_model.LinearRegression(fit_intercept=False)
    model.fit(x, y)
    observed_evolutionary_rate = model.coef_[0] # substitution rate per site per year
    return observed_evolutionary_rate, model

def factory_model(model_type: str):
    
    # Define the models
    def linear_model(x, A, B):
        return A*x + B

    def log_model(x, A, B, C):
        return A * np.log(B * x + C)

    def exp_model(x, A, B, C):
        return A * x**B + C

    def double_log_model(x, A, B, C, D, E, F):
        return A * np.log(B * x + C) + D * np.log(E * x + F)

    def tan_model(x, A, B, C, D):
        return A * np.tan(B * x - C) + D
    
    # select the model, assign parameters and return it
    if model_type == 'linear':
        model = Model(linear_model)
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=0)
        
        return model, params 
    
    elif model_type == 'log':
        model = Model(log_model)
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        
        return model, params 
    
    elif model_type == 'exp':
        model = Model(exp_model)
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0)
        
        return model, params 
    
    elif model_type == 'double_log':
        model = Model(double_log_model)
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=1, E=1, F=0)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        params['E'].set(min=0.000001)  
        params['F'].set(min=0.000001)
        
        return model, params 
    
    elif model_type == 'tan':
        model = Model(tan_model)
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=0)
        
        return model, params 
    
    if model_type == 'spline':
        # Define knot positions for the spline ensuring they are within range and well spaced
        knots = np.logspace(-3, 0, 5)  
        
        # Initialize the SplineModel with the defined knots
        model = SplineModel(prefix='spline_', xknots=knots)
        
        # Create parameters with initial guesses and set boundaries to avoid instability
        params = model.make_params()
        # for param in params:
        #     params[param].set(min=-2, max=2)  
        return model, params 
    
    else: raise ValueError('Invalid model selection')

def fit_observed_evolutionary_rate(experiment_name, model_type):
    # Read data from CSV
    data = om.read_u_e_values(experiment_name)
    x_data = data['evolutionary_rate'] 
    y_data = data['u']  
    
    weights = 1/y_data
    
    # Create the Model
    model, params = factory_model(model_type)
    
    # Fit the model to the data
    fit_result = model.fit(y_data, params, x=x_data, weights=weights)
    
    # Print the fit results
    print(fit_result.fit_report())
    return fit_result