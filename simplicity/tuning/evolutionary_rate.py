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
    tempest_regression = sklearn.linear_model.LinearRegression(fit_intercept=False)
    fitted_tempest_regression = tempest_regression.fit(x, y)
    
    return fitted_tempest_regression
    
def factory_model_func(model_type: str):
    # Define the models
    def linear_model(x, A, B):
        return A * x + B

    def log_model(x, A, B, C):
        return A * np.log(B * x + C)

    def exp_model(x, A, B, C):
        return A * x**B + C

    def double_log_model(x, A, B, C, D, E, F):
        return A * np.log(B * x + C) + D * np.log(E * x + F)

    def tan_model(x, A, B, C, D):
        return A * np.tan(B * x - C) + D
    
    # Define knot positions for the spline 
    knots = np.logspace(-3, 0, 5)
    spline_model = SplineModel(prefix='spline_', xknots=knots)
    
    # Model selection dictionary
    models = {
        "linear": linear_model,
        "log": log_model,
        "exp": exp_model,
        "double_log": double_log_model,
        "tan": tan_model,
        "spline": spline_model,
    }
    
    if model_type not in models:
        raise ValueError(f"Unknown model type: {model_type}")
    
    return models[model_type]
    
    return models[model_type]

def factory_model_lmfit(model_type: str):
    
    # select the model, assign parameters and return it
    if model_type == 'linear':
        model = Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=0)
        return model, params 
    
    elif model_type == 'log':
        model = Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        return model, params 
    
    elif model_type == 'exp':
        model = Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0)
        return model, params 
    
    elif model_type == 'double_log':
        model = Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=1, E=1, F=0)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        params['E'].set(min=0.000001)  
        params['F'].set(min=0.000001)
        return model, params 
    
    elif model_type == 'tan':
        model = Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=0)
        return model, params 
    
    if model_type == 'spline':
        model = factory_model_func(model_type)
        # Create parameters 
        params = model.make_params()
        return model, params 
    
    else: raise ValueError('Invalid model selection')

def fit_observed_evolutionary_rate_regressor(df, model_type, weights=None):
    
    x_data = df['evolutionary_rate'] 
    y_data = df['observed_evolutionary_rate']  
    
    # Create the Model
    model, params = factory_model_lmfit(model_type)
    
    # Fit the model to the data
    if weights is None:
        fit_result = model.fit(y_data, params, x=x_data)
    else:
        fit_result = model.fit(y_data, params, x=x_data, weights=weights)
    
    # Print the fit results
    print(fit_result.fit_report())
    return fit_result

def fit_weight(y_data):
    weights = 1/y_data
    return weights

def evaluate_model(model_type, params, x):
    model = factory_model_func(model_type)
    if isinstance(model, SplineModel):
        # Create initial parameters for lmfit
        param_obj = model.make_params()
        for key in param_obj.keys():
            if key in params:
                param_obj[key].set(value=params[key])
        return model.eval(params=param_obj, x=x)
    return model(x, **params)

