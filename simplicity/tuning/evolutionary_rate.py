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

'''
In this file we perform the TempEst linear regression to estimate the observed
evolutionary rate u from the simulated data of a SIMPLICITY run. There are also
the functions to plot E (model nucleotide substitution rate) vs u (observed substitution rate)
or vs any other simulation parameter.
'''

import sklearn.linear_model 
import numpy as np
import lmfit
import simplicity.output_manager as om
# from tqdm import tqdm

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
        observed substitution rate (OSR)
    model : func
        fitted model (sklearn linear regression).

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

    def tan_model(x, A, B, C, D):
        return A * np.tan(B * x - C) + D
    
    # Define knot positions for the spline 
    knots = np.logspace(-3, 0, 6)
    spline_model = lmfit.models.SplineModel(prefix='spline_', xknots=knots)
    
    # Model selection dictionary
    models = {
        "lin": linear_model,
        "log": log_model,
        "exp": exp_model,
        "tan": tan_model,
        "spline": spline_model,
    }
    
    if model_type not in models:
        raise ValueError(f"Unknown model type: {model_type}")
    
    return models[model_type]

def factory_model_lmfit(model_type: str):
    
    # select the model, assign parameters and return it
    if model_type == 'lin':
        model = lmfit.Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=0.01, B=0)
        return model, params 
    
    elif model_type == 'log':
        model = lmfit.Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=1)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        return model, params 
    
    elif model_type == 'exp':
        model = lmfit.Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=0.9, B=0.7, C=0.5)

        return model, params 
    
    elif model_type == 'tan':
        model = lmfit.Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=0)
        return model, params 
    
    if model_type == 'spline':
        model = factory_model_func(model_type)
        # Create parameters 
        params = model.make_params()
        return model, params 
    
    else: raise ValueError('Invalid model selection')

def fit_observed_substitution_rate_regressor(experiment_name, 
                                             df, model_type, weights=None,
                                             parameter_name='nucleotide_substitution_rate'):
    """
    Fits a model to the observed substitution rate against a varying parameter.
    
    Parameters:
        experiment_name (str): Name of the experiment
        df (pd.DataFrame): Data containing the parameter column and 'observed_substitution_rate'
        model_type (str): Type of model to fit ('lin', 'log', etc.)
        weights: Optional weights for fitting
        parameter_name (str): The column name of the independent variable (x-axis). 
                              Defaults to 'nucleotide_substitution_rate'.
    """
    
    if parameter_name not in df.columns:
        # Fallback for flexibility: try using the first column if specific name not found
        # This is useful if the DF structure varies slightly between steps
        print(f"[Fit] Warning: '{parameter_name}' not found in DF. Using column 0: {df.columns[0]}")
        x_data = df.iloc[:, 0]
    else:
        x_data = df[parameter_name] 
        
    y_data = df['observed_substitution_rate']  
    
    # Create the Model
    model, params = factory_model_lmfit(model_type)
    
    # Fit the model to the data
    try:
        if weights is None:
            fit_result = model.fit(y_data, params, x=x_data)
        else:
            fit_result = model.fit(y_data, params, x=x_data, weights=weights)
            
        # Print the fit results
        print(fit_result.fit_report())
        
        # save fit results
        om.write_fit_results_csv(experiment_name, model_type, fit_result)
        return fit_result
        
    except Exception as e:
        print(f"Fit failed for {model_type}: {e}")
        # Return a dummy object or None so the pipeline can handle it gracefully
        # If using lmfit, we might want to return None and check in the caller
        raise e

def fit_weight(df):
    x_data = df['nucleotide_substitution_rate']
    weights = 1/x_data 
    return weights

def fit_weight_time(df):
    x_data = df['nucleotide_substitution_rate']
    weights = 1/x_data * (df['simulation_final_time']/df['settings_final_time'])
    return weights

def evaluate_model(model_type, params, x):
    model = factory_model_func(model_type)
    if isinstance(model, lmfit.models.SplineModel):
        # Create initial parameters for lmfit
        param_obj = model.make_params()
        for key in param_obj.keys():
            if key in params:
                param_obj[key].set(value=params[key])
        return model.eval(params=param_obj, x=x)
    return model(x, **params)

# ==============================================================================
# INVERSE FUNCTIONS & BRIDGE
# ==============================================================================

def inverse_linear_regressor(OSR, params):
    """
    Computes inverse of: y = A * x + B
    Returns: x = (y - B) / A
    """
    A = params.get('A', 0)
    B = params.get('B', 0)
    
    if A == 0:
        raise ValueError("Parameter A is zero; cannot invert linear function.")
        
    return (OSR - B) / A

def inverse_log_regressor(OSR, params):
    """
    Computes inverse of: y = A * log(B * x + C)
    Returns: x = (exp(y / A) - C) / B
    """
    A = params.get('A', 1)
    B = params.get('B', 1)
    C = params.get('C', 0)
    
    if A == 0:
        raise ValueError("Parameter A is zero; cannot invert log function.")
    if B == 0:
        raise ValueError("Parameter B is zero; cannot invert log function.")
    
    # x = (exp(y/A) - C) / B
    val = (np.exp(OSR / A) - C) / B
    return val
    
def inverse_exp_regressor(OSR, params):
    """
    Computes inverse of: y = A * x**B + C
    Returns: x = ((y - C) / A) ** (1/B)
    """
    A = params.get('A', 1)
    B = params.get('B', 1)
    C = params.get('C', 0)
    
    if A == 0 or B == 0:
        raise ValueError("Parameters A and B must be non-zero.")
    
    # Ensure domain validity if possible
    # if np.any(OSR - C <= 0): ... (checked by caller or numpy raises warning)

    NSR = ((OSR - C) / A) ** (1 / B)
    return NSR

def inverse_tan_regressor(OSR, params):
    """
    Computes inverse of: y = A * tan(B * x - C) + D
    Returns: x = (arctan((y - D) / A) + C) / B
    """
    A = params.get('A', 1)
    B = params.get('B', 1)
    C = params.get('C', 0)
    D = params.get('D', 0)

    if A == 0 or B == 0:
        raise ValueError("Parameters A or B are zero.")

    # (y - D) / A = tan(Bx - C)
    # arctan(...) = Bx - C
    # x = (arctan(...) + C) / B
    val = (np.arctan((OSR - D) / A) + C) / B
    return val

def compute_calibrated_parameter(model_type, fit_result, target_osr):
    """
    The unified bridge function.
    Dispatches the calculation to the correct inverse function.

    Parameters:
        model_type (str): 'lin', 'log', 'exp', 'tan'.
        fit_result (lmfit.model.ModelResult): The result object from fitting.
        target_osr (float): The target Observed Substitution Rate.

    Returns:
        float: The calibrated parameter value (NSR or Factor).
    """
    dispatch_map = {
        'lin': inverse_linear_regressor,
        'log': inverse_log_regressor,
        'exp': inverse_exp_regressor,
        'tan': inverse_tan_regressor,
        # 'spline' is not supported analytically
    }

    if model_type not in dispatch_map:
        raise NotImplementedError(
            f"Inverse function for model '{model_type}' is not implemented."
        )

    # Extract parameters as a dictionary of values
    # fit_result.params is an lmfit Parameters object; values accessed via key lookups or .valuesdict()
    params_dict = fit_result.params.valuesdict()

    return dispatch_map[model_type](target_osr, params_dict)