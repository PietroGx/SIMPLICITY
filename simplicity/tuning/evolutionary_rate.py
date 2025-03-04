'''
In this file we perform the TempEst linear regression to estimate the observed
evolutionary rate u from the simulated data of a SIMPLICITY run. There are also
the functions to plot E (model evolutionary rate) vs u (observed evolutionary rate)
or vs any other simulation parameter.
'''

import sklearn.linear_model 
import numpy as np
import lmfit
import simplicity.output_manager as om

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
        observed evolutionary rate (OER)
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

    def double_log_model(x, A, B, C, D, E, F):
        return A * np.log(B * x + C) + D * np.log(E * x + F)

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
        params = model.make_params(A=1, B=1, C=0)
        return model, params 
    
    elif model_type == 'double_log':
        model = lmfit.Model(factory_model_func(model_type))
        # Set initial parameter guesses 
        params = model.make_params(A=1, B=1, C=0, D=1, E=1, F=0)
        # Set boundaries for parameters
        params['B'].set(min=0.000001)  
        params['C'].set(min=0.000001)
        params['E'].set(min=0.000001)  
        params['F'].set(min=0.000001)
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

def fit_observed_evolutionary_rate_regressor(experiment_name, 
                                             df, model_type, weights=None):
    
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
    # save fit results
    om.write_fit_results_csv(experiment_name, model_type, fit_result)
    return fit_result

def fit_weight(df):
    x_data = df['evolutionary_rate']
    weights = 1/x_data 
    return weights

def fit_weight_time(df):
    x_data = df['evolutionary_rate']
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

def inverse_log_regressor(OER, params):
    """
    Computes the inverse of the logarithmic regression function y = A * log(Bx + C)
    
    Parameters:
        OER:  desired observed evolutionary rate
        params (dict): A dictionary containing the parameters 'A', 'B', and 'C'.
    
    Returns:
        float or numpy array: The computed evolutionary rate value(s).
    """
    A = params.get('A', 0)
    B = params.get('B', 0)
    C = params.get('C', 0)
    
    # Compute the inverse function
    x = (np.exp(OER / A) - C) / B
    
    return x
    
    
    
    