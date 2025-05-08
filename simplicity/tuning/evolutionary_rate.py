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
from tqdm import tqdm

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

def fit_observed_substitution_rate_regressor(experiment_name, 
                                             df, model_type, weights=None):
    
    x_data = df['nucleotide_substitution_rate'] 
    y_data = df['observed_substitution_rate']  
    
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

def inverse_log_regressor(OSR, params):
    """
    Computes the inverse of the logarithmic regression function y = A * log(Bx + C)
    
    Parameters:
        OSR:  desired observed substitution rate
        params (dict): A dictionary containing the parameters 'A', 'B', and 'C'.
    
    Returns:
        float or numpy array: The computed nucleotide substitution rate value(s).
    """
    A = params.get('A', 0)
    B = params.get('B', 0)
    C = params.get('C', 0)
    
    # Compute the inverse function
    NSR = (np.exp(OSR / A) - C) / B
    
    return NSR
    

# def bootstrap_fit_ci(model_type, fit_result, x, y, num_bootstrap=1000, ci_percentile=95):
#     """
#     Compute the confidence intervals using bootstrapping.
    
#     Parameters:
#     - model_type: Type of model being used.
#     - fit_result: Fitted result from the model (used for initial parameters).
#     - x: Independent variable (nucleotide_substitution_rate).
#     - y: Dependent variable (observed_substitution_rate).
#     - num_bootstrap: Number of bootstrap resampling iterations (default 1000).
#     - ci_percentile: Percentile for the confidence interval (default 95%).
    
#     Returns:
#     - x: x values for plotting.
#     - lower_curve: Lower bound of the CI.
#     - upper_curve: Upper bound of the CI.
#     """
    
#     bootstrap_results = []

#     # resample the data and fit the model multiple times
#     for _ in tqdm(range(num_bootstrap), desc="Bootstrapping fit for CI", ncols=100):
#         # Resample data with replacement
#         resampled_indices = np.random.choice(range(len(x)), size=len(x), replace=True)
#         x_resampled = x[resampled_indices]
#         y_resampled = y[resampled_indices]
        
#         # Refit the model on the resampled data
#         model, params = factory_model_lmfit(model_type)
#         fit_result_resampled = model.fit(y_resampled, params, x=x_resampled)
        
#         bootstrap_results.append(fit_result_resampled.best_fit)
    
#     bootstrap_results = np.array(bootstrap_results)

#     # Calculate the lower and upper percentiles for the CI
#     lower_curve = np.percentile(bootstrap_results, (100 - ci_percentile) / 2, axis=0)
#     upper_curve = np.percentile(bootstrap_results, 100 - (100 - ci_percentile) / 2, axis=0)

#     return x, lower_curve, upper_curve

def bootstrap_fit_ci(model_type, fit_result, x, y, num_bootstrap=1000, ci_percentile=95, aic_threshold=-12000):
    """
    Compute the confidence intervals using bootstrapping with a progress bar, and reject bad fits.
    
    Parameters:
    - model_type: Type of model being used.
    - fit_result: Fitted result from the model (used for initial parameters).
    - x: Independent variable (nucleotide_substitution_rate).
    - y: Dependent variable (observed_substitution_rate).
    - num_bootstrap: Number of bootstrap resampling iterations (default 1000).
    - ci_percentile: Percentile for the confidence interval (default 95%).
    - aic_threshold: Optional AIC threshold to reject bad fits.
    
    Returns:
    - x: x values for plotting.
    - lower_curve: Lower bound of the CI.
    - upper_curve: Upper bound of the CI.
    """
    bootstrap_results = []
    accepted_aics = []  # List to store the AICs of accepted bootstrap fits
    bootstrap_size = int(0.9*len(x))
    
    for _ in tqdm(range(num_bootstrap), desc="Bootstrapping", ncols=100):
        # Resample data 
        resampled_indices = np.random.choice(len(x), size=bootstrap_size, replace=False)
        x_resampled = x[resampled_indices]
        y_resampled = y[resampled_indices]
        
        # Refit the model on the resampled data
        model, params = factory_model_lmfit(model_type)
        fit_result_resampled = model.fit(y_resampled, params, x=x_resampled)
        
        # If an AIC threshold is set, reject fits with poor AIC values
        if aic_threshold is not None and fit_result_resampled.aic > aic_threshold:
            continue  # Skip this iteration if AIC is too high

        # Store the fitted curve and accepted AIC
        bootstrap_results.append(fit_result_resampled.best_fit)
        accepted_aics.append(fit_result_resampled.aic)
    
    # Convert the list of bootstrap results into a numpy array
    bootstrap_results = np.array(bootstrap_results)

    # Calculate the lower and upper percentiles for the CI
    lower_curve = np.percentile(bootstrap_results, (100 - ci_percentile) / 2, axis=0)
    upper_curve = np.percentile(bootstrap_results, 100 - (100 - ci_percentile) / 2, axis=0)
    print('')
    print(f"Accepted {len(accepted_aics)} bootstrap samples based on AIC threshold")
    
    return x, lower_curve, upper_curve

    
    