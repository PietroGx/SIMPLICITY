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
        params = model.make_params(A=1, B=1, C=1)
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
    # print(fit_result.best_fit)
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
    


def bootstrap_fit_ci(model_type, fit_result, x, y, num_bootstrap=1000, ci_percentile=95, bootstrap_fraction=0.9):
    """
    Compute the confidence intervals using bootstrapping, ensuring correct handling of repeated x values.
    
    Parameters:
    - model_type: Type of model being used.
    - fit_result: Fitted result from the model (used for initial parameters).
    - x: Independent variable (nucleotide_substitution_rate).
    - y: Dependent variable (observed_substitution_rate).
    - num_bootstrap: Number of bootstrap resampling iterations (default 1000).
    - ci_percentile: Percentile for the confidence interval (default 95%).
    - bootstrap_fraction: Fraction of data to use in each bootstrap sample (default 0.9).
    
    Returns:
    - x: Unique x values for plotting.
    - lower_curve: Lower bound of the CI.
    - upper_curve: Upper bound of the CI.
    """
    bootstrap_results = []

    bootstrap_size = int(bootstrap_fraction * len(x))  # Calculate the size of each bootstrap sample

    # Perform bootstrapping: resample the data and fit the model multiple times
    for _ in tqdm(range(num_bootstrap), desc="Bootstrapping", ncols=100):
        # Resample data with replacement (but without replacement inside the bootstrap sample)
        resampled_indices = np.random.choice(len(x), size=bootstrap_size, replace=True)  # Bootstrap with replacement
        x_resampled = x[resampled_indices]
        y_resampled = y[resampled_indices]
        
        # Sort the resampled data by x values
        sorted_indices = np.argsort(x_resampled)  # Get the sorted indices for x_resampled
        x_resampled_sorted = x_resampled[sorted_indices]
        y_resampled_sorted = y_resampled[sorted_indices]

        # Refit the model on the sorted resampled data
        model, params = factory_model_lmfit(model_type)
        fit_result_resampled = model.fit(y_resampled_sorted, params, x=x_resampled_sorted)
        
        # Store the unique x values and the corresponding best_fit y values for this bootstrap iteration
        unique_x_resampled = np.unique(x_resampled_sorted)
        best_fit_unique = fit_result_resampled.best_fit[np.isin(x_resampled_sorted, unique_x_resampled)]

        
        # Store the unique x values and the best_fit corresponding to them
        bootstrap_results.append((unique_x_resampled, best_fit_unique))
    


    # Get unique x values for the plot and CI
    unique_x = np.unique(x)
    print('')
    print('unique_x')
    print(unique_x)
    print('')
    print('bootstrap_results')
    print(bootstrap_results[0])
    
    # Now, we will compute the lower and upper bounds for the confidence intervals
    lower_curve = np.zeros_like(unique_x_resampled)
    upper_curve = np.zeros_like(unique_x_resampled)

    # Calculate the lower and upper percentiles for the CI for each unique x value
    for i, x_val in enumerate(unique_x_resampled):
        # Get the corresponding y-values for this unique x[i] across the bootstrap iterations
        y_values_at_x = [best_fit[np.where(unique_x_resampled == x_val)[0][0]] for _, best_fit in bootstrap_results]

        # Calculate the percentiles for each unique x[i]
        lower_curve[i] = np.percentile(y_values_at_x, (100 - ci_percentile) / 2)
        upper_curve[i] = np.percentile(y_values_at_x, 100 - (100 - ci_percentile) / 2)
   

    print(f"Completed bootstrapping with {num_bootstrap} samples")
    
    return unique_x, lower_curve, upper_curve

