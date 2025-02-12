#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:20:25 2025

@author: pietro
"""

import numpy as np
from lmfit import Model
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import argparse 

import matplotlib
matplotlib.use("Agg")

# Define the model
def log_model(x, A, B, C):
    return A * np.log(B * x + C)

def fit_log_u(experiment_name):
    # Read data from CSV
    data = om.read_u_e_values(experiment_name)
    x_data = data['evolutionary_rate'] 
    y_data = data['u']  
    weights = 1 / y_data
    # Create the Model
    model = Model(log_model)
    
    # Set initial parameter guesses 
    params = model.make_params(A=1, B=1, C=1)
    
    # Set boundaries for parameters
    params['A'].set(min=0)  
    params['B'].set(min=0.000001)  
    params['C'].set(min=0.000001)
    
    # Fit the model to the data
    fit_result = model.fit(y_data, params, x=x_data)#, weights=weights)
    
    # Print the fit results
    print(fit_result.fit_report())
    return fit_result

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit curve to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    # Run the script with the provided parameter
    fit_result = fit_log_u(args.experiment_name)
    pm.plot_u_fit(args.experiment_name,fit_result,scale='semilog')
    pm.plot_u_fit(args.experiment_name,fit_result,scale='loglog')
    
if __name__ == "__main__":
    main()