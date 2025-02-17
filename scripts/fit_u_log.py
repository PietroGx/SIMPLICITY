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


def u_estimator(fit_result,e):
    A = fit_result.params['A'].value
    B = fit_result.params['B'].value
    C = fit_result.params['C'].value
    
    return A * np.log(B * e + C)

# Define the model
def log_model(x, A, B, C):
    return A * np.log(B * x + C)

def exp_model(x, A, B, C):
    return A * x**B + C

def double_log_model(x, A, B, C, D, E, F):
    return A * np.log(B * x + C) + D * np.log(E * x + F)

def tan_model(x, A, B, C, D):
    return A * np.tan(B * x - C) + D

def fit_log_u(experiment_name):
    # Read data from CSV
    data = om.read_u_e_values(experiment_name)
    x_data = data['evolutionary_rate'] 
    y_data = data['u']  
    # weights = 1 / y_data
    # Create the Model
    model = Model(double_log_model)
    
    # Set initial parameter guesses 
    params = model.make_params(A=1, B=1, C=1)
    
    # Set boundaries for parameters
    params['A'].set(min=0)  
    params['B'].set(min=0.000001)  
    params['C'].set(min=0.000001)
    params['D'].set(min=0)  
    params['E'].set(min=0.000001)  
    params['F'].set(min=0.000001)
    
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
    pm.plot_u_fit(args.experiment_name,fit_result,scale='lin')
    
if __name__ == "__main__":
    main()