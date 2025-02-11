#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:20:25 2025

@author: pietro
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import simplicity.output_manager as om
import argparse 
import os 
import simplicity.config as c 

# Define the model
def log_model(x, A, B, C, D):
    return A * np.log(B * x + C) + D

def fit_log_u(experiment_name):
    # Read data from CSV
    data = om.read_u_e_values(experiment_name)
    x_data = data['evolutionary_rate'] 
    y_data = data['u']  

    # Create the Model
    model = Model(log_model)
    
    # Set initial parameter guesses 
    params = model.make_params(A=3, B=2, C=1, D=2)
    
    # Set boundaries for parameters
    params['A'].set(min=0)  
    params['B'].set(min=0.000001)  
    params['C'].set(min=0.000001)
    
    # Fit the model to the data
    result = model.fit(y_data, params, x=x_data)
    
    # Print the fit results
    print(result.fit_report())
    
    # Plot the original data and the fitted curve
    plt.scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
    plt.plot(x_data, result.best_fit, label='Fitted curve', color='red', linewidth=2)
    plt.xlabel('Evolutionary Rate')
    plt.xlim(0)
    plt.ylim(0)
    plt.ylabel('u')
    plt.legend()
    plt.title('Logarithmic Fit to Data')
    file_path = os.path.join(c.get_experiment_dir(experiment_name),'ue_fitting.png')
    plt.savefig(file_path)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit curve to u data")
    parser.add_argument('experiment_name', type=str, required=True, help="experiment name")
    args = parser.parse_args()
    # Run the script with the provided parameter
    fit_log_u(args)
    
if __name__ == "__main__":
    main()