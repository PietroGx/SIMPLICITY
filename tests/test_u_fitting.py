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

# Read data from CSV
experiment_name = 'test_param_space_3'
data = om.read_u_e_values(experiment_name)
x_data = data['evolutionary_rate'] 
y_data = data['u']  

# Define the logarithmic model
def log_model(x, A, B, C, D):
    return A * np.log(B * x + C) + D

# Create the Model
model = Model(log_model)

# Set initial parameter guesses (A, B, C, D)
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
plt.ylabel('u')
plt.legend()
plt.title('Logarithmic Fit to Data')
plt.show()
