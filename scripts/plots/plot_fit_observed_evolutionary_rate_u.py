#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:20:25 2025

@author: pietro
"""
import simplicity.plots_manager as pm 
import simplicity.tuning.evolutionary_rate as er
import argparse 

# def u_estimator(fit_result,e):
#     A = fit_result.params['A'].value
#     B = fit_result.params['B'].value
#     C = fit_result.params['C'].value
    
#     return A * np.log(B * e + C)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    model_types = ['linear',
                   'log',
                   'exp',
                   'double_log',
                   'tan',
                   'spline']
    
    aic_models = {}
    for model_type in model_types:
        print('')
        print('###############################################################')
        print('')
        print(f'Fitting {model_type} to evolutionary rate data')
        print('')
        print('###############################################################')
        print('')
        
        try:
            fit_result = er.fit_observed_evolutionary_rate(args.experiment_name, model_type)
            aic_models[model_type] = fit_result.aic
            print(f'saving plot in {args.experiment_name}/04_Output/')
            pm.plot_observed_evolutionary_rate_fit(args.experiment_name, 
                                                   fit_result, 
                                                   model_type)
        except Exception as e:
            print(e)
            # print(f'model {model_type} could not be fit, too few data points!')
            
        print('')
        print('###############################################################')
        print('')
        
    # print AIC for each model fit
    sorted_aics = sorted(aic_models.items(), key=lambda item: item[1])
    print('Comparison of fits with Akaike Information Criterion')
    print("Key    Value")
    for key, value in sorted_aics:
        print(f"{key}      {value}")
    
if __name__ == "__main__":
    main()