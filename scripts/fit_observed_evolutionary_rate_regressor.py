#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:20:25 2025

@author: pietro
"""
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import argparse 

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    parameter = 'evolutionary_rate'
    # build the dataframe needed for the fit
    om.build_combined_observed_evolutionary_rate_vs_parameter_df(args.experiment_name, 
                                                               parameter, 
                                                               min_sim_lenght=0)
    # import the df needed for the fit
    df = om.read_combined_observed_evolutionary_rate_csv(args.experiment_name, parameter)
    
    # define model types for the fit
    model_types = ['linear',
                   'log',
                   'exp',
                   'double_log',
                   'tan',
                   'spline']
    # dic to store the aic of each model for fit report
    aic_models = {}
    # fit the models to the generated data
    for model_type in model_types:
        print('')
        print('###############################################################')
        print('')
        print(f'Fitting {model_type} to evolutionary rate data')
        print('')
        print('###############################################################')
        print('')
        try:
            fit_result = er.fit_observed_evolutionary_rate_regressor(df, model_type, weights=None)
            aic_models[model_type] = fit_result.aic
            print(f'saving plot in {args.experiment_name}/04_Output/')
            pm.plot_observed_evolutionary_rate_fit(args.experiment_name, 
                                                   fit_result, 
                                                   model_type)
        except Exception as e:
            print(e)
            
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