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

def fit_models(experiment_name, model_types, data_type):
    
    parameter = 'evolutionary_rate'
    min_seq_number = 500
    min_sim_lenght = 0
    
    if data_type == 'combined_rate':
        # build the dataframe needed for the fit
        om.write_combined_observed_evolutionary_rate_vs_parameter_csv(experiment_name, 
                                                                    parameter,
                                                                    min_seq_number,
                                                                    min_sim_lenght)
         # import the df needed for the fit
        df = om.read_combined_observed_evolutionary_rate_csv(experiment_name, parameter)
        weights = None
        # select plot function
        plot_fit = pm.plot_combined_observed_evolutionary_rate_fit
        kwargs = {"min_sim_lenght": min_sim_lenght}
        
    elif data_type == 'single_rates':
        # build the dataframe needed for the fit
        om.write_combined_observed_evolutionary_rate_vs_parameter_csv(experiment_name, 
                                                                    parameter, 
                                                                    min_seq_number,
                                                                    min_sim_lenght)
        om.write_observed_evolutionary_rates_vs_parameter_csv(experiment_name, 
                                                                    parameter, 
                                                                    min_seq_number,
                                                                    min_sim_lenght)
         # import the df needed for the fit
        df = om.read_observed_evolutionary_rates_csv(experiment_name, 
                                                     parameter,
                                                     min_sim_lenght)
        weights = None
        # select plot function
        plot_fit = pm.plot_observed_evolutionary_rates_fit
        kwargs = {"min_seq_number": min_seq_number, "min_sim_lenght": min_sim_lenght}
        
    else:
        raise ValueError('invalid data_type, check data_type sintax!')
        
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
            fit_result = er.fit_observed_evolutionary_rate_regressor(df, model_type, weights)
            aic_models[model_type] = fit_result.aic
            print(f'saving plot in {experiment_name}/.')
            plot_fit(experiment_name, fit_result, model_type, **kwargs)
        except Exception as e:
            print(e)
            
        print('')
        print('###############################################################')
        print('')
    
    return aic_models, len(df)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    # define model types for the fit
    model_types = ['linear']
                   # 'log',
                   # 'exp',
                   # 'double_log',
                   # 'tan',
                   # 'spline'] 
    
    # aic_models_combined = fit_models(args.experiment_name, model_types, 'combined_rate')
    aic_models, df_len = fit_models(args.experiment_name, model_types, 'single_rates')
    
    # # print AIC for each model fit
    # sorted_aics_combined = sorted(aic_models_combined.items(), key=lambda item: item[1])
    # print('Comparison of combined rates fits with Akaike Information Criterion')
    # print("Key    Value")
    # for key, value in sorted_aics_combined:
    #     print(f"{key}      {value}")
    
    # print AIC for each model fit
    sorted_aics = sorted(aic_models.items(), key=lambda item: item[1])
    print('Comparison of rates fits with Akaike Information Criterion')
    print("Key    Value")
    for key, value in sorted_aics:
        print(f"{key}      {value}")
    
    print('')
    print(df_len)
    
    
    
if __name__ == "__main__":
    main()