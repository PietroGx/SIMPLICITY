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

def plot_fit_log_model(experiment_name):
    
    parameter = 'nucleotide_substitution_rate'
    min_seq_number = 30
    min_sim_lenght = 100
    model_type = 'exp'
    weights = None
    
    # build the dataframe needed for the fit
    om.write_combined_OSR_vs_parameter_csv(experiment_name, 
                                            parameter, 
                                            min_seq_number,
                                            min_sim_lenght)
    om.write_OSR_vs_parameter_csv(experiment_name, 
                                    parameter, 
                                    min_seq_number,
                                    min_sim_lenght)
     # import the df needed for the fit
    df = om.read_OSR_vs_parameter_csv(experiment_name, 
                                         parameter,
                                         min_seq_number,
                                         min_sim_lenght)

    # fit log model to the generated data
    
    print('')
    print('###############################################################')
    print('')
    print('Fitting log model to OSR data')
    print('')
    print('###############################################################')
    print('')
    try:
        fit_result = er.fit_observed_substitution_rate_regressor(experiment_name,
                                                                 df, model_type, weights)
        print(f'saving plot in {experiment_name}/.')
        pm.plot_OSR_fit_figure(experiment_name, 
                        fit_result, 
                        model_type,
                        min_seq_number,
                        min_sim_lenght)
    except Exception as e:
        print(e)
        
    print('')
    print('###############################################################')
    print('')

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    plot_fit_log_model(args.experiment_name)
    
if __name__ == "__main__":
    main()