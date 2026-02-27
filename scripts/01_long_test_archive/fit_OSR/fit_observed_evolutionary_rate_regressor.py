# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:20:25 2025

@author: pietro
"""
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.dir_manager as dm
import simplicity.tuning.evolutionary_rate as er
import argparse 
import pandas as pd 
import os

def fit_models(experiment_name, model_types, data_type):
    
    parameter = 'nucleotide_substitution_rate'
    min_seq_number = 30
    min_sim_lenght = 100
    
    if data_type == 'combined_rate':
        # build the dataframe needed for the fit
        om.write_combined_OSR_vs_parameter_csv(experiment_name, 
                                                parameter,
                                                min_seq_number,
                                                min_sim_lenght)
         # import the df needed for the fit
        df = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                                   parameter,
                                                   min_seq_number,
                                                   min_sim_lenght)
        weights = None
        # select plot function
        plot_fit = pm.plot_combined_OSR_fit
        kwargs = {"min_seq_number": min_seq_number, "min_sim_lenght": min_sim_lenght}
        
    elif data_type == 'single_rates':
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
        weights = None
        # select plot function
        plot_fit = pm.plot_OSR_fit
        kwargs = {"min_seq_number": min_seq_number, "min_sim_lenght": min_sim_lenght}
        
    else:
        raise ValueError('invalid data_type, check data_type sintax!')
        
    # dic to store the aic of each model for fit report
    model_metrics = {}
    # fit the models to the generated data
    for model_type in model_types:
        print('')
        print('###############################################################')
        print('')
        print(f'Fitting {model_type} to OSR data')
        print('')
        print('###############################################################')
        print('')
        try:
            fit_result = er.fit_observed_substitution_rate_regressor(experiment_name,
                                                             df, model_type, weights)
                                                             
            metrics = {
                'AIC': fit_result.aic,
                'BIC': fit_result.bic,
                'R2': fit_result.rsquared 
            }
            model_metrics[model_type] = metrics
            
            plot_fit(experiment_name, fit_result, model_type, **kwargs)
            print(f'saved plot in {experiment_name}/.')
        except Exception as e:
            print(e)
            raise ValueError(e)
            
        print('')
        print('###############################################################')
        print('')
    
    return model_metrics, df

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    # define model types for the fit
    model_types = ['log',
                   'exp',
                   'lin',
                   'tan',
                   'spline'] 
    
    # aic_models_combined, _ = fit_models(args.experiment_name, model_types, 'combined_rate')
    model_metrics, df = fit_models(args.experiment_name, model_types, 'single_rates')
    
    # save fit scores
    fit_comparison_df = pd.DataFrame.from_dict(model_metrics, orient='index')
    
    # sort the DataFrame by AIC 
    fit_comparison_df = fit_comparison_df.sort_values(by='AIC')
    output_filename = f"{args.experiment_name}_goodness_of_fit_comparison.csv"
    output_filepath = os.path.join(dm.get_experiment_fit_result_dir(args.experiment_name), output_filename)
    fit_comparison_df.to_csv(output_filepath)
    
    print(f"\n Goodness of fit comparison table saved to: {output_filepath}\n")
    
    # # print AIC for each model fit
    # sorted_aics_combined = sorted(aic_models_combined.items(), key=lambda item: item[1])
    # print('Comparison of combined rates fits with Akaike Information Criterion')
    # print("Key    Value")
    # for key, value in sorted_aics_combined:
    #     print(f"{key}      {value}")
    
    # print AIC for each model fit
    sorted_aics = sorted(model_metrics.items(), key=lambda item: item[1]['AIC'])
    print('Comparison of rates fits with Akaike Information Criterion')
    print("Key    Value")
    for key, value in sorted_aics:
        print(f"{key:<10} {value['AIC']}") # Access the AIC value
    
    print('')
    print('df NSR / OSR:')
    df_counts = df['nucleotide_substitution_rate'].value_counts().reset_index()
    df_counts.columns = ['nucleotide_substitution_rate_value', 'count']
    df_sorted = df_counts.sort_values(by='nucleotide_substitution_rate_value')
    print(df_sorted)
    
if __name__ == "__main__":
    main()