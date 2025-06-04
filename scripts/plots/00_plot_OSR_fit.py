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
import simplicity.tuning.evolutionary_rate as er
import argparse 

def plot_fit_exp_model(experiment_name):
    
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
    OSR_single_sim_data = om.read_OSR_vs_parameter_csv(experiment_name, 
                                         parameter,
                                         min_seq_number,
                                         min_sim_lenght)

    # fit exp model to the generated data    
    print('')
    print('###############################################################')
    print('')
    print('Fitting log model to OSR data')
    print('')
    print('###############################################################')
    print('')
    try:
        fit_result = er.fit_observed_substitution_rate_regressor(experiment_name,
                                                                 OSR_single_sim_data, 
                                                                 model_type, weights)
    except Exception as e:
        print(e)
        
    # import combined regression data
    OSR_combined_sim_data = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                                          parameter,
                                                          min_seq_number,
                                                          min_sim_lenght)
    
    # Group by nucleotide_substitution_rate and compute mean and standard deviation for OSR
    OSR_mean_data = om.get_mean_std_OSR(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    print(f'saving plot in {experiment_name}/.')
    pm.plot_OSR_fit_figure(experiment_name, 
                    fit_result, 
                    OSR_single_sim_data,
                    OSR_combined_sim_data,
                    OSR_mean_data,
                    model_type,
                    min_seq_number,
                    min_sim_lenght)
    
    print('')
    print('###############################################################')
    print('')
    

def main():
  
    parser = argparse.ArgumentParser(description="Run script to fit different curves to u data")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    experiment_name= args.experiment_name
    
    plot_fit_exp_model(experiment_name)
    
if __name__ == "__main__":
    main()