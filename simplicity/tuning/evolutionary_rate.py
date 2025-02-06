'''
In this file we perform the TempEst linear regression to estimate the observed
evolutionary rate u from the simulated data of a SIMPLICITY run. There are also
the functions to plot E (model evolutionary rate) vs u (observed evolutionary rate)
or vs any other simulation parameter.
'''

import pandas as pd
import sklearn.linear_model 
import os
import glob


def create_joint_sequencing_df(seeeded_simulations_output_directory):
    '''
    seeeded_simulations_output_directory ==> path to subfolder of 
                                             experiment_name/04_Output/
    
    Join all sequencing_data_regression.csv files of different 
    seeded simulation runs (SAME PARAMETERS, different seeds) 
    in a single df and returns it (for tempest regression).
    '''
    target_folder = seeeded_simulations_output_directory
    csv_files = glob.glob(os.path.join(target_folder,'**',
                                       'sequencing_data_regression.csv'),
                                        recursive=True)
    # List to store individual DataFrames
    data_frames = []
    for csv_file in csv_files:
        # Read each CSV file into a DataFrame
        try:
            df = pd.read_csv(csv_file)
            data_frames.append(df)
        except:
            pass
    # Concatenate all DataFrames into one
    try:
        combined_df = pd.concat(data_frames, ignore_index=True)
    except:
        raise ValueError('No sequencing data available to plot!')
    return combined_df, csv_files

def tempest_regression(df):
    '''
    perform TempEst regression on joint dataframe of sequencing data

    Parameters
    ----------
    df : pandas df 
        output of create_joint_sequencing_df.

    Returns
    -------
    u : TYPE
        observed evolutionary rate.
    model : func
        fitted model (sklearn linear regression(.

    '''
    x = df['Sequencing_time'].values.reshape(-1, 1)
    y = df['Distance_from_root'].values
    model = sklearn.linear_model.LinearRegression(fit_intercept=False)
    model.fit(x, y)
    u = model.coef_[0]
    return u, model