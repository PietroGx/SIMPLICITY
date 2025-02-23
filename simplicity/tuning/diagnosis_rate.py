#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:26:33 2022

@author: pi
"""
import numpy as np
import scipy 
import simplicity.dir_manager as dm

def get_B(k_i,tau_3 = 7.5):
        '''
        Generate the matrix that defines the intra-host model of SARS-CoV-2 
        pathogenesis with diagnosis
    
        Returns
        -------
        B : matrix
            21x21 matrix that defines the intra-host model
    
        '''
        # Model parameters
        
        # subphases number for each phase 
        n_1 = 5   # pre-detection
        n_2 = 1   # pre-symptomatic
        n_3 = 13  # infectious
        n_4 = 1   # post-infectious
        # last state is recovered
        
        # parameters for each sub-phase
        tau_1 = 2.86 # pre-detection
        tau_2 = 3.91 # pre-symptomatic
        # tau_3 = 7.5  # infectious
        tau_4 = 8    # post-infectious
        
        compartments = [[n_1,tau_1],
                        [n_2,tau_2],
                        [n_3,tau_3],
                        [n_4,tau_4]]
        
        # create empty matrix to be filled
        dim = np.sum([i[0] for i in compartments])+1 # matrix dimensions
        B = np.zeros((dim,dim))
        # fill the matrix with the corresponding compartment parameters
        start = 0
        comp = 0
        for c in compartments:
            comp = comp + c[0]
            r = c[0]/c[1]
            for i in range(start,comp+1):
                if 5 <= i <=19:
                    B[i][i] = - (r+k_i)
                else:
                    B[i][i] = -r
                
                if B[i][i-1] == 0: 
                    B[i][i-1]= r
                B[i][-1] = 0
            start = start+c[0]
        return B

def get_diagnosis_rate_in_percent(k_d,tau_3 = 7.5):
    '''
    refer to method paper to read the math behind this.
    '''
    y = np.zeros((1,21))
    y[0][5:20] = k_d
    B = get_B(k_d,tau_3)
    B_aug = np.concatenate((B,y))
    z = np.zeros((22,1))
    B_aug = np.concatenate((B_aug,z),axis=1)
    B_ex = scipy.linalg.expm(B_aug)
    Bt = scipy.linalg.fractional_matrix_power(B_ex,1000)
    p_t0 = np.zeros((1,22))[0]
    p_t0[0] = 1
    prob_t = np.matmul(Bt,p_t0) 
    return prob_t[-1]    

def get_k_d_from_diagnosis_rate(target_diagnosis_rate_in_percent, tau_3 = 7.5):
    """
    linear search to find k_d value that correspond to desired diagnosis rate 
    in percent (0.00-1)

    :param diagnosis_rate_percent: Target value
    :return: value of k_d corresponding to diagnosis_rate_in_percent
    """
    max_iter=10000
    step = 0.0001
    k_d = 0.0001
    diagnosis_rate_in_percent = get_diagnosis_rate_in_percent(k_d)
    iter_count = 0

    # Loop until output is close to the target or max iterations reached
    while abs(diagnosis_rate_in_percent - target_diagnosis_rate_in_percent) > 0.001 and iter_count < max_iter:  # 0.001 is the tolerance
        # Adjust input_value depending on whether output is greater or smaller than the target
        if diagnosis_rate_in_percent < target_diagnosis_rate_in_percent:
            k_d += step
        else:
            k_d -= step
        
        diagnosis_rate_in_percent = get_diagnosis_rate_in_percent(k_d,tau_3)
        iter_count += 1

    return round(k_d,4)

def diagnosis_rate_table(diagnosis_rates):
    '''
    Solve IH augmented B matrix to find k_d values corresponging to the desired
    diagnosis rates (in decimal % of diagnosed individuals)

    Parameters
    ----------
    diagnosis_rates : lst
         human readable diagnosis rates (0.00-1) for which to find k_d values.

    Returns
    -------
    None.

    '''
    import pandas as pd
    
    dic = {}
    for diagnosis_rate in diagnosis_rates:
        
        k_d = get_k_d_from_diagnosis_rate(diagnosis_rate)
    
        dic[diagnosis_rate] = k_d
    
    df = pd.DataFrame(list(dic.values()),index=dic.keys(),columns=['k_d value'])
    
    return df

def get_effective_diagnosis_rate(simulation_output_dir):
    import os 
    import pandas as pd 
    import simplicity.dir_manager as dm
    effective_diagnosis_rates = []
    for seeded_simulation_output_dir in dm.get_seeded_simulation_output_dirs(simulation_output_dir):
        simulation_trajectory = os.path.join(seeded_simulation_output_dir, 'simulation_trajectory.csv')
    
        trajectory_data = pd.read_csv(simulation_trajectory)
        # Get the last entry 
        last_entry = trajectory_data.iloc[-1]
    
        diagnosed = last_entry['diagnosed']
        deceased = last_entry['deceased']
        recovered = last_entry['recovered']
    
        # Calculate the effective diagnosis rate 
        total = recovered + diagnosed
        if total > 0:
            effective_rate = (diagnosed + deceased) / total
            effective_diagnosis_rates.append(effective_rate)
    
    return np.mean(effective_diagnosis_rates), np.std(effective_diagnosis_rates)

def get_diagnosis_rates(simulation_output_dir):
    
    effective_rate = get_effective_diagnosis_rate(simulation_output_dir)
    theoretical_rate = dm.get_simulation_parameters_of_simulation_output_dir(
                       simulation_output_dir)['diagnosis_rate']
    
    return theoretical_rate, effective_rate



























