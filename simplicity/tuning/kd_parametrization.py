#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:26:33 2022

@author: pi
"""
import numpy as np
import scipy 

def get_A(tau_3 = 7.5):
        '''
        Generate the matrix that defines the intra-host model of SARS-CoV-2 
        pathogenesis
    
        Returns
        -------
        A : matrix
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
        tau_3 = 7.5  # infectious
        tau_4 = 8    # post-infectious
        
        compartments = [[n_1,tau_1],
                        [n_2,tau_2],
                        [n_3,tau_3],
                        [n_4,tau_4]]
        
        # create empty matrix to be filled
        dim = np.sum([i[0] for i in compartments])+1 # matrix dimensions
        A = np.zeros((dim,dim))
        # fill the matrix with the corresponding compartment parameters
        start = 0
        comp = 0
        for c in compartments:
            comp = comp + c[0]
            r = c[0]/c[1]
            for i in range(start,comp+1):
                A[i][i] = -r
                if A[i][i-1] == 0: 
                    A[i][i-1]= r
                A[i][-1] = 0
            start = start+c[0]
        return A
    
def integrate_p_inf(M,t):
    # returns integral of p(infectious)
    # used to obtain R (reproduction number) 
    y = np.zeros((1,21))
    y[0][5:19] = 1
    M_aug = np.concatenate((M,y))
    z = np.zeros((22,1))
    M_aug = np.concatenate((M_aug,z),axis=1)
    M_ex = scipy.linalg.expm(M_aug)
    Mt = scipy.linalg.fractional_matrix_power(M_ex,t)
    p_t0 = np.zeros((1,22))[0]
    p_t0[0] = 1
    prob_t = np.matmul(Mt,p_t0) 
    return prob_t

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

def detectrate(k_i,tau_3 = 7.5):
    y = np.zeros((1,21))
    y[0][5:20] = k_i
    B = get_B(k_i,tau_3)
    B_aug = np.concatenate((B,y))
    z = np.zeros((22,1))
    B_aug = np.concatenate((B_aug,z),axis=1)
    B_ex = scipy.linalg.expm(B_aug)
    Bt = scipy.linalg.fractional_matrix_power(B_ex,1000)
    p_t0 = np.zeros((1,22))[0]
    p_t0[0] = 1
    prob_t = np.matmul(Bt,p_t0) 
    return prob_t[-1]    

def get_k_i(diagnosis_rate, start, step, tau_3 = 7.5, max_iter=10000):
    """
    Adjusts the input to the function `func` until the output reaches the `target` value.

    :param func: Function to be adjusted
    :param target: Target output value
    :param start: Starting input value
    :param step: Step size to adjust the input
    :param max_iter: Maximum number of iterations to prevent infinite loops
    :return: The input value that makes the function output close to the target, and the actual output
    """
    input_value = start
    output = detectrate(input_value)
    iter_count = 0

    # Loop until output is close to the target or max iterations reached
    while abs(output - diagnosis_rate) > 0.001 and iter_count < max_iter:  # 0.001 is the tolerance
        # Adjust input_value depending on whether output is greater or less than the target
        if output < diagnosis_rate:
            input_value += step
        else:
            input_value -= step
        
        output = detectrate(input_value,tau_3)
        iter_count += 1

    return input_value

def k_d_to_rate(diagnosis_rates):
    '''
    Solve IH augmented B matrix to find ki values corresponging to the desired
    diagnosis rates and saves them to a file.

    Parameters
    ----------
    diagnosis_rates : lst
        rate of diagnosis for which to find k values.

    Returns
    -------
    None.

    '''
    import pandas as pd
    
    dic = {}
    for dr in diagnosis_rates:
        
        k_i = get_k_i(dr, 0.0001, 0.0001)
    
        dic[dr] = k_i
    
    df = pd.DataFrame(list(dic.values()),index=dic.keys(),columns=['k_d value'])
    
    return df



diagnosis_rates = [0.01,0.05,0.1]

k_d_to_rate(diagnosis_rates)
