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
Created on Fri Jul 29 13:26:33 2022

@author: pi
"""
import numpy as np
import scipy 

def get_B(tau_1, tau_2, tau_3, tau_4, k_d):
    '''
    Generate the matrix that defines the intra-host model of SARS-CoV-2 
    pathogenesis with diagnosis
    '''
    n_1, n_2, n_3, n_4 = 5, 1, 13, 1
    
    compartments = [[n_1,tau_1],
                    [n_2,tau_2],
                    [n_3,tau_3],
                    [n_4,tau_4]]
    
    dim = np.sum([i[0] for i in compartments])+1 
    B = np.zeros((dim,dim))
    
    start = 0
    comp = 0
    for c in compartments:
        comp = comp + c[0]
        r = c[0]/c[1]
        for i in range(start,comp+1):
            if 5 <= i <=19:
                B[i][i] = - (r+k_d)
            else:
                B[i][i] = -r
            
            if B[i][i-1] == 0: 
                B[i][i-1]= r
            B[i][-1] = 0
        start = start+c[0]
    return B

def get_diagnosis_rate_in_percent(tau_1, tau_2, tau_3, tau_4, k_d):
    '''
    refer to method paper to read the math behind this.
    '''
    y = np.zeros((1,21))
    y[0][5:20] = k_d
    B = get_B(tau_1, tau_2, tau_3, tau_4, k_d)
    B_aug = np.concatenate((B,y))
    z = np.zeros((22,1))
    B_aug = np.concatenate((B_aug,z),axis=1)
    # (e^A)^n = e^(nA) for integer n, so one expm replaces expm + a matrix power
    Bt = scipy.linalg.expm(1000 * B_aug)
    p_t0 = np.zeros((1,22))[0]
    p_t0[0] = 1
    prob_t = np.matmul(Bt,p_t0) 
    return prob_t[-1]    

def get_k_d_from_diagnosis_rate(target_diagnosis_rate_in_percent, tau_1, tau_2, tau_3, tau_4):
    """
    linear search to find k_d value that correspond to desired diagnosis rate 
    in percent (0.00-1)

    :param diagnosis_rate_percent: Target value
    :return: value of k_d corresponding to diagnosis_rate_in_percent
    """
    max_iter=10000
    step = 0.0001
    k_d = 0.0001
    diagnosis_rate_in_percent = get_diagnosis_rate_in_percent(tau_1, tau_2, tau_3, tau_4, k_d)
    iter_count = 0

    # Loop until output is close to the target or max iterations reached
    while abs(diagnosis_rate_in_percent - target_diagnosis_rate_in_percent) > 0.001 and iter_count < max_iter:  # 0.001 is the tolerance
        # Adjust input_value depending on whether output is greater or smaller than the target
        if diagnosis_rate_in_percent < target_diagnosis_rate_in_percent:
            k_d += step
        else:
            k_d -= step
        
        diagnosis_rate_in_percent = get_diagnosis_rate_in_percent(tau_1, tau_2, tau_3, tau_4, k_d)
        iter_count += 1

    return round(k_d,4)






























