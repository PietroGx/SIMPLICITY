#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intra host transient state model of SARS-COV-2 infection 

@author: Pietro Gerletti
"""
# imports
import numpy as np
import scipy
import scipy.linalg 

np.set_printoptions(precision=3, suppress=True)

def method_wrapper(func):
    def wrapper(self, *args):
        return func(self, *args)
    return wrapper

def _get_A_matrix(tau_3):
    '''
    Generate the matrix that defines the intra-host model of SARS-CoV-2 
    pathogenesis

    Parameters
    ----------
    tau_3 : float
        the model parameter that regulates the transition rate of individuals
        in the infected compartment. The standard value is, like the others,
        taken from the model of Van der Toorn et al. 
        It can be modified to obtain a model of immunocompromised individuals 
        that stay infected for a long time.

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
    tau_2 = 3.91 # pre-symptomatic (infectious)
    #tau_3 = 7.5  # infectious
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

class Host:
    '''
    The class defines the intra-host model and the functions needed to solve it.
    '''
    def __init__(self,tau_3=7.5):
        
        # setup intra-host model matrix and calculate matrix exponential
        self.A = self._get_A_matrix(tau_3)
        self.A_ex = scipy.linalg.expm(self.A)
        
        # possible states of individual in the system
        self.states = np.arange(0,21,1) # state 20 is "healed"
        
    # setup ODE system matrix A 
    @method_wrapper
    def _get_A_matrix(self, tau_3):
        return _get_A_matrix(tau_3)
    
    def _A_t(self,t):
        # the method returns A**t
        return(scipy.linalg.fractional_matrix_power(self.A_ex,t))
    
    def _p_t(self,A_t,state):
        '''
        Calculate p_t taking the system matrix exponential A_t 
        and the individual disease progression state as inputs
        
        Parameters
        ----------
        A_t: matrix
            matrix exponential
        state : int
            individual disease progression state
    
        Returns
        -------
        prob_t: list
        probability of infected individual to be in each of the disease
        progression states
        '''
        #create initial condition vector based on state subject is in
        p_t0 = np.zeros((1,21))[0]
        p_t0[state] = 1
        #solve the system using precomputed matrix
        prob_t = np.matmul(A_t,p_t0) # vector of prob. of being in a state at time t
        # returns p_t
        return prob_t
    
    def update_state(self,p_t,state,tau):
        '''
        Updates the disease state of a modelled individual
        taking p_t and the previous state as input
    
        Parameters
        ----------
        p_t : list
            vector of probabilities for the individual to be in a given disease
            progression state
        state : int
            starting state of the individual
        tau: float
            random variable used for rejection sampling
    
        Returns
        -------
        ns : int
            new disease progression state for the individual modelled
    
        '''
        p_cum = np.cumsum(p_t) # cumulative prob. of being in a state
        ns = np.where(tau<=p_cum)[0][0]
        # returns new state for individual
        return ns
       
    def probabilities_t(self,delta_t):
        '''
        Compute and store p_t for each state on the intra-host model.

        Parameters
        ----------
        delta_t : float
            time interval for which the probabilities are computed.

        Returns
        -------
        Update self.probabilities

        '''
        A_t = self._A_t(delta_t)
        probabilities = []
        # stores p_t for each state
        for i in self.states:
            probabilities.append(self._p_t(A_t,i))
        self.probabilities = probabilities

    def _data_plot_model(self,state,time,step):
        '''
        Compute:
            p_inf - probability of being infectious after a time t 
            p_dia - probability of being diagnosed  after a time t 
            p_rec - probability of being recovered  after a time t 

        Parameters
        ----------
        state : int
            Intra-host model starting state
        time : float
            Time for the intra-host model solution

        Returns
        -------
        Either p_inf, p_dia or p_red. The output is used to plot the 
        intra-host model results.
        '''
        t = np.arange(0,time,step)
    
        p_inf = []
        p_dia = []
        p_rec = []
        
        for time_point in t:
            p_i = self._p_t(self._A_t(time_point),state)
            p_inf.append(np.sum(p_i[5:19]))
            p_dia.append(np.sum(p_i[5:20]))
            p_rec.append(p_i[20])
        
        return p_inf
        

        

        















 


















