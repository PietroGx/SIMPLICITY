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

class Host:
    '''
    The class defines the intra-host model and the functions needed to solve it.
    '''
    def __init__(self,tau_3=7.5):
        
        # setup intra-host model matrix and calculate matrix exponential
        self.A = self._get_A_matrix(tau_3)
        # self.A_ex = scipy.linalg.expm(self.A)
        
        # possible states of individual in the system
        self.states = np.arange(0,21,1) # state 20 is "healed"
        
    # setup matrix A 
    def _get_A_matrix(self,tau_3):
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
    
    # def _A_t(self,t):
    #     # the method returns A**t
    #     return(scipy.linalg.fractional_matrix_power(self.A_ex,t))
    
    # def _A_t(self, t):
    #     return scipy.linalg.expm(self.A * t)
    
    def _A_t(self, t):
        if not hasattr(self, "_A_cache"):
            self._A_cache = {}
        if t not in self._A_cache:
            self._A_cache[t] = scipy.linalg.expm(self.A * t)
        return self._A_cache[t]

    
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
    
    def simulate_trajectory(self, delta_t, tau_sampler=None):
        """
        Simulate the trajectory of a single individual over time,
        using the intra-host model and rejection sampling.
    
        Parameters
        ----------
        delta_t : float
            Time step between state updates
        tau_sampler : callable
            Function to sample random tau values from [0, 1). Defaults to np.random.rand.
    
        Returns
        -------
        trajectory : list of int
            List of states visited at each time step
        time_points : list of float
            Corresponding time points
        info : dict
            Contains total duration, and durations in each phase category
        """
        if tau_sampler is None:
            tau_sampler = np.random.rand
    
        state = 0
        t = 0
        trajectory = [state]
        time_points = [t]
    
        while state < 20:
            self.probabilities_t(delta_t)
            p_t = self.probabilities[state]
            tau = tau_sampler()
            state = self.update_state(p_t, state, tau)
            t += delta_t
            trajectory.append(state)
            time_points.append(t)
    
        # Analyze trajectory phases
        total_time = time_points[-1]
        detect_start = 5   # starts after pre-detection (n1 = 5)
        infect_start = 5   # starts at pre-symptomatic (n1 + 1)
        infect_end = 18    # end of infectious phase (n1 + n2 + n3)
        
        info = {
            "total_duration": total_time,
            "detectable_duration": sum(
                delta_t for s in trajectory if s >= detect_start
            ),
            "infectious_duration": sum(
                delta_t for s in trajectory if infect_start <= s < infect_end
            ),
            "pre_infectious_duration": sum(
                delta_t for s in trajectory if s < infect_start
            )
        }
    
        return trajectory, time_points, info

import numpy as np
import matplotlib.pyplot as plt

def run_simulations(delta_t, n_runs=1000):
    """
    Run multiple simulations of the intra-host model and aggregate durations.

    Parameters
    ----------
    delta_t : float
        Time step size for simulation
    n_runs : int
        Number of individuals to simulate

    Returns
    -------
    durations : dict
        Lists of durations for each category across simulations
    """
    host = Host()
    durations = {
        "total_duration": [],
        "detectable_duration": [],
        "infectious_duration": [],
        "pre_infectious_duration": []
    }

    for _ in range(n_runs):
        _, _, info = host.simulate_trajectory(delta_t)
        for key in durations:
            durations[key].append(info[key])

    return durations

from tqdm import tqdm


from multiprocessing import Pool

def _simulate_single(args):
    delta_t = args
    host = Host()
    _, _, info = host.simulate_trajectory(delta_t)
    return info

def run_simulations_parallel(delta_t, n_runs=1000):
    durations = {
        "total_duration": [],
        "detectable_duration": [],
        "infectious_duration": [],
        "pre_infectious_duration": []
    }

    with Pool() as pool:
        results = list(tqdm(pool.imap(_simulate_single, [delta_t]*n_runs), total=n_runs))

    for info in results:
        for key in durations:
            durations[key].append(info[key])

    return durations

def plot_duration_distributions(durations, delta_t):
    """
    Plot histograms for duration types from the simulations.

    Parameters
    ----------
    durations : dict
        Dictionary containing lists of durations per category
    delta_t : float
        Time step used in simulation (for title)
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.flatten()
    keys = list(durations.keys())
    titles = ["Total Infection", "Detectability", "Infectiousness", "Pre-infectious"]

    for i, key in enumerate(keys):
        axs[i].hist(durations[key], bins=30, alpha=0.7, edgecolor='black')
        axs[i].set_title(f"{titles[i]} Duration (Δt = {delta_t})")
        axs[i].set_xlabel("Days")
        axs[i].set_ylabel("Count")
        axs[i].grid(True)

    plt.tight_layout()
    plt.show()

def collect_durations_for_deltas(delta_ts, n_runs=1000):
    """
    Run simulations for each delta_t and return all duration results.

    Returns
    -------
    results : dict
        { delta_t: {duration_key: list of values} }
    """
    results = {}
    for dt in delta_ts:
        durations = run_simulations(delta_t=dt, n_runs=n_runs)
        results[dt] = durations
    return results

def plot_durations_from_results(results):
    """
    Plot mean durations for each delta_t using precomputed results.

    Parameters
    ----------
    results : dict
        Output from collect_durations_for_deltas
    """
    duration_keys = list(next(iter(results.values())).keys())
    delta_ts = sorted(results.keys())

    # Prepare data
    means_by_category = {key: [] for key in duration_keys}

    for dt in delta_ts:
        for key in duration_keys:
            means_by_category[key].append(np.mean(results[dt][key]))

    # Plot
    plt.figure(figsize=(10, 6))
    for key, means in means_by_category.items():
        plt.plot(delta_ts, means, marker='o', label=key.replace('_', ' ').title())

    plt.xlabel("Δt (Time step in days)")
    plt.ylabel("Mean Duration (days)")
    plt.title("Effect of Δt on Estimated Durations")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

 
def main():    
    delta_ts = [0.001, 0.01,
                0.1]
    results = {}

    for dt in delta_ts:
        results[dt] = run_simulations_parallel(dt, n_runs=1000)

    # inspect results
    for dt in delta_ts:
        print(f"\nΔt = {dt}")
        for key, vals in results[dt].items():
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")

    # plot aggregate duration trends
    plot_durations_from_results(results)
    
    # Plot distributions per delta_t
    for dt in delta_ts:
        plot_duration_distributions(results[dt], delta_t=dt)
    plot_durations_from_results(results)
    

if __name__ == '__main__':
    import multiprocessing 
    multiprocessing.set_start_method("spawn", force=True)
    main()




        















 


















        