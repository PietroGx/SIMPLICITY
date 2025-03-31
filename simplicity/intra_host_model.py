#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intra host transient state model of SARS-COV-2 infection 

@author: Pietro Gerletti
"""
import numpy as np
import scipy.linalg
from multiprocessing import Pool
from tqdm import tqdm
import matplotlib.pyplot as plt

class Host:
    '''
    This class defines the intra-host model of SARS-CoV-2 pathogenesis.
    '''

    def __init__(self, tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8.0):
        self.tau_1 = tau_1
        self.tau_2 = tau_2
        self.tau_3 = tau_3
        self.tau_4 = tau_4
        self.A = self._get_A_matrix(tau_1, tau_2, tau_3, tau_4)
        self._A_cache = {}  # cache for expm(A * t)
        self.states = np.arange(0, 21)

    @staticmethod
    def _get_A_matrix(tau_1, tau_2, tau_3, tau_4):
        '''
        Generate the intra-host transition matrix A.
        '''
        n_1, n_2, n_3, n_4 = 5, 1, 13, 1
        compartments = [[n_1, tau_1], [n_2, tau_2], [n_3, tau_3], [n_4, tau_4]]
        dim = sum([n for n, _ in compartments]) + 1
        A = np.zeros((dim, dim))

        start = 0
        comp = 0
        for n, tau in compartments:
            comp += n
            r = n / tau
            for i in range(start, comp + 1):
                A[i][i] = -r
                if A[i][i - 1] == 0:
                    A[i][i - 1] = r
                A[i][-1] = 0
            start += n
        return A

    def _A_t(self, t, use_cache=True):
        '''
        Compute or retrieve matrix exponential expm(A * t).
        Disable cache when use_cache=False.
        '''
        if not use_cache:
            return scipy.linalg.expm(self.A * t)
        if t not in self._A_cache:
            self._A_cache[t] = scipy.linalg.expm(self.A * t)
        return self._A_cache[t]

    @staticmethod
    def _p_t(A_t, state):
        '''
        Compute state probability vector p(t) for given A^t and state.
        '''
        p0 = np.zeros(21)
        p0[state] = 1
        return np.matmul(A_t, p0)

    @staticmethod
    def update_state(p_t, tau):
        '''
        Sample next state based on rejection sampling.
        '''
        p_cum = np.cumsum(p_t)
        return np.where(tau <= p_cum)[0][0]

    def compute_all_probabilities(self, delta_t, use_cache=True):
        '''
        Compute probability distributions for all states.
        '''
        A_t = self._A_t(delta_t, use_cache=use_cache)
        return [self._p_t(A_t, i) for i in self.states]

    def simulate_trajectory(self, delta_t, rng=None, exponential_dt=False):
        '''
        Simulate disease progression trajectory.

        Parameters
        ----------
        delta_t : float
            Mean time step or fixed step size
        rng : np.random.Generator
            Optional random generator for reproducibility
        exponential_dt : bool
            If True, draw step from Exp(delta_t)

        Returns
        -------
        trajectory, time_points, info : tuple
            Full simulation result
        '''
        if rng is None:
            rng = np.random.default_rng()

        state = 0
        t = 0
        trajectory = [state]
        time_points = [t]

        while state < 20:
            dt = rng.exponential(delta_t) if exponential_dt else delta_t
            use_cache = not exponential_dt
            probabilities = self.compute_all_probabilities(dt, use_cache=use_cache)
            p_t = probabilities[state]
            tau = rng.random()
            state = self.update_state(p_t, tau)
            t += dt
            trajectory.append(state)
            time_points.append(t)

        detect_start, infect_start, infect_end = 5, 5, 18

        info = {
            "pre_infectious_duration": sum(dt for s in trajectory if s < infect_start),
            "infectious_duration": sum(dt for s in trajectory if infect_start <= s < infect_end),
            "detectable_duration": sum(dt for s in trajectory if s >= detect_start and s < 20),
            "total_duration": t
        }
        return trajectory, time_points, info

def _simulate_worker(args):
    delta_t, tau_1, tau_2, tau_3, tau_4, exponential_dt, seed = args
    rng = np.random.default_rng(seed)
    host = Host(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3, tau_4=tau_4)
    _, _, info = host.simulate_trajectory(delta_t, rng=rng, exponential_dt=exponential_dt)
    return info

def run_parallel_simulations(delta_t, n_runs=100, tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8.0, exponential_dt=False, base_seed=None):
    args_list = [
        (delta_t, tau_1, tau_2, tau_3, tau_4, exponential_dt, None if base_seed is None else base_seed + i)
        for i in range(n_runs)
    ]
    durations = {
        "pre_infectious_duration": [],
        "infectious_duration": [],
        "detectable_duration": [],
        "total_duration": []
    }
    with Pool() as pool:
        results = list(tqdm(pool.imap(_simulate_worker, args_list), total=n_runs))
        for info in results:
            for key in durations:
                durations[key].append(info[key])
    return durations

def plot_duration_vs_deltat(results):
    '''
    Plot average durations vs. delta_t for each category.
    '''
    delta_ts = sorted(results.keys())
    categories = list(next(iter(results.values())).keys())
    
    for category in categories:
        means = [np.mean(results[dt][category]) for dt in delta_ts]
        plt.plot(delta_ts, means, marker='o', label=category.replace('_', ' ').title())

    plt.xlabel("Delta t")
    plt.ylabel("Average Duration")
    plt.title("Average Durations vs. Delta t")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    delta_ts = [0.001, 0.01, 0.1, 1.0]
    n_runs = 100
    results = {}

    for dt in delta_ts:
        print(f"\n Running simulations for delta_t = {dt}")
        results[dt] = run_parallel_simulations(delta_t=dt, n_runs=n_runs, exponential_dt=True ,base_seed=13)
        
    for dt in delta_ts:
        print(f"\nΔt = {dt}")
        for key, vals in results[dt].items():
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")
   
    # Print results 
    ordered_keys = ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]
    for dt in delta_ts:
        print(f"Δt = {dt}")
        for key in ordered_keys:
            vals = results[dt][key]
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")
    plot_duration_vs_deltat(results)

if __name__ == '__main__':
    main()




        















 


















        