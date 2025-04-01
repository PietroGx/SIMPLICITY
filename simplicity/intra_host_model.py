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

    def __init__(self, tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8.0, update_mode = 'jump'):
        self.tau_1 = tau_1
        self.tau_2 = tau_2
        self.tau_3 = tau_3
        self.tau_4 = tau_4
        self.A = self._get_A_matrix(tau_1, tau_2, tau_3, tau_4)
        self._A_cache = {}  # cache for expm(A * t)
        self.states = np.arange(0, 21)
        self.update_mode = update_mode # either jump or matrix
    
    def get_update_mode(self):
        return self.update_mode
    
    def get_jump_rate(self, state):
        return -self.A[state][state]

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

    def _A_t(self, t, use_cache=False):
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

    def compute_all_probabilities(self, delta_t, use_cache=False):
        '''
        Compute probability vectors for all initial  states.
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
            use_cache = exponential_dt
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


import pickle
import os

def save_results(filename, data):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)

def load_results(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def _simulate_worker(args):
    delta_t, tau_1, tau_2, tau_3, tau_4, exponential_dt, seed = args
    rng = np.random.default_rng(seed)
    host = Host(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3, tau_4=tau_4)
    return host.simulate_trajectory(delta_t, rng=rng, exponential_dt=exponential_dt)

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
    trajectories = []
    with Pool() as pool:
        results = list(tqdm(pool.imap(_simulate_worker, args_list), total=n_runs))
        for trajectory, time_points, info in results:
            trajectories.append((time_points, trajectory))
            for key in durations:
                durations[key].append(info[key])
    return durations, trajectories

def compute_state_durations(trajectories, max_state=20):
    """
    For each trajectory (t_points, states), compute how long was spent in each state.
    Returns: dict of state -> list of durations
    """
    state_durations = {s: [] for s in range(max_state + 1)}
    
    for t_points, states in trajectories:
        for i in range(len(states) - 1):
            state = int(states[i])
            duration = t_points[i + 1] - t_points[i]
            state_durations[state].append(duration)
    
    return state_durations


def plot_state_duration_stats_grid(all_durations, keys, title_prefix):
    """
    Create a 2x2 grid of bar charts showing durations in each intra-host state for each Δt or λ value.
    """
    fig, axs = plt.subplots(2, 2, figsize=(14, 8), sharey=True)
    axs = axs.flatten()

    for idx, key in enumerate(keys):
        durations = compute_state_durations(all_durations[key])
        states = sorted(durations.keys())
        means = [np.mean(durations[s]) for s in states]
        stds = [np.std(durations[s]) for s in states]

        ax = axs[idx]
        ax.bar(states, means, yerr=stds, capsize=4, color='skyblue', edgecolor='black')
        ax.set_title(f"{title_prefix} {key}")
        ax.set_xlabel("State")
        ax.set_ylabel("Duration (days)")
        ax.set_xticks(states)
        ax.grid(True, axis='y')

    plt.tight_layout()
    plt.show()

def plot_duration_summary_scatter(fixed_results, exp_results, fixed_dts, exp_lambdas):
    """
    Compare durations vs. Δt (fixed) and λ (exp) using scatter plots with error bars.
    """
    categories = ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.flatten()

    for i, cat in enumerate(categories):
        ax = axs[i]

        fixed_means = [np.mean(fixed_results[dt][cat]) for dt in fixed_dts]
        fixed_stds  = [np.std(fixed_results[dt][cat]) for dt in fixed_dts]
        ax.errorbar(fixed_dts, fixed_means, yerr=fixed_stds, fmt='o', label='Fixed Δt', color='blue', capsize=4)

        exp_means = [np.mean(exp_results[lmbda][cat]) for lmbda in exp_lambdas]
        exp_stds  = [np.std(exp_results[lmbda][cat]) for lmbda in exp_lambdas]
        ax.errorbar(exp_lambdas, exp_means, yerr=exp_stds, fmt='s', label='Exp(λ)', color='green', capsize=4)

        ax.set_title(cat.replace('_', ' ').title())
        ax.set_xlabel("Δt / λ")
        ax.set_ylabel("Duration (days)")
        ax.grid(True)
        ax.legend()

    plt.suptitle("Phase Durations: Fixed Δt vs Exp(λ)", fontsize=14)
    plt.tight_layout()
    plt.show()




def plot_longest_shortest_trajectories(trajectories_dict, durations_dict, title_prefix=""):
    """
    For each Δt or λ, plot:
      - shortest and longest trajectory as stair plots
      - average residence time per state as horizontal bars
    """

    fig, axs = plt.subplots(2, 2, figsize=(14, 8), sharey=True)
    axs = axs.flatten()
    keys = list(trajectories_dict.keys())

    for idx, key in enumerate(keys):
        ax = axs[idx]
        trajectories = trajectories_dict[key]
        durations = durations_dict[key]

        # --- Find shortest and longest trajectory ---
        lengths = [tp[-1] for tp, _ in trajectories]
        shortest_idx = np.argmin(lengths)
        longest_idx = np.argmax(lengths)

        for label, (tp, states), color in zip(
            ["Shortest", "Longest"],
            [trajectories[shortest_idx], trajectories[longest_idx]],
            ["tab:blue", "tab:red"]
        ):
            ax.step(tp, states, where="post", label=label, color=color, lw=1.5)

        # --- Plot average residence time bars ---
        states = sorted(durations.keys())
        means = {s: np.mean(durations[s]) for s in states}
        stds = {s: np.std(durations[s]) for s in states}

        cumulative_time = 0
        for state in states:
            duration = means[state]
            ax.hlines(y=state, xmin=cumulative_time, xmax=cumulative_time + duration, color='black', lw=2)
            ax.errorbar(
                x=cumulative_time + duration / 2,
                y=state,
                xerr=stds[state] / 2,
                fmt='none',
                color='black',
                capsize=3
            )
            cumulative_time += duration

        ax.set_title(f"{title_prefix} Δt = {key}")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("State")
        ax.grid(True)
        ax.legend()

    plt.suptitle("Shortest & Longest Trajectories with Average Residence Times")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


def main():
    fixed_dts = [0.001, 0.01, 0.1, 1.0]
    exp_lambdas = [0.001, 0.01, 0.1, 1.0]
    n_runs = 1000

    # File paths
    fixed_results_file = 'fixed_results.pkl'
    exp_results_file = 'exp_results.pkl'

    # Load or compute fixed Δt results
    if os.path.exists(fixed_results_file):
        print("Loading fixed Δt results from file...")
        fixed_results, fixed_trajectories = load_results(fixed_results_file)
    else:
        fixed_results, fixed_trajectories = {}, {}
        for dt in fixed_dts:
            print(f"\nRunning fixed Δt = {dt}")
            res, traj = run_parallel_simulations(delta_t=dt, n_runs=n_runs, exponential_dt=False, base_seed=13)
            fixed_results[dt] = res
            fixed_trajectories[dt] = traj
        save_results(fixed_results_file, (fixed_results, fixed_trajectories))

    # Load or compute exponential Δt results
    if os.path.exists(exp_results_file):
        print("Loading exponential Δt results from file...")
        exp_results, exp_trajectories = load_results(exp_results_file)
    else:
        exp_results, exp_trajectories = {}, {}
        for lmbda in exp_lambdas:
            print(f"\nRunning exp(λ) Δt ~ Exp({lmbda})")
            res, traj = run_parallel_simulations(delta_t=lmbda, n_runs=n_runs, exponential_dt=True, base_seed=42)
            exp_results[lmbda] = res
            exp_trajectories[lmbda] = traj
        save_results(exp_results_file, (exp_results, exp_trajectories))

    # -------------------------------------
    # 🧮 Compute durations for each setting
    # -------------------------------------
    fixed_durations = [compute_state_durations(fixed_trajectories[dt]) for dt in fixed_dts]
    exp_durations = [compute_state_durations(exp_trajectories[lmbda]) for lmbda in exp_lambdas]

    # ----------------------------
    # 📋 Summary Stats Printout
    # ----------------------------
    print("\nFixed Δt Results:")
    for dt in fixed_dts:
        print(f"Δt = {dt}")
        for key in ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]:
            vals = fixed_results[dt][key]
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")

    print("\nExponential Δt Results:")
    for lmbda in exp_lambdas:
        print(f"λ = {lmbda}")
        for key in ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]:
            vals = exp_results[lmbda][key]
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")

    # ----------------------------
    # 📈 Plotting
    # ----------------------------

    print("\nPlotting: State Duration Stats (grid)...")
    plot_state_duration_stats_grid(fixed_trajectories, fixed_dts, title_prefix="Fixed")
    plot_state_duration_stats_grid(exp_trajectories, exp_lambdas, title_prefix="Exp")

    print("\nPlotting: Duration Summary (scatter)...")
    plot_duration_summary_scatter(fixed_results, exp_results, fixed_dts, exp_lambdas)

    print("\nPlotting: Shortest & Longest Trajectories vs Avg...")
    plot_longest_shortest_trajectories(
        trajectories_dict=fixed_trajectories,
        durations_dict=dict(zip(fixed_dts, fixed_durations)),
        title_prefix="Fixed"
    )   

    plot_longest_shortest_trajectories(
        trajectories_dict=exp_trajectories,
        durations_dict=dict(zip(exp_lambdas, exp_durations)),
        title_prefix="Exp"
    )

if __name__ == '__main__':
    main()






        















 


















        