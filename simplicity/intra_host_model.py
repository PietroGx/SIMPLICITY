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

"""
Intra host transient state model of SARS-COV-2 infection 

@author: Pietro Gerletti
"""
import numpy as np
import scipy.linalg
from multiprocessing import Pool
from tqdm import tqdm
import matplotlib.pyplot as plt
import simplicity.output_manager as om
import os
import pickle
from types import MethodType

class Host:
    '''
    This class defines the intra-host model of SARS-CoV-2 pathogenesis.
    '''

    def __init__(self, tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8.0, update_mode = 'matrix'):
        
        # set up the model matrix
        self.tau_1 = tau_1
        self.tau_2 = tau_2
        self.tau_3 = tau_3
        self.tau_4 = tau_4
        self.A = self._get_A_matrix(tau_1, tau_2, tau_3, tau_4)
        
        # attributes for model solution
        self.n_states = self.A.shape[0]
        self.states = np.arange(0, self.n_states)
        self.update_mode = update_mode # either jump or matrix
        self.use_precomputed_matrix = True
        if self.use_precomputed_matrix and update_mode == 'matrix':
            self.exp_table = self._load_or_precompute_exponentials()
            self.delta_t_not_in_table = []
       
        self.get_A_t = self.factory_get_A_t()

    def get_update_mode(self):
        return self.update_mode
    
    def get_jump_rate(self, state):
        return -self.A[state][state]

    @staticmethod
    def _get_A_matrix(tau_1, tau_2, tau_3, tau_4):
        '''
        Generate the intra-host transition matrix A.
        '''
        # subphases number for each phase 
        n_1 = 5   # pre-detection
        n_2 = 1   # pre-symptomatic
        n_3 = 13  # infectious
        n_4 = 1   # post-infectious
        # last state is recovered
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
    
    def factory_get_A_t(self):
        '''
        Compute or retrieve matrix exponential expm(A * t) from the precomputed table.
        '''
       
        if self.use_precomputed_matrix:
            def get_A_t(self, delta_t):
                key = round(delta_t, 8)
                try:
                    return self.exp_table[key]
                except KeyError:
                    self.delta_t_not_in_table.append(key)
                    return scipy.linalg.expm(self.A * delta_t)
        else:
            def get_A_t(self, delta_t):
                return scipy.linalg.expm(self.A * delta_t)
        
        return MethodType(get_A_t, self)
    
    def _load_or_precompute_exponentials(self):
        file_path = om.get_procomputed_matrix_table_filepath(self.tau_1,self.tau_2,self.tau_3,self.tau_4)
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                # print(f"Loaded matrix exponential table from {file_path}")
                return pickle.load(f)
        else:
            print('Precomputing matrix exponentials...')
            dts = self._generate_dts()
            exp_table = {round(dt, 8): scipy.linalg.expm(self.A * dt) for dt in dts}
            with open(file_path, "wb") as f:
                pickle.dump(exp_table, f)
                print(f"Saved {len(exp_table)} matrix exponentials to {file_path}")
            return exp_table

    @staticmethod
    def _generate_dts():
        dts_small = np.logspace(-5, -2, 200, endpoint=False)
        dts_large = np.linspace(0.01, 10, 100)
        return np.concatenate([dts_small, dts_large])
    
    def get_p_t(self, A_t, state):
        '''
        Compute state probability vector p(t) for given A^t and state.
        '''
        p0 = np.zeros(self.n_states)
        p0[state] = 1
        return np.matmul(A_t, p0)

    @staticmethod
    def update_state(p_t, tau):
        '''
        Sample next state based on rejection sampling.
        '''
        p_cum = np.cumsum(p_t)
        new_state = np.where(tau <= p_cum)[0][0]
        return new_state

    def compute_all_probabilities(self, delta_t):
        '''
        Compute probability vectors for all initial  states.
        '''
        A_t = self.get_A_t(delta_t)
        return [self.get_p_t(A_t, i) for i in self.states]
    
    def data_plot_ih_solution(self,state,time,step):
        '''
        Compute:
            p_inf - probability of being infectious after a time t 
            p_det - probability of being detectable  after a time t 
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
        p_det = []
        p_rec = []
        
        for time_point in t:
            A_t = self.get_A_t(time_point)
            p_i = self.get_p_t(A_t,state)
            p_inf.append(np.sum(p_i[5:19]))
            p_det.append(np.sum(p_i[5:20]))
            p_rec.append(p_i[20])
        
        return p_inf,p_det,p_rec

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
        trajectory = [(t, state)]  # initial state
        
        while state < 20:
            dt = rng.exponential(delta_t) if exponential_dt else delta_t
            probabilities = self.compute_all_probabilities(dt)
            p_t = probabilities[state]
            tau = rng.random()
            new_state = self.update_state(p_t, tau)
            t += dt
            if new_state != state:
                trajectory.append((t, new_state))

            state = new_state
            
        return trajectory

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

    trajectories = []
    with Pool() as pool:
        results = list(tqdm(pool.imap(_simulate_worker, args_list), total=n_runs))
        for trajectory in results:
            trajectories.append(trajectory)
            
    return trajectories

def compute_state_durations(trajectories, max_state=20):
    """
    Correctly compute residence time in each state, accounting for repeated states.
    Returns: dict of state -> list of durations
    """
    state_durations = {s: [] for s in range(max_state + 1)}

    for traj in trajectories:
        if not traj:
            continue
        prev_time, prev_state = traj[0]
        for curr_time, curr_state in traj[1:]:
            dt = curr_time - prev_time
            state_durations[prev_state].append(dt)
            prev_time, prev_state = curr_time, curr_state

    return state_durations

def compute_phase_durations(trajectories):
    """
    Compute durations in major infection phases for each trajectory.
    Returns: list of dicts with phase durations per individual
    """
    results = []
    infect_start=5
    infect_end=18
    detect_start=5
    detect_end = 19
    
    for traj in trajectories:
        if not traj or len(traj) < 2:
            continue

        phase_info = {
            "pre_infectious_duration": 0.0,
            "infectious_duration": 0.0,
            "detectable_duration": 0.0,
            "total_duration": traj[-1][0]  # time of final transition
        }

        for i in range(len(traj) - 1):
            t0, s = traj[i]
            t1, _ = traj[i + 1]
            dt = t1 - t0

            if s < infect_start:
                phase_info["pre_infectious_duration"] += dt
            if infect_start <= s <=infect_end:
                phase_info["infectious_duration"] += dt
            if detect_start <= s <= detect_end:
                phase_info["detectable_duration"] += dt

        results.append(phase_info)

    return results


def plot_state_duration_stats_grid(trajectories_dict, keys, title_prefix):
    """
    Create a 2x2 grid of bar charts showing durations in each intra-host state (0–20)
    for each Δt or λ value.
    
    Parameters
    ----------
    trajectories_dict : dict
        Dict of { Δt or λ : list of trajectories ([(t, s), ...]) }
    keys : list
        List of Δt or λ values to plot
    title_prefix : str
        Title prefix for each subplot (e.g., 'Fixed', 'Exp')
    """
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
    axs = axs.flatten()

    for idx, key in enumerate(keys):
        ax = axs[idx]
        durations = compute_state_durations(trajectories_dict[key])
        states = sorted(durations.keys())

        means = [np.mean(durations[s]) if durations[s] else 0 for s in states]
        stds  = [np.std(durations[s]) if durations[s] else 0 for s in states]

        ax.bar(states, means, yerr=stds, capsize=4, color='skyblue', edgecolor='black')
        ax.set_title(f"{title_prefix} {key}")
        ax.set_xlabel("State")
        ax.set_ylabel("Mean Duration (days)")
        ax.set_xticks(states)
        ax.grid(True, axis='y')

    plt.tight_layout()
    plt.show()


def plot_duration_summary_scatter(fixed_results, exp_results, fixed_dts, exp_lambdas, x_shift=0.015):
    """
    Compare phase durations vs. Δt (fixed) and λ (exp) using scatter plots with error bars.
    Points are shifted slightly for clarity.
    """
    categories = ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]

    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    axs = axs.flatten()

    for i, cat in enumerate(categories):
        ax = axs[i]

        # Shift fixed Δt slightly to the left
        fixed_xs = [dt - x_shift*dt for dt in fixed_dts]
        fixed_means = [np.mean(fixed_results[dt][cat]) for dt in fixed_dts]
        fixed_stds  = [np.std(fixed_results[dt][cat]) for dt in fixed_dts]
        ax.errorbar(fixed_xs, fixed_means, yerr=fixed_stds, fmt='o', label='Fixed Δt', color='blue', capsize=4)

        # Shift exp λ slightly to the right
        exp_xs = [lmbda + x_shift*lmbda for lmbda in exp_lambdas]
        exp_means = [np.mean(exp_results[lmbda][cat]) for lmbda in exp_lambdas]
        exp_stds  = [np.std(exp_results[lmbda][cat]) for lmbda in exp_lambdas]
        ax.errorbar(exp_xs, exp_means, yerr=exp_stds, fmt='s', label='Exp(λ)', color='green', capsize=4)

        ax.set_title(cat.replace('_', ' ').title())
        ax.set_xlabel("Δt / λ")
        ax.set_ylabel("Duration (days)")
        ax.grid(True)
        ax.legend()
        ax.set_xscale('log')

    plt.suptitle("Phase Durations: Fixed Δt vs Exp(λ)", fontsize=14)
    plt.tight_layout()
    plt.show()

def plot_infectious_duration_vs_step(fixed_phase_durations, exp_phase_durations, fixed_dts, exp_lambdas):
    """
    Plot infectious duration vs Δt or λ on log scale.

    Parameters
    ----------
    fixed_phase_durations : dict
        { Δt: list of dicts with phase durations }
    exp_phase_durations : dict
        { λ: list of dicts with phase durations }
    fixed_dts : list of floats
        Fixed step sizes
    exp_lambdas : list of floats
        Exponential step scales
    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(8, 5))

    # Fixed Δt
    x_fixed = fixed_dts
    y_fixed = [np.mean([d["infectious_duration"] for d in fixed_phase_durations[dt]]) for dt in fixed_dts]
    yerr_fixed = [np.std([d["infectious_duration"] for d in fixed_phase_durations[dt]]) for dt in fixed_dts]
    ax.errorbar(x_fixed, y_fixed, yerr=yerr_fixed, fmt='o', color='blue', label='Fixed Δt', capsize=4)

    # Exp(λ)
    x_exp = exp_lambdas
    y_exp = [np.mean([d["infectious_duration"] for d in exp_phase_durations[lmbda]]) for lmbda in exp_lambdas]
    yerr_exp = [np.std([d["infectious_duration"] for d in exp_phase_durations[lmbda]]) for lmbda in exp_lambdas]
    ax.errorbar(x_exp, y_exp, yerr=yerr_exp, fmt='s', color='green', label='Exp(λ)', capsize=4)

    ax.set_xscale('log')
    ax.set_xlabel("Δt / λ (log scale)")
    ax.set_ylabel("Infectious Duration (days)")
    ax.set_title("Infectious Duration vs Δt / λ")
    ax.grid(True, which='both', axis='both')
    ax.legend()
    plt.tight_layout()
    plt.show()


def plot_state_timeline_summary(state_durations_dict, phase_durations_dict, title_prefix=""):
    """
    Visualize average residence times for each state as a timeline-style plot (one per Δt or λ).
    Each state's duration is shown as a horizontal line, placed sequentially on the time axis.
    States are color-coded by the infection phase they belong to.

    Parameters
    ----------
    state_durations_dict : dict
        Dictionary of {Δt or λ: state_durations}, where each state_durations is a dict:
        { state_index -> list of durations }
        (e.g. output of compute_state_durations)
    title_prefix : str
        Prefix to add to each subplot title, e.g. "Fixed" or "Exp"
    """

    # Define which states belong to which biological phase (color-coded)
    def get_phase_color(state):
        if state < 5:
            return "Pre-infectious", "#1f77b4"  # blue
        elif 5 <= state <= 18:
            return "Infectious", "#d62728"      # red
        elif state == 19:
            return "Detectable", "#2ca02c"      # green
        else:
            return "Final", "gray"              # absorbing state 

    fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
    axs = axs.flatten()
    keys = list(state_durations_dict.keys())
    
    # compute global max time
    max_total_time = 0
    for durations in state_durations_dict.values():
        total_time = sum(np.mean(durations[s]) for s in range(20) if durations.get(s))
        max_total_time = max(max_total_time, total_time)

    for idx, key in enumerate(keys):
        ax = axs[idx]
        state_durations = state_durations_dict[key]

        center_time = 0  # center of first state bar

        prev_mean_dur = 0
        for state in range(20):
            durations = state_durations.get(state, [])
            if not durations:
                continue
            mean_dur = np.mean(durations)
            std_dur = np.std(durations)
            phase_label, color = get_phase_color(state)
        
            start = center_time - mean_dur / 2
            end = center_time + mean_dur / 2
        
            ax.hlines(
                y=state,
                xmin=start,
                xmax=end,
                color=color,
                linewidth=3,
                label=phase_label if state in [0, 5, 19] else ""
            )
        
            ax.errorbar(
                x=center_time,
                y=state,
                xerr=std_dur / 2,
                fmt='none',
                ecolor='black',
                capsize=3
            )
        
            # advance center to the next state center
            center_time += 0.5 * (prev_mean_dur + mean_dur)
            prev_mean_dur = mean_dur
            
            # === Compute mean durations for each phase and total ===
            phases = phase_durations_dict.get(key, [])
            if phases:
                pre_mean   = np.mean([p["pre_infectious_duration"] for p in phases])
                inf_mean   = np.mean([p["infectious_duration"] for p in phases])
                det_mean   = np.mean([p["detectable_duration"] for p in phases])
                total_mean = np.mean([p["total_duration"] for p in phases])
            
                # === Plot vertical lines for mean durations ===
                # Pre-infectious duration
                ax.axvline(pre_mean, color='blue', linestyle='--')
                ax.text(pre_mean, 11, f"Pre: {pre_mean:.1f}", color='blue', ha='left', fontsize=8)
            
                # Infectious duration (cumulative from pre)
                ax.axvline(inf_mean, color='red', linestyle='--')
                ax.text(inf_mean, 9, f"Inf: {inf_mean:.1f}", color='red', ha='left', fontsize=8)
            
                # Detectable duration (cumulative from pre + inf)
                ax.axvline(det_mean, color='green', linestyle='--')
                ax.text(det_mean, 7, f"Det: {det_mean:.1f}", color='green', ha='left', fontsize=8)
            
                # Total duration (may differ slightly due to overlap)
                ax.axvline(total_mean, color='black', linestyle=':', linewidth=1.5)
                ax.text(total_mean, 5, f"Total: {total_mean:.1f}", color='black', ha='left', fontsize=8)
                
            
        ax.set_title(f"{title_prefix} Δt = {key}")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("State")
        ax.set_yticks(range(0, 20))
        ax.grid(True, axis='x')
        ax.set_xlim(0, max_total_time)
        ax.legend()

    plt.suptitle("Timeline of Average State Residence Times (Phase Colored)", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    



def main():
    fixed_dts   = [0.0001,0.001, 0.01, 0.1, 1.0]#, 10]
    exp_lambdas = [0.0001,0.001, 0.01, 0.1, 1.0]#, 10]
    # fixed_dts   =  [1,2,4,8,10]
    # exp_lambdas =  [1,2,4,8,10]
    n_runs = 100

    # File paths
    fixed_results_file = 'fixed_results.pkl'
    exp_results_file = 'exp_results.pkl'

    ## -----------------------------
    # Load or run simulations
    # -----------------------------
    if os.path.exists(fixed_results_file):
        print("Loading fixed Δt results...")
        fixed_trajectories = load_results(fixed_results_file)
    else:
        fixed_trajectories = {}
        for dt in fixed_dts:
            print(f"Running fixed Δt = {dt}")
            fixed_trajectories[dt] = run_parallel_simulations(
                delta_t=dt, n_runs=n_runs, exponential_dt=False, base_seed=13
            )
        save_results(fixed_results_file, fixed_trajectories)
    
    if os.path.exists(exp_results_file):
        print("Loading exp(λ) results...")
        exp_trajectories = load_results(exp_results_file)
    else:
        exp_trajectories = {}
        for lmbda in exp_lambdas:
            print(f"Running exp(λ) = {lmbda}")
            exp_trajectories[lmbda] = run_parallel_simulations(
                delta_t=lmbda, n_runs=n_runs, exponential_dt=True, base_seed=42
            )
        save_results(exp_results_file, exp_trajectories)

    # -----------------------------
    # Compute state & phase durations
    # -----------------------------
    fixed_state_durations = {
        dt: compute_state_durations(fixed_trajectories[dt])
        for dt in fixed_dts
    }

    exp_state_durations = {
        lmbda: compute_state_durations(exp_trajectories[lmbda])
        for lmbda in exp_lambdas
    }

    fixed_phase_durations = {
        dt: compute_phase_durations(fixed_trajectories[dt])
        for dt in fixed_dts
    }

    exp_phase_durations = {
        lmbda: compute_phase_durations(exp_trajectories[lmbda])
        for lmbda in exp_lambdas
    }
    
    # -----------------------------
    # Print summary stats
    # -----------------------------
    print("\nFixed Δt Phase Durations:")
    for dt in fixed_dts:
        print(f"Δt = {dt}")
        for key in ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]:
            vals = [d[key] for d in fixed_phase_durations[dt]]
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")

    print("\nExp(λ) Phase Durations:")
    for lmbda in exp_lambdas:
        print(f"λ = {lmbda}")
        for key in ["pre_infectious_duration", "infectious_duration", "detectable_duration", "total_duration"]:
            vals = [d[key] for d in exp_phase_durations[lmbda]]
            print(f"  {key}: mean = {np.mean(vals):.2f}, std = {np.std(vals):.2f}")

    # ----------------------------
    # Plotting
    # ----------------------------

    print("\nPlotting: State Duration Stats (grid)...")
    plot_state_duration_stats_grid(fixed_trajectories, fixed_dts, title_prefix="Fixed")
    plot_state_duration_stats_grid(exp_trajectories, exp_lambdas, title_prefix="Exp")

    print("\nPlotting: Phase Duration Summary (scatter)...")
    fixed_phase_for_plot = {
        dt: {k: [d[k] for d in fixed_phase_durations[dt]] for k in fixed_phase_durations[dt][0]}
        for dt in fixed_dts
    }
    exp_phase_for_plot = {
        lmbda: {k: [d[k] for d in exp_phase_durations[lmbda]] for k in exp_phase_durations[lmbda][0]}
        for lmbda in exp_lambdas
    }
    plot_duration_summary_scatter(fixed_phase_for_plot, exp_phase_for_plot, fixed_dts, exp_lambdas)

    print("\nPlotting: Average State Residence Time Timelines...")
    plot_state_timeline_summary(fixed_state_durations, fixed_phase_durations, title_prefix="Fixed")
    plot_state_timeline_summary(exp_state_durations, exp_phase_durations, title_prefix="Exp")
    
    plot_infectious_duration_vs_step(
    fixed_phase_durations=fixed_phase_durations,
    exp_phase_durations=exp_phase_durations,
    fixed_dts=fixed_dts,
    exp_lambdas=exp_lambdas
    )
if __name__ == '__main__':
    main() 
    
# import cProfile
# import pstats

# if __name__ == '__main__':
#     with cProfile.Profile() as pr:
#         main()  # or whatever your entry function is

#     stats = pstats.Stats(pr)
#     stats.strip_dirs()
#     stats.sort_stats("cumtime").print_stats(20)  # top 20 slowest calls