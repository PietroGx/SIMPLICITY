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

# -*- coding: utf-8 -*-
import simplicity.phenotype.weight as w 
import simplicity.evolution.reference as ref 
import numpy as np

reference = ref.get_reference()

def get_seq_weights(data,t_sim):
    """
    preprocess lineage data to calculate weighted consensus. 
    calculate weight for sequence at time t_sim
    
    Parameters
    ----------
    data : list 
        list of list containing lineages data.
        data contains lineage_name, sequence, n_infected_t, t
    t_sim : float
        time of the simulation at which we evalute the weights.

    Returns
    -------
    seq_for_consensus : list
        [sequence, #_infected_t, w_t(t_sim)].

    """
    # get parameters for w_t(t_sim)
    params = w.w_t_params()

    k_e = params[0]
    k_a = params[1]
    
    seq_for_consensus = []
    # loop over data, compute weight and add entry to list
    for lst in data:
        seq_for_consensus.append([
            lst[0],
            lst[1],
            w.weights(lst[2],t_sim, k_e, k_a, w.t_max)
            ])
   
    return seq_for_consensus

def build_weighted_consensus_matrix(data):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    matrix : TYPE
        DESCRIPTION.
    bases : TYPE
        DESCRIPTION.
    unique_positions : TYPE
        DESCRIPTION.

    '''
    # Extract unique positions
    unique_positions = sorted({num for entry in data for num in entry[0]})
    num_index = {num: i for i, num in enumerate(unique_positions)}
    
    bases = ['A', 'T', 'C', 'G']  
    letter_index = {letter: i for i, letter in enumerate(bases)}
    
    # Initialize the matrix
    matrix = np.zeros((len(bases), len(unique_positions)), dtype=float)
    if not unique_positions:
        return matrix, bases, unique_positions

    # Every genome votes in every column: its own base where it stores one, the
    # reference otherwise. Since a genome differs from the reference at only a
    # handful of the columns, seed every column with the full weight on the
    # reference base and then correct where a genome differs. Identical result,
    # inner loop over mutations rather than over all columns.
    total = sum(n_infected * weight for _, n_infected, weight in data)
    ref_rows = np.fromiter((letter_index[reference[p]] for p in unique_positions),
                           dtype=np.intp, count=len(unique_positions))
    matrix[ref_rows, np.arange(len(unique_positions))] = total
    for seq, n_infected, weight in data:
        v = n_infected * weight
        for num, char in seq.items():
            col = num_index[num]
            matrix[ref_rows[col], col] -= v
            matrix[letter_index[char], col] += v
    
    return matrix, bases, unique_positions

def weighted_consensus(matrix, positions):
    '''
    calculate the weighted consensus sequence between individuals in the population
    '''
    bases = ['A', 'T', 'C', 'G']
    
    consensus = {}
    
    # Iterate through each column in the matrix to find the base with the highest count
    for col_idx, pos in enumerate(positions):
        # Use argmax to find the index of the maximum count in the column
        max_base_idx = np.argmax(matrix[:, col_idx])
        
        # Use the index to get the corresponding letter
        max_base= bases[max_base_idx]
        # if position is diffrent from wt, append it. (we encode sequences as only positions that differ from wt)
        if max_base != reference[pos]:
            # same encoding as a lineage genome: only positions differing from wt
            consensus[pos] = max_base
    
    return consensus

def consensus_distribution(matrix, bases, positions):
    '''Per-column distributions P(p,b), normalised column by column, and

        delta = sum over columns of [ P(p, consensus base) - P(p, reference) ]

    the offset used by the centred distributional distance.
    '''
    P = {}
    delta = 0.0
    for col, pos in enumerate(positions):
        total = matrix[:, col].sum()
        if not total:
            continue
        column = matrix[:, col] / total
        P[pos] = {b: float(column[i]) for i, b in enumerate(bases)}
        delta += float(column.max()) - P[pos][reference[pos]]
    return P, delta

def get_consensus(data,t):
    '''(consensus sequence, column distributions, delta). The argmax distance
    reads the first only; the distributional distance reads the other two.'''
    data = get_seq_weights(data,t)
    matrix, bases, positions= build_weighted_consensus_matrix(data)
    return (weighted_consensus(matrix, positions),
            *consensus_distribution(matrix, bases, positions))


# ## example use 
# data = [
#         [[[10,'A'],[45,'G']],10,10],
#         [[[13,'A']],       5,10],
#         [[[10,'C'],[45,'G']],10,10]
#         ]
# # data contains sequence, n_infected_t, t

# data = get_seq_weights(data,100)

# # Build the ndarray
# matrix, bases, positions= build_weighted_consensus_matrix(data)


# # Call the function with the matrix and unique numbers
# positions_with_max_letter = weighted_consensus(matrix, positions)


# # # To display the result similarly to a pandas DataFrame for visualization:
# # import pandas as pd
# # df = pd.DataFrame(matrix, index=bases, columns=positions)
# # print(df)

# # Display the result
# print(positions_with_max_letter)

# import simplicity.evolution.decoder as dec 

# example_consensus_decoded = dec.decode_genome(positions_with_max_letter)

# print(example_consensus_decoded)


# =============================================================================
# Incremental consensus
# -----------------------------------------------------------------------------
# get_consensus above rescans the whole snapshot on every rebuild, so its cost
# grows with elapsed days x circulating lineages and the snapshot itself grows
# without bound -- millions of genome dicts by day 1095 at N=5000.
#
# The Bateman kernel is a difference of two pure exponentials in dt, so the
# exponent separates:
#
#     e^(k(t_i - t)) = e^(k*t_i) * e^(-k*t)
#     sum_i v_i e^(k(t_i-t)) = e^(-k*t) * sum_i v_i e^(k*t_i)
#                                        \___ depends only on history ___/
#
# Keep that inner sum per (position, base) for each of the two exponentials and
# the matrix is recoverable at any t in O(positions). Exact, not approximate.
#
# Accumulators are held relative to a moving reference time t0, rescaled at
# every evaluation, so e^(k(t_i-t0)) never exceeds e^(k_a * rebuild interval)
# -- about 1.5 -- instead of e^(k_a * 1095) ~ 3e40.
#
# The layout mirrors build_weighted_consensus_matrix exactly rather than
# approximately: every stored genome entry adds its weight to its own base row
# AND subtracts it from the reference row, and the reference row is seeded with
# the running total at evaluation time. A genome holding a reference-valued
# base therefore cancels to zero in both implementations.
# =============================================================================

class ConsensusAccumulator:
    '''Running weighted consensus. add() per snapshot entry, at() to evaluate.

    Replaces keeping the snapshot list: nothing persists the list (only
    consensus_sequences_t is written out), so entries can be folded in and
    dropped.
    '''

    BASES = ['A', 'T', 'C', 'G']

    def __init__(self):
        k_e, k_a = w.w_t_params()
        self.k_e = k_e
        self.k_a = k_a
        self.norm = np.exp(-k_e * w.t_max) - np.exp(-k_a * w.t_max)
        self._base_index = {b: i for i, b in enumerate(self.BASES)}
        self.t0 = 0.0
        # position -> length-4 vector, one per exponential
        self._s_e = {}
        self._s_a = {}
        self._t_e = 0.0
        self._t_a = 0.0

    def add(self, genome, frequency, t):
        '''Fold one snapshot entry in. O(mutations in genome).

        `genome` is {position: base}; `frequency` is the entry's weight before
        the kernel; `t` is the time it was recorded.
        '''
        dt = t - self.t0
        v_e = frequency * np.exp(self.k_e * dt)
        v_a = frequency * np.exp(self.k_a * dt)
        self._t_e += v_e
        self._t_a += v_a
        for position, base in genome.items():
            col_e = self._s_e.get(position)
            if col_e is None:
                col_e = np.zeros(4)
                col_a = np.zeros(4)
                self._s_e[position] = col_e
                self._s_a[position] = col_a
            else:
                col_a = self._s_a[position]
            own = self._base_index[base]
            ref_row = self._base_index[reference[position]]
            col_e[own] += v_e
            col_e[ref_row] -= v_e
            col_a[own] += v_a
            col_a[ref_row] -= v_a

    def at(self, t):
        '''(matrix, bases, positions), the same triple
        build_weighted_consensus_matrix returns. Rescales to t afterwards, so
        the stored exponents stay bounded.'''
        dt = t - self.t0
        f_e = np.exp(-self.k_e * dt)
        f_a = np.exp(-self.k_a * dt)
        total = (f_e * self._t_e - f_a * self._t_a) / self.norm

        positions = sorted(self._s_e)
        matrix = np.zeros((len(self.BASES), len(positions)), dtype=float)
        if not positions:
            self._rescale(f_e, f_a, t)
            return matrix, list(self.BASES), positions

        for col, position in enumerate(positions):
            column = (f_e * self._s_e[position]
                      - f_a * self._s_a[position]) / self.norm
            column[self._base_index[reference[position]]] += total
            matrix[:, col] = column

        self._rescale(f_e, f_a, t)
        return matrix, list(self.BASES), positions

    def _rescale(self, f_e, f_a, t):
        '''Move the reference time to t, so later adds carry small exponents.'''
        for column in self._s_e.values():
            column *= f_e
        for column in self._s_a.values():
            column *= f_a
        self._t_e *= f_e
        self._t_a *= f_a
        self.t0 = t

    def consensus(self, t):
        '''(consensus sequence, column distributions, delta) -- the same triple
        get_consensus returns, so this is a drop-in for it.'''
        matrix, bases, positions = self.at(t)
        return (weighted_consensus(matrix, positions),
                *consensus_distribution(matrix, bases, positions))
