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
'''
Here there are the functions we use to calculate the hamming distance for the 
phenotype model. Please consider that they are adapted to the data sctructure we
use to store genomic data (we only store the positions and mutations that are 
                           different from the reference genome.)
'''
from simplicity.evolution import reference as ref
reference = ref.get_reference()


# def hamming(lineage):
#     # compute the hamming distance of a lineage from reference genome
#     distance = 1
    
#     for mutation in lineage:
#         if reference[mutation[0]] != mutation[1]:
#             distance += 1
#     return distance

def hamming(lineage):
    # distance of a lineage from the reference genome
    return sum(1 for pos, base in lineage.items() if reference[pos] != base)

def distributional(lineage, column_dist, delta):
    # centred distributional distance: the expected hamming distance to a
    # randomly drawn circulating genome, minus the consensus's own. Reduces to
    # hamming_iw exactly when every column is pure.
    distance = delta
    for position, base in lineage.items():
        column = column_dist.get(position)
        if column is None:
            distance += 1.0   # no snapshot carries this position: certain disagreement
        else:
            distance += column.get(reference[position], 0.0) - column.get(base, 0.0)
    # >= 0 by construction: the consensus base is each column's argmax, so the
    # sum above cannot fall below -delta. A lineage that IS the consensus hits
    # exactly 0 by near-total cancellation, which in doubles lands ~1e-17 below
    # it -- and that sign reaches rng4.choice as a weight once phi saturates.
    return distance if distance > 0.0 else 0.0

def hamming_iw(lineage,lineage2):
    # distance between two genomes. A position absent from a genome carries the
    # reference base, so an entry storing the reference is not a difference --
    # same rule as hamming, whatever is upstream.
    distance = 0
    for position, base in lineage.items():
        if lineage2.get(position, reference[position]) != base:
            distance += 1
    for position, base in lineage2.items():
        if position not in lineage and base != reference[position]:
            distance += 1
    return distance

    
