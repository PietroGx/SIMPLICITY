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
    # compute the hamming distance of a lineage from reference genome
    distance = 0
    
    for mutation in lineage:
        if reference[mutation[0]] != mutation[1]:
            distance += 1
            
    #epsilon = 1e-6  # fitness floor
    #distance = max(distance, epsilon)
    
    return distance 

def hamming_iw(lineage,lineage2):
    # compute the hamming distance of a lineage from consensus sequence
    
    # Convert sequences into dictionaries
    dict1 = {position: base for position, base in lineage}
    dict2 = {position: base for position, base in lineage2}
    
    # Get the set of all unique positions in either sequence
    all_positions = set(dict1.keys()) | set(dict2.keys())
    
    # Initialize distance
    distance = 0
    
    # Iterate through all positions and compare bases
    for position in all_positions:
        base1 = dict1.get(position)
        base2 = dict2.get(position)
        
        # If bases are different or only one sequence contains the position, increment distance
        if base1 != base2:
            distance += 1
    
    return distance 

    
