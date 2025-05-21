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

import simplicity.evolution.reference as ref

def decode_genome(single_encoded_genome):
    decoded_genome = ref.get_reference() # fetch reference genome
    if not isinstance(single_encoded_genome, list):
        raise ValueError(f"Expected a list, but got {type(single_encoded_genome)} instead.")
    
    if single_encoded_genome == []:
        return decoded_genome
    else:
        for substitution_coordinate in single_encoded_genome:            
            mutation_location = substitution_coordinate[0]
            mutated_Nbase     = substitution_coordinate[1]
            next_location_rejoin = substitution_coordinate[0]+1
            decoded_genome = decoded_genome[:mutation_location] + mutated_Nbase + decoded_genome[next_location_rejoin:]
        return decoded_genome