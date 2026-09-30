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
    if not isinstance(single_encoded_genome, dict):
        raise ValueError(f"Expected a dict, but got {type(single_encoded_genome)} instead.")

    for position, base in single_encoded_genome.items():
        decoded_genome = (decoded_genome[:position] + base
                          + decoded_genome[position + 1:])
    return decoded_genome