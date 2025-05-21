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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pietro
"""
from simplicity.tree.tree_builder import get_tree
from simplicity.tests.common_test_fixtures import fixture_tree_builder_type
from simplicity.tests.common_test_fixtures import fixture_infection_tree_subtype
from simplicity.tests.common_test_fixtures import fixture_phylogenetic_tree_subtype

output_directory = ''

##### run all tests
if __name__ == "__main__":
    
    # build each tree and each tree type for single simulation run
    for tree_type in fixture_tree_builder_type():
        
        if tree_type == 'infection':  
            for tree_subtype in fixture_infection_tree_subtype():
                get_tree(output_directory,
                         tree_type,
                         tree_subtype,
                         save_plot=True,
                         export='nexus',
                         dashplot=False)
                
        if tree_type == 'phylogenetic':  
            for tree_subtype in fixture_phylogenetic_tree_subtype():
                get_tree(output_directory,
                         tree_type,
                         tree_subtype,
                         save_plot=True,
                         export='nexus',
                         dashplot=False)
           