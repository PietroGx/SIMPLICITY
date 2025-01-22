#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pietro
"""
from simplicity.tree.tree_builder import get_tree
from tests.common_test_fixtures import fixture_tree_builder_type
from tests.common_test_fixtures import fixture_infection_tree_subtype
from tests.common_test_fixtures import fixture_phylogenetic_tree_subtype

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
           