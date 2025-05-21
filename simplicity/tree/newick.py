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
Created on Tue Jan 18 19:57:13 2022

@author: pietro
"""
import anytree
import collections

def get_newick_str_from_root(node_to_children) -> str:
    '''
    Generate newick string from python AnyTree tree. 

    Parameters
    ----------
    node_to_children : instance of AnyTree.Node

    Returns
    -------
    str
        Tree in Newick format.

    '''
    def newick_render_node(node_to_children) -> str:
        
        name = node_to_children['label']
        distance = node_to_children['distance']
        
        # Leaves
        if 'children' not in node_to_children.keys():
            return F'{name}:{distance}'
        # Nodes
        else:
                       
            children = node_to_children['children']
            children_strings = [newick_render_node(child) for child in children]
            
            children_strings = ",".join(children_strings)
            
            return F'({children_strings}):{distance}'

    newick_string = newick_render_node(node_to_children) + ';'

    return newick_string

def export_newick(root):
    # tree to ordered dictionary
    exporter = anytree.exporter.DictExporter(dictcls= collections.OrderedDict, attriter=sorted)
    dic = exporter.export(root)
    # ordered dictionary to newick format
    newick_tree = get_newick_str_from_root(dic)
    return newick_tree

def write_newick_file(root, newick_filepath):
    newick_tree = export_newick(root)
    with open(newick_filepath, 'w') as f:
        f.write(newick_tree)
        f.write('\n')





