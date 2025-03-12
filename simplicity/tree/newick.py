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





