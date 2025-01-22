#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 01:35:14 2022

@author: pietro
"""
import os

def export_nexus(newick_tree,tree,output_directory,filename):
    '''
    Generate nexus file from AnyTree tree and newick file

    Parameters
    ----------
    newick_tree : newick file
    tree : AnyTree tree
    output_directory : output_directory 
    -------
    Save tree in Nexus format.

    '''
    file_path = os.path.join(output_directory,filename)
    with open(file_path, 'w') as f:
        f.write('#NEXUS')
        f.write('\n')
        f.write('begin data;')
        f.write('\n')
        f.write('  Matrix')
        f.write('\n')
        for node in tree:
            if node.leaf or node.label=='root':
                name = node.name
                # sequence = node.sequence if we want sequence stored in tree file, 
                # we can use the node name (individual index), individuals_data.csv
                # and the decoder to save the full sequence data here
                f.write(F'{name}')# {sequence}')
                f.write('\n')
            
        f.write('    ;')        
        f.write('\n')
        f.write('end;')
        f.write('\n')
        f.write('\n')
        f.write('begin trees;')
        f.write('\n')
        f.write('    tree Tree = ('+newick_tree+');')
        f.write('\n')
        f.write('end;')