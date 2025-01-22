#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:51:20 2025

@author: pietro
"""
import anytree
import matplotlib
import pandas as pd
import os
import ast
import simplicity.evolution.decoder as decoder
from anytree.exporter import DotExporter
import numpy as np
import collections
import simplicity.tree.newick as nwk
import simplicity.tree.nexus as nx

def fitness_color(fitness,data):
    cmap = matplotlib.pyplot.get_cmap('cool')
    
    # Normalize the value
    norm = matplotlib.colors.Normalize(vmin=data.fitness.min(), vmax=data.fitness.max())
    normalized_value = norm(fitness)
    hexcolor = matplotlib.colors.rgb2hex(cmap(normalized_value))
    
    return hexcolor

def getcolor(state):
    if state == 'infected':
        return 'red'
    if state == 'recovered':
        return 'green'
    if state == 'diagnosed':
        return 'orange'
    if state == 'deceased':
        return 'black'

def infection_tree(output_directory):
    
    # import individuals data
    individuals_data_file_path = os.path.join(output_directory, 'individuals_data.csv')
    data = pd.read_csv(individuals_data_file_path)
    
    tree = []
    # tree root
    root = anytree.Node('root',label='root',distance=0,
                sequence='',leaf=False,
                color = 'white',fitness_color = 'white',infection_type ='normal')
    tree.append(root)
    # patients 0 nodes
    for row in data[data['parent'] =='root'].itertuples():
        
        tree.append(anytree.Node(
                         name           = str(row.Index), 
                         parent         = root,
                         distance       = 0,   
                         label          = str(row.Index),
                         leaf           = True,
                         
                         t_infection      = row.t_infection,
                         t_final          = row.t_final,
                         
                         viral_genomes    = row.viral_genomes,
                         sequence         = '',
                         
                         state            = row.state,
                         infection_type   = row.type,
                         fitness          = row.fitness,
                         IH_virus_fitness = row.IH_virus_fitness,
                         IH_virus_max     = row.IH_virus_max,
                         
                         color = 'black',
                         fitness_color = fitness_color(row.fitness,data)
                         ))
        data = data.drop(row.Index)
        
    # build the infection tree
    sorted_data = data.sort_values(by='t_infection')
    
    for row in sorted_data.itertuples():
        #
        if row.parent != None:
            try:
                parent_node = [node for node in tree if node.name == str(row.parent)
                               and node.leaf == True][0]
                # for node in tree: print(node.name, node.leaf)
                # print('')
            except: 
                # print('pass')
                pass
               
            # add new infection
            tree.append(anytree.Node(
                        name           = str(row.Index), 
                        parent         = parent_node,
                        distance       = row.t_infection - parent_node.t_infection,
                        label          = str(row.Index),
                        leaf           = True,
                        
                        t_infection    = row.t_infection,
                        t_final        = row.t_final,
                        
                        viral_genomes  = row.viral_genomes,
                        
                        state            = row.state,
                        infection_type   = row.type,
                        fitness          = row.fitness,
                        IH_virus_fitness = row.IH_virus_fitness,
                        IH_virus_max     = row.IH_virus_max,
                        
                        color = getcolor(row.state),
                        fitness_color = fitness_color(row.fitness,data)
                         ))
            
            # extend parent node
            tree.append(anytree.Node(
                        name           = parent_node.name, 
                        parent         = parent_node,
                        distance       = 0,
                        label          = str(parent_node.label),
                        leaf           = True,
                        
                        t_infection    = parent_node.t_infection,
                        t_final        = parent_node.t_final,
                        
                        viral_genomes  = parent_node.viral_genomes,
                        
                        state            = parent_node.state,
                        infection_type   = parent_node.infection_type,
                        fitness          = parent_node.fitness,
                        IH_virus_fitness = parent_node.IH_virus_fitness,
                        IH_virus_max     = parent_node.IH_virus_max,
                        
                        color          = getcolor(parent_node.state),
                        fitness_color = fitness_color(parent_node.fitness,data)
                         ))
            
            
            # update parent node
            parent_node.leaf=False # labels internal nodes (inactive)
            # parent_node.color = 'white'
            parent_node.label =  parent_node.name + '=>'+ str(row.Index) #+'}' '{'+
            
    return(tree)

def phylogenetic_tree(output_directory):
    
    # import phylogenetic_data 
    phylogenetic_data_file_path = os.path.join(output_directory, 'phylogenetic_data.csv')
    phylogenetic_data = pd.read_csv(phylogenetic_data_file_path)
    phylogenetic_data['genome'] = phylogenetic_data['genome'].apply(ast.literal_eval)
    if phylogenetic_data.empty:
        print('No evolution took place during this simulation')
        return 
   
    tree = []
    # tree root
    root = anytree.Node(name='wt',label='wt',distance=0,
                leaf=True,
                color = 'white',
                sequence  = None,
                genome = [],
                individual = None,
                host_type      = None,
                fitness        = 1,
                fitness_color = 'black')
    tree.append(root)
        
    # build the phylogenetic tree
    for row in phylogenetic_data.itertuples():
            try:
                parent_node = [node for node in tree if node.name == row.parent
                               and node.leaf == True][0]
                # for node in tree: print(node.name, node.leaf)
                # print('')
            except: 
                print('pass')
                pass
    
            # add new variant
            tree.append(anytree.Node(
                        name           = str(row.name), 
                        parent         = parent_node,
                        label          = str(row.name),
                        distance       = 1,   
                        leaf           = True,
                        
                        time           = row.time,
                        
                        genome         = row.genome,
                        sequence       = None,
                        individual     = row.individual,
                        
                        host_type      = row.host_type,
                        fitness        = row.fitness,
                        fitness_color = fitness_color(row.fitness,phylogenetic_data)
                         ))
            
            # extend parent node
            tree.append(anytree.Node(
                
                        name           = str(parent_node.name), 
                        parent         = parent_node,
                        label          = str(parent_node.label)+ '.m',
                        distance       = 0,   
                        leaf           = True,
                        
                        time           = row.time,
                        
                        genome         = parent_node.genome,
                        sequence       = None,
                        individual     = parent_node.individual,
                        
                        host_type      = parent_node.host_type,
                        fitness        = parent_node.fitness,
                        fitness_color = fitness_color(parent_node.fitness,phylogenetic_data)
                        ))
            
            # update patent node
            parent_node.leaf=False # labels internal nodes (inactive)
            # parent_node.color = 'white'
            # parent_node.label += '.m'
             
    # save all the mutations positions in the leaves nodes
    for node in tree[1:]:
        node.sequence = decoder.decode_genome(node.genome)
            
    return(tree)

def tree_fitness_legend(tree_data, tree_type, output_directory): 
    '''
    Create legend for fitness color in plotted trees 
    
    tree_data: 
        phylogenetic data OR individuals data
    tree: str 
        'infection' or 'phylogenetic' 
    '''
    # Create a dummy invisible image.
    d = np.linspace(0, tree_data.fitness.max(),
                    int(tree_data.fitness.max())).reshape(1, -1)
    d = np.vstack((d, d))

    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 2))
    fig.subplots_adjust(bottom=0.5, top=0.99, left=0.01, right=0.8)

    # Set the extent to cover the range of data
    extent = [0, tree_data.fitness.max(), 0, 1]

    # The imshow plots the dummy data.
    cax = ax.imshow(d, aspect='auto',
                    cmap=matplotlib.pyplot.get_cmap('cool'),
                    extent=extent)

    # Set the ticks at the beginning and end of the color bar
    ax.set_xticks([0, tree_data.fitness.max()])
    # ax.set_xticks(np.arange(0,data.fitness.max(),10))
    # Set the labels "low" and "high" for the ticks
    ax.set_xticklabels(["low", "high"])
    
    # Remove y-ticks and their labels
    ax.set_yticklabels([])
    ax.set_yticks([])

    # Set the title for the plot
    ax.set_title("Relative Fitness")
    legend_path = os.path.join(output_directory,
                               tree_type+"_tree_fitness_legend.png")
    matplotlib.pyplot.savefig(legend_path,
                  dpi=600, bbox_inches='tight')
    matplotlib.pyplot.close()
    
def get_tree(output_directory,
             tree_type,
             tree_subtype='binary',
             save_plot=True,
             export='nexus',
             dashplot=False):
    '''
    Build the infection tree or the phylogenetic tree of the simulation.

    Parameters
    ----------
    tree_type: str
        'infection' or 'phylogenetic'. Selects the type of tree to render.
   
    tree_subtype :  str
    
    __________________________FOR INFECTION TREE___________________________
    binary - binary infection tree where each internal node is an infection 
             event that has as offspring the newly infected individual and 
             the parent that can continue to infect more individuals; the
             tree is colored based on the compartment each individual is in.
    
    compact - infection tree, colored as above, but not binary. Each node 
                is an individual connected with all the people they infected
    
    fitness - infection tree colored with heatmap of individuals fitness value
    _________________________FOR PHYLOGENETIC TREE__________________________
    binary - binary phylogenetic tree where each internal node is 
             substitution event happening in the simulation
    
    non-binary - lineages tree where each edge connects parent and 
                offspring variants
    ________________________________________________________________________
    save_plot : bool
    
    export : str
       file type to export tree (newick or nexus). The default is 'nexus'.
    '''
    img_filename    = tree_type + '_tree_'+tree_subtype+'.png'
    newick_filename = tree_type + '_tree_newick.txt'
    nexus_filename  = tree_type + '_tree_nexus.txt'
    # path to save the plot
    img_path = os.path.join(output_directory,img_filename) 
    
    if tree_type == 'infection':
        # import individuals data
        individuals_data_file_path = os.path.join(output_directory, 'individuals_data.csv')
        tree_data = pd.read_csv(individuals_data_file_path)
        # infection tree
        tree = infection_tree(output_directory)
        root = tree[0]
        
        # Visualize the tree
        if dashplot:
            for pre, fill, node in anytree.RenderTree(root):
                print("%s%s" % (pre, node.label))
        
        if save_plot:
            # path to save the plot
            tree_img_filepath = os.path.join(output_directory,'infection_tree_'+tree_subtype+'.png')
            
            if tree_subtype == 'binary':
               
                def nodenamefunc(node):
                    return node.label 
                
                def nodeattrfunc(node):
                    if node.infection_type == 'normal':
                        return 'color="{}", label="{}"'.format(node.color, 
                                                               node.label)
                    else:
                        return 'color="{}", label="{}", shape=diamond, style=filled'.format(
                            node.color, node.label)
                    
                DotExporter(root,
                            nodeattrfunc=nodeattrfunc,
                            nodenamefunc=nodenamefunc,
                            ).to_picture(tree_img_filepath)
                
            elif tree_subtype == 'compact':
                def nodenamefunc(node):
                    return node.name 
                
                def edgeattrfunc(node, child):
                    if node.name == child.name:
                        return 'color=transparent'
                
                def nodeattrfunc(node):
                    if node.infection_type == 'normal':
                        return 'color="{}", label="{}"'.format(node.color, 
                                                               node.label)
                    else:
                        return 'color="{}", label="{}", shape=diamond, style=filled'.format(
                            node.color, node.label)
            
                DotExporter(root,
                            nodeattrfunc=nodeattrfunc,
                            nodenamefunc=nodenamefunc,
                            edgeattrfunc=edgeattrfunc,
                            ).to_picture(tree_img_filepath)
             
            elif tree_subtype == 'fitness':
                def nodenamefunc(node):
                    return node.name 
                
                def edgeattrfunc(node, child):
                    if node.name == child.name:
                        return 'color=transparent'
                    
                def nodeattrfunc(node):
                    if node.infection_type == 'normal':
                        return 'color="{}", label="{}"'.format(node.fitness_color, 
                                                               node.label)
                    else:
                        return 'color="{}", label="{}", shape=diamond, style=filled'.format(
                            node.fitness_color, node.label)
              
                    
                tree_fitness_legend(tree_data, tree_type, output_directory)
                DotExporter(root,
                            nodeattrfunc=nodeattrfunc,
                            nodenamefunc=nodenamefunc,
                            edgeattrfunc=edgeattrfunc,
                            ).to_picture(tree_img_filepath)
                
    elif tree_type == 'phylogenetic':
        # import phylogenetic_data 
        phylogenetic_data_file_path = os.path.join(output_directory, 'phylogenetic_data.csv')
        phylogenetic_data = pd.read_csv(phylogenetic_data_file_path)
        tree_data = phylogenetic_data
        
        if phylogenetic_data.empty:
            print('No evolution took place during this simulation')
            return 
        
        # phylogenetic tree
        tree = phylogenetic_tree(output_directory)
        root = tree[0]
        
        # Visualize the tree
        if dashplot:
            for pre, fill, node in anytree.RenderTree(root):
                print("%s%s" % (pre, node.name))
        
        if save_plot:
            # define nodeattrfunc for DotExporter to format the graph picture
            def nodeattrfunc(node):
                    return 'color="{}", label="{}"'.format(node.fitness_color, 
                                                           node.name)
            # plot and save legend for fitness color scale
            tree_fitness_legend(tree_data, tree_type, output_directory)
            # format tree for binary tree
            if tree_subtype == 'binary':
                # define nodenamefunc for DotExporter so that it uses labels 
                # instead of names to build the tree picture
                def nodenamefunc(node):
                    return node.label 
                # save tree picture
                anytree.exporter.DotExporter(root,
                            nodeattrfunc=nodeattrfunc,
                            nodenamefunc=nodenamefunc,
                            ).to_picture(img_path)
                
            # format tree for non-binary tree
            elif tree_subtype == 'non-binary':
                
                # not needed but to add transparency to the code I explicitly
                # left the nodenamefunc doing what it already does
                def nodenamefunc(node):
                    return node.name 
                
                # define nodeattrfunc for DotExporter for it to not display edges 
                # that connect collapsed nodes
                def edgeattrfunc(node, child):
                    if node.name == child.name:
                        return 'color=transparent'
                # save tree picture
                anytree.exporter.DotExporter(root,
                            nodeattrfunc=nodeattrfunc,
                            nodenamefunc=nodenamefunc,
                            edgeattrfunc=edgeattrfunc,
                            ).to_picture(img_path)
        
    # tree to ordered dictionary
    exporter = anytree.exporter.DictExporter(dictcls= collections.OrderedDict, attriter=sorted)
    dic = exporter.export(root)
    # ordered dictionary to newick format
    newick_tree = nwk.newick(dic)
        
    if export =='nexus':
        nx.export_nexus(newick_tree, tree, output_directory, nexus_filename) 
       
    if export == 'newick':
        newick_filepath = os.path.join(output_directory, newick_filename) 
        with open(newick_filepath, 'w') as f:
            f.write(newick_tree)
            f.write('\n')