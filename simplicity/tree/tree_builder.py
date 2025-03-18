#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:51:20 2025

@author: pietro
"""
import anytree
from ete3 import  Tree
# import simplicity.evolution.decoder as decoder
import simplicity.output_manager as om
import simplicity.plots_manager as pm

def infection_tree(seeded_simulation_output_dir):
    
    # import individuals data
    data = om.read_individuals_data(seeded_simulation_output_dir)
    tree = []
    # tree root
    root = anytree.Node('root',label='root',distance=0,
                sequence='',leaf=False,t_infection = 0,
                fitness=0,infection_type ='normal',
                lineage = 'wt')
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
                         
                         sequence         = '',
                         
                         state            = row.state,
                         infection_type   = row.type,
                         fitness          = row.fitness,
                         lineage          = row.IH_lineages[0]
                         
                         ))
        data = data.drop(row.Index)
        
    # build the infection tree
    sorted_data = data.sort_values(by='t_infection')
    
    for row in sorted_data.itertuples():
        #
        if row.parent != None:
            try:
                parent_node = [node for node in tree if node.label == str(row.parent)
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
                        
                        state            = row.state,
                        infection_type   = row.type,
                        fitness          = row.fitness,
                        lineage          = row.IH_lineages[0]
                         ))
            
            # extend parent node
            tree.append(anytree.Node(
                        name           = parent_node.name, 
                        parent         = parent_node,
                        distance       = 0,
                        label          = parent_node.label,
                        leaf           = True,
                        
                        t_infection    = parent_node.t_infection,
                        t_final        = parent_node.t_final,
                        
                        state            = parent_node.state,
                        infection_type   = parent_node.infection_type,
                        fitness          = parent_node.fitness,
                        lineage          = parent_node.lineage
                         ))
            
            
            # update parent node
            parent_node.leaf=False # labels internal nodes (inactive)
            # parent_node.color = 'white'
            parent_node.name =  parent_node.label + '=>'+ str(row.Index) #+'}' '{'+
            
    return(tree)

def phylogenetic_tree(seeded_simulation_output_dir):
    
    # import phylogenetic_data 
    phylogenetic_data = om.read_phylogenetic_data(seeded_simulation_output_dir)
    first_row = phylogenetic_data.iloc[0]
    tree = []
    # tree root
    root = anytree.Node(
                name           = str(first_row.lineage_name), 
                parent         = None,
                lineage        = str(first_row.lineage_name),
                label          = str(first_row.lineage_name),
                distance       = 0,   
                leaf           = True,
                
                time_emergence = 0,
                
                sequence       = None,
                individual     = None,
                
                host_type      = first_row.host_type,
                fitness        = first_row.fitness)
    
    tree.append(root)
    phylogenetic_data.drop(phylogenetic_data.index[0], inplace=True)
    # build the phylogenetic tree
    internal_node_names = []
    for row in phylogenetic_data.itertuples():
            
            parent_node = [node for node in tree if node.label == row.parent
                           and node.leaf == True][0]
           
            time_distance = row.time_emergence - parent_node.time_emergence
            # add new variant
            tree.append(anytree.Node(
                        name           = str(row.lineage_name), 
                        parent         = parent_node,
                        lineage        = str(row.lineage_name),
                        label          = str(row.lineage_name),
                        distance       = time_distance,   
                        leaf           = True,
                        
                        time_emergence = row.time_emergence,
                        
                        sequence       = None,
                        individual     = row.individual,
                        
                        host_type      = row.host_type,
                        fitness        = row.fitness
                         ))
            
            # extend parent node
            tree.append(anytree.Node(
                
                        name           = parent_node.name,
                        parent         = parent_node,
                        lineage        = parent_node.lineage,
                        label          = parent_node.label,
                        distance       = time_distance,   
                        leaf           = True,
                        
                        time_emergence = row.time_emergence,
                        
                        genome         = parent_node.genome,
                        sequence       = None,
                        individual     = parent_node.individual,
                        
                        host_type      = parent_node.host_type,
                        fitness        = parent_node.fitness
                        ))
            
            # update parent node
            parent_node.leaf=False # labels internal nodes (inactive)
            internal_node_name = parent_node.name+ f'_time:{row.time_emergence:.2f}'
            parent_node.name=internal_node_name
            if internal_node_name in internal_node_names:
                internal_node_name += '_' 
                parent_node.name=internal_node_name
            
            internal_node_names.append(internal_node_name)
    return(tree)
   
def get_tree(experiment_name,
             seeded_simulation_output_dir,
             tree_type,
             tree_subtype='binary',
             coloring = 'lineage',
             save_plot=True,
             export_filetype='json',
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
             the parent that can continue to infect more individuals;
    
    compact - infection tree. Each node is an individual 
              connected with all the people they infected
    
    _________________________FOR PHYLOGENETIC TREE__________________________
    binary - binary phylogenetic tree where each internal node is 
             substitution event happening in the simulation
    
    compact - lineages tree where each edge connects parent and 
              offspring lineages
    ________________________________________________________________________
    save_plot : bool
    
    export : str
       file type to export tree (newick or nexus). The default is 'nexus'.
    '''
    # import individuals data
    infection_tree_data = om.read_individuals_data(seeded_simulation_output_dir)
    # import phylogenetic_data 
    phylogenetic_data = om.read_phylogenetic_data(seeded_simulation_output_dir)
    # stop execution if tree cannot be created
    if phylogenetic_data.empty:
        print('No evolution took place during this simulation')
        
        if tree_type == 'phylogenetic':
            print('Cannot create phylogenetic tree, no evolution happened!')
            return
            
        elif coloring == 'lineage':
            print('Cannot color by lineage, no evolution happened!')
            return
    # get lineages colormap
    colormap_df = pm.make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')
   
    # infection tree
    if tree_type == 'infection':
        # infection tree
        tree = infection_tree(seeded_simulation_output_dir)
        root = tree[0]
        
        # Visualize the tree
        if dashplot:
            for pre, fill, node in anytree.RenderTree(root):
                print("%s%s" % (pre, node.label))
        
        if save_plot:
            
            tree_plot_filepath = om.get_tree_plot_filepath(experiment_name,
                                  seeded_simulation_output_dir,
                                  tree_type,
                                  tree_subtype)
            
            pm.plot_infection_tree(root,
                                       infection_tree_data,
                                       tree_subtype,
                                       coloring,
                                       colormap_df,
                                       tree_plot_filepath)
            
            if coloring == 'fitness':
                pm.tree_fitness_legend(infection_tree_data, tree_type, tree_plot_filepath)
       
    # phylogenetic tree            
    elif tree_type == 'phylogenetic':
        tree = phylogenetic_tree(seeded_simulation_output_dir)
        root = tree[0]
        
        # Visualize the tree
        if dashplot:
            for pre, fill, node in anytree.RenderTree(root):
                print("%s%s" % (pre, node.name))
        
        if save_plot:
            tree_plot_filepath = om.get_tree_plot_filepath(experiment_name,
                                  seeded_simulation_output_dir,
                                  tree_type,
                                  tree_subtype)
            
            pm.plot_phylogenetic_tree(root,
                                       phylogenetic_data,
                                       tree_subtype,
                                       coloring,
                                       colormap_df,
                                       tree_plot_filepath)
            if coloring == 'fitness':
                # plot and save legend for fitness color scale
                pm.tree_fitness_legend(phylogenetic_data,tree_type, tree_plot_filepath)
       
    # ------------------- export tree to file -----------------------------        
    om.export_tree(tree,
                    experiment_name,
                    seeded_simulation_output_dir,
                    tree_type,
                    tree_subtype,
                    export_filetype)
    
# Convert anytree node -> ETE node
def build_ete_from_anytree(any_node):
    ete_node = Tree()
    ete_node.name = any_node.name
    # Copy all attributes from any_node
    for key, val in any_node.__dict__.items():
        if key not in ("children", "parent", "name"):
            ete_node.add_features(**{key: val})
    # Recursively add children
    for child in any_node.children:
        ete_child = build_ete_from_anytree(child)
        ete_node.add_child(ete_child)
    return ete_node