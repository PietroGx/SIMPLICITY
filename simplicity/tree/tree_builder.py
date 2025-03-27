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

def create_infection_node(row, parent_node):
    """Creates a new infection node as a child of the given parent node."""
    return anytree.Node(
        
        name=str(row.Index),
        parent=parent_node,
        distance=row.t_infection - parent_node.t_infection,
        label=str(row.Index),
        leaf=True,
        
        t_infection=row.t_infection,
        t_infectious=row.t_infectious,
        t_not_infectious=row.t_not_infectious,
        
        state=row.state,
        infection_type=row.type,
        fitness=row.fitness,
        lineage=row.IH_lineages[0]
    )

def extend_parent_node(parent_node):
    """Creates a duplicate of the parent node as a leaf child to prolong the branch."""
    return anytree.Node(
        name=parent_node.name,
        parent=parent_node,
        distance=0,
        label=parent_node.label,
        leaf=True,
        
        t_infection=parent_node.t_infection,
        t_infectious=parent_node.t_infectious,
        t_not_infectious=parent_node.t_not_infectious,
        
        state=parent_node.state,
        infection_type=parent_node.infection_type,
        fitness=parent_node.fitness,
        lineage=parent_node.lineage
    )

def infection_tree(seeded_simulation_output_dir):
    # Import individuals data
    data = om.read_individuals_data(seeded_simulation_output_dir)
    tree = []

    # Initialize root node
    root = anytree.Node(
        name='root',
        parent=None,
        distance=0,
        label='root',
        leaf=False,
        
        t_infection=0,
        t_infectious=0,
        t_not_infectious=0,
        
        fitness=0,
        infection_type='normal',
        lineage='wt'
    )
    tree.append(root)

    # Lookup for current leaves by label
    leaf_lookup = {}

    # Handle patient zero(s)
    patient_zero_data = data[data['parent'] == 'root']
    for row in patient_zero_data.itertuples():
        node = create_infection_node(row, root)
        tree.append(node)
        leaf_lookup[node.label] = node

    # Drop patient zero entries from data
    data = data.drop(patient_zero_data.index)

    # Build the infection tree from remaining data
    sorted_data = data.sort_values(by='t_infection')

    for row in sorted_data.itertuples():
        if row.parent is not None:
            parent_label = str(row.parent)
            parent_node = leaf_lookup.get(parent_label)

            if parent_node is None:
                print(f"Warning: parent node '{parent_label}' not found. Skipping.")
                continue

            # Mark the parent node as internal
            parent_node.leaf = False
            parent_node.name = f"{parent_node.label}=>{row.Index}"

            # Create and add new infection node
            child_node = create_infection_node(row, parent_node)
            tree.append(child_node)
            leaf_lookup[child_node.label] = child_node

            # Create and add a "prolonged" leaf version of the parent
            prolonged_parent = extend_parent_node(parent_node)
            tree.append(prolonged_parent)
            leaf_lookup[parent_label] = prolonged_parent

    return tree

def phylogenetic_tree(seeded_simulation_output_dir):
    
    # import phylogenetic_data 
    phylogenetic_data = om.read_phylogenetic_data(seeded_simulation_output_dir)
    first_row = phylogenetic_data.iloc[0]
    tree = []
    # tree root
    root = anytree.Node(
                name           = str(first_row.Lineage_name), 
                parent         = None,
                lineage        = str(first_row.Lineage_name),
                label          = str(first_row.Lineage_name),
                distance       = 0,   
                leaf           = True,
                
                time_emergence = 0,
                
                host_type      = first_row.Host_type,
                )
    
    tree.append(root)
    phylogenetic_data.drop(phylogenetic_data.index[0], inplace=True)
    # build the phylogenetic tree
    internal_node_names = []
    for row in phylogenetic_data.itertuples():
            
            parent_node = [node for node in tree if node.label == row.Lineage_parent
                           and node.leaf == True][0]
           
            time_distance = row.Time_emergence - parent_node.time_emergence
            # add new variant
            tree.append(anytree.Node(
                        name           = row.Lineage_name, 
                        parent         = parent_node,
                        lineage        = row.Lineage_name,
                        label          = row.Lineage_name,
                        distance       = time_distance,   
                        leaf           = True,
                        
                        time_emergence = row.Time_emergence,
                        
                        host_type      = row.Host_type
                         ))
            
            # extend parent node
            tree.append(anytree.Node(
                
                        name           = parent_node.name,
                        parent         = parent_node,
                        lineage        = parent_node.lineage,
                        label          = parent_node.label,
                        distance       = time_distance,   
                        leaf           = True,
                        
                        time_emergence = row.Time_emergence,
                        
                        host_type      = parent_node.host_type
                        ))
            
            # update parent node
            parent_node.leaf=False # labels internal nodes (inactive)
            internal_node_name = parent_node.name+ f'_time:{row.Time_emergence:.2f}'
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