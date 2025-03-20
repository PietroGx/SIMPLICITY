#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import simplicity.tree.json_tree as jt
import simplicity.tree.tree_builder as tb 
import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm

def generate_and_plot_trees(experiment_name,seeded_simulation_output_dir):
    ''' Generate and plot infection and phylogenetic tree of seeded simulation.
    ''' 
    print(seeded_simulation_output_dir)
    # import individuals data
    individuals_data = om.read_individuals_data(seeded_simulation_output_dir)
    individuals_lineages = [lin[0] for lin in individuals_data.IH_lineages.to_list()]
    colormap_df = pm.make_lineages_colormap(seeded_simulation_output_dir, cmap_name='gist_rainbow')
    print('Importing trees...')
    # get infection tree filepath
    inf_tree_filepath = om.get_tree_file_filepath(experiment_name,
                          seeded_simulation_output_dir,
                          tree_type = 'infection',
                          tree_subtype = 'binary',
                          file_type='json')
    # get phylogenetic tree filepath
    phylo_tree_filepath = om.get_tree_file_filepath(experiment_name,
                          seeded_simulation_output_dir,
                          tree_type = 'phylogenetic',
                          tree_subtype = 'binary',
                          file_type='json')
    # import infection and phylogenetic trees from json
    try:
        inf_tree = jt.import_tree(inf_tree_filepath)
    except:
        print('')
        print('Infection tree file not found, building it from data...')
        # build infection tree
        tb.get_tree(experiment_name,
                     seeded_simulation_output_dir,
                     tree_type='infection',
                     tree_subtype='binary',
                     coloring = 'lineage',
                     save_plot=False,
                     export_filetype='json',
                     dashplot=False)
        inf_tree = jt.import_tree(inf_tree_filepath)
    try:
        phylo_tree = jt.import_tree(phylo_tree_filepath)
    except:
        print('')
        print('Phylogenetic tree file not found, building it from data...')
        # build phylogenetic tree    
        tb.get_tree(experiment_name,
                     seeded_simulation_output_dir,
                     tree_type='phylogenetic',
                     tree_subtype='binary',
                     coloring = 'lineage',
                     save_plot=False,
                     export_filetype='json',
                     dashplot=False)
        phylo_tree = jt.import_tree(phylo_tree_filepath)
    print('')
    print('Converting trees into ete3 format...')
    # convert anytree trees into ete3 trees
    ete_inf_tree = tb.build_ete_from_anytree(inf_tree)
    ete_phylo_tree = tb.build_ete_from_anytree(phylo_tree)
    # get ete (circular) trees filepaths to save them
    ete_inf_tree_filepath =om.get_tree_plot_filepath(experiment_name,
                          seeded_simulation_output_dir,
                          tree_type = 'infection',
                          tree_subtype = 'circular')
    ete_phylo_tree_filepath =om.get_tree_plot_filepath(experiment_name,
                          seeded_simulation_output_dir,
                          tree_type = 'phylogenetic',
                          tree_subtype = 'circular')
    print('')
    print('Plotting circular trees...')
    pm.plot_circular_tree(ete_inf_tree,
                          'infection',
                          colormap_df,
                          individuals_lineages,
                          ete_inf_tree_filepath)
    pm.plot_circular_tree(ete_phylo_tree,
                          'phylogenetic',
                          colormap_df,
                          individuals_lineages,
                          ete_phylo_tree_filepath)

def select_seeded_simulations_to_plot(experiment_name, time_threshold, lineage_number_threshold):
    """
    Filters simulation output directories by:
    1. Final time.
    2. Lineages in the simulation.
    3. Common seed numbers across all filtered directories.

    Parameters:
        experiment_name (str): Name of the experiment.
        time_threshold (float): Minimum final time.
        lineage_number_threshold (int): Minimum lineage number.

    Returns:
        list of str: Flattened list of directories that meet all criteria.
    """
   
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    all_ssods = [dm.get_seeded_simulation_output_dirs(dir_) for dir_ in simulation_output_dirs]

    # Flatten the list and apply filters
    filtered_dirs = []
    seed_pattern = re.compile(r"seed_(\d+)")
    seeds_per_dir = {}

    for ssod_list in all_ssods:
        for ssod in ssod_list:
            final_time = om.read_final_time(ssod)
            if final_time > time_threshold:
                # Read phylogenetic data
                phylogenetic_data = om.read_phylogenetic_data(ssod)
                lineages_number = len(phylogenetic_data['Lineage_name'])

                if lineages_number > lineage_number_threshold:
                    filtered_dirs.append(ssod)
                    
                    # Extract seed number
                    match = seed_pattern.search(ssod)
                    
                    if match:
                        seed = match.group(1)
                        
                        if seed in seeds_per_dir:
                            seeds_per_dir[seed].append(ssod)
                        else:
                            seeds_per_dir[seed] = [ssod]
    
    # Keep only directories with seeds that appear in every filtered list
    common_seeds = [key for key in seeds_per_dir if len(seeds_per_dir[key]) == len(simulation_output_dirs)]
    final_filtered_dirs = [dir_ for seed, dirs in seeds_per_dir.items() if seed in common_seeds for dir_ in dirs]
    kept_seeds = sorted(list(common_seeds))
    return final_filtered_dirs, kept_seeds

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    # Run the script 
    time_threshold = 600
    lineage_number_threshold =  30
    selected_ssods, seeds = select_seeded_simulations_to_plot(args.experiment_name, 
                                                       time_threshold, 
                                                       lineage_number_threshold)
    for selected_ssod in selected_ssods:
        generate_and_plot_trees(args.experiment_name, selected_ssod)
    print('Tree plotting completed.')
    print(f'Selected seeds: {seeds}')
        
if __name__ == "__main__":
    main()
    