#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import simplicity.tree.json_tree as jt
import simplicity.tree.tree_builder as tb 
import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm

def generate_and_plot_trees(experiment_name,time_threshold,lineage_threshold):
    SSODs = []
    for simulation_output_dir in dm.get_simulation_output_dirs(experiment_name):
        SSODs += (dm.get_seeded_simulation_output_dirs(simulation_output_dir))
    
    
    for ssod in SSODs:
        seeded_simulation_output_dir = ssod
        final_time = om.read_final_time(seeded_simulation_output_dir)
        if final_time > time_threshold:
            # import phylogenetic_data 
            phylogenetic_data = om.read_phylogenetic_data(seeded_simulation_output_dir)
            # get lineages data
            lineages = phylogenetic_data['lineage_name'].tolist()   
            lineages_number = len(lineages)
            if lineages_number > lineage_threshold:
                print(ssod)
                print(lineages_number)
                # import individuals data
                individuals_data = om.read_individuals_data(seeded_simulation_output_dir)
                individuals_lineages = [lin[0] for lin in individuals_data.IH_virus_names.to_list()]
                # build infection tree
                tb.get_tree(experiment_name,
                             seeded_simulation_output_dir,
                             tree_type='infection',
                             tree_subtype='binary',
                             coloring = 'lineage',
                             save_plot=False,
                             export_filetype='json',
                             dashplot=False)
                # build phylogenetic tree    
                tb.get_tree(experiment_name,
                             seeded_simulation_output_dir,
                             tree_type='phylogenetic',
                             tree_subtype='binary',
                             coloring = 'lineage',
                             save_plot=False,
                             export_filetype='json',
                             dashplot=False)
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
                inf_tree = jt.import_tree(inf_tree_filepath)
                phylo_tree = jt.import_tree(phylo_tree_filepath)
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
                # plot circular trees
                pm.plot_circular_tree(ete_inf_tree,
                                      lineages,
                                      'infection',
                                      individuals_lineages,
                                      ete_inf_tree_filepath)
                pm.plot_circular_tree(ete_phylo_tree,
                                      'phylogenetic',
                                      lineages,
                                      individuals_lineages,
                                      ete_phylo_tree_filepath)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    # Run the script 
    time_threshold = 300
    lineage_threshold =  20
    generate_and_plot_trees(args.experiment_name,time_threshold,lineage_threshold)
    
if __name__ == "__main__":
    main()
    