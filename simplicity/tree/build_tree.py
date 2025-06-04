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
import simplicity.tree.json_tree as jt
import simplicity.tree.tree_builder as tb
import simplicity.plots_manager as pm
import simplicity.output_manager as om
    
def draw_tree(ssod, experiment_name, tree_type):
    """
    Builds and plots a circular tree of the given type ('phylogenetic' or 'infection')
    for the specified seeded simulation output directory.
    """

    assert tree_type in {"phylogenetic", "infection"}, f"Unsupported tree_type: {tree_type}"

    # Load individual data for lineage assignment
    individuals_data = om.read_individuals_data(ssod)
    individuals_lineages = [lin[0] for lin in individuals_data.IH_lineages.to_list()]
    
    # Build colormap
    colormap_df = pm.make_lineages_colormap(ssod, cmap_name='gist_rainbow')

    # Locate or generate tree file
    tree_file_path = om.get_tree_file_filepath(
        experiment_name,
        ssod,
        tree_type=tree_type,
        tree_subtype='binary',
        file_type='json'
    )

    try:
        tree = jt.import_tree(tree_file_path)
    except FileNotFoundError:
        print(f"[INFO] {tree_type.capitalize()} tree not found — building for {ssod}")
        tb.get_tree(
            experiment_name,
            ssod,
            tree_type=tree_type,
            tree_subtype='binary',
            coloring='lineage',
            save_plot=False,
            export_filetype='json',
            dashplot=False
        )
        tree = jt.import_tree(tree_file_path)

    ete_tree = tb.build_ete_from_anytree(tree)

    # Determine path for circular tree image
    tree_image_path = om.get_tree_plot_filepath(
        experiment_name,
        ssod,
        tree_type=tree_type,
        tree_subtype='circular'
    )

    # Plot and save the tree image
    pm.plot_circular_tree(
        ete_tree,
        tree_type,
        colormap_df,
        individuals_lineages,
        tree_image_path
    )
