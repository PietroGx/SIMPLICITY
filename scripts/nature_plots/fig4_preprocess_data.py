import os
import sys

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tree.tree_builder as tb
import baltic as bt

# get_lineage_labels/_clade_analysis already work against plain sod/ssod
# paths (no dependency on the old grid) -- reuse rather than duplicate.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_preprocess import get_lineage_labels

# _experiment_sod is the same "resolve one scenario's sod, or None if that
# scenario hasn't been run yet" helper already used for Figure 3 -- same
# directory, plain sibling import.
from fig3_preprocess_data import _experiment_sod

EXP_NAME = "impact_long_shedders"


def _get_name(k):
    """
    baltic stores node identity as a plain `.name` attribute for leaves, but
    (after simplicity.tree.newick's label_internal_nodes fix) as
    `.traits['label']` for internal nodes -- normalize the difference here.
    """
    return getattr(k, 'name', None) or k.traits.get('label')


def _get_newick_path(ssod):
    """Build-if-missing helper for the phylogenetic Newick file of one ssod
    (mirrors scripts/nextstrain.py's get_json_tree, for Newick instead of JSON)."""
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    path = om.get_tree_file_filepath(
        experiment_name, ssod, tree_type="phylogenetic", tree_subtype="binary",
        file_type="newick")

    if not os.path.isfile(path):
        tb.get_tree(
            experiment_name=experiment_name,
            seeded_simulation_output_dir=ssod,
            tree_type="phylogenetic",
            tree_subtype="binary",
            coloring="lineage",
            save_plot=False,
            export_filetype="newick",
            dashplot=False,
        )

    return path


def get_tree_for_scenario(exp_num, scenario, seed, cluster_threshold=5):
    """
    Loads scenario's phylogenetic tree (for `seed`) into a baltic tree object,
    with every branch's `.traits['cluster']` set to its long/standard/mixed/
    founder label (same color language as Figure 3). Returns None if the
    scenario hasn't been run yet.
    """
    sod = _experiment_sod(exp_num, scenario, exp_name=EXP_NAME)
    if sod is None:
        return None

    ssod = dm.get_ssod(sod, int(seed))
    newick_path = _get_newick_path(ssod)

    tre = bt.loadNewick(newick_path, tip_regex=None, sortBranches=False, absoluteTime=False)
    lin2label = get_lineage_labels(ssod, cluster_threshold)
    for k in tre.Objects:
        k.traits['cluster'] = lin2label.get(_get_name(k), 'founder')

    return tre
