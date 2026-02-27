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

import simplicity.output_manager as om
from anytree import Node
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import copy
import simplicity.plots_manager as pm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# ---------------------------
# 1) Build the Infection Tree Using Transmitted Lineage Info
# ---------------------------
# For each infection event, we record the parent's transmitted lineage (inherited_lineage)
# by looking into the parent's new_infections list. For patient zero, inherited_lineage is None.

def create_individual_node(row, parent):
    """
    Create a new node representing an individual's infection.
    
    - Branch length is computed as the difference between t_infection times.
    - The node's final (current) lineage is row.IH_lineages[0].

    """
    branch_length = row.t_infection - (parent.t_infection if parent is not None else 0)
    if branch_length < 0:
        raise ValueError(f"Negative branch length for individual {row.Index} (t_infection={row.t_infection}) vs parent's t_infection={parent.t_infection}")
    
    if parent is not None and hasattr(parent, "new_infections"):
        parent_new_infections = parent.new_infections
    else:
        parent_new_infections = []
    
    inherited = None
    for event in parent_new_infections:
        # Compare individual_infected field with current row.Index as strings.
        if str(event.get('individual_infected')) == str(row.Index):
            inherited = event.get('transmitted_lineage')
            break

    current_lineage = row.IH_lineages[0]  # need to update for multiple lineages
    
    node = Node(
        name=str(row.Index),         # internal name for bookkeeping
        parent=parent,
        distance=branch_length,
        label=str(row.Index),        # display label: individual ID
        leaf=True,
        t_infection=row.t_infection,
        t_infectious=row.t_infectious,
        t_not_infectious=row.t_not_infectious,
        state=row.state,
        infection_type=row.type,
        fitness_score=row.fitness_score,
        lineage=current_lineage       # final (current) lineage
    )
    node.inherited_lineage = inherited  # store parent's transmitted lineage
    # Directly store new_infections if present.
    node.new_infections = row.new_infections if hasattr(row, "new_infections") else []
    return node

def infection_tree(data):
    """
    Constructs an infection tree from CSV data.
    
    Each CSV row represents one infection event.
    The 'parent' column is used to attach nodes.

    """
    # Create the root node.
    root = Node(
        name='root',
        parent=None,
        distance=0,
        label='root',
        leaf=False,
        t_infection=0,
        t_infectious=0,
        t_not_infectious=0,
        fitness_score=0,
        infection_type='standard',
        lineage='wt'
    )
    root.inherited_lineage = None  # Patient zero did not inherit from anyone.
    # If the root row has new_infections, store it.
    if 'new_infections' in data.columns:
        root.new_infections = data.iloc[0].new_infections
    else:
        root.new_infections = []
    
    nodes_lookup = {}
    # Process patient zero entries (rows where parent == 'root').
    patient_zero_data = data[data['parent'] == 'root']
    for row in patient_zero_data.itertuples():
        node = create_individual_node(row, root)
        nodes_lookup[node.label] = node
    data = data.drop(patient_zero_data.index)
    
    # Process remaining rows in order of t_infection.
    sorted_data = data.sort_values(by='t_infection')
    for row in sorted_data.itertuples():
        parent_label = str(row.parent)
        parent_node = nodes_lookup.get(parent_label)
        if parent_node is None:
            print(f"Warning: Parent '{parent_label}' not found for row {row.Index}. Skipping.")
            continue
        node = create_individual_node(row, parent_node)
        parent_node.leaf = False
        nodes_lookup[node.label] = node
        
    return root

# ---------------------------
# 2) Conversion from anytree to ete3 Tree
# ---------------------------
def anytree_to_ete(node):
    """
    Recursively converts an anytree Node into an ete3 Tree.
    """
    t = Tree()
    t.name = node.label  # display label (individual ID)
    t.dist = getattr(node, "distance", 0)
    t.add_feature("t_infection", getattr(node, "t_infection", 0))
    t.add_feature("t_not_infectious", getattr(node, "t_not_infectious", 0))
    t.add_feature("lineage", getattr(node, "lineage", ""))
    inherited = getattr(node, "inherited_lineage", None)
    t.add_feature("inherited_lineage", inherited if inherited is not None else "")
    for child in node.children:
        t.add_child(anytree_to_ete(child))
    return t

# ---------------------------
# 3) Filtering Function (by transmitted lineage)
# ---------------------------
def filter_by_transmitted_lineage(ete_tree, target_lineage):
    """
    Returns a deep-copied, pruned version of ete_tree such that:
      - We keep the parent node if it has at least one child whose inherited_lineage == target_lineage,
      - We also keep those child nodes themselves,
      - And we keep all ancestors of those parent nodes.
    
    This ensures the entire parent–child transmission link for the target_lineage is in the final subtree,
    preserving a connected chain of transmissions for that lineage.
    """
    # Collect the nodes (parents) that transmitted the target lineage
    # AND also collect the child nodes that inherited the target lineage.
    nodes_to_keep = set()

    for nd in ete_tree.traverse():
        for child in nd.get_children():
            if child.inherited_lineage == target_lineage:
                # Keep the parent
                nodes_to_keep.add(nd.name)
                # Keep the child
                nodes_to_keep.add(child.name)
                # Keep only those ancestors that have inherited_lineage equal to target_lineage.
                for anc in nd.get_ancestors():
                    # if getattr(anc, "inherited_lineage", None) == target_lineage:
                    nodes_to_keep.add(anc.name)

    if not nodes_to_keep:
        return None
    
    # Prune the deep copy of the original tree
    subtree = copy.deepcopy(ete_tree)
    subtree.prune(list(nodes_to_keep), preserve_branch_length=True)
    
    return subtree

# ---------------------------
# 4) Plotting with ete3
# ---------------------------
def plot_infection_tree_ete(anytree_root, target_lineage=None, layout="r", 
                            color=None, outfile=None):
    """
    Converts the anytree infection tree to an ete3 Tree and plots it.
    
    Parameters:
      anytree_root: the root of the anytree-based infection tree.
      target_lineage: if provided, the tree is pruned to include only transmission events
                      where the child's inherited_lineage equals target_lineage (plus ancestors).
                      Also used in the title.
      layout: "r" for rectangular (time-oriented) or "c" for circular (radial).
      outfile: if provided, the tree is rendered to this file (e.g., "tree.png")
               instead of being displayed interactively.
    
    In rectangular mode, branch lengths (t_infection differences) appear horizontally.
    Each node displays its name (individual ID) and its final lineage.
    """
    ete_tree = anytree_to_ete(anytree_root)
    
    if target_lineage is not None:
        filtered_subtree = filter_by_transmitted_lineage(ete_tree, target_lineage)
        if filtered_subtree is not None:
            ete_tree = filtered_subtree
        else:
            print(f"Warning: No transmission events found with inherited_lineage '{target_lineage}'. Plotting full tree.")
    else:
        target_lineage = "all"
    
    ts = TreeStyle()
    if layout == "r":
        ts.mode = "r"
        # ts.title.add_face(TextFace(f"Infection Tree (Rectangular) - Transmissions of '{target_lineage}'", fsize=12), column=0)
    elif layout == "c":
        ts.mode = "c"
        # ts.title.add_face(TextFace(f"Infection Tree (Radial) - Transmissions of '{target_lineage}'", fsize=12), column=0)
    else:
        raise ValueError("layout must be 'r' or 'c'.")
    
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False  # disable automatic leaf labels
    ts.show_scale = False       # disable the scale bar (black bar)

    def layout_fn(node):
        ns = NodeStyle()
        ns["shape"] = "circle"
        ns["size"] = 5
        if color is not None:
            ns["fgcolor"] = color
            
        node.set_style(ns)
        # --- Display the node's name (ID).
        # node.add_face(TextFace(node.name, fsize=8), column=0, position="branch-right")
        # --- Display the node's final lineage in red.
        # node.add_face(TextFace(" [" + node.lineage + "]", fsize=8, fgcolor="red"), column=1, position="branch-right")
        # --- inherited lineage label
        # inherited_label = f" {{{node.inherited_lineage}}}" if node.inherited_lineage else ""
        # node.add_face(TextFace(inherited_label, fsize=8, fgcolor="green"), column=1, position="branch-right")
        # --- Display the node's time of infection
        node.add_face(TextFace(f"t:{node.t_infection:.1f}", fsize=8, fgcolor="black"), column=1, position="branch-right")
    ts.layout_fn = layout_fn

    if outfile:
        ete_tree.render(outfile, w=1600, dpi=600, tree_style=ts)
        print(f"Tree rendered to {outfile}")
    else:
        ete_tree.show(tree_style=ts)

# ---------------------------------------------------------------------------------
# Build the infection tree from your imported dataframe.
ssod = 'Data/Tree_data/04_Output/T_1095/seed_0007'
df = om.read_individuals_data(ssod)
anytree_root = infection_tree(df)

# Specify target lineage for filtering: we want the full transmission chain
# where a child inherited the specified lineage from its parent.
lineage_name = "wt.7"

colors_df = pm.make_lineages_colormap(ssod)
color = colors_df[colors_df['Lineage_name'] == lineage_name]['Color'].iloc[0]


# 1. Plot the full tree in rectangular mode.
# plot_infection_tree_ete(anytree_root, target_lineage=None, layout="r", outfile='tree_linear.png')

# 2. Plot only the subtree corresponding to transmission events where the child's inherited_lineage equals "wt".
plot_infection_tree_ete(anytree_root, target_lineage=lineage_name, layout="r", color=color,
                        outfile=f"infection_tree_filtered_{lineage_name}.png")

# --------------------------------------------------------------------------------
import numpy as np
# Prepare the histogram data
filtered_df = df[df['inherited_lineage'] == lineage_name]

# Determine the x-axis limits from the data
xmin = filtered_df['t_infection'].min()
xmax = filtered_df['t_infection'].max()

# Create combined plot with two subplots (histogram above, tree image below)
fig, (ax_hist, ax_tree) = plt.subplots(2, 1, figsize=(10, 8))

# Plot the histogram on the top subplot
ax_hist.hist(filtered_df['t_infection'], bins=20, color=color, edgecolor='black')
ax_hist.set_ylabel("New Infections")
ax_hist.set_title("Infection Histogram and Transmission Tree")
ax_hist.set_xlim(0,xmax)
ax_hist.set_xticks(np.arange(0, xmax + 1, 7))


# Load the tree image and display it on the bottom subplot
outfile=f"infection_tree_filtered_{lineage_name}.png"
tree_img = mpimg.imread(outfile)
# Here we set the extent so that the image's x-axis corresponds to [xmin, xmax]
ax_tree.imshow(tree_img, extent=[xmin, xmax, 0, 1], aspect='auto')
ax_tree.axis('off')  # Turn off axis labels and ticks for the tree image
ax_tree.set_xlabel("Time of Infection")

plt.tight_layout()
plt.show()

