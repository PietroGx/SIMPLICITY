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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 19:57:13 2022

@author: pietro
"""
import anytree
import collections

def get_newick_str_from_root(node_to_children, label_internal_nodes=False) -> str:
    '''
    Generate newick string from python AnyTree tree.

    label_internal_nodes: internal nodes get no label by default (e.g.
    infection-tree internal nodes are transmission events -- the lineage an
    individual happens to carry at that moment isn't "the lineage that
    emerged here", so labeling them would misrepresent the tree). Pass True
    only where every internal node's 'lineage'/'label' genuinely identifies
    that node (e.g. phylogenetic trees, where internal nodes ARE lineage
    emergence points).
    '''

    def _quote_if_needed(s: str) -> str:
        # Newick-safe quoting for spaces/specials
        if any(ch in s for ch in " \t():,;[]'\""):
            return "'" + s.replace("'", "''") + "'"
        return s

    def newick_render_node(node_to_children) -> str:
        # branch length = time delta
        distance = node_to_children.get('distance', 0.0)
        time_val = node_to_children.get('time_emergence', None)

        # leaf label: lineage, fallback to label
        leaf_label = node_to_children.get('lineage', node_to_children.get('label', ''))

        if 'children' not in node_to_children:
            # leaves
            lbl = _quote_if_needed(str(leaf_label)) if leaf_label is not None else ''
            if time_val is not None:
                lbl = f"{lbl}[&time={time_val}]"
            return f"{lbl}:{distance}"
        else:
            # internal
            children = node_to_children['children']
            children_strings = [newick_render_node(child) for child in children]
            children_strings = ",".join(children_strings)

            if label_internal_nodes and leaf_label is not None:
                lbl = _quote_if_needed(str(leaf_label))
            else:
                lbl = ''
            if time_val is not None:
                lbl = f"{lbl}[&time={time_val}]"

            return f"({children_strings}){lbl}:{distance}"

    newick_string = newick_render_node(node_to_children) + ';'
    return newick_string


def export_newick(root, label_internal_nodes=False):
    # tree to ordered dictionary
    exporter = anytree.exporter.DictExporter(dictcls= collections.OrderedDict, attriter=sorted)
    dic = exporter.export(root)
    # ordered dictionary to newick format
    newick_tree = get_newick_str_from_root(dic, label_internal_nodes=label_internal_nodes)
    return newick_tree

def write_newick_file(root, newick_filepath, label_internal_nodes=False):
    newick_tree = export_newick(root, label_internal_nodes=label_internal_nodes)
    with open(newick_filepath, 'w') as f:
        f.write(newick_tree)
        f.write('\n')





