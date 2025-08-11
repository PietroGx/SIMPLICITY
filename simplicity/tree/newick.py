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

def get_newick_str_from_root(node_to_children) -> str:
    '''
    Generate newick string from python AnyTree tree. 
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
        internal_label = ''

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

            lbl = internal_label
            if time_val is not None:
                lbl = f"{lbl}[&time={time_val}]"

            return f"({children_strings}){lbl}:{distance}"

    newick_string = newick_render_node(node_to_children) + ';'
    return newick_string


def export_newick(root):
    # tree to ordered dictionary
    exporter = anytree.exporter.DictExporter(dictcls= collections.OrderedDict, attriter=sorted)
    dic = exporter.export(root)
    # ordered dictionary to newick format
    newick_tree = get_newick_str_from_root(dic)
    return newick_tree

def write_newick_file(root, newick_filepath):
    newick_tree = export_newick(root)
    with open(newick_filepath, 'w') as f:
        f.write(newick_tree)
        f.write('\n')





