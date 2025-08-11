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
Created on Wed Jan 19 01:35:14 2022

@author: pietro
"""

import simplicity.tree.newick as nwk

def write_nexus_file(tree, nexus_filepath):
    '''
    Generate NEXUS file from AnyTree tree, including branch lengths in time units
    and time_emergence annotations 
    '''
    root = tree[0]

    newick_tree = nwk.export_newick(root)  

    with open(nexus_filepath, 'w') as f:
        f.write('#NEXUS\n')
        
        # --- TAXA block (lists tip names) ---
        tips = [n for n in tree if getattr(n, 'leaf', False)]
        if tips:
            f.write('begin taxa;\n')
            f.write(f'  dimensions ntax={len(tips)};\n')
            f.write('  taxlabels\n')
            for n in tips:
                f.write(f'    {getattr(n, "lineage", n.name)}\n')
            f.write('  ;\nend;\n\n')

        # # --- DATA block ---
        # f.write('begin data;\n')
        # f.write('  Matrix\n')
        # for node in tree:
        #     if getattr(node, 'leaf', False) or getattr(node, 'label', '') == 'root':
        #         f.write(f'    {getattr(node, "lineage", node.name)}\n')
        # f.write('  ;\nend;\n\n')

        # --- Trees block with time-scaled Newick ---
        f.write('begin trees;\n')
        # newick_tree already has trailing semicolon, so no wrapping or extra ";"
        f.write(f'  tree TimeScaled = {newick_tree}\n')
        f.write('end;\n')
