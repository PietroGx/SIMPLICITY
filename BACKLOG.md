#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 22:01:47 2025

@author: pietro
"""

## Backlog
- [ ] write unit_runs for running with/wo plotting (technical debt item 1)

### Bug Fixes
- [ ] slurm runner gets stuck if there are abrupt job terminations (like OOM or OOT errors): the job is not reported as failed and the runner does not submit a new job to slurm

### Features
- [ ] implement long-shedders subpopulations

### Technical Debt
- [ ] plot_trajectory should be passed to runners with partial function (not changing the function's signature)
- [ ] memory and runtime requests to slurm hardcoded into runner
- [ ] env variables set in dir_manager (should be set somewhere else or rename dir_manager)
- [ ] tau_2 value hard coded into extrande. it should be passed and read from parameters, as maybe should all the other tau's values 
- [ ] in tree builder line 43,78: node of tree will only store the first of the IH lineages (for correct tree lineage coloring).  
      If multiple lineages are stored, decide on coloring strategy.
- [ ] We are not explicitely handling double substitutions in same spot in genome