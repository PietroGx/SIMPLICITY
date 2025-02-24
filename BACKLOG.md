#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 22:01:47 2025

@author: pietro
"""

## Backlog
- [x] function in dir_manager to list all simulation output folders (hardcoded into check_simulations.py)

### Bug Fixes
- [x] memory and runtime options need to be added to slurm runner to avoid Runtime and OOM errors

### Features
- [ ] 

### Technical Debt
- [ ] plot_trajectory should be passed to runners with partial function (not changing the function's signature)
- [ ] memory and runtime requests to slurm hardcoded into runner
- [ ] env variables set in dir_manager (should be set somewhere else or rename dir_manager)
- [ ] tau_2 value hard coded into extrande. it should be passed and read from parameters, as maybe should all the other tau's values 