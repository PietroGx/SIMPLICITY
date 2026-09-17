#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ============================================================================
# Shared scenario resolution + output paths for the paper figures.
# ----------------------------------------------------------------------------
# The four figures must agree on which scenarios exist, so the list comes from
# the pipeline's OWN config rather than being repeated in each script:
#
#   impact_long_shedders          -> 4 scenarios
#   impact_long_shedders_unbound  -> 5 scenarios, including edge_case
#
# resolve_scenarios also reports which of those have no output on disk. A
# scenario is DROPPED but WARNED about -- never silently omitted. Discovering
# the list from disk instead would turn a half-finished run into a smaller
# figure with no indication anything was missing, which is the failure mode
# this repo has been bitten by before (see BACKLOG).
# ============================================================================

import os
import sys

import simplicity.dir_manager as dm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'experiments'))

BOUND_EXP_NAME = "impact_long_shedders"
UNBOUND_EXP_NAME = "impact_long_shedders_unbound"

# exp_name -> the module holding that pipeline's SCENARIOS
_CONFIG_MODULES = {
    BOUND_EXP_NAME: "impact_long_shedders_config",
    UNBOUND_EXP_NAME: "impact_long_shedders_unbound_config",
}

# exp_name -> short label for output filenames
_ARM_LABELS = {
    BOUND_EXP_NAME: "bound",
    UNBOUND_EXP_NAME: "unbound",
}


def scenario_names(exp_name=BOUND_EXP_NAME):
    """Every scenario this pipeline defines, in config order."""
    module_name = _CONFIG_MODULES.get(exp_name)
    if module_name is None:
        raise SystemExit(
            f"No scenario config known for --exp-name {exp_name!r}. "
            f"Known: {sorted(_CONFIG_MODULES)}.")
    module = __import__(module_name)
    return [s["name"] for s in module.SCENARIOS]


def has_output(exp_name, scenario, exp_num):
    """True if this scenario's experiment produced output."""
    try:
        return bool(dm.get_simulation_output_dirs(
            f"{exp_name}_{scenario}_#{exp_num}"))
    except Exception:
        # dir_manager raises rather than returning [] for an absent experiment
        return False


def resolve_scenarios(exp_name, exp_num, required=None, verbose=True):
    """(present, missing) for this pipeline at this experiment number.

    `required` optionally restricts to a subset of the config's scenarios --
    used by panels that deliberately plot only some of them.
    """
    names = scenario_names(exp_name)
    if required is not None:
        names = [n for n in names if n in set(required)]

    present = [n for n in names if has_output(exp_name, n, exp_num)]
    missing = [n for n in names if n not in present]

    if verbose:
        print(f"[scenarios] {exp_name} #{exp_num}: "
              f"{len(present)}/{len(names)} with output -> {present}")
        for n in missing:
            print(f"[scenarios][warn] '{n}' is defined in the config but has no "
                  f"output for #{exp_num} -- it will be left out of the figure.")
    return present, missing


def arm_label(exp_name):
    """Short label for the pipeline arm, for output filenames."""
    return _ARM_LABELS.get(exp_name, exp_name)


def figures_dir():
    """Data/figures -- created on demand."""
    path = os.path.join(dm.get_data_dir(), "figures")
    os.makedirs(path, exist_ok=True)
    return path


def figure_path(number, exp_name, exp_num, fmt):
    """Data/figures/Figure_<n>_<arm>_#<exp_num>.<fmt>.

    The arm and run number are in the name so the two pipelines cannot
    overwrite each other -- all four figures previously wrote a bare
    Figure_N.<fmt> into the working directory.
    """
    return os.path.join(
        figures_dir(),
        f"Figure_{number}_{arm_label(exp_name)}_#{exp_num}.{fmt}")
