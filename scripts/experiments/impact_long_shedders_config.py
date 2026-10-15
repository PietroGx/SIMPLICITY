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
# impact_long_shedders -- PIPELINE CONFIG
# ----------------------------------------------------------------------------
# The single place that knows HOW TO BUILD the parameter settings for every
# stage. impact_long_shedders_cal_1.py / _cal_2.py / _exp.py only RUN things:
# parse CLI args, call a build_*_settings() function below to get a
# ready-to-submit settings callable, dispatch it via run_experiment_script,
# and handle stage-specific I/O. They never construct parameter dicts
# themselves -- if a parameter dict changes, this is the only file to touch.
# SCENARIOS and NSR_RANGES are plain code-level data, not external reference
# files: retuning them is a normal code edit.
#
# Each cohort is calibrated ALONE and they are then let loose together.
# NSR_long comes from a 100% long-shedder population, the standard NSR from a
# 0% long-shedder population, and the experiment phase applies both unchanged.
# The population-level clock is an OUTPUT -- free to rise above
# target_osr_std -- rather than a constraint, and the two calibrations are
# independent: stage 2 does not read stage 1.
#
# Until v2.4.90 a second, "bound" pipeline calibrated the standard NSR IN
# CONTEXT instead, per scenario inside the mixed population, so the standard
# rate was pushed down however far it took to hold the population clock at
# target. Run #7 drove HIV_high to 9.3e-7, 150x below its own control and
# below the swept floor. That result is the argument for the design kept
# here; the code that produced it is gone.
#
# Measurement code is deliberately NOT duplicated: stage 1 reuses
# long_nsr_calibration_plot.plot_and_fit_long_nsr_calibration and therefore
# evolutionary_rate.extract_ih_regression_data. Parameter construction lives
# here, in the config.
# ============================================================================

import os

import numpy as np
import pandas as pd

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm

# =============================================================================
# SHARED MATH -- tau_3_long derivation (R_long is a plain SCENARIOS value,
# nothing to derive)
# =============================================================================
def phase_offset(sp):
    """Sum of the non-long-shedding infection phases (tau_1+tau_2+tau_4)."""
    return sp["tau_1"] + sp["tau_2"] + sp["tau_4"]


def derive_tau_3_long(scenario, sp):
    """Corrected long-shedding duration: inf_duration_long - phase_offset(sp)
    for a long-shedder scenario, or sp's own tau_3_long default for control."""
    if scenario["long_shedders_ratio"] == 0.0:
        return float(sp["tau_3_long"])
    return float(scenario["inf_duration_long"]) - phase_offset(sp)


# =============================================================================
# SCENARIOS
# =============================================================================
# edge_case is back. It was dropped in v2.4.24 for raw-Hamming saturation at
# "~0.34 substitutions/site", but that figure dates from before the intra-host
# measurement was corrected, when the calibration was returning NSR_long ~0.116
# -- about 100x too high. At the rate the corrected pipeline actually
# calibrates (~1.1e-3), a 350-day infection accumulates 0.0020 subs/site, some
# 50x below the saturation threshold. Its 0.1% prevalence puts exactly one
# carrier in the initial 1000, which is the point: what does a single
# super-long shedder do to an outbreak?
SCENARIOS = [
    {"name": "control",   "inf_duration_long": None,  "long_shedders_ratio": 0.000, "R_long": None, "susceptibility_long": 1.0},
    {"name": "SOT",       "inf_duration_long": 63.0,  "long_shedders_ratio": 0.010, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "HIV_low",   "inf_duration_long": 109.0, "long_shedders_ratio": 0.010, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "HIV_high",  "inf_duration_long": 109.0, "long_shedders_ratio": 0.120, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "edge_case", "inf_duration_long": 365.0, "long_shedders_ratio": 0.001, "R_long": 1.1,  "susceptibility_long": 1.0},
]

TAU_ROUND = 3  # decimals for tau_3_long dict-key / group matching

# distribution first: it is the default and the mode the paper uses. argmax
# stays available for code-level investigations.
#
# Declared here, above every builder, because ALL THREE stages write it into
# fixed_params. They did not always: until v2.4.93 only the experiment phase
# set it and the two calibration stages inherited whatever
# standard_values.json happened to say. That file is machine-local and
# gitignored, and read_standard_parameters_values merges it OVER the code
# defaults, so flipping the default in simplicity/ does not reach a machine
# that already has the file. Stating it here is what makes a calibration's
# consensus mode a property of the pipeline instead of the host.
CONSENSUS_MODES = ("distribution", "argmax")
DEFAULT_CONSENSUS = "distribution"

# ONE experiment per run (stage 8e). Every phase is a set of named GROUPS
# inside impact_long_shedders_#N; there is no longer an experiment per
# calibration stage or per scenario.
EXP_NAME = "impact_long_shedders"

# setup_data keeps its own directory: it holds the frozen table, which is an
# input to the next phase rather than simulation output, and it is what the
# experiment phase reads to know what to run.
SETUP_DIR_TEMPLATE = "Data/impact_long_shedders_setup_data_#{exp_num}"
TABLE_FILENAME = "nsr_calibration_table.csv"

# ---------------------------------------------------------------------------
# GROUP NAMES -- the identity of a phase's slice of the experiment
# ---------------------------------------------------------------------------
# The group name IS the identity, in five places at once: the `groups` list in
# settings.json, the path in repeats.json, the state directory, 04_Output/
# <group>/, and the fit-results filename. So it must be unique across the
# WHOLE experiment now that phases share one, and it must be FLAT -- no "/" --
# because the same string is both a directory component and part of a
# filename.
#
# The prefix is what keeps a calibration cohort and an experiment-phase
# scenario of the same name apart: "SOT" is a cal_1 cohort AND a scenario, and
# before the prefixes they would have collided into one group mixing a sweep
# with a single run.
CAL2_GROUP_NAME = 'standard_only'

CALIB_LONG_PREFIX = "calib_long__"
CALIB_STD_PREFIX = "calib_std__"
EXP_PREFIX = "exp__"


def calib_long_group(cohort):
    """cal_1's group for one (tau_3_long, R_long) cohort."""
    return f"{CALIB_LONG_PREFIX}{cohort}"


def calib_std_group():
    """cal_2's single standard-only group."""
    return f"{CALIB_STD_PREFIX}{CAL2_GROUP_NAME}"


def exp_group(scenario_name):
    """The experiment phase's group for one scenario."""
    return f"{EXP_PREFIX}{scenario_name}"


def scenario_experiment(exp_name, exp_num):
    """The ONE experiment a run writes. Since 8e there is no experiment per
    scenario: a scenario is a GROUP inside this."""
    return f"{exp_name}_#{exp_num}"


def scenario_output_dirs(exp_name, scenario, exp_num):
    """Every simulation output directory of ONE scenario, and nothing else.

    THE call every figure and sanity plot makes. Before 8e each scenario was
    its own experiment, so "<exp_name>_<scenario>_#N" named it and an
    unscoped listing could not reach another scenario. Now all five
    scenarios, both calibration sweeps and every one of their grid points
    live in one experiment, and an unscoped listing returns the lot -- a
    panel would plot a calibration sweep as if it were its scenario. That is
    the exact shape CLAUDE.md records as having bitten the figures twice.

    Returns [] when the scenario has no group here, which is what the callers
    already handle; resolve_scenarios is what warns about it.
    """
    experiment = scenario_experiment(exp_name, exp_num)
    group = exp_group(scenario)
    try:
        sods = dm.get_simulation_output_dirs(experiment, group=group)
    except Exception:                                           # noqa: BLE001
        return []                      # no such experiment, or no such group
    # The group is a path component of everything it scopes. If that stops
    # holding, this series is reading another group's rows under this
    # scenario's label -- loudly, rather than as a plausible-looking panel.
    stray = [sod for sod in sods if group not in sod.split(os.sep)]
    if stray:
        raise RuntimeError(
            f"{experiment}: group {group!r} returned directories that are not "
            f"in it: {stray[:3]}. A panel labelled {scenario!r} would be "
            f"plotting another group's simulations.")
    return sods


def scenario_label(exp_name, scenario, exp_num):
    """What to call this series in a message: experiment AND group, because
    the experiment name alone no longer identifies the scenario."""
    return f"{scenario_experiment(exp_name, exp_num)}/{exp_group(scenario)}"


def phase_of(group):
    """Which phase a group belongs to, or None. The inverse of the three
    helpers above, for readers that are handed a bare group name."""
    for phase, prefix in (("calib_long", CALIB_LONG_PREFIX),
                          ("calib_std", CALIB_STD_PREFIX),
                          ("exp", EXP_PREFIX)):
        if group.startswith(prefix):
            return phase
    return None


def groups_of_phase(groups, phase):
    """Every group of `groups` belonging to `phase`, order preserved."""
    return [g for g in groups if phase_of(g) == phase]


def derive_scenario_params(scenario, sp):
    """Frozen long-side parameters for one production scenario. Mirrors the
    bound pipeline's version -- same tau derivation, same directly-specified
    R_long, same diagnosis-only sequencing (sequence_long_shedders False)."""
    name = scenario["name"]
    ratio = scenario["long_shedders_ratio"]
    tau_3_long = derive_tau_3_long(scenario, sp)
    is_long = ratio > 0.0

    return {
        "scenario_name": name,
        "long_shedders_ratio": float(ratio),
        "susceptibility_long": float(scenario["susceptibility_long"]),
        "tau_3_long": tau_3_long,
        "R_long": float(scenario["R_long"]) if is_long else 0.0,
        "sequence_long_shedders": False,
        "is_long": is_long,
    }


def unique_long_tau_r_long_pairs(sp):
    """Sorted unique (tau_3_long, R_long) pairs across long-shedder scenarios.
    HIV_low/HIV_high share one; SOT and edge_case each add their own, so
    stage 1 fits three groups here against the bound pipeline's two."""
    pairs = {
        (round(derive_tau_3_long(s, sp), TAU_ROUND),
         round(float(s["R_long"]), TAU_ROUND))
        for s in SCENARIOS if s["long_shedders_ratio"] > 0.0
    }
    return sorted(pairs)


def group_name_for_pair(pair, sp):
    """The name of the cal_1 group covering one (tau_3_long, R_long) pair,
    over THIS pipeline's SCENARIOS -- which include edge_case.

    Recorded on the group at definition. long_nsr_calibration_plot used to map a
    pair back to a name against the BOUND pipeline's scenario list by default,
    so edge_case fell through to a raw 'tau=350.23,R_long=1.1' fallback in
    unbound run #1 -- and that string landed in the fit-results filename.
    """
    names = [s["name"] for s in SCENARIOS
             if s["long_shedders_ratio"] > 0.0
             and (round(derive_tau_3_long(s, sp), TAU_ROUND),
                  round(float(s["R_long"]), TAU_ROUND)) == pair]
    return "+".join(names) if names else f'tau_{pair[0]}_Rl_{pair[1]}'


def lookup_long_nsr(long_nsr_by_group, tau_3_long, r_long):
    """Exact (rounded) match on (tau_3_long, R_long). Raises on miss rather
    than guessing."""
    key = (round(float(tau_3_long), TAU_ROUND), round(float(r_long), TAU_ROUND))
    if key not in long_nsr_by_group:
        raise KeyError(
            f"No calibrated long NSR for (tau_3_long, R_long)={key}. "
            f"Available: {sorted(long_nsr_by_group)}.")
    return long_nsr_by_group[key]


# =============================================================================
# STAGE 1 (cal_1) -- isolated long-shedder NSR calibration
# ----------------------------------------------------------------------------
# A 100% long-shedder population, used ONLY to calibrate the long-shedder NSR
# in isolation -- deliberately different from USER_FIXED_PARAMS (Stage 2/3's
# mixed production population), not a duplicate of it.
# =============================================================================
CAL1_ISOLATED_FIXED_PARAMS = {
    "long_shedders_ratio": 1.0,
    "R": 1.0,
    "population_size": 1000,
    "infected_individuals_at_start": 10,
    "sequence_long_shedders": True,
    "final_time": 365,
    # R_long absent: build_cal1_settings sets it per (tau_3_long, R_long)
    # pair straight from SCENARIOS.
}


# -----------------------------------------------------------------------------
# The simulated window is PER GROUP, not global.
# -----------------------------------------------------------------------------
# The window is PER GROUP, not global. A group needs room for a few complete
# infection cycles, and infection duration ranges from 63 days (SOT) to 365
# (edge_case). SOT and HIV fit several cycles into the standard 365-day window;
# edge_case does not fit even one, and its census only fires on RECOVERED long
# shedders, so a 365-day window would starve the min_seq gate -- HPC test #1 in
# BACKLOG found exactly that and confirmed it clears at 3 years. The rule below
# gives 365 to SOT and HIV and 3x365 to edge_case, without hardcoding a
# scenario name: any new duration gets a window that fits it.
CAL1_FINAL_TIME = 365
CAL1_MIN_CYCLES = 3


def cal1_final_time(tau_3_long, sp):
    """Simulated window for one cal_1 group: the standard one, unless the
    group's own infection duration needs more to complete
    CAL1_MIN_CYCLES cycles."""
    inf_duration = tau_3_long + phase_offset(sp)
    needed = CAL1_MIN_CYCLES * inf_duration
    return int(max(CAL1_FINAL_TIME,
                   np.ceil(needed / CAL1_FINAL_TIME) * CAL1_FINAL_TIME))


def build_cal1_settings(seeds, ranges, R=None, ih_virus_emergence_rate=None,
                        consensus=DEFAULT_CONSENSUS):
    """One group per (tau_3_long, R_long) pair, each sweeping the same
    NSR_long grid in a 100% long-shedder population, each with its own
    simulated window (see cal1_final_time). The consensus mode is stated,
    not inherited: it must match the one the experiment phase runs, or the
    rates calibrated here do not transfer."""
    sp = sm.read_standard_parameters_values()
    fixed_params = CAL1_ISOLATED_FIXED_PARAMS.copy()
    fixed_params["final_time"] = CAL1_FINAL_TIME
    fixed_params["consensus"] = consensus
    if R is not None:
        fixed_params["R"] = R
    fixed_params["IH_virus_emergence_rate"] = (
        USER_FIXED_PARAMS['IH_virus_emergence_rate']
        if ih_virus_emergence_rate is None else ih_virus_emergence_rate)

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        scenario_groups = [
            {
                'name': calib_long_group(group_name_for_pair((tau, r_long), sp)),
                'nucleotide_substitution_rate_long': nsr_values,
                'tau_3_long': tau,
                'R_long': r_long,
                # per-group scalar: generate_experiment_settings applies group
                # overrides after fixed_params, so this wins for this group only
                'final_time': cal1_final_time(tau, sp),
            }
            for tau, r_long in unique_long_tau_r_long_pairs(sp)
        ]
        return ({'_scenario_groups': scenario_groups}, fixed_params, seeds)

    return make_settings


# =============================================================================
# STAGE 2 -- standard individuals alone
# =============================================================================
# A 0%-long-shedder population at production settings: the same context every
# scenario runs in, minus the long side. One sweep, one calibrated standard
# NSR, shared by every scenario -- there is nothing scenario-specific left to
# fit once calibration stops compensating for the long side.
#
# final_time 365, matching the standard cal_1 window. Production runs 1095, so
# control will not sit exactly on target_osr_std -- but in this pipeline that
# is not what the standard rate is for: it defines "the rate at which standard
# individuals mutate", measured in a clean population, and the population clock
# is deliberately an output. Worth watching control's production clock all the
# same, since it is the reference every elevation is measured against.
#
# There is no sequencing_rate to set: a run writes the COMPLETE diagnosed
# record, which is what the old rate 1.0 produced and what both calibration
# stages always asked for (simplicity/sequencing.py). Any lower ratio is a
# post-hoc subsample via reconstruct_sequencing_data. The CAL2_SEQUENCING_RATE
# constant both configs carried was dead by then and is not folded in.
# Fixed experiment overrides shared by cal_2 and exp (mirrors the population
# context of the impact experiment). Single source of truth -- previously
# duplicated verbatim in both scripts.
USER_FIXED_PARAMS = {
    "population_size": 1000,
    "infected_individuals_at_start": 50,
    "R": 1.03,
    "final_time": 1095,
    "IH_virus_emergence_rate": 0.1,
}

# cal_2 calibrates over a shorter window than the experiment phase runs on
# (USER_FIXED_PARAMS['final_time']); build_cal2_settings applies it.
CAL2_FINAL_TIME = 365

# Shared color palette for the calibration-fit plots (long-NSR per tau group,
# standard-NSR per scenario), keyed by scenario name so both plots look
# consistent with each other.
DEFAULT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]


def build_cal2_settings(seeds, ranges, R=None, ih_virus_emergence_rate=None,
                        consensus=DEFAULT_CONSENSUS):
    """Single standard-only sweep: long_shedders_ratio 0, standard NSR swept.
    Same consensus rule as cal_1 -- stated, not inherited."""
    fixed_params = USER_FIXED_PARAMS.copy()
    fixed_params.update({
        "long_shedders_ratio": 0.0,
        "sequence_long_shedders": False,
        "final_time": CAL2_FINAL_TIME,
        "consensus": consensus,
    })
    if R is not None:
        fixed_params["R"] = R
    if ih_virus_emergence_rate is not None:
        fixed_params["IH_virus_emergence_rate"] = ih_virus_emergence_rate

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        # Named, so the group on disk and the 'standard_only' the fit used to
        # hardcode are the same string by construction rather than by habit.
        varying_params = {'_scenario_groups': [
            {'name': calib_std_group(),
             'nucleotide_substitution_rate': nsr_values}]}
        return (varying_params, fixed_params, seeds)

    return make_settings


# =============================================================================
# STAGE 3 -- the experiment phase
# =============================================================================
# The suffix keeps the two modes from colliding in Data/ at the same
# --exp-num. The DEFAULT mode owns the bare experiment name, so this flipped
# with the default in v2.4.92: distribution is now "", argmax "_amx". Runs
# made before that carry the old convention -- a bare name is argmax and a
# "_dist" name is distribution -- so a directory's name alone does not say
# which mode produced it across that boundary. Check the run's own
# simulation parameters, which record `consensus` either way.
CONSENSUS_SUFFIX = {"distribution": "", "argmax": "_amx"}


# Outputs a calibration run does not need. cal_1's fit and cal_2's fit both
# read only final_time.csv and sequencing_data_regression.csv, and cal_1's plot
# and the sanity regressions additionally go through
# evolutionary_rate.extract_ih_regression_data, which reads individuals_data.csv
# and phylogenetic_data.csv. Nothing in either path opens the other three.
#
# lineage_frequency.csv alone reached 100.6 GB on the calibrated grid -- a
# single simulation of the top NSR_long sweep point weighed 1.1 GB, and that
# point exists only to bracket the fit from above.
CALIBRATION_SKIPPED_OUTPUTS = ("lineage_frequency", "simulation_trajectory",
                               "sequencing_data")


def set_calibration_output_env():
    """Tell the calibration simulations which outputs to leave unwritten.

    An env var rather than a parameter because parameter_specs.json and
    standard_values.json are machine-local and gitignored, so a new parameter
    would not travel with the code. simplicity.runners.slurm forwards the whole
    environment through sbatch, so this reaches every task.
    """
    os.environ[dm.SKIP_OUTPUTS_ENV] = ",".join(CALIBRATION_SKIPPED_OUTPUTS)


def exp_name(consensus):
    """Experiment-phase name prefix for one consensus mode."""
    if consensus not in CONSENSUS_MODES:
        raise ValueError(f"consensus must be one of {CONSENSUS_MODES}, got {consensus!r}")
    return EXP_NAME + CONSENSUS_SUFFIX[consensus]


def build_exp_scenario_settings(row, n_seeds, consensus=DEFAULT_CONSENSUS):
    """One simulation point per frozen-table row. Both rates are applied as
    calibrated; nothing here compensates for their interaction.

    `consensus` selects how a lineage's distance from the population consensus
    is measured -- 'distribution' (the default) against the full column
    distribution, 'argmax' against the consensus sequence. Everything else is
    identical, which is what makes the two comparable.
    """
    if consensus not in CONSENSUS_MODES:
        raise ValueError(f"consensus must be one of {CONSENSUS_MODES}, got {consensus!r}")
    is_control = float(row["long_shedders_ratio"]) == 0.0

    fixed = USER_FIXED_PARAMS.copy()
    fixed.update({
        "R": float(row["R"]),
        "IH_virus_emergence_rate": float(row["IH_virus_emergence_rate"]),
        "long_shedders_ratio": float(row["long_shedders_ratio"]),
        "susceptibility_long": float(row["susceptibility_long"]),
        "tau_3_long": float(row["tau_3_long"]),
        "R_long": float(row["R_long"]),
        "nucleotide_substitution_rate": float(row["nucleotide_substitution_rate"]),
        "sequence_long_shedders": False,
        "consensus": consensus,
    })

    if not is_control:
        long_nsr = row["nucleotide_substitution_rate_long"]
        if pd.isna(long_nsr):
            raise ValueError(
                f"Scenario '{row['scenario_name']}' has long_shedders_ratio>0 "
                f"but no nucleotide_substitution_rate_long in the table.")
        fixed["nucleotide_substitution_rate_long"] = float(long_nsr)

    def make_settings():
        # One named group, so the output says which scenario it is
        # and the fit-results filename cannot disagree with the output tree.
        # An empty group dict is still one group: generate_experiment_settings
        # expands it to combinations = [()], exactly as empty varying_params did.
        return ({'_scenario_groups': [{'name': exp_group(row['scenario_name'])}]},
                fixed.copy(), n_seeds)

    return make_settings


def build_exp_phase_settings(rows, n_seeds, consensus=DEFAULT_CONSENSUS):
    """The WHOLE experiment phase as one phase: one named group per scenario.

    Each group carries its own frozen parameters as scalar overrides, which
    generate_experiment_settings applies after fixed_params, so a group's
    values win for that group only. The per-scenario dict is built by
    build_exp_scenario_settings, not rebuilt here -- the two must not be able
    to drift.

    Before 8e each scenario was its own EXPERIMENT, and the stage submitted
    five of them from five threads with the Slurm concurrency cap divided
    between them. One experiment means one submission: the whole phase is a
    single array and the cap applies to it once.
    """
    groups = []
    for row in rows:
        _, fixed, _ = build_exp_scenario_settings(row, n_seeds, consensus)()
        name = exp_group(row["scenario_name"])
        # every value here is a scalar; a list would be read as a sweep axis
        bad = sorted(k for k, v in fixed.items() if isinstance(v, (list, tuple)))
        if bad:
            raise ValueError(
                f"{name}: frozen parameters must be scalars, got lists for "
                f"{bad} -- those would become a sweep axis.")
        groups.append({"name": name, **fixed})

    names = [g["name"] for g in groups]
    if len(set(names)) != len(names):
        raise ValueError(
            f"the frozen table names a scenario twice: {sorted(names)}")

    def make_settings():
        return ({"_scenario_groups": [dict(g) for g in groups]},
                USER_FIXED_PARAMS.copy(), n_seeds)

    return make_settings


# =============================================================================
# NSR SWEEP RANGES
# =============================================================================
# cal_long reaches a decade lower than the bound pipeline's old 1e-4..1e-2:
# edge_case sheds for 350 days against HIV's 109, so it has far longer to
# accumulate and should calibrate below the others.
# cal_std brackets bound run #7's control, which landed at 1.38e-4 -- the
# closest thing to this stage that pipeline had.
_NSR_STEPS = 10

NSR_RANGES = {
    "cal_long": {"min": 1e-05, "max": 1e-02, "steps": _NSR_STEPS},
    "cal_std":  {"min": 1e-05, "max": 1e-03, "steps": _NSR_STEPS},
}

TABLE_COLUMNS = [
    "scenario_name", "tau_3_long", "long_shedders_ratio", "susceptibility_long",
    "R", "IH_virus_emergence_rate", "R_long",
    "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
    "target_osr_long", "model_type", "long_calib_source", "std_calib_source",
]


def print_fixed_params(fixed_params, label="Fixed parameters"):
    """Single shared format for echoing a stage's fixed_params dict to the
    console/log before it's submitted."""
    print(f"{label}:")
    for k, v in fixed_params.items():
        print(f"  {k}: {v}")


# =============================================================================
# SLURM per-task resource request -- shared by every stage + the orchestrator
# ----------------------------------------------------------------------------
# Read by simplicity.runners.slurm.submit_simulations from the
# SIMPLICITY_SLURM_MEM/SIMPLICITY_SLURM_TIME env vars (same pattern as
# SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM) rather than a function
# argument, so it never has to touch the run_seeded_simulations(...) interface
# shared uniformly by the serial/multiprocessing/slurm runner modules.
#
# The baseline default is set exactly once, in simplicity.dir_manager,
# alongside every other simplicity-wide env var -- not redefined here, so
# there's a single literal source for it. Originally 2G, based on real sacct
# data on the isolated long-NSR calibration grid (population_size=1000,
# final_time=1200): MaxRSS ranged 0.29-0.45G across all 36 grid points --
# real headroom over that, well below the old blanket 8G every job used to
# request regardless of actual need. Raised to 4G after real run #1's much
# wider NSR sweep (see BACKLOG) produced 284 OOM-killed tasks (~18% of the
# run) at 3G -- higher NSR means more cumulative mutation events, which
# means more memory (see the get_lineage_genome/phylogenetic_data growth
# analysis, same BACKLOG entry) -- 2G was tuned for a much narrower range
# and stopped being enough once the range widened. Time limit was bumped
# from 1 day to 2: some high-duration tau groups run well past 1 day of
# wall-clock time to reach final_time and were getting killed by the
# walltime mid-simulation -- see simplicity.runners.slurm.
# reconcile_terminated_tasks, which now detects and signals that case
# regardless, but giving cheap extra headroom avoids the kill (and the
# wasted compute) in the first place.
# =============================================================================
def add_slurm_resource_args(parser):
    """Add --slurm-mem/--slurm-time to an argparse parser (shared wiring so
    every stage script exposes the same two flags identically). Defaults come
    from the env vars dir_manager already set on import, not a second
    hardcoded copy of the same values."""
    parser.add_argument('--slurm-mem', type=str, default=os.environ["SIMPLICITY_SLURM_MEM"],
                        help="Per-task SLURM memory request "
                            f"(default: {os.environ['SIMPLICITY_SLURM_MEM']}).")
    parser.add_argument('--slurm-time', type=str, default=os.environ["SIMPLICITY_SLURM_TIME"],
                        help="Per-task SLURM time limit "
                            f"(default: {os.environ['SIMPLICITY_SLURM_TIME']}).")


def set_slurm_resource_env(mem, time):
    """Set the env vars simplicity.runners.slurm reads. Harmless no-op for
    the serial/multiprocessing runners, which never read them."""
    os.environ["SIMPLICITY_SLURM_MEM"] = mem
    os.environ["SIMPLICITY_SLURM_TIME"] = time
