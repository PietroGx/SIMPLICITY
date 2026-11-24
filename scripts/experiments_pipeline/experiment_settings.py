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
# impact_long_shedders -- EXPERIMENT SETTINGS FILES
# ----------------------------------------------------------------------------
# One JSON file in scripts/experiments_settings/ holds everything a pipeline
# run is launched with. impact_long_shedders_config reads it once, at import,
# from $SIMPLICITY_EXPERIMENT_SETTINGS or DEFAULT_SETTINGS, and every stage
# builds from those values. The pipeline writes the file it ran with, CLI
# overrides applied, into the run's setup_data/ and points the env var there,
# so every stage subprocess and the sanity job read exactly what was recorded.
#
# Imports nothing from the config: the pipeline must be able to set the env
# var before the config loads.
# ============================================================================

import json
import os
import re

SCHEMA = 2
SETTINGS_ENV = "SIMPLICITY_EXPERIMENT_SETTINGS"
REPO_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
SETTINGS_DIR = os.path.join(REPO_DIR, "scripts", "experiments_settings")
DEFAULT_SETTINGS = os.path.join(SETTINGS_DIR, "impact_long_shedders_N3000_T1095.json")
RUN_SETTINGS_FILENAME = "experiment_settings.json"
CONSENSUS_MODES = ("distribution", "argmax")

# top-level key -> nested keys, None for a scalar. Unknown keys are refused:
# a misspelt key would otherwise leave the value it meant to set at default.
KEYS = {
    "schema": None,
    "description": None,
    "population": ("population_size", "infected_individuals_at_start", "R",
                   "IH_virus_emergence_rate", "final_time"),
    "scenarios": None,
    "cal_1": ("R", "final_time", "min_cycles", "nsr_long"),
    "cal_2": ("nsr_std",),
    "targets": ("osr_std", "osr_long"),
    "seeds": ("cal", "exp"),
    "consensus": None,
    "slurm": ("mem", "mem_cal2", "time"),
}
SCENARIO_KEYS = ("name", "inf_duration_long", "long_shedders_ratio", "R_long",
                 "susceptibility_long")
RANGE_KEYS = ("min", "max", "steps")
_NAME = re.compile(r"^[A-Za-z0-9]+(_[A-Za-z0-9]+)*$")
_MEM = re.compile(r"^[1-9][0-9]*[KMGT]?$")
_TIME = re.compile(r"^([0-9]+-)?[0-9]{1,2}(:[0-9]{2}){0,2}$")


class SettingsError(ValueError):
    pass


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def validate(s):
    """Every problem in one pass, as a list of messages; [] when valid."""
    errs = []

    def obj(d, keys, where):
        if not isinstance(d, dict):
            errs.append(f"{where}: expected an object")
            return False
        for k in keys:
            if k not in d:
                errs.append(f"{where}: missing '{k}'")
        for k in d:
            if k not in keys:
                errs.append(f"{where}: unknown key '{k}'")
        return all(k in d for k in keys)

    def num(v, where, lo=None, hi=None, integer=False, strict=True):
        if not (_is_int(v) if integer else _is_num(v)):
            errs.append(f"{where}: expected {'an integer' if integer else 'a number'}, got {v!r}")
            return False
        if lo is not None and (v <= lo if strict else v < lo):
            errs.append(f"{where}: must be {'>' if strict else '>='} {lo}, got {v}")
            return False
        if hi is not None and v > hi:
            errs.append(f"{where}: must be <= {hi}, got {v}")
            return False
        return True

    def nsr_range(r, where):
        if obj(r, RANGE_KEYS, where):
            ok = num(r["min"], f"{where}.min", 0) & num(r["max"], f"{where}.max", 0)
            num(r["steps"], f"{where}.steps", 2, integer=True, strict=False)
            if ok and r["max"] <= r["min"]:
                errs.append(f"{where}: max ({r['max']}) must exceed min ({r['min']})")

    if not obj(s, tuple(KEYS), "settings"):
        if not isinstance(s, dict):
            return errs
    if s.get("schema") != SCHEMA:
        errs.append(f"schema: this code reads schema {SCHEMA}, got {s.get('schema')!r}")
    if "description" in s and not isinstance(s["description"], str):
        errs.append("description: expected a string")

    p = s.get("population")
    if obj(p, KEYS["population"], "population"):
        if num(p["population_size"], "population.population_size", 1, integer=True, strict=False) \
                and num(p["infected_individuals_at_start"],
                        "population.infected_individuals_at_start", 1, integer=True,
                        strict=False) \
                and p["infected_individuals_at_start"] >= p["population_size"]:
            errs.append("population.infected_individuals_at_start must be below population_size")
        num(p["R"], "population.R", 0)
        num(p["IH_virus_emergence_rate"], "population.IH_virus_emergence_rate", 0)
        num(p["final_time"], "population.final_time", 0, integer=True)

    sc = s.get("scenarios")
    if not isinstance(sc, list) or not sc:
        errs.append("scenarios: expected a non-empty list")
    else:
        names = []
        for i, x in enumerate(sc):
            where = f"scenarios[{i}]"
            if not obj(x, SCENARIO_KEYS, where):
                continue
            name = x["name"]
            where = f"scenarios[{i}] '{name}'"
            if not isinstance(name, str) or not _NAME.match(name):
                errs.append(f"{where}: name must be letters/digits joined by single "
                            f"underscores (it becomes a directory and part of a group name)")
            names.append(name)
            if not num(x["long_shedders_ratio"], f"{where}.long_shedders_ratio", 0, 1,
                       strict=False):
                continue
            num(x["susceptibility_long"], f"{where}.susceptibility_long", 0)
            if x["long_shedders_ratio"] == 0:
                for k in ("inf_duration_long", "R_long"):
                    if x[k] is not None:
                        errs.append(f"{where}: {k} must be null when long_shedders_ratio is 0")
            else:
                num(x["inf_duration_long"], f"{where}.inf_duration_long", 0)
                num(x["R_long"], f"{where}.R_long", 0)
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            errs.append(f"scenarios: duplicate names {dup}")

    c1 = s.get("cal_1")
    if obj(c1, KEYS["cal_1"], "cal_1"):
        num(c1["R"], "cal_1.R", 0)
        num(c1["final_time"], "cal_1.final_time", 0, integer=True)
        num(c1["min_cycles"], "cal_1.min_cycles", 1, integer=True, strict=False)
        nsr_range(c1["nsr_long"], "cal_1.nsr_long")
    c2 = s.get("cal_2")
    if obj(c2, KEYS["cal_2"], "cal_2"):
        nsr_range(c2["nsr_std"], "cal_2.nsr_std")
    t = s.get("targets")
    if obj(t, KEYS["targets"], "targets"):
        num(t["osr_std"], "targets.osr_std", 0)
        num(t["osr_long"], "targets.osr_long", 0)
    sd = s.get("seeds")
    if obj(sd, KEYS["seeds"], "seeds"):
        num(sd["cal"], "seeds.cal", 1, integer=True, strict=False)
        num(sd["exp"], "seeds.exp", 1, integer=True, strict=False)
    if "consensus" in s and s["consensus"] not in CONSENSUS_MODES:
        errs.append(f"consensus: one of {list(CONSENSUS_MODES)}, got {s['consensus']!r}")
    sl = s.get("slurm")
    if obj(sl, KEYS["slurm"], "slurm"):
        if not (isinstance(sl["mem"], str) and _MEM.match(sl["mem"])):
            errs.append(f"slurm.mem: a Slurm size such as '3G', got {sl['mem']!r}")
        if sl["mem_cal2"] is not None and not (isinstance(sl["mem_cal2"], str)
                                               and _MEM.match(sl["mem_cal2"])):
            errs.append(f"slurm.mem_cal2: a Slurm size or null, got {sl['mem_cal2']!r}")
        if sl["time"] is not None and not (isinstance(sl["time"], str)
                                           and _TIME.match(sl["time"])):
            errs.append(f"slurm.time: D-HH:MM:SS or null, got {sl['time']!r}")
    return errs


# real-valued fields, stored as floats whatever the JSON spelt: a browser
# writes 1.0 as 1, and an int R would name output directories differently
FLOATS = {"population": ("R", "IH_virus_emergence_rate"),
          "cal_1": ("R",), "targets": ("osr_std", "osr_long")}
SCENARIO_FLOATS = ("inf_duration_long", "long_shedders_ratio", "R_long",
                   "susceptibility_long")


def normalize(s):
    """A deep copy of a VALID settings dict with every real-valued field a float."""
    s = json.loads(json.dumps(s))
    for sec, keys in FLOATS.items():
        for k in keys:
            s[sec][k] = float(s[sec][k])
    for x in s["scenarios"]:
        for k in SCENARIO_FLOATS:
            if x[k] is not None:
                x[k] = float(x[k])
    for r in (s["cal_1"]["nsr_long"], s["cal_2"]["nsr_std"]):
        r["min"], r["max"] = float(r["min"]), float(r["max"])
    return s


def check(s, source="settings"):
    errs = validate(s)
    if errs:
        raise SettingsError(f"{source}: {len(errs)} problem(s)\n  " + "\n  ".join(errs))
    return normalize(s)


def load(path):
    with open(path) as fh:
        try:
            s = json.load(fh)
        except json.JSONDecodeError as exc:
            raise SettingsError(f"{path}: not valid JSON ({exc})") from None
    return check(s, path)


def dumps(s):
    return json.dumps(s, indent=2) + "\n"


def write(s, path):
    s = check(s, path)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(dumps(s))
    os.replace(tmp, path)
    return path


def settings_path():
    """The settings file this process runs with."""
    return os.path.abspath(os.environ.get(SETTINGS_ENV) or DEFAULT_SETTINGS)


def activate(path):
    """Point this process and its children at `path`. Must run before
    impact_long_shedders_config is first imported."""
    path = os.path.abspath(path)
    load(path)
    os.environ[SETTINGS_ENV] = path
    return path


def activate_from_argv(argv):
    """`--settings PATH` / `--settings=PATH` from a script's argv, before its
    argparse runs: the config reads the env var when it is imported."""
    for i, a in enumerate(argv):
        if a == "--settings" and i + 1 < len(argv):
            return activate(argv[i + 1])
        if a.startswith("--settings="):
            return activate(a.split("=", 1)[1])
    return None


def reservoir(population_size):
    """Individuals a simulation of this size allocates (simplicity/population.py)."""
    from simplicity.population import RESERVOIR_PER_CAPITA
    return RESERVOIR_PER_CAPITA * population_size
