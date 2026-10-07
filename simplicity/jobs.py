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

"""One repeat = one simulation plus a seed = one unit of work.

This module owns what a repeat IS and what state it is in, so the runners only
decide when to start things.

WHAT REPLACED WHAT

Each repeat used to cost up to seven files under 03_Seeded_simulation_parameters:
a seed_NNNN.json that was its simulation's parameters with `seed` overwritten,
plus .submitted/.released/.started/.completed|.failed/.progress/.launch_retries
beside it, each carrying one bit by existing. The json held nothing that was not
already in the simulation's own parameters file plus an integer, which is the
repeat's own position. Now there is one repeats.json per group and one state file
per repeat.

ORDER IS EXPLICIT, AND THAT IS THE POINT

A Slurm array task finds its repeat by POSITION. The old list came from an
unsorted os.walk over a directory that gains signal files while the run is in
flight, so "which simulation is task 7" depended on readdir order in two
different processes at two different times. Nothing would have reported a
mismatch. Here the order is written down once, in repeats.json, and every reader
uses that.

EVERYTHING INSIDE A GROUP IS NUMBERED FROM 0

Its simulations, its repeats, and each simulation's seeds. The seed is the RNG
seed, so a simulation must get 0..n_seeds-1 wherever it sits or the same
parameters stop reproducing byte-identically. A group's subtree is therefore
self-contained: appending one never reopens another's files.

A SUBMISSION CAN SPAN GROUPS

Dispatch granularity is deliberately NOT the group. impact_long_shedders_cal_2
submits all of its (cell, scenario) groups as one Slurm array, and since
run_seeded_simulations blocks, splitting that into one array per group would
serialise what is currently parallel. So a submission records its own ordered
map from array position to (group, index), and a task resolves itself by looking
its position up in it.

STATE IS ONE FILE PER REPEAT, WRITTEN WHOLE

Per-repeat files are kept deliberately: hundreds of tasks on different nodes
write their own state with no lock and no coordination, which is why the old
scheme survived every monitor bug. What changes is that one file now carries the
state AND the attempt count, instead of several files carrying one bit each.
Every write is temp-then-rename, atomic on POSIX, so a reader never sees half a
state.

PROGRESS IS A SEPARATE FILE, ON PURPOSE

A running simulation rewrites its progress every 30 seconds. Putting that in the
state file would make it a read-modify-write, and a progress write that read the
record before a terminal COMPLETED and landed after it would put STARTED back --
which is exactly the monitor hang that has been fixed twice. Separate files make
that impossible rather than merely unlikely.
"""
import json
import os
import time

import simplicity.dir_manager as dm

# A repeat is only ever in one of these.
SUBMITTED = "submitted"
RELEASED = "released"
STARTED = "started"
COMPLETED = "completed"
FAILED = "failed"

# Terminal: the polling loop stops waiting on a repeat in one of these.
RESOLVED = (COMPLETED, FAILED)

# The lifecycle is monotonic, and set_state enforces it. This also answers "has
# this repeat reached at least X?", which the old scheme answered by having five
# separate signal files that each kept existing.
#
# WHY IT IS ENFORCED RATHER THAN ASSUMED: release_simulations records RELEASED
# for a batch AFTER scontrol returns. With signal files that was harmless --
# touching .released next to an existing .completed changes nothing. With one
# state field, a task that managed to start and finish in that window would have
# its COMPLETED overwritten by RELEASED, and the polling loop would wait on it
# forever. Narrow in production, certain under a synchronous fake controller,
# and the same class as the monitor hangs already fixed twice.
STATE_RANK = {SUBMITTED: 1, RELEASED: 2, STARTED: 3, COMPLETED: 4, FAILED: 4}

REPEATS_FILENAME = "repeats.json"
SUBMISSIONS_DIRNAME = "submissions"


def _write_json(path, payload):
    """Temp-then-rename, so a reader never sees a half-written file."""
    tmp = f"{path}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w") as handle:
            json.dump(payload, handle, indent=1)
        os.replace(tmp, path)
    except OSError:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _read_json(path, default=None):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        # a truncated read means a writer is mid-rename; the next poll sees it
        return default


# ------------------------------------------------------------------- repeats

def repeats_path(experiment_name, group):
    return os.path.join(dm.get_repeats_dir(experiment_name, group),
                        REPEATS_FILENAME)


def write_repeats(experiment_name):
    """Record each group's ordered repeat list, from the experiment settings.

    Simulations in id order, seeds 0..n_seeds-1 within each. This order IS the
    Slurm array mapping, so it is written down rather than rediscovered.
    """
    import simplicity.settings_manager as sm
    simulations = sm.read_simulations(experiment_name)
    written = {}
    for group in sm.read_groups(experiment_name):
        name = group['name']
        records = []
        for simulation in simulations:
            if simulation['group'] != name:
                continue
            stem = sm.simulation_stem(simulation)
            for seed in range(group['n_seeds']):
                records.append({'index': len(records),
                                'simulation': simulation['id'],
                                'stem': stem,
                                'seed': seed})
        _write_json(repeats_path(experiment_name, name), records)
        written[name] = records
    return written


def read_repeats(experiment_name, group):
    """One group's ordered repeats. Index i is the i-th repeat, always."""
    return _read_json(repeats_path(experiment_name, group), default=[])


def get_repeat(experiment_name, group, index):
    return read_repeats(experiment_name, group)[index]


def all_repeats(experiment_name, groups=None):
    """[(group, record), ...] across groups, in group order then index order."""
    import simplicity.settings_manager as sm
    if groups is None:
        groups = [g['name'] for g in sm.read_groups(experiment_name)]
    return [(group, record)
            for group in groups
            for record in read_repeats(experiment_name, group)]


def repeat_label(record):
    """Short human name, also the output sub-directory for this repeat."""
    return f"seed_{record['seed']:04d}"


def repeat_description(group, record):
    return f"{group}/{record['stem']}/{repeat_label(record)}"


# --------------------------------------------------------------------- state

def state_dir(experiment_name, group):
    path = os.path.join(dm.get_repeats_dir(experiment_name, group), "state")
    os.makedirs(path, exist_ok=True)
    return path


def state_path(experiment_name, group, index):
    return os.path.join(state_dir(experiment_name, group), f"{index:06d}.json")


def set_state(experiment_name, group, index, state, **extra):
    """Write a repeat's whole state. Never partial: a reader sees the previous
    state or this one, never a mixture."""
    record = get_state(experiment_name, group, index) or {}
    record.update(extra)
    record["index"] = index
    # never backwards: a late RELEASED must not unseat a COMPLETED
    if STATE_RANK.get(state, 0) >= STATE_RANK.get(record.get("state"), 0):
        record["state"] = state
    record["updated"] = now = time.time()
    # Stamped once, kept by every later write (set_state merges the previous
    # record). The old scheme read a task's wall time from the .started and
    # .completed file mtimes; one file per repeat would otherwise lose the start
    # the moment the terminal state overwrote it.
    if state == STARTED and "started" not in record:
        record["started"] = now
    try:
        _write_json(state_path(experiment_name, group, index), record)
    except OSError:
        pass  # state reporting must never take a simulation down with it
    return record


def get_state(experiment_name, group, index):
    """A repeat's state record, or None if it has none yet."""
    return _read_json(state_path(experiment_name, group, index))


def state_of(experiment_name, group, index):
    """Just the state name, or None."""
    record = get_state(experiment_name, group, index)
    return record.get("state") if record else None


def is_resolved(experiment_name, group, index):
    return state_of(experiment_name, group, index) in RESOLVED


def wall_seconds(experiment_name, group, index):
    """How long this repeat ran, or None if it has not finished.

    updated - started, which is what the .started/.completed mtime difference
    used to give.
    """
    record = get_state(experiment_name, group, index)
    if not record or record.get("state") not in RESOLVED:
        return None
    started = record.get("started")
    updated = record.get("updated")
    if started is None or updated is None or updated < started:
        return None
    return updated - started


def clear_state(experiment_name, group, index):
    for path in (state_path(experiment_name, group, index),
                 progress_path(experiment_name, group, index)):
        try:
            os.unlink(path)
        except OSError:
            pass


def read_all_states(experiment_name, group):
    """{index: record} for every repeat of this group that has reported."""
    states = {}
    try:
        names = os.listdir(state_dir(experiment_name, group))
    except OSError:
        return states
    for name in names:
        if not name.endswith(".json") or ".tmp" in name:
            continue
        try:
            index = int(name[:-len(".json")])
        except ValueError:
            continue
        record = get_state(experiment_name, group, index)
        if record is not None:
            states[index] = record
    return states


def count_states(experiment_name, groups=None):
    """{state: count} across groups -- what the polling loop reports."""
    import simplicity.settings_manager as sm
    if groups is None:
        groups = [g['name'] for g in sm.read_groups(experiment_name)]
    counts = {}
    for group in groups:
        for record in read_all_states(experiment_name, group).values():
            state = record.get("state")
            if state:
                counts[state] = counts.get(state, 0) + 1
    return counts


# ------------------------------------------------------------------ progress

def progress_dir(experiment_name, group):
    path = os.path.join(dm.get_repeats_dir(experiment_name, group), "progress")
    os.makedirs(path, exist_ok=True)
    return path


def progress_path(experiment_name, group, index):
    return os.path.join(progress_dir(experiment_name, group), f"{index:06d}.json")


def get_progress(experiment_name, group, index):
    return _read_json(progress_path(experiment_name, group, index))


# --------------------------------------------------------------- submissions

def submissions_dir(experiment_name):
    path = os.path.join(dm.get_repeats_dir(experiment_name), SUBMISSIONS_DIRNAME)
    os.makedirs(path, exist_ok=True)
    return path


def submission_path(experiment_name, submission_id):
    return os.path.join(submissions_dir(experiment_name), f"{submission_id}.json")


def write_submission(experiment_name, groups=None, only=None):
    """Record one array's ordered map from position to (group, index).

    `only` restricts it to a set of (group, index) pairs -- a resume submits the
    same array shape, so this stays the whole list unless a caller narrows it.
    Returns (submission_id, [[group, index], ...]).
    """
    pairs = [[group, record['index']]
             for group, record in all_repeats(experiment_name, groups)]
    if only is not None:
        wanted = {(g, i) for g, i in only}
        pairs = [p for p in pairs if (p[0], p[1]) in wanted]

    existing = [name for name in os.listdir(submissions_dir(experiment_name))
                if name.endswith('.json')]
    submission_id = f'{len(existing) + 1:03d}'
    _write_json(submission_path(experiment_name, submission_id), pairs)
    return submission_id, pairs


def read_submission(experiment_name, submission_id):
    return _read_json(submission_path(experiment_name, submission_id), default=[])


def resolve_task(experiment_name, submission_id, position):
    """(group, record) for one array position. The whole point of the map: a
    task resolves itself by lookup, not by walking a directory."""
    pairs = read_submission(experiment_name, submission_id)
    group, index = pairs[position]
    return group, get_repeat(experiment_name, group, index)
