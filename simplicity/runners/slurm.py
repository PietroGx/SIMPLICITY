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
@author: jbescudie


Slurm will create new processes (on possibily other hosts).
We need to let these processes load user's "run_seeded_simulation" function.

Given the user's provided function
1. we determine a python importable path
2. in processes started by Slurm, the path is imported

Caveat:
   this imposes the following restriction: the user must import the run_seeded_simulation
   function before calling simplicity.runners.slurm.run_seeded_simulations(...).
   A test verifies this and raises before submitting jobs to Slurm.

Kown bugs:
1. if the run_seeded_function is importable when submitting, but not in the processes
   started by Slurm. They will fail just after the repeat is recorded STARTED.

2. Interrupting the program before completion will leave held jobs in Slurm's queue,
   and leaving the Data directory in a incomplete state (like for other
   simplicity.runners implementations).

   In this case, manually cleaning Slurm (hint: squeue) and the Data directory
   (removing or archiving the experiment folder) is highly recommanded before retrying.

3. A task killed directly by Slurm (walltime/--time or memory/--mem limit)
   never gets a chance to run job()'s except block and record FAILED itself
   -- the process is just terminated. Left unresolved, poll_simulations_status
   would count that task as permanently "left" and run_seeded_simulations'
   polling loop would never see status.left reach 0. Handled by
   reconcile_terminated_tasks, called periodically (RECONCILE_INTERVAL_S)
   from run_seeded_simulations: it cross-checks STARTED-but-unresolved repeats
   against `sacct` (which retains terminal-state history after a job leaves
   squeue) and records FAILED on Slurm's behalf for any task Slurm reports
   as TIMEOUT/OUT_OF_MEMORY/CANCELLED/etc.

"""
import typing, os, json, time, pathlib, subprocess, platform
import simplicity.dir_manager as dm
import simplicity.jobs as jobs

# How often (seconds) run_seeded_simulations' polling loop reports the
# internal state (time/final_time/infected) of simulations that have been
# running for a while -- separate from the SimulationsStatus line, which
# prints on change; this is meant as an occasional deeper look.
LONG_RUNNING_REPORT_INTERVAL_S = 600
# How long a simulation must have been running (since it was recorded STARTED)
# before it's included in that report. Was 3600, at which nothing ever
# qualified: the slowest of unbound #1's 1,400 tasks took 2,423s, and the
# v2.4.41-43 speedups cut runtimes a further ~6x. The report prints nothing
# at all when nothing qualifies, so that read as the feature being broken.
LONG_RUNNING_THRESHOLD_S = 300

# How often (seconds) run_seeded_simulations' polling loop reconciles
# STARTED-but-unresolved repeats against Slurm's own accounting (sacct). This
# is what unblocks a run stuck on an externally-killed task (see
# reconcile_terminated_tasks below) -- unlike the purely-informational
# LONG_RUNNING_REPORT above, so it runs on its own, shorter timer.
RECONCILE_INTERVAL_S = 900

# Seconds to wait on any Slurm command issued from the polling loop. Without
# one, an unresponsive slurmdbd makes subprocess.run block forever and the loop
# stops entirely -- no status line, no [long-running] report, no reconciliation,
# just silence. Seen on the calibrated grid: cell #17 went quiet at 01:46 with
# 95 started / 61 completed / 34 unresolved and an empty queue, and was still
# frozen 6.5 hours later. Every one of these is advisory; timing out and
# retrying on the next pass is always better than waiting forever.
SLURM_QUERY_TIMEOUT_S = 120

# sacct job states that mean "Slurm itself terminated this task" -- i.e. the
# task's own process never got a chance to record completed/failed (job()'s
# except block only runs on a catchable Python exception, not a SIGTERM/
# SIGKILL from Slurm hitting --time/--mem). RUNNING/PENDING/COMPLETING are
# deliberately excluded (still in flight); COMPLETED is excluded too (job()
# should have already recorded completed itself in that case).
SLURM_TERMINAL_FAILURE_STATES = {
    "TIMEOUT", "OUT_OF_MEMORY", "FAILED", "NODE_FAIL", "CANCELLED",
    "DEADLINE", "PREEMPTED", "BOOT_FAIL",
}

# The other way a task can be finished-but-unresolved. Slurm reports COMPLETED
# only once the job script has exited 0, and job() records completed before
# returning -- so COMPLETED with no completed state means the write did not
# land or is not yet visible on the shared filesystem, NOT that the simulation
# failed. Seen on profile_grid_#910: 116 completed signals for 120 tasks, all
# 120 reported COMPLETED by sacct and all 120 having written a full profile.
# Nothing matched these, so `left` never reached 0 and the polling loop ran for
# hours against a stale progress snapshot. Reconciled as completed rather than
# failed: marking an exit-0 simulation failed would under-report the run.
SLURM_TERMINAL_SUCCESS_STATES = {"COMPLETED"}

# Launch-failure retry/reconciliation -- distinct from reconcile_terminated_tasks
# above: applies to tasks that never reached job() at all (Slurm failed to
# launch them on their assigned node and auto-requeued them into a held
# state, e.g. squeue reason "(launch failed requeued held)"). Checked more
# often than the OOT/OOM reconciliation since a re-release is cheap and the
# whole point is prompt recovery from what's usually a transient node issue.
LAUNCH_FAILURE_RECONCILE_INTERVAL_S = 120
LAUNCH_FAILURE_MAX_RETRIES = 3
# Substring match against squeue's free-text Reason column (case-insensitive)
# -- Slurm doesn't expose this as a stable enum the way job State is.
LAUNCH_FAILURE_REASON_MARKERS = ("launch failed",)


def expand_array_task_ids(field):
    """Task ids from one squeue ArrayTaskID field.

    Slurm collapses an array's pending tasks that share a state and a reason
    into ONE squeue row, and writes the task field as a list expression:
    "5", "1-3", "1,3,5", "1-3,7", sometimes with a "%" throttle suffix.

    This was `int(task_id_str)` inside a try/except that skipped anything else.
    Measured on slurm 26.05.4 (tests/check_slurm_interface.py): a three-task
    held array reports a single row with ArrayTaskID "1-3", so every pending
    task was silently skipped. That is the only state a launch-failed task is
    ever in, and such a task never reaches STARTED either, so
    reconcile_terminated_tasks never sees it -- the polling loop would wait on
    it forever. The collapsing is not new in 26.05; the parser was written
    against a shape Slurm only produces when exactly one task carries a given
    reason.
    """
    for part in field.split(","):
        part = part.strip().partition("%")[0]      # drop any throttle suffix
        if not part:
            continue
        low, dash, high = part.partition("-")
        try:
            if dash:
                yield from range(int(low), int(high) + 1)
            else:
                yield int(low)
        except ValueError:
            continue

class SimulationsStatus(typing.NamedTuple):
    total    : int
    submitted: int
    released : int
    left     : int
    pending  : int
    started  : int
    running  : int
    completed: int
    failed   : int

    
def print_simulations_status(status, experiment_name=None):
    """One timestamped status line. Printed on every change of status and
    never otherwise: a fixed-interval reprint of an unchanged line buries
    the transitions that actually carry information. The experiment name is
    on the line because several experiments can now poll concurrently (see
    impact_long_shedders_unbound_exp.dispatch_all) and their output
    interleaves."""
    tag = f" {experiment_name}" if experiment_name else ""
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]{tag} {status}")


def get_platform_executable_extension():
    """file extension to use when calling Slurm command-line utilities (sbatch, squeue, scontrol)."""
    return ".exe" if platform.system() == "Windows" else  ""


# Set to "1" to resume an experiment instead of refusing to touch it: every
# repeat that already reached COMPLETED is kept and skipped, and every other one
# is cleared back to unstarted and re-run. Written by a --rerun flag rather than
# read from one, because the pipeline reaches this code through three layers of
# subprocess and an env var is the only thing that survives.
RESUME_ENV = "SIMPLICITY_RESUME"

# The submission whose task map a Slurm task resolves itself through. A
# submission may span several groups, so a task cannot infer its repeat from the
# group alone.
SUBMISSION_ENV = "SIMPLICITY_SUBMISSION"


def resuming():
    return os.environ.get(RESUME_ENV) == "1"


def prepare_resume(experiment_name, groups=None):
    """Keep every finished repeat, clear the rest so they run again.

    A FAILED repeat is cleared too: on a resume it is work still to do, not a
    settled answer. Returns (kept, to_run) for the caller to report.
    """
    kept = to_run = 0
    for group, record in jobs.all_repeats(experiment_name, groups):
        if jobs.state_of(experiment_name, group, record["index"]) == jobs.COMPLETED:
            kept += 1
            continue
        to_run += 1
        jobs.clear_state(experiment_name, group, record["index"])
    return kept, to_run


def find_stale_state(experiment_name, groups=None):
    """{state: count} of repeat states already on disk for this experiment.

    poll_simulations_status reads these, not Slurm. A set left behind by an
    earlier attempt at the same --exp-num therefore decides the run before it
    starts: with COMPLETED present, `left` is 0 on the first poll, the polling
    loop never runs, release_simulations is never called -- it is only reached
    from inside that loop -- and the held array sits in the queue forever while
    the pipeline walks on to a fit with no data behind it. Nothing deletes these
    files, and a Data tree copied between machines carries them along without
    the outputs.
    """
    return jobs.count_states(experiment_name, groups)


def raise_on_stale_state(experiment_name, groups=None):
    """Refuse to submit over an earlier attempt's state.

    Checked BEFORE sbatch, so a refused run leaves no held array behind. This
    never deletes anything: which of the two runs matters is not something this
    code can know.
    """
    if resuming():
        kept, to_run = prepare_resume(experiment_name, groups)
        print(f"[resume] {experiment_name}: keeping {kept} completed "
              f"repeat(s), running {to_run}")
        return
    counts = find_stale_state(experiment_name, groups)
    if not counts:
        return
    summary = ", ".join(f"{n}x {state}" for state, n in sorted(counts.items()))
    repeats_dir = dm.get_repeats_dir(experiment_name)
    raise RuntimeError(
        f"{experiment_name} already carries repeat state from an earlier "
        f"attempt: {summary}.\n"
        f"Status is read from these files, so submitting now would either "
        f"skip the run entirely (if any are completed) or miscount it. "
        f"Nothing has been submitted.\n"
        f"Either run at an --exp-num with nothing under it, or clear them:\n"
        f"  rm -rf {os.path.join(repeats_dir, '*', 'state')}")


def submit_simulations(experiment_name: str,
                       run_seeded_simulation: typing.Callable,
                       submission_id: str,
                       pairs: list):
    """sbatch one held array for this submission's (group, index) pairs."""
    # run_seeded_simulation to qualname
    fn_name = run_seeded_simulation.__name__
    run_seeded_simulation_module_qualname = run_seeded_simulation.__module__
    if run_seeded_simulation_module_qualname == "__main__":
        import simplicity.runners.slurm
        help(simplicity.runners.slurm)
        raise Exception("run_seeded_simulation must be imported.")
    run_seeded_simulation_qualname = run_seeded_simulation_module_qualname + "." + fn_name
    
    # job is to run the same python to call this module main() (which in turn will call user's run_seeded_simulation function)
    import sys
    stdin = "\n".join((
        # run the same python
        f"#!{sys.executable}",
        # call this module main()
        "import simplicity.runners.slurm",
        "simplicity.runners.slurm.job()",
    )).encode()
    
    # ! slurm array indexing is 1 based. Position p in the submission map is
    # array task p+1, which is the only mapping job() needs.
    batch_start = 1
    batch_end   = len(pairs)
   
    # Define the output and error file paths
    
    slurm_logs_dir = dm.get_slurm_logs_dir(experiment_name)
    
    output_file = f"{slurm_logs_dir}/{experiment_name}-%A_%a.out"  # %A = job ID, %a = array index
    error_file  = f"{slurm_logs_dir}/{experiment_name}-%A_%a.err"  # %A = job ID, %a = array index
    
    # send logs to /dev/null 
    # output_file = "/dev/null"
    # error_file  = "/dev/null"
        
    # Per-task resource request. Configurable via env var (set by the calling
    # script from a --slurm-mem/--slurm-time CLI flag) rather than a function
    # argument, so this doesn't touch the run_seeded_simulations(experiment_name,
    # run_seeded_simulation) interface shared uniformly by serial/multiprocessing/
    # slurm -- same pattern as SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM.
    # Baseline default is set once in dir_manager (imported above); a hard
    # subscript (not .get) matches how that other env var is read, and relies
    # on dir_manager always having run first to set it.
    max_runtime = os.environ["SIMPLICITY_SLURM_TIME"]
    mem_request = os.environ["SIMPLICITY_SLURM_MEM"]

    # submit the job
    slurm_process = subprocess.run((args:=[
        # calls sbatch
        "sbatch" + get_platform_executable_extension(),
        # to create the job array on hold
        f"--array={batch_start}-{batch_end}",
        "--hold",
        f"--time={max_runtime}",
        f"--mem={mem_request}",
        # with a name (used later for lookup)
        f"--job-name={experiment_name}",
        f"--output={output_file}",  # Specify the output file path
        f"--error={error_file}",  # Specify the error file path
    ]), env=(env:={
        **os.environ,
        "SIMPLICITY_EXPERIMENT_NAME": experiment_name,
        SUBMISSION_ENV: submission_id,
        "USER_RUN_SEEDED_SIMULATION": run_seeded_simulation_qualname
        # "PLOT_TRAJECTORY"           : str(plot_trajectory)
    }), input=stdin)
    assert slurm_process.returncode == 0, f"Slurm was called with the following arguments:\n{' '.join(args)}\n{env}\n=== stdin\n{stdin}\n=== /stdin"
    
    for group, index in pairs:
        jobs.set_state(experiment_name, group, index, jobs.SUBMITTED)


def poll_simulations_status(experiment_name, pairs):
    states = [jobs.state_of(experiment_name, group, index)
              for group, index in pairs]
    ranks = [jobs.STATE_RANK.get(state, 0) for state in states]
    reached = lambda state: sum(r >= jobs.STATE_RANK[state] for r in ranks)
    completed = states.count(jobs.COMPLETED)
    failed = states.count(jobs.FAILED)
    submitted = reached(jobs.SUBMITTED)
    released = reached(jobs.RELEASED)
    started = reached(jobs.STARTED)
    return SimulationsStatus(
        total    = len(pairs),
        submitted= submitted,
        released = released,
        left     = submitted - (completed + failed),
        pending  = released - started,
        started  = started,
        running  = started - (completed + failed),
        completed= completed,
        failed   = failed,
    )


def report_long_running_simulations(experiment_name, pairs,
                                   threshold_seconds=LONG_RUNNING_THRESHOLD_S):
    """Print the internal state (time/final_time/infected) of every repeat that
    has been running for more than threshold_seconds.

    State is read from the periodic progress snapshot
    simplicity.extrande.ProgressReporter writes to jobs.progress_path, so this
    works without a live terminal attached to the simulation process. Elapsed
    time comes from the state record's `updated`, which for a STARTED repeat is
    when it started -- progress is a separate file, so it does not move it.
    """
    now = time.time()
    for group, index in pairs:
        record = jobs.get_state(experiment_name, group, index)
        if not record or record.get("state") != jobs.STARTED:
            continue
        elapsed = now - record.get("updated", now)
        if elapsed < threshold_seconds:
            continue

        repeat = jobs.get_repeat(experiment_name, group, index)
        name = jobs.repeat_description(group, repeat)
        snapshot = jobs.get_progress(experiment_name, group, index)
        if not snapshot:
            print(f"[long-running] {name}: running {elapsed/3600:.1f}h, "
                  f"no progress snapshot yet")
            continue
        print(f"[long-running] {name}: running {elapsed/3600:.1f}h, "
              f"time={snapshot.get('time')}/{snapshot.get('final_time')}, "
              f"infected={snapshot.get('infected')}")


def _build_slurm_id_map_index(experiment_name):
    """Reverse-index job()'s slurm_id_map_dir CSV files (job()-written, never
    previously read): {(group, index): (job_id, task_id)}.
    Filenames are "{experiment_name}_{job_id}_{task_id}.csv"; content is
    "{group}/{index}"."""
    index = {}
    map_dir = pathlib.Path(dm.get_slurm_id_map_dir(experiment_name))
    prefix = f"{experiment_name}_"
    for map_file in map_dir.glob(f"{prefix}*.csv"):
        stem = map_file.stem[len(prefix):]
        job_id, _, task_id = stem.rpartition("_")
        if not job_id:
            continue
        try:
            content = map_file.read_text().strip()
        except OSError:
            continue
        group, _, index_str = content.rpartition("/")
        if not group:
            continue
        try:
            index[(group, int(index_str))] = (job_id, task_id)
        except ValueError:
            continue
    return index


def reconcile_terminated_tasks(experiment_name, pairs):
    """
    Find every repeat in STARTED that never reached COMPLETED/FAILED, and check
    Slurm's own accounting (sacct) for whether Slurm itself already killed that
    task (TIMEOUT/OUT_OF_MEMORY/CANCELLED/...).

    This exists because job()'s try/except (simplicity.runners.slurm.job,
    below) can only record FAILED on a catchable Python exception -- a Slurm
    walltime or memory kill terminates the process directly, so job() never
    gets to run its except block, and that repeat would otherwise stay STARTED
    forever. poll_simulations_status's `left` count then never reaches 0 for it,
    and run_seeded_simulations' polling loop can only exit through the unrelated
    "no held task" exception in release_simulations. Unlike a live task's own
    job() process, sacct retains terminal-state history after a job leaves
    squeue, so it's the only reliable place to check this from outside the
    (now-dead) task's own process.

    Any task confirmed terminated by Slurm is recorded FAILED on its behalf,
    with a note explaining it was an externally-detected kill (not a
    Python-level failure) -- this is what actually unblocks the polling loop.

    A task Slurm reports as COMPLETED is recorded COMPLETED instead: the job
    script exited 0, so the simulation ran and its output is on disk, and only
    the state write is missing. That case blocks the loop exactly as hard as a
    kill does, and used to match nothing here at all.
    """
    stuck = [(group, index) for group, index in pairs
             if jobs.state_of(experiment_name, group, index) == jobs.STARTED]
    if not stuck:
        return

    id_index = _build_slurm_id_map_index(experiment_name)

    # Group stuck repeats by Slurm array job id, so each distinct job id is
    # queried via sacct exactly once (normally there's only one, but nothing
    # here assumes that).
    by_job_id = {}
    for key in stuck:
        ids = id_index.get(key)
        if ids is None:
            # job() hasn't written its map file yet (narrow window right after
            # STARTED, before the map file write) -- not stuck, just not
            # checkable yet. Skip silently; found on a later pass.
            continue
        job_id, task_id = ids
        by_job_id.setdefault(job_id, {})[f"{job_id}_{task_id}"] = key

    for job_id, task_id_to_key in by_job_id.items():
        try:
            sacct_process = subprocess.run([
                "sacct" + get_platform_executable_extension(),
                    "-j", job_id, "--format=JobID,State", "--noheader", "--parsable2", "-X",
            ], stdout=subprocess.PIPE, timeout=SLURM_QUERY_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            print(f"[reconciled] sacct timed out after {SLURM_QUERY_TIMEOUT_S}s "
                  f"for job {job_id}; retrying on the next pass")
            continue
        if sacct_process.returncode != 0:
            # sacct itself unreachable/erroring -- try again on the next
            # reconcile pass rather than guessing at these repeats' state now.
            continue

        for line in sacct_process.stdout.decode().splitlines():
            if "|" not in line:
                continue
            sacct_job_id, state = line.split("|", maxsplit=1)
            key = task_id_to_key.get(sacct_job_id.strip())
            if key is None:
                continue
            # sacct states can carry a suffix, e.g. "CANCELLED by 12345".
            state = state.strip().split()[0] if state.strip() else state.strip()
            if state not in SLURM_TERMINAL_FAILURE_STATES \
                    and state not in SLURM_TERMINAL_SUCCESS_STATES:
                continue
            group, index = key
            # Re-check immediately before writing: `stuck` was built at the top
            # of this function, which on a large grid is thousands of reads ago,
            # and on a shared filesystem the task's own state may have become
            # visible in between. Without this a finished repeat could be
            # overwritten by a stale verdict.
            if jobs.is_resolved(experiment_name, group, index):
                continue
            name = jobs.repeat_description(
                group, jobs.get_repeat(experiment_name, group, index))
            if state in SLURM_TERMINAL_SUCCESS_STATES:
                jobs.set_state(experiment_name, group, index, jobs.COMPLETED,
                               reconciled=f"slurm:{state}")
                print(f"[reconciled] {name}: Slurm reports {state} -- "
                      f"marking completed (task finished but never signaled "
                      f"itself; its output is on disk)")
            else:
                jobs.set_state(experiment_name, group, index, jobs.FAILED,
                               reconciled=f"slurm:{state}")
                print(f"[reconciled] {name}: Slurm reports {state} -- "
                      f"marking failed (task never signaled itself)")


def reconcile_launch_failures(experiment_name, pairs):
    """
    Detect tasks Slurm auto-requeued into a held state after failing to launch
    them (e.g. squeue reason "(launch failed requeued held)") -- distinct from
    reconcile_terminated_tasks above: these tasks never reached job() at all
    (still RELEASED, never STARTED), so there is no job-ID mapping file to look
    them up by either (job() only writes that after STARTED). Array task IDs
    map directly to the submission map's positions instead (1-based), the same
    mapping submit_simulations and job() already rely on.

    Re-releases each affected task up to LAUNCH_FAILURE_MAX_RETRIES times
    (counted in the repeat's own state record); past that, records it FAILED so
    the polling loop stops waiting on it indefinitely.
    """
    # Candidates: already released once (Slurm tried to launch it) but never
    # started. Position in `pairs` is the array task id, minus one.
    candidates = {}  # 1-based array task id -> (group, index)
    for position, (group, index) in enumerate(pairs):
        if jobs.state_of(experiment_name, group, index) == jobs.RELEASED:
            candidates[position + 1] = (group, index)

    if not candidates:
        return

    try:
        squeue_process = subprocess.run([
            "squeue" + get_platform_executable_extension(),
                "--name", experiment_name,
                "--Format=ArrayJobID,ArrayTaskID,Reason", "--noheader",
        ], stdout=subprocess.PIPE, timeout=SLURM_QUERY_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f"[reconciled] squeue timed out after {SLURM_QUERY_TIMEOUT_S}s; "
              f"retrying on the next pass")
        return
    if squeue_process.returncode != 0:
        # squeue itself unreachable/erroring -- try again on the next pass.
        return

    to_release = {}  # "{job_id}_{task_id}" -> ((group, index), retry count)
    for line in squeue_process.stdout.decode().splitlines():
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        job_id, task_id_str, reason = parts
        if not any(marker in reason.lower() for marker in LAUNCH_FAILURE_REASON_MARKERS):
            continue
        # one row can cover several tasks -- see expand_array_task_ids
        for task_id in expand_array_task_ids(task_id_str):
            key = candidates.get(task_id)
            if key is None:
                continue

            group, index = key
            record = jobs.get_state(experiment_name, group, index) or {}
            retries = record.get("attempts", 0)
            name = jobs.repeat_description(
                group, jobs.get_repeat(experiment_name, group, index))

            if retries >= LAUNCH_FAILURE_MAX_RETRIES:
                jobs.set_state(experiment_name, group, index, jobs.FAILED,
                               reconciled="launch_failed")
                print(f"[reconciled] {name}: launch failed {retries} times "
                     f"(Slurm reason: {reason.strip()}) -- giving up, marking failed")
                continue

            to_release[f"{job_id}_{task_id}"] = (key, retries)

    if not to_release:
        return

    job_list = ",".join(to_release.keys())
    try:
        scontrol_process = subprocess.run([
            "scontrol" + get_platform_executable_extension(),
                "release", job_list
        ], timeout=SLURM_QUERY_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f"[reconciled] scontrol release timed out after "
              f"{SLURM_QUERY_TIMEOUT_S}s for {job_list}; retrying next pass")
        return
    if scontrol_process.returncode != 0:
        print(f"[reconciled] scontrol release failed for launch-failed tasks: {job_list}")
        return

    for (group, index), retries in to_release.values():
        jobs.set_state(experiment_name, group, index, jobs.RELEASED,
                       attempts=retries + 1)
        name = jobs.repeat_description(
            group, jobs.get_repeat(experiment_name, group, index))
        print(f"[reconciled] {name}: launch failed, re-released "
             f"(attempt {retries + 1}/{LAUNCH_FAILURE_MAX_RETRIES})")


def release_simulations(experiment_name, pairs, n: int):
    # up to n repeats that were submitted but never released. Position in
    # `pairs` is the array task id, minus one.
    to_release = {}  # 1-based array task id -> (group, index)
    for position, (group, index) in enumerate(pairs):
        if jobs.state_of(experiment_name, group, index) == jobs.SUBMITTED:
            to_release[position + 1] = (group, index)
        if len(to_release) >= n:
            break

    if not to_release:
        # Nothing actually needs releasing right now. This happens once every
        # task has already been released at least once, while `n` (the caller's
        # estimate of remaining capacity) can still come out > 0 because
        # status.left also counts repeats stuck on an externally-killed Slurm
        # task that hasn't been reconciled yet (see reconcile_terminated_tasks)
        # -- those inflate `left` without inflating `pending`. Querying squeue
        # here would be pointless (nothing to release) and, once the whole array
        # job has aged out of squeue's listing, would incorrectly raise "no held
        # task" for a totally benign state. Just return; the reconciler is what
        # actually clears this up.
        return

    # slurm find array job id from job name
    try:
        slurm_process = subprocess.run([
            "squeue" + get_platform_executable_extension(),
                "--Format=ArrayJobID", f"--name={experiment_name}" 
        ], stdout=subprocess.PIPE, timeout=SLURM_QUERY_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        # unlike the reconcilers, this one is on the release path: returning
        # leaves the tasks held, and the next poll tries again.
        print(f"[release] squeue timed out after {SLURM_QUERY_TIMEOUT_S}s; "
              f"leaving these tasks held for the next pass")
        return
    assert slurm_process.returncode == 0
    array_job_id_set = set(line.strip() for line in slurm_process.stdout.decode().splitlines(keepends=False)[1:])
    if len(array_job_id_set) == 0:
        raise Exception("no held task in the job array. Hint: check slurm log output for error occurring before the started state ")
    assert len(array_job_id_set) == 1, f"Expect exactly one job array with name {experiment_name}. If several slurm array job share the name {experiment_name}. please clean slurm's queue before continuing.\nKnown bug if Slurm's hold a job in CG state.\nHint: better Slurm's squeue parsing in 'simplicity.runners.slurm' may resolve CG status case.\n" + slurm_process.stdout.decode()
    SLURM_ARRAY_JOB_ID = next(iter(array_job_id_set))

    job_list = ",".join(f"{SLURM_ARRAY_JOB_ID}_{task_id}"
                        for task_id in to_release)
    print(job_list)

    # slurm release
    try:
        slurm_process = subprocess.run([
            "scontrol" + get_platform_executable_extension(),
                "release", job_list
        ], timeout=SLURM_QUERY_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f"[release] scontrol release timed out after "
              f"{SLURM_QUERY_TIMEOUT_S}s; leaving these tasks held")
        return
    assert slurm_process.returncode == 0

    for group, index in to_release.values():
        jobs.set_state(experiment_name, group, index, jobs.RELEASED)


def run_seeded_simulations(experiment_name, run_seeded_simulation, groups=None):
    """the simplicity.runner.run_seeded_simulations function

    `groups` restricts the run to some of the experiment's groups; None means
    all of them. One submission deliberately spans whatever it is given, because
    cal_2 runs all of its (cell, scenario) groups as a single Slurm array and
    this function blocks -- one array per group would serialise them.
    """
    SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM = int(os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"])

    # Before sbatch and before the submission is recorded, so a refused run
    # leaves neither a held array nor a half-written map behind.
    raise_on_stale_state(experiment_name, groups)

    submission_id, pairs = jobs.write_submission(experiment_name, groups)
    if not pairs:
        raise RuntimeError(f"{experiment_name} has no repeats to run"
                           f"{'' if groups is None else f' in {list(groups)}'}. "
                           f"Was jobs.write_repeats called?")

    submit_simulations(experiment_name, run_seeded_simulation,
                       submission_id, pairs)
    # "submitted" used to read as "N are running now", which it never was: the
    # whole array goes to Slurm at once, on hold, and release_simulations lets
    # them go a capped number at a time from the polling loop below.
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {experiment_name}: queued "
          f"{len(pairs)} repeats as one held Slurm array (submission "
          f"{submission_id}); at most "
          f"{SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM} "
          f"released to run at a time")

    # loop until no repeat left to release
    last_status  = None
    last_long_running_report = time.time()
    last_reconcile = time.time()
    last_launch_failure_reconcile = time.time()
    while (status := poll_simulations_status(experiment_name, pairs)).left > 0:
        # print only when a repeat actually changed status
        if last_status != status:
            print_simulations_status(status, experiment_name)
        last_status = status

        # release repeats (silently -- this can fire every poll cycle once jobs
        # start turning over, and the periodic status line above already conveys
        # progress without repeating this on every release)
        n = min(status.left - status.pending, SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM - status.pending - status.running)
        if n:
            release_simulations(experiment_name, pairs, n)

        # periodically reconcile started-but-unresolved repeats against Slurm's
        # own accounting -- this is what unblocks a run stuck on a task Slurm
        # killed externally (walltime/OOM) without ever getting a chance to
        # record completed/failed itself.
        if (time.time() - last_reconcile) >= RECONCILE_INTERVAL_S:
            reconcile_terminated_tasks(experiment_name, pairs)
            last_reconcile = time.time()

        # periodically reconcile tasks Slurm auto-requeued into a held state
        # after a launch failure -- these never reach STARTED at all, so
        # reconcile_terminated_tasks above never sees them.
        if (time.time() - last_launch_failure_reconcile) >= LAUNCH_FAILURE_RECONCILE_INTERVAL_S:
            reconcile_launch_failures(experiment_name, pairs)
            last_launch_failure_reconcile = time.time()

        # once an hour, report the internal state of any repeat that's been
        # running for more than an hour (time/final_time, infected count)
        if (time.time() - last_long_running_report) >= LONG_RUNNING_REPORT_INTERVAL_S:
            report_long_running_simulations(experiment_name, pairs)
            last_long_running_report = time.time()

        # sleep
        time.sleep(7.)

    # completed
    print_simulations_status(status, experiment_name)

    
def job():
    """This runs in a process started by Slurm"""
    import os, sys
    print("<simplicity.runners.slurm.job>")

    # print useful paths to debug import errors (hint: help(simplicity.runners.slurm))
    print(os.path.abspath(os.curdir))
    print(sys.path)

    # raise if not called by sbatch (see submit_simulations)
    if "SLURM_ARRAY_TASK_ID" not in os.environ:
        raise Exception("this code is meant to be executed as a Slurm task. Hint: Slurm jobs management is handled by 'simplicity.runners.slurm'.")

    # retrieve arguments value (set by submit_simulations) 
    experiment_name                = os.environ["SIMPLICITY_EXPERIMENT_NAME"]
    submission_id                  = os.environ[SUBMISSION_ENV]
    run_seeded_simulation_qualname = os.environ["USER_RUN_SEEDED_SIMULATION"]

    # Resolve this task's repeat by looking its array position up in the
    # submission's recorded map. One lookup: the order was written down once, by
    # the process that submitted the array, instead of being rediscovered here
    # by walking a directory that gains files while the run is in flight.
    position = int(os.environ["SLURM_ARRAY_TASK_ID"]) - int(os.environ["SLURM_ARRAY_TASK_MIN"])
    group, repeat = jobs.resolve_task(experiment_name, submission_id, position)
    index = repeat["index"]

    # save the mapping from slurm job to repeat
    slurm_id_map_dir  = dm.get_slurm_id_map_dir(experiment_name)
    slurm_array_job_id = os.getenv('SLURM_ARRAY_JOB_ID')
    slurm_array_task_id = os.getenv('SLURM_ARRAY_TASK_ID')
    map_file = f"{slurm_id_map_dir}/{experiment_name}_{slurm_array_job_id}_{slurm_array_task_id}.csv"  
    with open(map_file, mode='w', newline='') as file:
        file.write(f"{group}/{index}")

    # A resumed experiment submits the whole array again -- the task ids are
    # positional, so there is no way to submit a subset -- and the repeats that
    # already finished bail out here, before any state is written, so the counts
    # stay right and nothing recomputes.
    if resuming() and jobs.state_of(experiment_name, group, index) == jobs.COMPLETED:
        print(f"<skip> {jobs.repeat_description(group, repeat)} already completed")
        return

    try:
        jobs.set_state(experiment_name, group, index, jobs.STARTED)

        # import run_seeded_simulation function from its qualname
        import importlib
        run_seeded_simulation_module_qualname, fn_name = run_seeded_simulation_qualname.rsplit(".", maxsplit=1)
        run_seeded_simulation_module = importlib.import_module(run_seeded_simulation_module_qualname)
        run_seeded_simulation = getattr(run_seeded_simulation_module, fn_name)

        run_seeded_simulation(experiment_name, group, index)

    except Exception as exc:
        jobs.set_state(experiment_name, group, index, jobs.FAILED)
        raise exc
    else:
        jobs.set_state(experiment_name, group, index, jobs.COMPLETED)

    print("</simplicity.runners.slurm.job>")

    
if __name__ == "__main__":
    job()
