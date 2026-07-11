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
   started by Slurm. They will fail just after the  .started signal is emmitted.

2. Interrupting the program before completion will leave held jobs in Slurm's queue,
   and leaving the Data directory in a incomplete state (like for other
   simplicity.runners implementations).

   In this case, manually cleaning Slurm (hint: squeue) and the Data directory
   (removing or archiving the experiment folder) is highly recommanded before retrying.

3. A task killed directly by Slurm (walltime/--time or memory/--mem limit)
   never gets a chance to run job()'s except block and touch .failed itself
   -- the process is just terminated. Left unresolved, poll_simulations_status
   would count that task as permanently "left" and run_seeded_simulations'
   polling loop would never see status.left reach 0. Handled by
   reconcile_terminated_tasks, called periodically (RECONCILE_INTERVAL_S)
   from run_seeded_simulations: it cross-checks .started-but-unresolved tasks
   against `sacct` (which retains terminal-state history after a job leaves
   squeue) and touches .failed on Slurm's behalf for any task Slurm reports
   as TIMEOUT/OUT_OF_MEMORY/CANCELLED/etc.

"""
import typing, os, json, pathlib, subprocess, platform
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm

# How often (seconds) run_seeded_simulations' polling loop reports the
# internal state (time/final_time/infected) of simulations that have been
# running for a while -- separate from the ~17s SimulationsStatus printout,
# since this is meant as an occasional deeper look, not a repeated ping.
LONG_RUNNING_REPORT_INTERVAL_S = 3600
# How long a simulation must have been running (since its .started signal)
# before it's included in that report.
LONG_RUNNING_THRESHOLD_S = 3600

# How often (seconds) run_seeded_simulations' polling loop reconciles
# .started-but-unresolved tasks against Slurm's own accounting (sacct). This
# is what unblocks a run stuck on an externally-killed task (see
# reconcile_terminated_tasks below) -- unlike the purely-informational
# LONG_RUNNING_REPORT above, so it runs on its own, shorter timer.
RECONCILE_INTERVAL_S = 900

# sacct job states that mean "Slurm itself terminated this task" -- i.e. the
# task's own process never got a chance to touch .completed/.failed (job()'s
# except block only runs on a catchable Python exception, not a SIGTERM/
# SIGKILL from Slurm hitting --time/--mem). RUNNING/PENDING/COMPLETING are
# deliberately excluded (still in flight); COMPLETED is excluded too (job()
# should have already touched .completed itself in that case).
SLURM_TERMINAL_FAILURE_STATES = {
    "TIMEOUT", "OUT_OF_MEMORY", "FAILED", "NODE_FAIL", "CANCELLED",
    "DEADLINE", "PREEMPTED", "BOOT_FAIL",
}

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

    
def get_platform_executable_extension():
    """file extension to use when calling Slurm command-line utilities (sbatch, squeue, scontrol)."""
    return ".exe" if platform.system() == "Windows" else  ""


def submit_simulations(experiment_name: str, 
                       run_seeded_simulation: typing.Callable, 
                       n:int):
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
    
    # ! slurm array indexing is 1 based
    batch_start = 1
    batch_size  = n
    batch_end   = batch_start + batch_size - 1
   
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
        "USER_RUN_SEEDED_SIMULATION": run_seeded_simulation_qualname
        # "PLOT_TRAJECTORY"           : str(plot_trajectory)
    }), input=stdin)
    assert slurm_process.returncode == 0, f"Slurm was called with the following arguments:\n{' '.join(args)}\n{env}\n=== stdin\n{stdin}\n=== /stdin"
    
    # ! python __getitem__ indexing is 0 based
    seeded_simulation_parameters = sm.get_seeded_simulation_parameters_paths(experiment_name)
    for seeded_simulation_parameters_path in seeded_simulation_parameters[batch_start-1:batch_end]:
        signal_submitted_path   = seeded_simulation_parameters_path + ".submitted"
        pathlib.Path(signal_submitted_path).touch()
        
    
def poll_simulations_status(experiment_name):
    # get the seeded simulation parameters files paths
    seeded_simulation_parameters = sm.get_seeded_simulation_parameters_paths(experiment_name)
    total = len(seeded_simulation_parameters)
    # read signals written by this script
    submitted = 0
    released  = 0
    # read signals written by run_seeded_simulations.py
    started   = 0
    completed = 0
    failed    = 0
    for seeded_simulation_parameters_path in seeded_simulation_parameters:
        signal_submitted_path = seeded_simulation_parameters_path + ".submitted"
        signal_released_path  = seeded_simulation_parameters_path + ".released"
        signal_started_path   = seeded_simulation_parameters_path + ".started"
        signal_completed_path = seeded_simulation_parameters_path + ".completed"
        signal_failed_path    = seeded_simulation_parameters_path + ".failed"
        submitted += pathlib.Path(signal_submitted_path ).exists()
        released  += pathlib.Path(signal_released_path  ).exists()
        started   += pathlib.Path(signal_started_path   ).exists()
        completed += pathlib.Path(signal_completed_path ).exists()
        failed    += pathlib.Path(signal_failed_path    ).exists()
    return SimulationsStatus(
        total    = total,
        submitted= submitted,
        released = released,
        left     = submitted - (completed+failed),
        pending  = released - started,
        started  = started,
        running  = started -  (completed+failed),
        completed= completed,
        failed   = failed,
    )

def report_long_running_simulations(experiment_name, threshold_seconds=LONG_RUNNING_THRESHOLD_S):
    """Print the internal state (time/final_time/infected) of every seeded
    simulation that has been running for more than threshold_seconds.

    State is read from the periodic progress snapshot
    simplicity.extrande.ProgressReporter writes next to the seeded params
    file (same "<seeded_simulation_parameters_path>.xyz" convention as the
    .started/.completed/.failed signals below), so this works without a live
    terminal attached to the actual simulation process.
    """
    import time
    seeded_simulation_parameters = sm.get_seeded_simulation_parameters_paths(experiment_name)
    now = time.time()
    for seeded_simulation_parameters_path in seeded_simulation_parameters:
        signal_started_path   = seeded_simulation_parameters_path + ".started"
        signal_completed_path = seeded_simulation_parameters_path + ".completed"
        signal_failed_path    = seeded_simulation_parameters_path + ".failed"
        started_path = pathlib.Path(signal_started_path)
        if not started_path.exists():
            continue
        if pathlib.Path(signal_completed_path).exists() or pathlib.Path(signal_failed_path).exists():
            continue
        elapsed = now - started_path.stat().st_mtime
        if elapsed < threshold_seconds:
            continue

        name = os.path.basename(seeded_simulation_parameters_path)
        progress_path = pathlib.Path(seeded_simulation_parameters_path + ".progress")
        if not progress_path.exists():
            print(f"[long-running] {name}: running {elapsed/3600:.1f}h, no progress snapshot yet")
            continue
        try:
            with open(progress_path) as f:
                snapshot = json.load(f)
        except (OSError, json.JSONDecodeError):
            print(f"[long-running] {name}: running {elapsed/3600:.1f}h, progress snapshot unreadable")
            continue
        print(f"[long-running] {name}: running {elapsed/3600:.1f}h, "
              f"time={snapshot.get('time')}/{snapshot.get('final_time')}, "
              f"infected={snapshot.get('infected')}")


def _build_slurm_id_map_index(experiment_name):
    """Reverse-index job()'s slurm_id_map_dir CSV files (job()-written, never
    previously read): {seeded_simulation_parameters_path: (job_id, task_id)}.
    Filenames are "{experiment_name}_{job_id}_{task_id}.csv"; content is the
    seeded_simulation_parameters_path itself."""
    index = {}
    map_dir = pathlib.Path(dm.get_slurm_id_map_dir(experiment_name))
    prefix = f"{experiment_name}_"
    for map_file in map_dir.glob(f"{prefix}*.csv"):
        # "{experiment_name}_{job_id}_{task_id}" -> split off the prefix, then
        # the remaining "{job_id}_{task_id}" splits on the last underscore.
        stem = map_file.stem[len(prefix):]
        job_id, _, task_id = stem.rpartition("_")
        if not job_id:
            continue
        try:
            path = map_file.read_text().strip()
        except OSError:
            continue
        if path:
            index[path] = (job_id, task_id)
    return index


def reconcile_terminated_tasks(experiment_name):
    """
    Find every seeded simulation with a .started signal but no .completed/
    .failed, and check Slurm's own accounting (sacct) for whether Slurm
    itself already killed that task (TIMEOUT/OUT_OF_MEMORY/CANCELLED/...).

    This exists because job()'s try/except (simplicity.runners.slurm.job,
    above) can only touch .failed on a catchable Python exception -- a Slurm
    walltime or memory kill terminates the process directly, so job() never
    gets to run its except block, and that task would otherwise stay
    "started" forever. poll_simulations_status's `left` count then never
    reaches 0 for it, and run_seeded_simulations' polling loop can only exit
    through the unrelated (and previously silently-swallowed) "no held task"
    exception in release_simulations -- see the module Kown bugs note. Unlike
    a live task's own job() process, sacct retains terminal-state history
    after a job leaves squeue, so it's the only reliable place to check this
    from outside the (now-dead) task's own process.

    Any task confirmed terminated by Slurm gets .failed touched on its
    behalf, with a note explaining it was an externally-detected kill (not a
    Python-level failure) -- this is what actually unblocks the polling loop.
    """
    seeded_simulation_parameters = sm.get_seeded_simulation_parameters_paths(experiment_name)
    stuck_paths = []
    for seeded_simulation_parameters_path in seeded_simulation_parameters:
        signal_started_path   = seeded_simulation_parameters_path + ".started"
        signal_completed_path = seeded_simulation_parameters_path + ".completed"
        signal_failed_path    = seeded_simulation_parameters_path + ".failed"
        if not pathlib.Path(signal_started_path).exists():
            continue
        if pathlib.Path(signal_completed_path).exists() or pathlib.Path(signal_failed_path).exists():
            continue
        stuck_paths.append(seeded_simulation_parameters_path)

    if not stuck_paths:
        return

    id_index = _build_slurm_id_map_index(experiment_name)

    # Group stuck paths by Slurm array job id, so each distinct job id is
    # queried via sacct exactly once (normally there's only one, but nothing
    # here assumes that).
    by_job_id = {}
    unmapped = []
    for path in stuck_paths:
        ids = id_index.get(path)
        if ids is None:
            # job() hasn't written its map file yet (narrow window right
            # after .started, before the map file write) -- not stuck, just
            # not checkable yet. Skip silently; it'll be found on a later pass.
            unmapped.append(path)
            continue
        job_id, task_id = ids
        by_job_id.setdefault(job_id, {})[f"{job_id}_{task_id}"] = path

    for job_id, task_id_to_path in by_job_id.items():
        sacct_process = subprocess.run([
            "sacct" + get_platform_executable_extension(),
                "-j", job_id, "--format=JobID,State", "--noheader", "--parsable2", "-X",
        ], stdout=subprocess.PIPE)
        if sacct_process.returncode != 0:
            # sacct itself unreachable/erroring -- try again on the next
            # reconcile pass rather than guessing at these tasks' state now.
            continue

        for line in sacct_process.stdout.decode().splitlines():
            if "|" not in line:
                continue
            sacct_job_id, state = line.split("|", maxsplit=1)
            sacct_job_id = sacct_job_id.strip()
            path = task_id_to_path.get(sacct_job_id)
            if path is None:
                continue
            # sacct states can carry a suffix, e.g. "CANCELLED by 12345".
            state = state.strip().split()[0] if state.strip() else state.strip()
            if state in SLURM_TERMINAL_FAILURE_STATES:
                pathlib.Path(path + ".failed").touch()
                name = os.path.basename(path)
                print(f"[reconciled] {name}: Slurm reports {state} -- "
                     f"marking .failed (task never signaled itself)")


def release_simulations(experiment_name, n: int):
    # get the seeded simulation parameters files paths
    seeded_simulation_parameters_paths = sm.get_seeded_simulation_parameters_paths(experiment_name)
    
    # use signals to find up to n submitted but not released simulations
    i_th_seeds = {}
    for i_th_seed, seeded_simulation_parameters_path in enumerate(seeded_simulation_parameters_paths):
        signal_submitted_path = seeded_simulation_parameters_path + ".submitted"
        signal_released_path  = seeded_simulation_parameters_path + ".released"
        if pathlib.Path(signal_submitted_path).exists() and not pathlib.Path(signal_released_path).exists():
            i_th_seeds[i_th_seed] = seeded_simulation_parameters_path
        if len(i_th_seeds) >= n:
            break

    if not i_th_seeds:
        # Nothing actually needs releasing right now. This happens once
        # every task has already been released at least once, while `n` (the
        # caller's estimate of remaining capacity) can still come out > 0
        # because status.left also counts tasks stuck on an externally-killed
        # Slurm task that hasn't been reconciled yet (see
        # reconcile_terminated_tasks) -- those inflate `left` without
        # inflating `pending`. Querying squeue here would be pointless (nothing
        # to release) and, once the whole array job has aged out of squeue's
        # listing, would incorrectly raise "no held task" for a totally benign
        # state. Just return; the reconciler is what actually clears this up.
        return

    # slurm find array job id from job name
    slurm_process = subprocess.run([
        "squeue" + get_platform_executable_extension(),
            "--Format=ArrayJobID", f"--name={experiment_name}" 
    ], stdout=subprocess.PIPE)
    assert slurm_process.returncode == 0
    array_job_id_set = set(line.strip() for line in slurm_process.stdout.decode().splitlines(keepends=False)[1:])
    if len(array_job_id_set) == 0:
        raise Exception("no held task in the job array. Hint: check slurm log output for error occurring before .started signal ")
    assert len(array_job_id_set) == 1, f"Expect exactly one job array with name {experiment_name}. If several slurm array job share the name {experiment_name}. please clean slurm's queue before continuing.\nKnown bug if Slurm's hold a job in CG state.\nHint: better Slurm's squeue parsing in 'simplicity.runners.slurm' may resolve CG status case.\n" + slurm_process.stdout.decode()
    SLURM_ARRAY_JOB_ID = next(iter(array_job_id_set))
    
    # ! slurm task array indexing is 1 based
    SLURM_ARRAY_TASK_IDs = [i_th_seed + 1 for i_th_seed in i_th_seeds]
    job_list = ",".join(f"{SLURM_ARRAY_JOB_ID}_{SLURM_ARRAY_TASK_ID}" for SLURM_ARRAY_TASK_ID in SLURM_ARRAY_TASK_IDs)
    print(job_list)
    
    # slurm release
    slurm_process = subprocess.run([
        "scontrol" + get_platform_executable_extension(),
            "release", job_list
    ])
    assert slurm_process.returncode == 0
    
    # signal as released
    for i_th_seed, seeded_simulation_parameters_path in i_th_seeds.items():
        signal_released_path  = seeded_simulation_parameters_path + ".released"
        pathlib.Path(signal_released_path).touch()
        
def run_seeded_simulations(experiment_name, run_seeded_simulation):
    """the simplicity.runner.run_seeded_simulations function"""
    import time, os
    SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM = int(os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"])

    # retrieve status
    status = poll_simulations_status(experiment_name)

    # submit simulations
    submit_simulations(experiment_name, run_seeded_simulation, n=status.total)
    print(f"submitted {status.total} seeded simulations"); last_printed = time.time()
    
    # loop until no simulation left to release
    last_status  = None
    last_long_running_report = time.time()
    last_reconcile = time.time()
    while (status := poll_simulations_status(experiment_name)).left > 0:
        # print if status changed or after 17 seconds
        if last_status != status or (time.time() - last_printed) > 17.:
            print(status); last_printed = time.time()
        last_status = status

        # release simluations (silently -- this can fire every poll cycle
        # once jobs start turning over, and the periodic status line above
        # already conveys progress without repeating this on every release)
        n = min(status.left - status.pending, SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM - status.pending - status.running)
        if n:
            release_simulations(experiment_name, n)

        # periodically reconcile .started-but-unresolved tasks against
        # Slurm's own accounting -- this is what unblocks a run stuck on a
        # task Slurm killed externally (walltime/OOM) without ever getting a
        # chance to touch .completed/.failed itself.
        if (time.time() - last_reconcile) >= RECONCILE_INTERVAL_S:
            reconcile_terminated_tasks(experiment_name)
            last_reconcile = time.time()

        # once an hour, report the internal state of any simulation that's
        # been running for more than an hour (time/final_time, infected count)
        if (time.time() - last_long_running_report) >= LONG_RUNNING_REPORT_INTERVAL_S:
            report_long_running_simulations(experiment_name)
            last_long_running_report = time.time()

        # sleep
        time.sleep(7.)

    # completed
    print(status); last_printed = time.time()

    
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

    # retrieve arguments value (set by submit_simulation) 
    experiment_name                = os.environ["SIMPLICITY_EXPERIMENT_NAME"]
    run_seeded_simulation_qualname = os.environ["USER_RUN_SEEDED_SIMULATION"]
    # plot_trajectory = os.getenv("PLOT_TRAJECTORY", "False").lower() == "true"
    
    # resolve i_th seeded_simulation given Slurm's given rank
    i_th_seeded_simulation = int(os.environ["SLURM_ARRAY_TASK_ID"]) - int(os.environ["SLURM_ARRAY_TASK_MIN"])
    
    # resolve seeded_simulation_parameters_path
    seeded_simulation_parameters_paths = sm.get_seeded_simulation_parameters_paths(experiment_name)
    seeded_simulation_parameters_path  = seeded_simulation_parameters_paths[i_th_seeded_simulation]
    
    # save the mapping from slurm job to seeded simulation number
    slurm_id_map_dir  = dm.get_slurm_id_map_dir(experiment_name)
    slurm_array_job_id = os.getenv('SLURM_ARRAY_JOB_ID')
    slurm_array_task_id = os.getenv('SLURM_ARRAY_TASK_ID')
    map_file = f"{slurm_id_map_dir}/{experiment_name}_{slurm_array_job_id}_{slurm_array_task_id}.csv"  
    with open(map_file, mode='w', newline='') as file:
        file.write(seeded_simulation_parameters_path)
        
    # define signals 
    signal_started_path   = seeded_simulation_parameters_path + ".started"
    signal_failed_path    = seeded_simulation_parameters_path + ".failed"
    signal_completed_path = seeded_simulation_parameters_path + ".completed"

    # call and set signals
    try:
        # signal seeded simulation started
        pathlib.Path(signal_started_path).touch()
        
        # import run_seeded_simulation function from its qualname
        import importlib
        run_seeded_simulation_module_qualname, fn_name = run_seeded_simulation_qualname.rsplit(".", maxsplit=1)
        run_seeded_simulation_module = importlib.import_module(run_seeded_simulation_module_qualname)
        run_seeded_simulation = getattr(run_seeded_simulation_module, fn_name)
    
        # call run_seeded_simulation
        # print(type(plot_trajectory))
        run_seeded_simulation(seeded_simulation_parameters_path, experiment_name)
        
    except Exception as exc:
        # signal seeded simulation failed
        pathlib.Path(signal_failed_path).touch()
        # raise
        raise exc
    else:
        # signal seeded simulation completed
        pathlib.Path(signal_completed_path).touch()

    print("</simplicity.runners.slurm.job>")

    
if __name__ == "__main__":
    job()
