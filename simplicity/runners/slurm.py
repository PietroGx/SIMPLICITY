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

"""
import typing, os, pathlib, subprocess, platform


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
                       plot_trajectory,
                       n:int):
    # run_seeded_simulation to qualname
    fn_name = run_seeded_simulation.__name__
    run_seeded_simulation_module_qualname = run_seeded_simulation.__module__
    if run_seeded_simulation_module_qualname == "__main__":
        import simplicity.runners.slurm
        help(simplicity.runners.slurm)
        raise Exception(f"run_seeded_simulation must be imported.")
    run_seeded_simulation_qualname = run_seeded_simulation_module_qualname + "." + fn_name
    
    # job is to run the same python to call this module main() (which in turn will call user's run_seeded_simulation function)
    import sys
    stdin = "\n".join((
        # run the same python
        f"#!{sys.executable}",
        # call this module main()
        f"import simplicity.runners.slurm",
        f"simplicity.runners.slurm.job()",
    )).encode()
    
    # ! slurm array indexing is 1 based
    batch_start = 1
    batch_size  = n
    batch_end   = batch_start + batch_size - 1
    
    # submit the job
    slurm_process = subprocess.run((args:=[
        # calls sbatch
        "sbatch" + get_platform_executable_extension(),
        # to create the job array on hold 
        f"--array={batch_start}-{batch_end}", "--hold",
        # with a name (used later for lookup)
        f"--job-name={experiment_name}",
    ]), env=(env:={
        **os.environ,
        "SIMPLICITY_EXPERIMENT_NAME": experiment_name,
        "USER_RUN_SEEDED_SIMULATION": run_seeded_simulation_qualname,
        "PLOT_TRAJECTORY"           : plot_trajectory
    }), input=stdin)
    assert slurm_process.returncode == 0, f"Slurm was called with the following arguments:\n{' '.join(args)}\n{env}\n=== stdin\n{stdin}\n=== /stdin"
    
    # ! python __getitem__ indexing is 0 based
    import simplicity.settings_manager as sm
    seeded_simulation_parameters = sm.get_seeded_simulation_parameters_paths(experiment_name)
    for seeded_simulation_parameters_path in seeded_simulation_parameters[batch_start-1:batch_end]:
        signal_submitted_path   = seeded_simulation_parameters_path + ".submitted"
        pathlib.Path(signal_submitted_path).touch()
        
    
def poll_simulations_status(experiment_name):
    # get the seeded simulation parameters files paths
    import simplicity.settings_manager as sm
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

def release_simulations(experiment_name, n: int):
    # get the seeded simulation parameters files paths
    import simplicity.settings_manager as sm
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
        
def run_seeded_simulations(experiment_name, run_seeded_simulation, plot_trajectory):
    """the simplicity.runner.run_seeded_simulations function"""
    import time, os
    SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM = int(os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"])

    # retrieve status
    status = poll_simulations_status(experiment_name)

    # submit simulations
    submit_simulations(experiment_name, run_seeded_simulation, plot_trajectory, n=status.total)
    print(f"submitted {status.total} seeded simulations"); last_printed = time.time()
    
    # loop until no simulation left to release
    last_status  = None
    while (status := poll_simulations_status(experiment_name)).left > 0:
        # print if status changed or after 17 seconds
        if last_status != status or (time.time() - last_printed) > 17.:
            print(status); last_printed = time.time()
        last_status = status
        
        # release simluations
        n = min(status.left - status.pending, SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM - status.pending - status.running)
        if n:
            print(f"release up to {n} seeded simulations"); last_printed = time.time()
            release_simulations(experiment_name, n)
            
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
    plot_trajectory                = os.environ["PLOT_TRAJECTORY"]
    
    # resolve i_th seeded_simulation given Slurm's given rank
    i_th_seeded_simulation = int(os.environ["SLURM_ARRAY_TASK_ID"]) - int(os.environ["SLURM_ARRAY_TASK_MIN"])
    
    # resolve seeded_simulation_parameters_path
    import simplicity.settings_manager as sm
    seeded_simulation_parameters_paths = sm.get_seeded_simulation_parameters_paths(experiment_name)
    seeded_simulation_parameters_path  = seeded_simulation_parameters_paths[i_th_seeded_simulation]

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
        run_seeded_simulation(seeded_simulation_parameters_path, experiment_name, plot_trajectory)
        
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
