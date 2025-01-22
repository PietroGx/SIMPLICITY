#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jbescudie

"""
import simplicity.settings_manager as sm
import concurrent.futures


def run_seeded_simulations(experiment_name, run_seeded_simulation):
    """Implements run_seeded_simulations (see simplicity.runners).

    Uses CPython built-in concurrent.futures.ProcessPoolExecutor.
    """
    import time, os
    SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS = int(os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS"])
    
    seeded_simulation_parameters_paths = sm.get_seeded_simulation_parameters_paths(experiment_name)
    # run a Simplicity simulation for each seeded parameters
    with concurrent.futures.ProcessPoolExecutor(max_workers=SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS) as pool:
        # track status
        futures_to_args    = {}
        args_to_exceptions = {}
        
        # stats
        submitted = completed = failed = 0

        # submit
        for seeded_simulation_parameters_path in seeded_simulation_parameters_paths:
            args = (run_seeded_simulation, seeded_simulation_parameters_path, experiment_name)
            futures_to_args[pool.submit(*args)] = args
            submitted += 1
        
        # join
        futures = tuple(futures_to_args)
        for future_i, future in enumerate(concurrent.futures.as_completed(futures)):
            args = futures_to_args[future]
            try:
                future.result()
                completed += 1
            except Exception as exc:
                failed    += 1
                args_to_exceptions[args] = exc
                
            print(f"{future_i}/{submitted} simulation finished. Completed: {completed} Failed: {failed}")

    # raise if any exception occured
    if args_to_exceptions:
        exceptions_fmt = [f"{args}\n{exc}" for args,exc in args_to_exceptions.items()]
        raise Exception(f"{failed} failed simulations.\n\nExceptions:\n" + "\n".join(exceptions_fmt))

    # success
    print(f"Finished running seeded simulations. Submitted: {submitted} Completed: {completed} Failed: {failed}")
    
