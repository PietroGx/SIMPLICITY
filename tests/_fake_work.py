#!/usr/bin/env python3
'''Stand-in for a simulation, for tests of the Slurm control plane.

Records that it ran -- which is how a test asserts that every array position
resolved to a distinct repeat, exactly once -- and then does whatever
$FAKEWORK_BEHAVIOUR says. It never runs a simulation: this is a test of the
control plane, and an intra-host matrix exponential has nothing to say about
whether array task 7 found the right repeat.
'''
import os
import time


def _ran_dir():
    path = os.environ['FAKEWORK_RAN_DIR']
    os.makedirs(path, exist_ok=True)
    return path


def run_seeded_simulation(experiment_name, group, index):
    """The run_seeded_simulation contract: (experiment, group, index)."""
    with open(os.path.join(_ran_dir(), f'{group}_{index:06d}'), 'a') as handle:
        handle.write(f'{time.time()}\n')

    behaviour = os.environ.get('FAKEWORK_BEHAVIOUR', 'ok')
    if behaviour == 'raise':
        raise RuntimeError('fake work failed on purpose')
    if behaviour == 'hang':
        time.sleep(600)     # the fake controller kills this
