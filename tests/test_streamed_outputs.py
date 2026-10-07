#!/usr/bin/env python3
'''lineage_frequency and simulation_trajectory stream to disk; prove both paths
write the same bytes.

Population now writes these two files row by row instead of accumulating them
(555 MB and 47 MB in the worst N=5000 task of profile_grid_#910, measured). With
an output directory the rows go to "<name>.partial" and output_manager renames
it; without one they are buffered and pandas writes the file as before. The
buffered path is no longer exercised by a real run -- simulation.py always
passes a directory -- so it is only this test that keeps the two in step.

What is checked, and why each one is here rather than assumed:

  identical bytes   the streamed writer has no dtype inference. pandas read the
                    t=0 row's int 0 as part of a float64 column and wrote "0.0",
                    and csv.writer defaults to \\r\\n line endings; either alone
                    changes every line of the file. Both are pinned in code, and
                    this is what catches it if one is edited away.

  the rename        a run killed mid-flight must not leave a file under its real
                    name: scripts/check_completed_simulations.py:32 reads the
                    presence of lineage_frequency.csv as "this simulation
                    finished". Streaming without the rename would make every
                    killed run look complete.

    python tests/test_streamed_outputs.py
'''
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import pandas as pd

import simplicity.population as pop
import simplicity.output_manager as om

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def make(output_directory):
    """A Population with only what the two write paths touch.

    Built without __init__ on purpose: constructing a real one precomputes the
    intra-host matrix exponentials, which has nothing to do with writing csv
    rows and would make this test minutes long.
    """
    population = object.__new__(pop.Population)
    # __init__ is bypassed, so set what the write path reads: an empty skip set
    # means "write everything", which is what this test is about.
    population._skip = frozenset()
    population._streams = {}
    population.lineage_frequency = []
    population.trajectory = []
    if output_directory is not None:
        population._open_stream(output_directory, pop.LINEAGE_FREQUENCY_FILE,
                                pop.LINEAGE_FREQUENCY_COLUMNS)
        population._open_stream(output_directory, pop.TRAJECTORY_FILE,
                                pop.TRAJECTORY_COLUMNS)
    return population


def drive(population):
    """The same rows a short run produces, including the integer t=0 rows that
    pandas widens to float."""
    population._record(pop.LINEAGE_FREQUENCY_FILE, ['wt', 0.0, 1.0, 50])
    population._record(pop.TRAJECTORY_FILE,
                       [0.0, 50, 0, 0, 0, 0, 1450, 1])
    for step in range(1, 40):
        t = 27.0188943573 + step
        population._record(pop.LINEAGE_FREQUENCY_FILE,
                           [f'wt.{step}.{step % 7}', float(t),
                            0.6766917293233082 / (1 + step % 13), 91 + step])
        population._record(pop.TRAJECTORY_FILE,
                           [float(t), 50 + step, step, step * 2, step,
                            step, 1450 - step, 1])


def main():
    root = tempfile.mkdtemp(prefix='stream_test_')
    try:
        streamed_dir = os.path.join(root, 'streamed')
        buffered_dir = os.path.join(root, 'buffered')
        os.makedirs(streamed_dir)
        os.makedirs(buffered_dir)

        streamed = make(streamed_dir)
        drive(streamed)
        buffered = make(None)
        drive(buffered)

        print('\nwhile the run is in flight')
        check('streamed rows are on disk, under .partial only',
              sorted(os.listdir(streamed_dir)),
              ['lineage_frequency.csv.partial',
               'simulation_trajectory.csv.partial'])
        check('nothing is held in memory when streaming',
              (len(streamed.lineage_frequency), len(streamed.trajectory)),
              (0, 0))
        check('buffered run holds every row instead',
              (len(buffered.lineage_frequency), len(buffered.trajectory)),
              (40, 40))

        # the real save path, both ways
        for population, directory in ((streamed, streamed_dir),
                                      (buffered, buffered_dir)):
            om.save_lineage_frequency(population, directory)
            om.save_simulation_trajectory(population, directory)

        print('\nafter saving')
        check('the streamed files are published under their real names',
              sorted(os.listdir(streamed_dir)),
              ['lineage_frequency.csv', 'simulation_trajectory.csv'])

        for name in (pop.LINEAGE_FREQUENCY_FILE, pop.TRAJECTORY_FILE):
            with open(os.path.join(streamed_dir, name), 'rb') as handle:
                streamed_bytes = handle.read()
            with open(os.path.join(buffered_dir, name), 'rb') as handle:
                buffered_bytes = handle.read()
            check(f'{name}: streamed and buffered agree byte for byte',
                  streamed_bytes == buffered_bytes, True)
            check(f'{name}: no CR in the line endings',
                  b'\r' not in streamed_bytes, True)
            # the t=0 row is the one pandas widens to float
            first = streamed_bytes.decode().splitlines()[1]
            check(f'{name}: t=0 row keeps the float form', first.split(',')[1],
                  '0.0' if name == pop.LINEAGE_FREQUENCY_FILE else '50')

        print('\nstill readable by the pipeline')
        frame = om.read_lineage_frequency(streamed_dir)
        check('read_lineage_frequency returns every row', len(frame), 40)
        check('columns are unchanged', list(frame.columns),
              pop.LINEAGE_FREQUENCY_COLUMNS)
        check('Frequency_at_t is still numeric',
              pd.api.types.is_float_dtype(frame['Frequency_at_t']), True)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
