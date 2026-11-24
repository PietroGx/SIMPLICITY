#!/usr/bin/env python3
'''A calibration run must write exactly the files its own fits and plots read.

Calibration reads four outputs and used to write seven. The three extras were
not small: on the calibrated grid lineage_frequency.csv alone reached 100.6 GB
across Data/, with a single simulation of the top NSR_long sweep point weighing
1.1 GB -- a grid point swept only so the fit can bracket it from above and
reject it. /home hit its quota mid-run and git could no longer pull.

CALIBRATION_SKIPPED_OUTPUTS now names the three. The danger in that list is
obvious: put a file in it that something downstream opens, and the calibration
fails hours in, or worse, fits on whatever it can still find. So this does not
check the list against a copy of itself -- it reads the calibration and sanity
scripts and checks that

  nothing skipped is read       no reader for a skipped file appears anywhere
                                in the calibration or sanity path. This is what
                                fails if someone later adds one.

  everything kept is read       each file still written is actually opened by
                                that path. A file nobody reads should be on the
                                skip list, not written "just in case".

  the switch is declared        the same call the calibration stages make
                                (set_calibration_output_env, at cal_1.py:118
                                and cal_2.py:279) is what reaches dir_manager,
                                and dm.skipping -- the function om.save_*
                                consults -- agrees with the list. Undeclared
                                raises rather than defaulting to "skip
                                nothing", which is the default that let a
                                stage write the 100.6 GB in the first place.

    python tests/test_calibration_outputs.py
'''
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments_pipeline'))

import simplicity.dir_manager as dm
from impact_long_shedders_config import (
    CALIBRATION_SKIPPED_OUTPUTS, set_calibration_output_env,
)

# Everything a simulation can write, and the output_manager accessors that open
# it -- alternatives separated by "|". extract_ih_regression_data reads
# individuals_data and phylogenetic_data directly, so those two are matched on
# their file names as well. lineage_frequency is matched on the stem
# "lineage_frequency.csv", which catches the gzipped name too, and on
# lineage_frequency_path, which reaches the file without naming it.
READERS = {
    'final_time': ('read_final_time', 'final_time.csv'),
    'sequencing_data_regression': ('read_sequencing_data_regression',
                                   'sequencing_data_regression.csv'),
    'individuals_data': ('read_individuals_data', 'individuals_data.csv'),
    'phylogenetic_data': ('read_phylogenetic_data', 'phylogenetic_data.csv'),
    'lineage_frequency': ('read_lineage_frequency|lineage_frequency_path',
                          'lineage_frequency.csv'),
    'simulation_trajectory': ('read_simulation_trajectory',
                              'simulation_trajectory.csv'),
    'sequencing_data': ('read_sequencing_data', 'sequencing_data.csv'),
}

# The stages a calibration run actually goes through: both fits, and the sanity
# regressions the pipeline dispatches afterwards.
CALIBRATION_PATH = [
    'scripts/experiments_pipeline/impact_long_shedders_cal_1.py',
    'scripts/experiments_pipeline/impact_long_shedders_cal_2.py',
    'scripts/experiments_pipeline/long_nsr_calibration_plot.py',
    'scripts/experiments_pipeline/plot_sot_sanity_regressions.py',
    # cal_1's plot and the sanity regressions both reach the intra-host clock
    # through this, so what IT reads is part of the calibration path
    'simplicity/tuning/evolutionary_rate.py',
]

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def is_read(name, source):
    """Does the calibration path open this output, by reader call or by name?"""
    reader, filename = READERS[name]
    return bool(re.search(rf'\b(?:{reader})\s*\(', source)) or filename in source


def undeclared():
    """True when dir_manager refuses to answer instead of defaulting to "skip
    nothing" -- the default a stage could forget its way past."""
    try:
        dm.skipped_outputs()
    except RuntimeError:
        return True
    return False


def read_path():
    text = ''
    for rel in CALIBRATION_PATH:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            failures.append(f'missing from the calibration path: {rel}')
            continue
        with open(path) as handle:
            # comments describe what a stage does NOT read as often as what it
            # does; matching them would make this test unfailable
            text += '\n'.join(line.split('#')[0] for line in handle)
    return text


def main():
    source = read_path()
    skipped = set(CALIBRATION_SKIPPED_OUTPUTS)
    kept = set(READERS) - skipped

    print('\nwhat the calibration path opens')
    for name in sorted(READERS):
        where = 'skipped' if name in skipped else 'written'
        print(f'    {name:<28} {where:<8} read={is_read(name, source)}')

    print('\nnothing skipped is read')
    for name in sorted(skipped):
        check(f'{name} is not read by the calibration path',
              is_read(name, source), False)

    print('\neverything kept is read')
    for name in sorted(kept):
        check(f'{name} is read, so writing it is justified',
              is_read(name, source), True)

    print('\nthe switch reaches dir_manager')
    os.environ.pop(dm.SKIP_OUTPUTS_ENV, None)
    check('an undeclared run refuses to guess', undeclared(), True)
    # the one call impact_long_shedders_cal_1.py:118 and cal_2.py:279 make
    set_calibration_output_env()
    check('the calibration set is applied',
          sorted(dm.skipped_outputs()), sorted(skipped))
    check('dm.skipping agrees, which is what om.save_* consults',
          sorted(n for n in dm.ALL_OUTPUTS if dm.skipping(n)), sorted(skipped))
    os.environ.pop(dm.SKIP_OUTPUTS_ENV, None)
    check('and clears again', undeclared(), True)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
