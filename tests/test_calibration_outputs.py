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

    python tests/test_calibration_outputs.py
'''
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
from impact_long_shedders_unbound_config import (
    CALIBRATION_SKIPPED_OUTPUTS, set_calibration_output_env,
)

# Everything a simulation can write, and the output_manager reader that opens
# it. extract_ih_regression_data reads individuals_data and phylogenetic_data
# directly, so those two are matched on their file names as well.
READERS = {
    'final_time': ('read_final_time', 'final_time.csv'),
    'sequencing_data_regression': ('read_sequencing_data_regression',
                                   'sequencing_data_regression.csv'),
    'individuals_data': ('read_individuals_data', 'individuals_data.csv'),
    'phylogenetic_data': ('read_phylogenetic_data', 'phylogenetic_data.csv'),
    'lineage_frequency': ('read_lineage_frequency', 'lineage_frequency.csv'),
    'simulation_trajectory': ('read_simulation_trajectory',
                              'simulation_trajectory.csv'),
    'sequencing_data': ('read_sequencing_data', 'sequencing_data.csv'),
}

# The stages a calibration run actually goes through: both fits, and the sanity
# regressions the pipeline dispatches afterwards.
CALIBRATION_PATH = [
    'scripts/experiments/impact_long_shedders_unbound_cal_1.py',
    'scripts/experiments/impact_long_shedders_unbound_cal_2.py',
    'scripts/experiments/long_nsr_calibration_plot.py',
    'scripts/experiments/plot_sot_sanity_regressions.py',
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
        reader, filename = READERS[name]
        used = bool(re.search(rf'\b{reader}\s*\(', source)) or filename in source
        where = 'skipped' if name in skipped else 'written'
        print(f'    {name:<28} {where:<8} read={used}')

    print('\nnothing skipped is read')
    for name in sorted(skipped):
        reader, filename = READERS[name]
        used = bool(re.search(rf'\b{reader}\s*\(', source)) or filename in source
        check(f'{name} is not read by the calibration path', used, False)

    print('\neverything kept is read')
    for name in sorted(kept):
        reader, filename = READERS[name]
        used = bool(re.search(rf'\b{reader}\s*\(', source)) or filename in source
        check(f'{name} is read, so writing it is justified', used, True)

    print('\nthe switch reaches dir_manager')
    before = dm.skipped_outputs()
    check('nothing skipped by default', set(before), set())
    set_calibration_output_env()
    check('the calibration set is applied',
          sorted(dm.skipped_outputs()), sorted(skipped))
    os.environ.pop(dm.SKIP_OUTPUTS_ENV, None)
    check('and clears again', set(dm.skipped_outputs()), set())

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
