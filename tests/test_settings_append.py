#!/usr/bin/env python3
'''Stage 8a: one experiment must be able to gain a phase's groups later.

    python tests/test_settings_append.py

WHY

Stage 8 folds the pipeline into ONE experiment whose phases are groups. The
phases are sequential and data-dependent -- the experiment phase's parameters
come from the calibration phase's fit -- so its groups cannot be written until
the calibration has finished. The experiment must therefore be APPENDED to,
not rewritten.

The state layer already allows this (jobs.write_repeats guards per group,
write_submission takes a `groups` filter, raise_on_stale_state takes `groups`).
Two things in settings_manager did not:

  _as_simulations assigned id = enumerate(...), a GLOBAL position restarting at
  0 on every call, and simulation_stem puts that id in the directory name. A
  second call reissued ids 0..n and collided with the first phase's output
  directories -- the exact "two parameter sets sharing one filename and
  silently overwriting each other" failure its docstring says is impossible.

  write_experiment_settings opened 'w' (truncating the previous phase away) and
  stamped ONE scalar n_seeds on every group, so --cal-seeds and --exp-seeds
  could not coexist in one experiment.

Every check below is on the real settings/jobs path, with real files.
'''
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm

os.environ.setdefault('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM', '4')

failures = []


def check(label, ok, detail=''):
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}" + (f'   {detail}' if detail else ''))
    if not ok:
        failures.append(f'{label}: {detail}')


def raises(label, call, *, naming=()):
    try:
        call()
    except Exception as exc:                                    # noqa: BLE001
        text = str(exc)
        missing = [token for token in naming if token not in text]
        check(label, not missing,
              f'{type(exc).__name__}' +
              (f' -- but did not name {missing}' if missing else ''))
        return
    check(label, False, 'did not raise')


FIXED = {'population_size': 20, 'final_time': 2}

CAL_GROUPS = [
    {'name': 'calib_long__SOT',
     'nucleotide_substitution_rate_long': [1e-4, 2e-4, 3e-4]},
    {'name': 'calib_long__edge_case',
     'nucleotide_substitution_rate_long': [1e-4, 2e-4]},
]
EXP_GROUPS = [
    {'name': 'exp__SOT', 'nucleotide_substitution_rate': 5e-5},
    {'name': 'exp__control', 'nucleotide_substitution_rate': 6e-5},
]
CAL_SEEDS, EXP_SEEDS = 10, 3


class Scratch:
    def __enter__(self):
        self.root = tempfile.mkdtemp(prefix='settings_append_')
        dm.set_data_dir(self.root)
        sm._data_dir = self.root          # a copy taken at import time
        return self

    def __exit__(self, *exc):
        shutil.rmtree(self.root, ignore_errors=True)


def write_phase(name, groups, n_seeds, append):
    settings = sm.generate_experiment_settings(
        {'_scenario_groups': groups}, FIXED,
        start_id=sm.next_simulation_id(name) if append else 0)
    sm.write_experiment_settings(name, settings, n_seeds, append=append)
    sm.read_settings_and_write_simulation_parameters(name)
    return settings


# ===========================================================================
def test_append_keeps_the_first_phase():
    name = 'merged'
    dm.create_directories(name)
    write_phase(name, CAL_GROUPS, CAL_SEEDS, append=False)
    before = sm.read_simulations(name)
    check('calibration phase wrote its simulations', len(before) == 5,
          f'{len(before)} (3 + 2)')

    write_phase(name, EXP_GROUPS, EXP_SEEDS, append=True)
    after = sm.read_simulations(name)

    check('appending kept every earlier simulation, unchanged',
          after[:len(before)] == before, f'{len(before)} -> {len(after)}')
    check('and added the experiment phase', len(after) == 7, f'{len(after)}')

    ids = [record['id'] for record in after]
    check('every simulation id is distinct', len(set(ids)) == len(ids), f'{ids}')
    check('ids continue rather than restarting', ids == list(range(7)), f'{ids}')

    standard = sm.read_standard_parameters_values()
    stems = [sm.simulation_stem(record, standard) for record in after]
    check('every output directory stem is distinct -- THE failure 8a fixes',
          len(set(stems)) == len(stems),
          f'{len(set(stems))} unique of {len(stems)}')

    groups = [g['name'] for g in sm.read_groups(name)]
    check('groups are both phases, in order',
          groups == ['calib_long__SOT', 'calib_long__edge_case',
                     'exp__SOT', 'exp__control'], f'{groups}')
    return name


# ===========================================================================
def test_per_group_seeds(name):
    seeds = {g['name']: g['n_seeds'] for g in sm.read_groups(name)}
    check('calibration groups kept --cal-seeds',
          seeds['calib_long__SOT'] == CAL_SEEDS
          and seeds['calib_long__edge_case'] == CAL_SEEDS, f'{seeds}')
    check('experiment groups carry --exp-seeds',
          seeds['exp__SOT'] == EXP_SEEDS and seeds['exp__control'] == EXP_SEEDS,
          f'{seeds}')
    check('read_group_n_seeds agrees',
          sm.read_group_n_seeds(name, 'exp__SOT') == EXP_SEEDS)

    # the legacy scalar cannot answer this, so it must refuse rather than
    # silently hand back one phase's number for the whole experiment
    raises('read_n_seeds_file refuses a mixed experiment, naming the groups',
           lambda: sm.read_n_seeds_file(name),
           naming=('calib_long__SOT', 'exp__SOT', str(CAL_SEEDS),
                   str(EXP_SEEDS)))


def test_scalar_still_works_when_groups_agree():
    name = 'uniform'
    dm.create_directories(name)
    write_phase(name, CAL_GROUPS, CAL_SEEDS, append=False)
    write_phase(name, EXP_GROUPS, CAL_SEEDS, append=True)
    check('read_n_seeds_file unchanged when every group agrees',
          sm.read_n_seeds_file(name) == {'n_seeds': CAL_SEEDS})


# ===========================================================================
def test_duplicate_group_is_refused():
    name = 'dup'
    dm.create_directories(name)
    write_phase(name, CAL_GROUPS, CAL_SEEDS, append=False)
    raises('appending a group name that already exists is refused',
           lambda: write_phase(name, [{'name': 'calib_long__SOT',
                                       'nucleotide_substitution_rate': 1e-5}],
                               EXP_SEEDS, append=True),
           naming=('calib_long__SOT',))


# ===========================================================================
def test_repeats_and_dispatch(name):
    '''The state layer was already group-scoped. Prove it end to end.'''
    first_group = 'calib_long__SOT'
    new_group = 'exp__SOT'

    jobs.write_repeats(name)
    cal_repeats = jobs.read_repeats(name, first_group)
    check('calibration group has n_simulations x cal_seeds repeats',
          len(cal_repeats) == 3 * CAL_SEEDS, f'{len(cal_repeats)}')
    check('experiment group has its own, at exp_seeds',
          len(jobs.read_repeats(name, new_group)) == 1 * EXP_SEEDS)

    # mark the calibration phase finished, as a real run would
    for record in cal_repeats:
        jobs.set_state(name, first_group, record['index'], 'completed')

    identity_before = jobs.repeat_identity(jobs.read_repeats(name, first_group))
    jobs.write_repeats(name)
    check('rewriting repeats does not renumber a group that carries state',
          jobs.repeat_identity(jobs.read_repeats(name, first_group))
          == identity_before)

    import simplicity.runners.slurm as slurm
    try:
        slurm.raise_on_stale_state(name, [new_group])
        check("the stale guard ignores another phase's completed state", True)
    except Exception as exc:                                    # noqa: BLE001
        check("the stale guard ignores another phase's completed state",
              False, f'{type(exc).__name__}: {exc}')
    raises('but still refuses when that phase is included',
           lambda: slurm.raise_on_stale_state(name, None))

    submission_id, pairs = jobs.write_submission(name, [new_group])
    check('a submission can carry only the new phase',
          {group for group, _ in pairs} == {new_group}, f'{len(pairs)} pair(s)')
    check('submissions keep incrementing', submission_id == '001', submission_id)
    second, _ = jobs.write_submission(name, [new_group])
    check('and the next one is 002', second == '002', second)


# ===========================================================================
def test_negative_control():
    '''Without the offset, the append collides. Proves the check has teeth.

    The collision needs the id AND the label to coincide, because the stem is
    both. Two phases that sweep the same parameter over overlapping ranges --
    a calibration and a later refinement of it -- produce exactly that: the
    shared grid point has identical parameters in identical declaration order,
    so only the id keeps the two output directories apart.

    NOTE, found writing this: the label is built by iterating
    record['parameters'] in INSERTION order, and the grouped path inserts a
    group's list-valued keys before fixed_params and its scalar keys after. So
    the same parameter set declared as a sweep and as a scalar yields
    NSR_0p0001_N_20_T_2 and N_20_T_2_NSR_0p0001 -- different stems for
    identical parameters. That accident hides the collision in the
    scalar-vs-list case; it is not protection, and it is why this control uses
    matching declaration shapes. Reported, not fixed here: out of 8a's scope.
    '''
    swept = [{'name': 'calib_std__standard_only',
              'nucleotide_substitution_rate': [1e-4, 2e-4]}]
    calibrated_onto_a_grid_point = [
        {'name': 'calib_std__refine', 'nucleotide_substitution_rate': [1e-4]}]

    name = 'collide'
    dm.create_directories(name)
    write_phase(name, swept, CAL_SEEDS, append=False)
    before = sm.read_simulations(name)
    standard = sm.read_standard_parameters_values()
    old_stems = {sm.simulation_stem(r, standard) for r in before}

    naive = sm.generate_experiment_settings(
        {'_scenario_groups': calibrated_onto_a_grid_point}, FIXED)  # start_id 0
    naive_stems = {sm.simulation_stem(r, standard) for r in naive}
    check('restarting ids at 0 DOES collide -- the output would be overwritten',
          bool(old_stems & naive_stems),
          f'shared: {sorted(old_stems & naive_stems)}')

    offset = sm.generate_experiment_settings(
        {'_scenario_groups': calibrated_onto_a_grid_point}, FIXED,
        start_id=sm.next_simulation_id(name))
    offset_stems = {sm.simulation_stem(r, standard) for r in offset}
    check('the offset separates them although the parameters are identical',
          not (old_stems & offset_stems),
          f'{sorted(offset_stems)} vs {sorted(old_stems)}')

    raises('and write_experiment_settings refuses the colliding ids outright',
           lambda: sm.write_experiment_settings(name, naive, EXP_SEEDS,
                                                append=True),
           naming=('already used',))


def main():
    print(__doc__.strip().splitlines()[0])
    with Scratch():
        print('\n--- append')
        name = test_append_keeps_the_first_phase()
        print('\n--- per-group seeds')
        test_per_group_seeds(name)
        test_scalar_still_works_when_groups_agree()
        print('\n--- refusals')
        test_duplicate_group_is_refused()
        print('\n--- repeats, stale guard, dispatch')
        test_repeats_and_dispatch(name)
        print('\n--- negative control')
        test_negative_control()

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILED:')
        for line in failures:
            print(f'  {line}')
    else:
        print('one experiment can gain a phase without disturbing the last one')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
