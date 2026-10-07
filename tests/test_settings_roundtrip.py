#!/usr/bin/env python3
'''Stage 0 of the settings data flow refactor: the oracle.

Builds a real experiment through runme.run_experiment with a stubbed runner and
asks one question of every output directory it produces -- do the parameters you
can read back from here match the parameters the experiment was defined with?

This must pass on the CURRENT code before any stage of the refactor moves, and
on every stage after. It is the only thing standing between the refactor and a
silent renaming bug, which is the failure mode this repo keeps hitting: four
production runs were burned on a measurement that described something other
than what its name claimed.

    python tests/test_settings_roundtrip.py

WHAT EACH CHECK IS FOR

  round trip        04_Output/<name>/ resolves to 02_Parameters/<name>.json by
                    STRING today. get_parameter_value_from_simulation_output_dir
                    has 43 call sites and keeps its signature through the
                    refactor, so this pins its meaning rather than its code.

  agrees with the   the value the resolver reports must equal the value the
  task's own input  task actually ran with -- i.e. what is in the seeded
                    parameters file. Two sources for one number is exactly how
                    they drift apart.

  counts            one output dir per parameter combination, n_seeds repeats
                    in each. Catches a combination that silently vanished.

  ordering          an unsorted os.walk currently decides which repeat is Slurm
                    array task 7, in two processes at two different times.
                    Nothing would report a mismatch, so this does.
                    Since stage 4 this is an equality against repeats.json
                    AND against the settings order, so it holds by
                    construction. It used to pass by filesystem coincidence:
                    both builds ran in one process, on one filesystem, creating
                    directories in the same sequence, so readdir agreed.

  group scoping     get_simulation_output_dirs(exp, group=...) must partition
                    the experiment: disjoint per group, union equal to the
                    unscoped listing. That is what lets a read say which of an
                    experiment's independent sweeps it means, instead of
                    silently mixing them into one plot series.

  path robustness   get_experiment_foldername_from_SSOD counted back a fixed
                    number of path components. Checked here against a path with
                    an extra level, which is what the old form got wrong.

  groups            all 8 production configs go through _scenario_groups, so
                    the grouped path is the one that matters. Checks that a
                    named group keeps its name, an unnamed one gets a
                    positional name, group order is declared order, ids are
                    dense, and 'name' never reaches 'parameters' -- which
                    check_parameters_names would reject.

  repeats           repeats.json is the Slurm array mapping, so it must be
                    derived from settings order rather than rediscovered: one
                    list per group, numbered from 0, simulations in id order,
                    seeds 0..n_seeds-1 within each. Also that the submission
                    map can span groups while each group stays 0-based, which
                    is what lets cal_2 keep submitting all 18 of its groups as
                    one array.

  collision         KNOWN BAD until stage 2, a real check since.
                    generate_filename_from_params encodes floats lossily, so
                    two swept values could round to ONE filename and the second
                    write_simulation_parameters call silently overwrote the
                    first -- half the sweep never ran, and the surviving
                    folder's name claimed a value it did not hold.
                    simulation_stem now prefixes the id, so colliding LABELS no
                    longer mean colliding simulations. The probe asks the real
                    namer for a collision rather than hardcoding a format,
                    because a diagnostic that assumes its stage's internals
                    stops describing that stage the moment they change.
'''
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

failures = []
known_bad = []
fixed = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')
    return ok


def check_known(label, got, want, why):
    '''A check expected to fail until a named stage fixes it.'''
    if got == want:
        print(f"  [FIXED] {label}: {got!r}")
        fixed.append(label)
    else:
        print(f"  [known] {label}: {got!r}, want {want!r}  <- {why}")
        known_bad.append(label)


# ---------------------------------------------------------------- environment

def redirect_data_dir(path):
    '''Point the whole package at a scratch Data/.

    Both modules, deliberately: dir_manager.set_data_dir mutates its own global,
    but settings_manager._data_dir (settings_manager.py:32) is a copy taken at
    import time and does not follow. Setting only the first gets you an
    experiment whose settings land in the real Data/ while its output lands in
    the temp one.
    '''
    import simplicity.dir_manager as dm
    import simplicity.settings_manager as sm
    dm.set_data_dir(path)
    sm._data_dir = path
    return dm, sm


def build(experiment_name, varying, fixed_params, n_seeds):
    '''Run the real setup path; stub out only the simulating.

    The stub stands in for unit_run.run_seeded_simulation: it creates the output
    directory through the real output_manager and stops there. That is faithful
    to how directories come to exist, without spending minutes on intra-host
    matrix exponentials that have nothing to do with this test.
    '''
    import simplicity.runme as runme
    import simplicity.output_manager as om
    import simplicity.settings_manager as sm

    import simplicity.jobs as jobs
    created = []

    class StubRunner:
        @staticmethod
        def run_seeded_simulations(exp, _run_seeded_simulation, groups=None):
            for group, record in jobs.all_repeats(exp, groups):
                created.append((group, record,
                                om.setup_output_directory(exp, group, record)))

    runme.run_experiment(experiment_name,
                         lambda: (varying, fixed_params, n_seeds),
                         simplicity_runner=StubRunner,
                         archive_experiment=False)
    return created


# --------------------------------------------------------------------- checks

def roundtrip(dm, sm, experiment_name, created, varying, n_seeds):
    settings = sm.read_experiment_settings(experiment_name)
    swept = sorted(varying)

    print('\nthe experiment was defined with')
    for key in swept:
        print(f'    {key} = {varying[key]}')

    sods = dm.get_simulation_output_dirs(experiment_name)
    check('one output dir per parameter combination',
          len(sods), len(settings))

    for sod in sorted(sods):
        reps = dm.get_seeded_simulation_output_dirs(sod)
        check(f'{os.path.basename(sod)[:34]}: {n_seeds} repeats',
              len(reps), n_seeds)

    # Every value the resolver reports must be a value the experiment actually
    # asked for -- compared exactly, because the whole question is whether full
    # precision survives the trip through a filename.
    print('\nround trip, exactly')
    for key in swept:
        wanted = set(varying[key])
        got = set()
        for sod in sods:
            got.add(sm.get_parameter_value_from_simulation_output_dir(sod, key))
        check(f'{key}: every swept value is readable back', got, wanted)

    # The resolver and the task read two different files for one number. This is
    # the check that notices if they ever stop agreeing.
    print("\nresolver agrees with the task's own input")
    # Still two files for one number: the task reads 02_Simulations/<stem>.json,
    # the resolver reads 04_Output/<group>/<stem>/parameters.json. This is the
    # check that notices if they ever stop agreeing.
    disagreements = []
    for group, record, output_dir in created:
        task_input = sm.read_simulation_parameters(experiment_name, record['stem'])
        sod = os.path.dirname(output_dir)
        for key in swept:
            resolved = sm.get_parameter_value_from_simulation_output_dir(sod, key)
            if resolved != task_input[key]:
                disagreements.append(
                    f'{os.path.basename(sod)}/{os.path.basename(output_dir)} '
                    f'{key}: resolver {resolved!r} vs task input {task_input[key]!r}')
    check('no output dir disagrees with what its task was handed',
          disagreements, [])


def ordering(sm, varying, fixed_params, n_seeds):
    """The order repeats are listed in IS the Slurm array mapping."""
    print('\nordering is recorded, not rediscovered')
    import simplicity.jobs as jobs
    first = build('order_probe_a', varying, fixed_params, n_seeds)
    build('order_probe_b', varying, fixed_params, n_seeds)

    def shape(name):
        return [f"{group}/{record['stem']}/{jobs.repeat_label(record)}"
                for group, record in jobs.all_repeats(name)]

    check('two identical experiments list their repeats in the same order',
          shape('order_probe_a'), shape('order_probe_b'))

    # The check above used to pass by filesystem coincidence: both builds ran in
    # one process on one filesystem and created their directories in the same
    # sequence, so readdir happened to agree. Now the order comes from
    # repeats.json, so it is checkable against the definition itself.
    expected = [f"{r['group']}/{sm.simulation_stem(r)}/seed_{seed:04d}"
                for r in sm.read_simulations('order_probe_a')
                for seed in range(n_seeds)]
    check('and that order is the settings order, simulations then seeds',
          shape('order_probe_a'), expected)

    check('the stub created one output directory per repeat',
          len(first), len(expected))


def grouped(dm, sm, fixed_params):
    """The _scenario_groups path: one sweep per group, each with its own axis."""
    print('\ngroups survive to disk')
    groups = [{'name': 'cal_long', 'tau_3_long': 133.5,
               'nucleotide_substitution_rate_long': [1e-4, 2e-4]},
              {'name': 'production', 'tau_3_long': 60.0, 'R_long': 1.06},
              {'tau_3_long': 200.0, 'R_long': 1.07}]          # deliberately unnamed

    build('grouped_probe', {'_scenario_groups': groups}, fixed_params, 2)

    check('group order is the declared order, unnamed named by position',
          dm.get_groups('grouped_probe'),
          ['cal_long', 'production', 'group_02'])
    check('every group carries an n_seeds',
          [g['n_seeds'] for g in sm.read_groups('grouped_probe')], [2, 2, 2])

    records = sm.read_simulations('grouped_probe')
    check('ids are dense and in order',
          [r['id'] for r in records], list(range(len(records))))
    check("each group's own sweep length is preserved",
          [r['group'] for r in records],
          ['cal_long', 'cal_long', 'production', 'group_02'])

    # check_parameters_names validates every key against STANDARD_VALUES, so an
    # identity key inside 'parameters' is not a style problem -- it raises.
    check("'name' never reaches parameters",
          [k for r in records for k in r['parameters'] if k == 'name'], [])
    check('the per-group fixed override still lands on its own group only',
          [r['parameters']['tau_3_long'] for r in records],
          [133.5, 133.5, 60.0, 200.0])
    check('the old flat contract still answers',
          len(sm.read_experiment_settings('grouped_probe')), len(records))

    print('\ngroup-scoped listing partitions the experiment')
    everything = dm.get_simulation_output_dirs('grouped_probe')
    check('unscoped listing is in definition order',
          [os.path.basename(d).split('__')[0] for d in everything],
          ['sim_000', 'sim_001', 'sim_002', 'sim_003'])

    per_group = {name: dm.get_simulation_output_dirs('grouped_probe', group=name)
                 for name in dm.get_groups('grouped_probe')}
    for name, dirs in per_group.items():
        print(f'    {name:<12} {[os.path.basename(d).split("__")[0] for d in dirs]}')
    check('the groups cover every simulation',
          sorted(d for dirs in per_group.values() for d in dirs),
          sorted(everything))
    check('and they are disjoint',
          sum(len(dirs) for dirs in per_group.values()), len(everything))
    check('cal_long holds exactly its two sweep points',
          len(per_group['cal_long']), 2)

    try:
        dm.get_simulation_output_dirs('grouped_probe', group='no_such_group')
        raised = False
    except ValueError:
        raised = True
    check('an unknown group raises instead of returning nothing', raised, True)

    print('\nthe experiment name survives a deeper tree')
    # the shape 04_Output would have if the group became a path component --
    # exactly the case path_parts[-4] got wrong, returning the group name
    check('with an extra level between 04_Output and the simulation',
          dm.get_experiment_foldername_from_SSOD(
              os.path.join('x', 'Data', 'myexp', dm.OUTPUT_DIRNAME,
                           'some_group', 'sim_000__a', 'seed_0000')),
          'myexp')
    check('and with the current shape',
          dm.get_experiment_foldername_from_SSOD(
              os.path.join('x', 'Data', 'myexp', dm.OUTPUT_DIRNAME,
                           'sim_000__a', 'seed_0000')),
          'myexp')

    sods = dm.get_simulation_output_dirs('grouped_probe')
    check('repeats come back sorted',
          [os.path.basename(r) for r in dm.get_seeded_simulation_output_dirs(sods[0])],
          ['seed_0000', 'seed_0001'])


def repeats(sm, jobs, experiment_name, n_seeds):
    """The repeat list and its state, per group, numbered from 0."""
    print('\nrepeats are explicitly ordered')
    groups = dm_groups = [g['name'] for g in sm.read_groups(experiment_name)]

    for group in groups:
        records = jobs.read_repeats(experiment_name, group)
        sims = [r['id'] for r in sm.read_simulations(experiment_name)
                if r['group'] == group]
        check(f'{group}: indices are dense and start at 0',
              [r['index'] for r in records], list(range(len(records))))
        check(f'{group}: simulations in id order, seeds within each',
              [(r['simulation'], r['seed']) for r in records],
              [(sid, seed) for sid in sims for seed in range(n_seeds)])
        check(f'{group}: every repeat carries its simulation stem',
              all(r['stem'] for r in records), True)

    print('\na submission spans groups; groups stay 0-based')
    submission_id, pairs = jobs.write_submission(experiment_name)
    check('the map covers every repeat of every group',
          len(pairs), sum(len(jobs.read_repeats(experiment_name, g))
                          for g in groups))
    check('group order then index order',
          pairs[:3], [[groups[0], 0], [groups[0], 1], [groups[0], 2]][:len(pairs)])
    check('each group restarts at index 0',
          sorted({g for g, _ in pairs}) == sorted(groups)
          and all(0 in [i for g, i in pairs if g == group] for group in groups),
          True)

    group, record = jobs.resolve_task(experiment_name, submission_id, 1)
    check('a task resolves its array position to one repeat',
          (group, record['index']), (pairs[1][0], pairs[1][1]))

    print('\nstate and progress are separate files')
    jobs.set_state(experiment_name, groups[0], 0, jobs.STARTED)
    check('state round-trips', jobs.state_of(experiment_name, groups[0], 0),
          jobs.STARTED)
    check('not resolved while started',
          jobs.is_resolved(experiment_name, groups[0], 0), False)
    jobs.set_state(experiment_name, groups[0], 0, jobs.COMPLETED, attempts=2)
    done = jobs.get_state(experiment_name, groups[0], 0)
    check('a later write keeps earlier fields', done.get('attempts'), 2)
    check('resolved when completed',
          jobs.is_resolved(experiment_name, groups[0], 0), True)
    check('counts are what the monitor reports',
          jobs.count_states(experiment_name).get(jobs.COMPLETED), 1)
    # the reason progress is its OWN file: a 30s read-modify-write of the state
    # record could land after a terminal COMPLETED and put STARTED back
    check('progress lives beside state, not inside it',
          os.path.dirname(jobs.progress_path(experiment_name, groups[0], 0))
          != os.path.dirname(jobs.state_path(experiment_name, groups[0], 0)),
          True)
    jobs.clear_state(experiment_name, groups[0], 0)
    check('clearing state removes it',
          jobs.state_of(experiment_name, groups[0], 0), None)


def collision(dm, sm, fixed_params):
    """Two swept values whose LABELS coincide must still be two simulations."""
    print('\ncolliding labels, distinct simulations')
    a = 0.0001041
    b = a * (1 + 1e-9)          # differs in the 10th figure: collides at any
                                # precision the label could reasonably use
    varying = {'nucleotide_substitution_rate': [a, b]}
    labels = [sm.generate_filename_from_params(
                  {'nucleotide_substitution_rate': v}) for v in (a, b)]
    print(f'    {a!r} -> {labels[0]}')
    print(f'    {b!r} -> {labels[1]}')
    check('the probe is valid: two distinct values, one label',
          (a != b, labels[0] == labels[1]), (True, True))

    build('collision_probe', varying, fixed_params, 1)
    sods = dm.get_simulation_output_dirs('collision_probe')
    check('a 2-point sweep produces 2 output dirs', len(sods), 2)

    got = {sm.get_parameter_value_from_simulation_output_dir(
               sod, 'nucleotide_substitution_rate') for sod in sods}
    check('both swept values are readable back, exactly', got, {a, b})

    check('their stems differ only in the id',
          sorted(os.path.basename(sod).split('__')[0] for sod in sods),
          ['sim_000', 'sim_001'])


# ----------------------------------------------------------------------- main

def main():
    root = tempfile.mkdtemp(prefix='roundtrip_')
    try:
        dm, sm = redirect_data_dir(root)
        print(f'scratch Data/ : {root}')

        # Small and cheap: nothing here runs a simulation, so only the setup
        # path's cost matters. The sweep values are deliberately far enough
        # apart to survive the 2-sig-fig name -- the colliding case is tested
        # separately, below.
        fixed_params = {'population_size': 50, 'final_time': 5}
        varying = {'nucleotide_substitution_rate': [0.0001, 0.0004, 0.0009],
                   'R': [1.05, 1.5]}
        n_seeds = 2

        created = build('roundtrip_probe', varying, fixed_params, n_seeds)
        check('the stub ran every repeat',
              len(created), len(varying['nucleotide_substitution_rate'])
              * len(varying['R']) * n_seeds)

        roundtrip(dm, sm, 'roundtrip_probe', created, varying, n_seeds)
        ordering(sm, varying, fixed_params, n_seeds)
        grouped(dm, sm, fixed_params)

        import simplicity.jobs as jobs
        repeats(sm, jobs, 'grouped_probe', 2)
        repeats(sm, jobs, 'roundtrip_probe', n_seeds)
        collision(dm, sm, fixed_params)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if fixed:
        print('NEWLY FIXED (a known-bad check started passing -- update this '
              'file and say which stage did it):')
        for label in fixed:
            print(f'  {label}')
    if known_bad:
        print(f'{len(known_bad)} known-bad check(s), expected until stage 2 '
              f'retires the filename as a lookup key:')
        for label in known_bad:
            print(f'  {label}')
    if failures:
        print(f'\n{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
        print('\nThe round trip is broken. Nothing in the refactor should '
              'proceed past this.')
    else:
        print('\nall non-known checks passed')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
