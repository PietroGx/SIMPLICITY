#!/usr/bin/env python3
"""The experiment settings file: what the config loads from it, what the
pipeline records, and the reservoir it implies.

    python tests/test_experiment_settings.py

No simulation runs. Scratch files are named test_experiment_settings_* and
removed at exit.
"""
import argparse
import atexit
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
PIPE = os.path.join(REPO, 'scripts', 'experiments_pipeline')
sys.path[:0] = [REPO, PIPE]

import simplicity.dir_manager as dm                               # noqa: E402
import simplicity.population as pop                               # noqa: E402
import simplicity.settings_manager as sm                          # noqa: E402

import experiment_settings as es                                  # noqa: E402
import impact_long_shedders_config as cfg                         # noqa: E402
import run_impact_long_shedders_pipeline as runner                # noqa: E402

SCRATCH = tempfile.mkdtemp(prefix='test_experiment_settings_')
shutil.copytree(os.path.join(dm.get_data_dir(), '00_Reference_parameters'),
                os.path.join(SCRATCH, 'Data', '00_Reference_parameters'))
dm.set_data_dir(os.path.join(SCRATCH, 'Data'))
atexit.register(shutil.rmtree, SCRATCH, True)

fails = []


def check(label, got, want):
    ok = got == want
    if not ok:
        fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:60s} got={got!r}")


def default():
    return json.loads(json.dumps(es.load(es.DEFAULT_SETTINGS)))


def scratch_file(name, s):
    path = os.path.join(SCRATCH, f'test_experiment_settings_{name}.json')
    with open(path, 'w') as fh:
        json.dump(s, fh)
    return path


def test_config_reads_the_file():
    """Every value the config builds from comes out of the file it loaded."""
    s = es.load(cfg.SETTINGS_PATH)
    check('without the env var the config loads the default file',
          cfg.SETTINGS_PATH, os.path.abspath(es.DEFAULT_SETTINGS))
    p = s['population']
    check('population_size', cfg.POPULATION_SIZE, p['population_size'])
    check('infected_individuals_at_start', cfg.INFECTED_AT_START,
          p['infected_individuals_at_start'])
    check('USER_FIXED_PARAMS', cfg.USER_FIXED_PARAMS,
          {'population_size': p['population_size'],
           'infected_individuals_at_start': p['infected_individuals_at_start'],
           'R': p['R'], 'final_time': p['final_time'],
           'IH_virus_emergence_rate': p['IH_virus_emergence_rate']})
    check('cal_1 R and window', (cfg.CAL1_ISOLATED_FIXED_PARAMS['R'], cfg.CAL1_FINAL_TIME,
                                 cfg.CAL1_MIN_CYCLES),
          (s['cal_1']['R'], s['cal_1']['final_time'], s['cal_1']['min_cycles']))
    check('cal_2 window is the experiment phase\'s', cfg.CAL2_FINAL_TIME, p['final_time'])
    check('NSR ranges', cfg.NSR_RANGES,
          {'cal_long': s['cal_1']['nsr_long'], 'cal_std': s['cal_2']['nsr_std']})
    check('scenarios', cfg.SCENARIOS, s['scenarios'])
    check('the code-level default consensus does not follow the file',
          cfg.DEFAULT_CONSENSUS, 'distribution')


def test_validation():
    """A file the pipeline cannot use is refused, naming the field."""
    cases = {
        'unknown top-level key': (lambda s: s.update(extra=1), "unknown key 'extra'"),
        'missing section': (lambda s: s.pop('seeds'), "missing 'seeds'"),
        'misspelt nested key': (lambda s: s['population'].update(Rr=1.0), "unknown key 'Rr'"),
        'I0 not below N': (lambda s: s['population'].update(infected_individuals_at_start=3000),
                           'must be below population_size'),
        'bool is not an integer': (lambda s: s['seeds'].update(cal=True), 'seeds.cal'),
        'control with a shedding duration': (
            lambda s: s['scenarios'][0].update(inf_duration_long=63.0), 'must be null'),
        'long scenario without R_long': (lambda s: s['scenarios'][1].update(R_long=None),
                                         'R_long'),
        'duplicate scenario names': (lambda s: s['scenarios'][2].update(name='SOT'),
                                     'duplicate names'),
        'a name that would nest a directory': (
            lambda s: s['scenarios'][1].update(name='SOT/x'), 'name must be'),
        'prevalence above 1': (lambda s: s['scenarios'][1].update(long_shedders_ratio=1.5),
                               'long_shedders_ratio'),
        'sweep max not above min': (lambda s: s['cal_2']['nsr_std'].update(max=1e-6),
                                    'must exceed min'),
        'a one-point sweep': (lambda s: s['cal_1']['nsr_long'].update(steps=1), 'steps'),
        'bad Slurm memory': (lambda s: s['slurm'].update(mem='3 GB'), 'slurm.mem'),
        'bad Slurm time': (lambda s: s['slurm'].update(time='5h'), 'slurm.time'),
        'unknown consensus': (lambda s: s.update(consensus='mode'), 'consensus'),
        'the old schema, with exp_num': (lambda s: s.update(schema=1, exp_num=2), 'schema'),
    }
    check('the default file is valid', es.validate(default()), [])
    for label, (mutate, needle) in cases.items():
        s = default()
        mutate(s)
        errs = es.validate(s)
        check(f'refused: {label}', any(needle in e for e in errs), True)
    s = default()
    s['population']['R'] = 1
    check('an integer R is stored as a float', type(es.check(s)['population']['R']), float)


def test_env_var_reaches_a_fresh_interpreter():
    """What every stage subprocess does: import the config under the env var."""
    s = default()
    s['population'].update(population_size=10000, infected_individuals_at_start=300,
                           final_time=1825)
    s['scenarios'] = s['scenarios'][:2]
    path = scratch_file('n10000', s)
    code = ('import sys; sys.path[:0] = %r\n'
            'import json, impact_long_shedders_config as c\n'
            'v, f, _ = c.build_cal2_settings(1, c.NSR_RANGES["cal_std"])()\n'
            'g, f1, _ = c.build_cal1_settings(1, c.NSR_RANGES["cal_long"])()\n'
            'print(json.dumps([c.SETTINGS_PATH, f["population_size"], f["final_time"],\n'
            '  f1["population_size"], [x["name"] for x in c.SCENARIOS]]))\n' % [REPO, PIPE])
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                       env=dict(os.environ, **{es.SETTINGS_ENV: path}))
    out = json.loads(r.stdout.strip().splitlines()[-1]) if r.returncode == 0 else r.stderr[-300:]
    check('the subprocess loaded that file', out[0] if isinstance(out, list) else out, path)
    check('cal_2 builds at its population and window',
          out[1:3] if isinstance(out, list) else out, [10000, 1825])
    check('cal_1 builds at its population', out[3] if isinstance(out, list) else out, 10000)
    check('and its scenarios', out[4] if isinstance(out, list) else out, ['control', 'SOT'])

    s['population']['R'] = -1
    bad = scratch_file('invalid', s)
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                       env=dict(os.environ, **{es.SETTINGS_ENV: bad}))
    check('an invalid file stops the import', r.returncode != 0 and 'population.R' in r.stderr,
          True)


def test_pipeline_reads_settings_from_argv():
    """--settings is applied before the config loads: the orchestrator's own
    --help names, as the --settings default, the file the config read."""
    path = scratch_file('argv', default())
    script = os.path.join(PIPE, 'run_impact_long_shedders_pipeline.py')
    r = subprocess.run([sys.executable, script, '--settings', path, '--help'],
                       capture_output=True, text=True, cwd=REPO)
    check('--help runs', r.returncode, 0)
    check('the config loaded the file from argv',
          os.path.basename(path) in r.stdout.replace('\n', '').replace(' ', ''), True)
    check('activate_from_argv takes both spellings',
          (es.activate_from_argv(['--settings', es.DEFAULT_SETTINGS]),
           es.activate_from_argv([f'--settings={es.DEFAULT_SETTINGS}']),
           es.activate_from_argv(['--stamp', '3'])),
          (os.path.abspath(es.DEFAULT_SETTINGS),) * 2 + (None,))
    os.environ.pop(es.SETTINGS_ENV, None)


def test_stamps():
    """GLOSSARY.md: stamp. A new one is the current minute; an old number is
    still one; nothing that could nest or split a name is."""
    import re
    check('a new stamp is YYYYMMDD-HHMM',
          bool(re.fullmatch(r'\d{8}-\d{4}', cfg.new_stamp())), True)
    for ok in ('20261122-1533', '2', '910'):
        check(f'accepted: {ok}', cfg.check_stamp(ok), ok)
    for bad in ('a/b', 'x_y', '', '-1', '2026 11'):
        try:
            cfg.check_stamp(bad)
            got = 'accepted'
        except ValueError:
            got = 'refused'
        check(f'refused: {bad!r}', got, 'refused')

    # latest = most recently STARTED, not the highest or newest-named
    for stamp, started in (('20261122-1533', '2026-11-22T15:33:00'),
                           ('3', '2026-11-25T09:00:00'),
                           ('20261124-0800', '2026-11-24T08:00:00')):
        name = cfg.run_exp_name(stamp)
        os.makedirs(os.path.join(dm.get_data_dir(), name))
        with open(dm.get_experiment_state_path(name), 'w') as fh:
            json.dump({'experiment': name, 'started': started, 'completed': None}, fh)
    os.makedirs(os.path.join(dm.get_data_dir(), f'{cfg.EXP_NAME}_amx_#20991231-2359'))
    check('stamps on disk, oldest first', cfg.stamps_on_disk(),
          ['20261122-1533', '20261124-0800', '3'])
    with contextlib.redirect_stderr(io.StringIO()) as err:
        got = cfg.resolve_stamp('latest')
    check('latest is the last one started', got, '3')
    check('and it says which', 'impact_long_shedders_#3' in err.getvalue(), True)
    check('another mode has its own latest',
          cfg.resolve_stamp('latest', f'{cfg.EXP_NAME}_amx'), '20991231-2359')
    check('an explicit stamp is passed through', cfg.resolve_stamp('20261122-1533'),
          '20261122-1533')
    try:
        cfg.resolve_stamp('latest', 'test_experiment_settings_nothing')
        got = 'resolved'
    except SystemExit:
        got = 'refused'
    check('latest with nothing on disk is refused', got, 'refused')
    for d in os.listdir(dm.get_data_dir()):
        if d.startswith(cfg.EXP_NAME):
            shutil.rmtree(os.path.join(dm.get_data_dir(), d))


def _args(**kw):
    base = dict(stamp='991', r_cal1=1.0, r_cal2=1.05, ih_virus_emergence_rate=0.1,
                target_osr_std=0.0013, target_osr_long=0.00205, cal_seeds=30, exp_seeds=50,
                slurm_mem='3G', slurm_mem_cal2='1G', slurm_time='05:00:00')
    base.update(kw)
    return argparse.Namespace(**base)


def test_the_run_records_what_it_ran_with():
    eff = runner.effective_settings(_args(exp_seeds=7, r_cal2=1.1), 'distribution')
    check('a CLI override reaches the recorded settings', eff['seeds']['exp'], 7)
    check('so does R for cal_2 and the experiment phase', eff['population']['R'], 1.1)
    check('the stamp is not a setting', 'stamp' in eff, False)
    check('the rest is the file', eff['scenarios'], cfg.SCENARIOS)
    saved = os.environ.get(es.SETTINGS_ENV)
    try:
        path = runner.record_settings(eff, '991', 'distribution')
        check('written into the run\'s setup_data',
              path, os.path.join(cfg.setup_dir('991'), es.RUN_SETTINGS_FILENAME))
        check('and read back unchanged', es.load(path), eff)
        check('stages are pointed at the recorded copy', os.environ.get(es.SETTINGS_ENV), path)
        check('the same settings again are accepted',
              runner.record_settings(eff, '991', 'distribution'), path)
        other = runner.effective_settings(_args(exp_seeds=8), 'distribution')
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                runner.record_settings(other, '991', 'distribution')
            outcome = 'accepted'
        except SystemExit as exc:
            outcome = 'refused' if 'seeds' in str(exc) else str(exc)
        check('different settings into the same run are refused, naming the field',
              outcome, 'refused')
    finally:
        if saved is None:
            os.environ.pop(es.SETTINGS_ENV, None)
        else:
            os.environ[es.SETTINGS_ENV] = saved


def test_reservoir_scales_with_population():
    dm.declare_outputs()
    sp = sm.read_standard_parameters_values()
    for n in (50, 300):
        p = pop.create_population(dict(sp, population_size=n,
                                       infected_individuals_at_start=5, seed=0))
        check(f'N={n}: reservoir is {pop.RESERVOIR_PER_CAPITA} x N', p.reservoir,
              pop.RESERVOIR_PER_CAPITA * n)
        check(f'N={n}: every slot is an individual', len(p.individuals), p.reservoir)
        check(f'N={n}: the rest wait in the reservoir', len(p.reservoir_i), p.reservoir - n)
    check('the settings module reports the same number', es.reservoir(3000),
          pop.RESERVOIR_PER_CAPITA * 3000)


if __name__ == '__main__':
    print('-- the config reads the file --');      test_config_reads_the_file()
    print('-- validation --');                     test_validation()
    print('-- env var, fresh interpreter --');     test_env_var_reaches_a_fresh_interpreter()
    print('-- --settings on the orchestrator --'); test_pipeline_reads_settings_from_argv()
    print('-- stamps --');                         test_stamps()
    print('-- the run records its settings --');   test_the_run_records_what_it_ran_with()
    print('-- reservoir --');                      test_reservoir_scales_with_population()
    print('\nEXPERIMENT SETTINGS PASSED' if not fails else f'\nFAILED: {fails}')
    raise SystemExit(1 if fails else 0)
