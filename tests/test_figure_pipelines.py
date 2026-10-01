'''
Verifies _scenarios resolves all three figure pipelines -- bound, unbound
(argmax consensus) and unbound_dist (distributional consensus) -- and that the
three write to distinct figure paths. The dist pipeline previously raised
SystemExit: no config module was registered for it.
'''
import os, sys, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PLOTS = os.path.join(REPO, 'scripts', 'nature_plots')
sys.path.insert(0, REPO)
sys.path.insert(0, PLOTS)

import _scenarios as sc
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))
import impact_long_shedders_unbound_config as ucfg

EXPECTED_LABELS = {
    sc.BOUND_EXP_NAME:         'bound',
    sc.UNBOUND_EXP_NAME:       'unbound',
    sc.UNBOUND_DIST_EXP_NAME:  'unbound_distribution',
}

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f"{label}: got {got!r}, want {want!r}")


def test_dist_name_matches_the_pipeline_config():
    '''The figure code must not hardcode a name the pipeline doesn't produce.'''
    built = ucfg.prod_exp_name('distribution')
    check('dist experiment name agrees with the config',
          sc.UNBOUND_DIST_EXP_NAME, built)


def test_every_pipeline_resolves():
    for name in EXPECTED_LABELS:
        got = sc.scenario_names(name)          # raised SystemExit for dist
        check(f'{name} scenarios', bool(got), True)
    check('both unbound pipelines define the same scenarios',
          sc.scenario_names(sc.UNBOUND_EXP_NAME),
          sc.scenario_names(sc.UNBOUND_DIST_EXP_NAME))
    check('unbound includes edge_case',
          'edge_case' in sc.scenario_names(sc.UNBOUND_EXP_NAME), True)


def test_labels():
    for name, want in EXPECTED_LABELS.items():
        check(f'label for {name}', sc.pipeline_label(name), want)
    check('arm_label is gone', hasattr(sc, 'arm_label'), False)
    check('_ARM_LABELS is gone', hasattr(sc, '_ARM_LABELS'), False)


def test_figure_paths_are_distinct():
    paths = {n: sc.figure_path(3, n, 2, 'pdf') for n in EXPECTED_LABELS}
    for n, p in paths.items():
        print(f'    {n} -> {os.path.basename(p)}')
    check('three pipelines, three filenames', len(set(paths.values())), 3)
    check('dist filename', os.path.basename(paths[sc.UNBOUND_DIST_EXP_NAME]),
          'Figure_3_unbound_distribution_#2.pdf')
    check('argmax filename unchanged',
          os.path.basename(paths[sc.UNBOUND_EXP_NAME]),
          'Figure_3_unbound_#2.pdf')
    seeded = sc.figure_path(4, sc.UNBOUND_DIST_EXP_NAME, 2, 'png', seed='3')
    check('seed still carried', os.path.basename(seeded),
          'Figure_4_unbound_distribution_#2_seed3.png')


def test_save_scripts_still_parse():
    '''Catches an import or f-string break in the renamed help text.'''
    for script in ['fig1_save', 'fig2_save', 'fig3_save', 'fig3S_save',
                   'fig4_save']:
        r = subprocess.run([sys.executable, os.path.join(PLOTS, script + '.py'),
                            '--help'], capture_output=True, text=True)
        check(f'{script} --help', r.returncode, 0)
        if r.returncode != 0:
            print(r.stderr.strip()[-400:])


if __name__ == '__main__':
    for t in [test_dist_name_matches_the_pipeline_config,
              test_every_pipeline_resolves, test_labels,
              test_figure_paths_are_distinct, test_save_scripts_still_parse]:
        print(f'\n{t.__name__}')
        t()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
