'''
Verifies _scenarios resolves both consensus modes -- impact_long_shedders
(distribution) and impact_long_shedders_amx (argmax) -- and that the two
write to distinct figure paths. The non-default mode once raised SystemExit:
no config module was registered for it.

v2.4.91 deleted the third, bound, entry with the pipeline behind it, and the
labels stopped naming an arm ("bound"/"unbound") and started naming the only
thing that still differs: the consensus mode. v2.4.92 then made distribution
the default and gave it the bare name, so the suffix moved to argmax.
'''
import os, sys, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PLOTS = os.path.join(REPO, 'scripts', 'nature_plots')
sys.path.insert(0, REPO)
sys.path.insert(0, PLOTS)

import _scenarios as sc
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))
import impact_long_shedders_config as ucfg

EXPECTED_LABELS = {
    sc.EXP_NAME:         'distribution',
    sc.ARGMAX_EXP_NAME:  'argmax',
}

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f"{label}: got {got!r}, want {want!r}")


def test_names_match_the_pipeline_config():
    '''The figure code must not hardcode a name the pipeline doesn't produce.'''
    check('argmax experiment name agrees with the config',
          sc.ARGMAX_EXP_NAME, ucfg.exp_name('argmax'))
    check('distribution experiment name agrees with the config',
          sc.EXP_NAME, ucfg.exp_name('distribution'))
    check('the default mode owns the bare name',
          ucfg.exp_name(ucfg.DEFAULT_CONSENSUS), ucfg.EXP_NAME)


def test_every_pipeline_resolves():
    for name in EXPECTED_LABELS:
        got = sc.scenario_names(name)          # raised SystemExit for dist
        check(f'{name} scenarios', bool(got), True)
    check('both consensus modes define the same scenarios',
          sc.scenario_names(sc.EXP_NAME),
          sc.scenario_names(sc.ARGMAX_EXP_NAME))
    check('and they include edge_case',
          'edge_case' in sc.scenario_names(sc.EXP_NAME), True)
    check('the deleted bound pipeline is not resolvable',
          'impact_long_shedders_bound' in sc._CONFIG_MODULES, False)


def test_labels():
    for name, want in EXPECTED_LABELS.items():
        check(f'label for {name}', sc.pipeline_label(name), want)
    check('arm_label is gone', hasattr(sc, 'arm_label'), False)
    check('_ARM_LABELS is gone', hasattr(sc, '_ARM_LABELS'), False)


def test_figure_paths_are_distinct():
    paths = {n: sc.figure_path(3, n, 2, 'pdf') for n in EXPECTED_LABELS}
    for n, p in paths.items():
        print(f'    {n} -> {os.path.basename(p)}')
    check('two modes, two filenames', len(set(paths.values())), 2)
    check('distribution filename', os.path.basename(paths[sc.EXP_NAME]),
          'Figure_3_distribution_#2.pdf')
    check('argmax filename',
          os.path.basename(paths[sc.ARGMAX_EXP_NAME]),
          'Figure_3_argmax_#2.pdf')
    seeded = sc.figure_path(4, sc.EXP_NAME, 2, 'png', seed='3')
    check('seed still carried', os.path.basename(seeded),
          'Figure_4_distribution_#2_seed3.png')


def test_save_scripts_still_parse():
    '''Catches an import or f-string break in the renamed help text.'''
    for script in ['fig1_save', 'fig2_save', 'fig3_save', 'fig2S_save',
                   'fig4_save']:
        r = subprocess.run([sys.executable, os.path.join(PLOTS, script + '.py'),
                            '--help'], capture_output=True, text=True)
        check(f'{script} --help', r.returncode, 0)
        if r.returncode != 0:
            print(r.stderr.strip()[-400:])


if __name__ == '__main__':
    for t in [test_names_match_the_pipeline_config,
              test_every_pipeline_resolves, test_labels,
              test_figure_paths_are_distinct, test_save_scripts_still_parse]:
        print(f'\n{t.__name__}')
        t()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
