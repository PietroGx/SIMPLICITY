#!/usr/bin/env python3
"""Verification for the two consensus modes, argmax and distribution.

Exercises the actual dispatch with run_experiment_phase_script stubbed,
per the repo rule that "the config is correct" is not verification.

Run from the repo root:
    python tests/test_consensus_pipelines.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'scripts', 'experiments'))

import impact_long_shedders_config as cfg
import impact_long_shedders_exp as exp

_fails = []


def check(label, got, want):
    ok = got == want
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:52s} got={got!r}")


ROW = {
    "scenario_name": "HIV_high", "tau_3_long": 109.0, "long_shedders_ratio": 0.12,
    "susceptibility_long": 1.0, "R": 1.03, "IH_virus_emergence_rate": 0.0,
    "R_long": 3.0, "nucleotide_substitution_rate_long": 8.8e-05,
    "nucleotide_substitution_rate": 1.4e-04,
}
CONTROL = dict(ROW, scenario_name="control", long_shedders_ratio=0.0,
               nucleotide_substitution_rate_long=float("nan"))


def test_builder():
    for mode in cfg.CONSENSUS_MODES:
        _, fixed, seeds = cfg.build_exp_scenario_settings(pd.Series(ROW), 30, mode)()
        check(f"builder carries consensus={mode}", fixed.get("consensus"), mode)
        check(f"builder seeds ({mode})", seeds, 30)
    _, fixed, _ = cfg.build_exp_scenario_settings(pd.Series(ROW), 30)()
    check("builder defaults to distribution",
          fixed.get("consensus"), "distribution")
    check("and DEFAULT_CONSENSUS says so", cfg.DEFAULT_CONSENSUS,
          "distribution")
    try:
        cfg.build_exp_scenario_settings(pd.Series(ROW), 30, "nonsense")
        check("bad consensus rejected", False, True)
    except ValueError:
        check("bad consensus rejected", True, True)
    # The DEFAULT mode owns the bare name; the other one is suffixed. Tied
    # together on purpose: flipping one without the other would silently put
    # default runs under a suffix and the investigation mode under the bare
    # name, and every reader keyed on the name would follow it.
    check("the default mode owns the bare experiment name",
          cfg.exp_name(cfg.DEFAULT_CONSENSUS), cfg.EXP_NAME)
    check("argmax is the suffixed one",
          cfg.exp_name("argmax"), cfg.EXP_NAME + "_amx")
    check("exactly one mode is unsuffixed",
          sum(1 for m in cfg.CONSENSUS_MODES
              if cfg.exp_name(m) == cfg.EXP_NAME), 1)
    check("every mode has a distinct name",
          len({cfg.exp_name(m) for m in cfg.CONSENSUS_MODES}),
          len(cfg.CONSENSUS_MODES))


def test_calibration_states_its_consensus():
    """Both calibration stages must WRITE consensus into fixed_params.

    They used to inherit it from standard_values.json, which is machine-local
    and gitignored -- and read_standard_parameters_values merges that file
    OVER the code defaults, so the mode a calibration ran under was a
    property of the host, not of the pipeline. Two machines could calibrate
    the same --exp-num differently and nothing would say so.
    """
    ranges = {'min': 1e-4, 'max': 1e-2, 'steps': 3}
    for label, build in (('cal_1', cfg.build_cal1_settings),
                         ('cal_2', cfg.build_cal2_settings)):
        _, fixed, _ = build(2, ranges)()
        check(f"{label} states a consensus mode",
              fixed.get("consensus"), cfg.DEFAULT_CONSENSUS)
        _, fixed, _ = build(2, ranges, consensus="argmax")()
        check(f"{label} honours an explicit override",
              fixed.get("consensus"), "argmax")

    # the whole point: calibration and the experiment phase agree by default
    _, cal1, _ = cfg.build_cal1_settings(2, ranges)()
    _, cal2, _ = cfg.build_cal2_settings(2, ranges)()
    _, exp3, _ = cfg.build_exp_scenario_settings(pd.Series(ROW), 2)()
    check("all three stages default to the same consensus mode",
          {cal1["consensus"], cal2["consensus"], exp3["consensus"]},
          {cfg.DEFAULT_CONSENSUS})


def test_orchestrator_refuses_two_modes():
    """One calibration cannot serve two consensus modes, so the orchestrator
    must refuse rather than silently calibrate under one and run under both."""
    import run_impact_long_shedders_pipeline as r
    argv = sys.argv
    sys.argv = ['run', '--exp-num', '999',
                '--consensus', 'distribution', 'argmax']
    try:
        r.main()
        check("two modes refused", False, True)
    except SystemExit as exc:
        text = str(exc)
        check("two modes refused, naming both",
              'argmax' in text and 'distribution' in text, True)
    except Exception as exc:                                    # noqa: BLE001
        check("two modes refused", f'{type(exc).__name__}: {exc}', 'SystemExit')
    finally:
        sys.argv = argv
    check("no stray pipeline log was written",
          os.path.isfile(os.path.join('Data', 'pipeline_logs',
                                      'impact_long_shedders_pipeline_#999.log')),
          False)


def test_dispatch():
    """Run exp.main() for both pipelines with the dispatcher stubbed."""
    table = pd.DataFrame([CONTROL, ROW])
    exp.load_calibration_table = lambda path: table
    exp.set_slurm_resource_env = lambda *a, **k: None

    assert hasattr(exp, 'run_experiment_phase_script'), (
        'the stage was rewired; this stub would not intercept and the test '
        'would run real simulations')

    for mode in cfg.CONSENSUS_MODES:
        seen = []

        def record(runner, exp_num, settings_func, name, append=True, _s=seen):
            varying, fixed, n_seeds = settings_func()
            groups = [g['name'] for g in varying['_scenario_groups']]
            _s.append((name, groups, [g.get('consensus') for g in
                                      varying['_scenario_groups']], n_seeds,
                       append))
            return groups

        exp.run_experiment_phase_script = record
        sys.argv = ['exp', '--exp-num', '99', '--runner', 'serial',
                    '--seeds', '30', '--consensus', mode]
        exp.main()

        # ONE dispatch for the whole phase since 8e, not one per scenario
        check(f"{mode}: one dispatch for the phase", len(seen), 1)
        name, groups, consensus, n_seeds, append = seen[0]
        check(f"{mode}: into the run's experiment", name, cfg.exp_name(mode))
        check(f"{mode}: both scenarios as groups", groups,
              [cfg.exp_group('control'), cfg.exp_group('HIV_high')])
        check(f"{mode}: consensus reached every group",
              sorted(set(consensus)), [mode])
        check(f"{mode}: seeds reached the phase", n_seeds, 30)
        check(f"{mode}: appended, never overwriting the calibration",
              append, True)

    # the two modes must not share an experiment name
    check("the two modes' experiment names do not collide",
          cfg.exp_name('argmax') == cfg.exp_name('distribution'), False)


# ---------------------------------------------------------------------------
# sanity plots: one grid per consensus pipeline
# ---------------------------------------------------------------------------
import io
import contextlib
import subprocess as _sp

import run_impact_long_shedders_pipeline as runner


def test_sanity_plots():
    """submit_sanity_plots must issue one sbatch per pipeline, each naming that
    pipeline's experiment, with the dispatcher stubbed."""
    calls = []

    class FakeResult:
        returncode = 0
        stdout = "Submitted batch job 12345"
        stderr = ""

    real_run = _sp.run
    runner.subprocess.run = lambda cmd, **k: (calls.append(cmd), FakeResult())[1]
    try:
        log = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()):
            got = runner.submit_sanity_plots(7, 0.0013, 0.00205, log,
                                             list(cfg.CONSENSUS_MODES))
    finally:
        runner.subprocess.run = real_run

    check("one sbatch per consensus pipeline", len(calls), len(cfg.CONSENSUS_MODES))
    for mode, cmd in zip(cfg.CONSENSUS_MODES, calls):
        check(f"{mode}: sbatch carries its own experiment name",
              cmd[-1], cfg.exp_name(mode))
        check(f"{mode}: exp_num and targets passed",
              cmd[2:5], ['7', '0.0013', '0.00205'])
    check("returned pipelines are tagged by mode",
          [m for m, _ in got], list(cfg.CONSENSUS_MODES))
    check("job ids returned for the waiter",
          all('12345' in s for _, s in got), True)


def test_archive_lists_both():
    """write_artifacts_archive must look for a sanity plot per pipeline."""
    import tempfile
    buf = io.StringIO()
    # In a temp cwd: the candidate paths are relative to Data/, and after the
    # v2.4.91 rename the argmax name collides with the bound run #7 still on
    # this machine -- a found file prints nothing, so the check would pass or
    # fail on what happens to be on disk.
    saved = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        arc = os.path.join(tmp, 'out', 'a.zip')
        os.chdir(tmp)
        try:
            with contextlib.redirect_stdout(buf):
                runner.write_artifacts_archive(arc, 7,
                                               os.path.join(tmp, 'no.log'))
        finally:
            os.chdir(saved)
    looked_for = buf.getvalue()
    for mode in cfg.CONSENSUS_MODES:
        want = f"{cfg.exp_name(mode)}_sanity_#7"
        check(f"archive looks for the {mode} sanity plot",
              want in looked_for, True)


def test_shell_passes_exp_name():
    """The sbatch wrapper must forward a 4th argument and default it."""
    sh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                           'scripts', 'experiments',
                           'submit_sanity_plot.sh')).read()
    check("shell reads a 4th argument with a default",
          'EXP_NAME="${4:-impact_long_shedders}"' in sh, True)
    check("shell forwards it to the plotting script",
          '--exp-name "$EXP_NAME"' in sh, True)
    check("hardcoded experiment name is gone",
          '--exp-name impact_long_shedders \\' in sh, False)


if __name__ == '__main__':
    print("-- builder --");            test_builder()
    print("-- calibration consensus --"); test_calibration_states_its_consensus()
    print("-- orchestrator guard --");  test_orchestrator_refuses_two_modes()
    print("-- dispatch (stubbed) --"); test_dispatch()
    print("-- sanity plots --");       test_sanity_plots()
    print("-- artifacts archive --");  test_archive_lists_both()
    print("-- sbatch wrapper --");     test_shell_passes_exp_name()
    print("\nCONSENSUS PIPELINES PASSED" if not _fails else f"\nFAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
