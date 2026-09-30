#!/usr/bin/env python3
"""Verification for the two-arm unbound pipeline.

Exercises the actual dispatch with run_experiment_script stubbed, per the
repo rule that "the config is correct" is not verification.

Run from the repo root:
    python tests/test_consensus_arms.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'scripts', 'experiments'))

import impact_long_shedders_unbound_config as cfg
import impact_long_shedders_unbound_exp as exp

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
    check("builder defaults to argmax", fixed.get("consensus"), "argmax")
    try:
        cfg.build_exp_scenario_settings(pd.Series(ROW), 30, "nonsense")
        check("bad consensus rejected", False, True)
    except ValueError:
        check("bad consensus rejected", True, True)
    check("argmax keeps the existing experiment prefix",
          cfg.prod_exp_name("argmax"), cfg.PROD_EXP_NAME)
    check("distribution arm is suffixed",
          cfg.prod_exp_name("distribution"), cfg.PROD_EXP_NAME + "_dist")


def test_dispatch():
    """Run exp.main() for both arms with the dispatcher stubbed."""
    table = pd.DataFrame([CONTROL, ROW])
    exp.load_calibration_table = lambda path: table
    exp.set_slurm_resource_env = lambda *a, **k: None

    for mode in cfg.CONSENSUS_MODES:
        seen = []
        exp.run_experiment_script = (
            lambda runner, exp_num, settings_func, name, _s=seen:
                _s.append((name, settings_func()[1].get("consensus"),
                           settings_func()[2])))
        sys.argv = ['exp', '--exp-num', '99', '--runner', 'serial',
                    '--seeds', '30', '--consensus', mode]
        exp.main()
        prefix = cfg.prod_exp_name(mode)
        check(f"{mode}: dispatched both scenarios",
              [n for n, _, _ in seen],
              [f"{prefix}_control", f"{prefix}_HIV_high"])
        check(f"{mode}: consensus reached every settings dict",
              sorted({c for _, c, _ in seen}), [mode])
        check(f"{mode}: seeds reached every settings dict",
              sorted({s for _, _, s in seen}), [30])

    # the two arms must not share experiment names
    a = {f"{cfg.prod_exp_name('argmax')}_{s}" for s in ('control', 'HIV_high')}
    b = {f"{cfg.prod_exp_name('distribution')}_{s}" for s in ('control', 'HIV_high')}
    check("the arms' experiment names do not collide", bool(a & b), False)


def test_bound_untouched():
    import impact_long_shedders_config as bound
    import inspect
    sig = inspect.signature(bound.build_exp_scenario_settings)
    check("bound pipeline's builder is unchanged",
          list(sig.parameters), ['row', 'n_seeds'])


# ---------------------------------------------------------------------------
# sanity plots: one grid per consensus arm
# ---------------------------------------------------------------------------
import io
import contextlib
import subprocess as _sp

import run_impact_long_shedders_unbound_pipeline as runner


def test_sanity_plots():
    """submit_sanity_plots must issue one sbatch per arm, each naming that
    arm's experiment, with the dispatcher stubbed."""
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

    check("one sbatch per consensus arm", len(calls), len(cfg.CONSENSUS_MODES))
    for mode, cmd in zip(cfg.CONSENSUS_MODES, calls):
        check(f"{mode}: sbatch carries its own experiment name",
              cmd[-1], cfg.prod_exp_name(mode))
        check(f"{mode}: exp_num and targets passed",
              cmd[2:5], ['7', '0.0013', '0.00205'])
    check("returned arms are tagged by mode",
          [m for m, _ in got], list(cfg.CONSENSUS_MODES))
    check("job ids returned for the waiter",
          all('12345' in s for _, s in got), True)


def test_archive_lists_both():
    """write_artifacts_archive must look for a sanity plot per arm."""
    import tempfile
    buf = io.StringIO()
    with tempfile.TemporaryDirectory() as tmp:
        arc = os.path.join(tmp, 'a.zip')
        with contextlib.redirect_stdout(buf):
            runner.write_artifacts_archive(arc, 7, os.path.join(tmp, 'no.log'))
    looked_for = buf.getvalue()
    for mode in cfg.CONSENSUS_MODES:
        want = f"{cfg.prod_exp_name(mode)}_sanity_#7"
        check(f"archive looks for the {mode} sanity plot",
              want in looked_for, True)


def test_shell_passes_exp_name():
    """The sbatch wrapper must forward a 4th argument and default it."""
    sh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                           'scripts', 'experiments',
                           'submit_sanity_plot_unbound.sh')).read()
    check("shell reads a 4th argument with a default",
          'EXP_NAME="${4:-impact_long_shedders_unbound}"' in sh, True)
    check("shell forwards it to the plotting script",
          '--exp-name "$EXP_NAME"' in sh, True)
    check("hardcoded experiment name is gone",
          '--exp-name impact_long_shedders_unbound \\' in sh, False)


if __name__ == '__main__':
    print("-- builder --");            test_builder()
    print("-- dispatch (stubbed) --"); test_dispatch()
    print("-- bound pipeline --");     test_bound_untouched()
    print("-- sanity plots --");       test_sanity_plots()
    print("-- artifacts archive --");  test_archive_lists_both()
    print("-- sbatch wrapper --");     test_shell_passes_exp_name()
    print("\nCONSENSUS ARMS PASSED" if not _fails else f"\nFAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
