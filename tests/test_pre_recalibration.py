#!/usr/bin/env python3
'''Everything that has to be settled before recalibrating, in one run.

    cd ~/SIMPLICITY && git pull
    python tests/test_pre_recalibration.py --exp-num 909

Two independent stages. Either can be run alone with --stages, and a failure in
one does not stop the other -- the report says what ran.

  grid       the one untested configuration point: 2% seeding at N=5000,
             R=1.06. 1% and 3% are both measured; 1% was chosen on evenness
             (scenario spread 10 against 60) and 3% buys 190 days of window at
             32 points of saturation. 2% is the only cell that could change
             that call, and the window is the single axis where the chosen
             configuration is worse than run #3 (595 days against 786).
             250 tasks, about an hour. Expect it to lose -- this is insurance,
             not a contender.

  consensus  settles the one claim about v2.4.57 that rests on an uncontrolled
             test. Its before/after captures were taken in separate processes
             before PYTHONHASHSEED was dealt with (v2.4.60), so string-hash
             order alone could have produced the 15 differing genealogy files
             that were attributed to the accumulator. This re-runs both sides
             at 6f5c185 (v2.4.56, the snapshot list) and 037e427 (v2.4.57, the
             accumulator) with the hash seed pinned, in throwaway git
             worktrees, so neither your working tree nor Data/ is touched.

Not a stage, deliberately: whether a 595-day window is long enough to fit the
intra-host clock. It cannot be answered by truncating run #3's output --
extract_ih_regression_data measures each lineage against its OWN host's founder
genome and time since that host's infection, so there is no absolute-time axis
to truncate without re-deriving the clock, which is how a diagnostic drifts
from the stage it diagnoses. Calibration answers it directly: if cal_1 fits
cleanly at the chosen configuration the window was enough, and if it does not,
the lever is the horizon (T 1095 -> 1460 turns a 595-day window into 960), not
the population.
'''
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

BEFORE = ('6f5c185', 'v2.4.56 snapshot list')
AFTER = ('037e427', 'v2.4.57 accumulator')
CONSENSUS_MODES = ('argmax', 'distribution')
# Written out but never read back into the model under the distributional
# consensus, so a flip here is cosmetic -- same carve-out test_consensus_wiring
# makes, kept identical on purpose.
ARGMAX_ARTEFACT = 'consensus_sequences_t.csv'
# Driven by compartment counts, which no fitness value can move. If these ever
# differ, something changed the dynamics rather than the rounding.
EPIDEMIC_FILES = ('final_time.csv', 'simulation_trajectory.csv')

W = 92
results = {}


def rule(char='-'):
    print(char * W)


def heading(title):
    rule()
    print(title)
    rule()


def hms(seconds):
    seconds = int(round(seconds))
    if seconds < 60:
        return f'{seconds}s'
    if seconds < 3600:
        return f'{seconds // 60}m {seconds % 60:02d}s'
    return f'{seconds // 3600}h {(seconds % 3600) // 60:02d}m'


# ------------------------------------------------------------------ stage: grid

def stage_grid(args):
    heading('STAGE grid      2% seeding at N=5000, R=1.06')
    command = [sys.executable, os.path.join(HERE, 'test_config_grid.py'),
               '--runner', args.runner, '--seeds', str(args.seeds),
               '--exp-num', str(args.exp_num),
               '--populations', '5000', '--i0-fractions', '0.02',
               '--r-values', '1.06']
    if args.slurm_mem:
        command += ['--slurm-mem', args.slurm_mem]
    if args.slurm_time:
        command += ['--slurm-time', args.slurm_time]
    print('  ' + ' '.join(command) + '\n')
    start = time.monotonic()
    code = subprocess.call(command, cwd=REPO)
    elapsed = time.monotonic() - start
    results['grid'] = ('ok' if code == 0 else f'FAILED (exit {code})', elapsed)
    print(f'\n  stage grid: exit {code} in {hms(elapsed)}')
    print(f'  report: Data/config_grid_report_#{args.exp_num}.txt')
    return code == 0


# ------------------------------------------------------- stage: consensus

def add_worktree(root, commit, label, slot):
    # named by SLOT, not by commit: a control run puts the same commit on both
    # sides, and naming by commit made the second checkout collide with the
    # first instead of comparing anything.
    path = os.path.join(root, f'{slot}_{commit}')
    print(f'  worktree {commit} ({label}) -> {path}')
    subprocess.check_call(['git', 'worktree', 'add', '--detach', path, commit],
                          cwd=REPO, stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    return path


def seed_worktree(worktree, template, experiment, seeds, consensus, horizon,
                  population):
    """Write the parameter set this comparison runs, as a plain json file.

    Deliberately NOT an on-disk experiment tree: this test runs one side inside
    a worktree at an older commit, and the two checkouts do not necessarily lay
    Data/ out the same way. Each side builds its own tree with its own code (see
    run_in_worktree), from this one file.
    """
    import json
    if template:
        with open(template) as handle:
            params = json.load(handle)
    else:
        import simplicity.settings_manager as sm
        params = dict(sm.read_standard_parameters_values())
    params['final_time'] = horizon
    params['population_size'] = population
    params['infected_individuals_at_start'] = max(1, population // 20)
    params['consensus'] = consensus
    params.pop('seed', None)            # the repeat supplies it
    path = os.path.join(worktree, 'comparison_parameters.json')
    with open(path, 'w') as handle:
        json.dump(params, handle, indent=1)
    return path


def run_in_worktree(worktree, experiment, params_path, seeds):
    """Build and run the experiment with the WORKTREE's own code.

    Through run_experiment, so this works on either side of a commit that
    changed how an experiment is laid out or how a repeat is addressed. The test
    compares outputs; it does not need to know either shape.
    """
    env = dict(os.environ)
    # pin the hash seed: without it two runs of IDENTICAL code disagree, which
    # is what made the original comparison unreadable. Drop anything that would
    # make the worktree import the main checkout instead of its own.
    env['PYTHONHASHSEED'] = '0'
    env.pop('PYTHONPATH', None)
    env.pop('SIMPLICITY_PROFILE_DIR', None)
    code = subprocess.call(
        [sys.executable, '-c',
         'import json, sys\n'
         'import simplicity.runme as runme\n'
         'import simplicity.runners.serial as serial\n'
         'name, path, seeds = sys.argv[1], sys.argv[2], int(sys.argv[3])\n'
         'fixed = json.load(open(path))\n'
         'runme.run_experiment(name, lambda: ({}, fixed, seeds),\n'
         '                     simplicity_runner=serial,\n'
         '                     archive_experiment=False)\n',
         experiment, os.path.relpath(params_path, worktree), str(seeds)],
        cwd=worktree, env=env,
        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return code == 0, None if code == 0 else params_path


def output_dirs(worktree, experiment):
    root = os.path.join(worktree, 'Data', experiment, '04_Output')
    found = {}
    for dirpath, _, files in os.walk(root):
        if 'final_time.csv' in files:
            found[os.path.basename(dirpath)] = dirpath
    return found


def compare(before_dir, after_dir, mode):
    '''Returns (identical, genealogy_differing, epidemic_differing, artefact).'''
    before, after = output_dirs(before_dir, mode), output_dirs(after_dir, mode)
    keys = sorted(set(before) & set(after))
    identical = genealogy = epidemic = artefact = 0
    for key in keys:
        names = sorted(set(os.listdir(before[key])) | set(os.listdir(after[key])))
        for name in names:
            if not name.endswith('.csv'):
                continue
            a = os.path.join(before[key], name)
            b = os.path.join(after[key], name)
            if not (os.path.isfile(a) and os.path.isfile(b)):
                genealogy += 1
                continue
            with open(a, 'rb') as f1, open(b, 'rb') as f2:
                if f1.read() == f2.read():
                    identical += 1
                elif name in EPIDEMIC_FILES:
                    epidemic += 1
                elif name == ARGMAX_ARTEFACT:
                    artefact += 1
                else:
                    genealogy += 1
    return identical, genealogy, epidemic, artefact, len(keys)


def stage_consensus(args):
    heading('STAGE consensus      v2.4.56 vs v2.4.57, hash seed pinned')
    if not os.path.isfile(args.params):
        print(f'  [skip] no template parameters file at {args.params}')
        print('         pass --params <a seeded simulation .json>')
        results['consensus'] = ('skipped (no template params)', 0.0)
        return False
    root = tempfile.mkdtemp(prefix='consensus_recheck_')
    start = time.monotonic()
    worktrees = []
    try:
        for slot, (commit, label) in enumerate((BEFORE, AFTER)):
            worktrees.append((commit, label,
                              add_worktree(root, commit, label, slot)))
        print()
        for commit, label, worktree in worktrees:
            for mode in CONSENSUS_MODES:
                params_path = seed_worktree(
                    worktree, args.params, mode, args.cseeds, mode,
                    args.chorizon, args.cpopulation)
                print(f'  running {commit} {mode:<13} '
                      f'{args.cseeds} seeds ...', end='', flush=True)
                ok, failed = run_in_worktree(worktree, mode, params_path,
                                             args.cseeds)
                print(' ok' if ok else f' FAILED on {failed}')
                if not ok:
                    print('    (a simulation failed in this worktree; the '
                          'comparison below will be incomplete)')
        print()
        before_dir = worktrees[0][2]
        after_dir = worktrees[1][2]
        verdict = []
        for mode in CONSENSUS_MODES:
            identical, genealogy, epidemic, artefact, seeds = compare(
                before_dir, after_dir, mode)
            print(f'  {mode:<13} {seeds} seeds compared: {identical} identical, '
                  f'{genealogy} genealogy differ, {epidemic} epidemic differ'
                  + (f', {artefact} {ARGMAX_ARTEFACT}' if artefact else ''))
            verdict.append((mode, identical, genealogy, epidemic))
        print()
        print('  How to read it, now that the hash seed is pinned and both')
        print('  sides are therefore deterministic:')
        print('    genealogy differ = 0  -> the accumulator was byte-neutral,')
        print('                             and the 15 differing files in the')
        print('                             original capture were hash order.')
        print('                             test_consensus_wiring\'s docstring')
        print('                             then needs correcting.')
        print('    genealogy differ > 0  -> the accumulator really does reorder')
        print('                             the genealogy; the original reading')
        print('                             stands and nothing changes.')
        print('    epidemic differ  > 0  -> neither explanation holds and the')
        print('                             dynamics moved. That would be the')
        print('                             serious case.')
        total_gen = sum(v[2] for v in verdict)
        total_epi = sum(v[3] for v in verdict)
        results['consensus'] = (
            f'{total_gen} genealogy / {total_epi} epidemic files differ',
            time.monotonic() - start)
        return total_epi == 0
    finally:
        for _, _, worktree in worktrees:
            subprocess.call(['git', 'worktree', 'remove', '--force', worktree],
                            cwd=REPO, stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)
        shutil.rmtree(root, ignore_errors=True)
        subprocess.call(['git', 'worktree', 'prune'], cwd=REPO,
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


# ---------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--stages', nargs='+', default=['grid', 'consensus'],
                        choices=['grid', 'consensus'])
    parser.add_argument('--exp-num', type=int, default=909)
    parser.add_argument('--seeds', type=int, default=50)
    parser.add_argument('--runner', default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--slurm-mem', default=None)
    parser.add_argument('--slurm-time', default=None)
    parser.add_argument('--params', default=None,
        help='template parameters file for the consensus stage; '
             'omitted means standard values plus the overrides below')
    parser.add_argument('--cseeds', type=int, default=3,
                        help='seeds per consensus mode (default %(default)s)')
    parser.add_argument('--chorizon', type=float, default=400.0)
    parser.add_argument('--cpopulation', type=int, default=1000)
    args = parser.parse_args()

    rule('=')
    print('PRE-RECALIBRATION CHECKS')
    rule('=')
    print(f'stages : {args.stages}')
    print()

    start = time.monotonic()
    if 'grid' in args.stages:
        stage_grid(args)
        print()
    if 'consensus' in args.stages:
        try:
            stage_consensus(args)
        except subprocess.CalledProcessError as exc:
            print(f'  stage consensus could not run: {exc}')
            results['consensus'] = (f'FAILED ({exc})', 0.0)
        print()

    heading('SUMMARY')
    for stage in args.stages:
        outcome, elapsed = results.get(stage, ('did not run', 0.0))
        print(f'  {stage:<12}{outcome:<52}{hms(elapsed):>12}')
    print(f'\n  total {hms(time.monotonic() - start)}')
    print()
    print('  Next, whatever these say: recalibrate at the chosen configuration')
    print('  (N=5000, I0=50, R=1.06). Every NSR in the frozen table was fitted')
    print('  at R=1.03, N=1000 and is wrong for it.')
    rule('=')


if __name__ == '__main__':
    main()
