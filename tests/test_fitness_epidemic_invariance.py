#!/usr/bin/env python3
'''Guard: a one-ulp fitness change must not move the epidemic.

WHAT THIS DOES AND DOES NOT ESTABLISH -- read before trusting it.

It perturbs fitness by one ulp at fitness_from_distance (the chokepoint every
phenotype model routes through) and requires the epidemic (final_time,
simulation_trajectory) byte-identical while the genealogy is free to differ.
It passes.

It does NOT currently discriminate. Run it against code where the
fitness-weighted parent draw shares rng4 with the recipient draw and it passes
too, because numpy consumes a FIXED amount of stream for both of the draws
involved:

    choice(seq, p=w)    same stream state afterwards for any w
    integers(0, n)      same stream state afterwards for any n

both verified directly. So fitness cannot displace rng4, and the coupling this
test was written to catch does not exist. It is kept as a regression guard: it
would fail if someone later introduced a draw whose stream consumption depends
on fitness (rejection sampling on a fitness-dependent bound, a variable-length
shuffle), which is a real way to reintroduce the problem silently.

STILL OPEN. Comparing 6f5c185 (v2.4.56) with 037e427 (v2.4.57) under a pinned
PYTHONHASHSEED, 4 of 12 epidemic files differ -- controlled against the same
commit on both sides, which gives 42/42 identical, so the difference is real.
Fitness displacing rng4 was the proposed explanation and it is now ruled out.
The mechanism is unknown. Until it is found, do not assume the epidemic is
stable across code changes; test_consensus_wiring reports rather than asserts
it for that reason.

    python tests/test_fitness_epidemic_invariance.py
'''
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

EPIDEMIC_FILES = ('final_time.csv', 'simulation_trajectory.csv')
# A parameters file to start from. None means "standard values plus the
# overrides below", which is what makes this runnable without a Data/ tree.
DEFAULT_PARAMS = None

# Patches fitness_from_distance -- the single chokepoint every phenotype model
# routes through -- to return the next representable double. math.nextafter is
# exactly one ulp, so this is the smallest perturbation that exists; anything
# the epidemic does with it, it would do with a reordered sum.
RUNNER = '''
import math, sys
import simplicity.phenotype.update as U
if "{perturb}" == "yes":
    original = U.fitness_from_distance
    def nudged(population, distance):
        return math.nextafter(original(population, distance), math.inf)
    U.fitness_from_distance = nudged
# Through the public entry point, so this test does not encode any on-disk
# layout: run_experiment writes whatever shape this checkout uses and the
# serial runner invokes it with whatever contract this checkout has.
import json
import simplicity.runme as runme
import simplicity.runners.serial as serial
name, params_path, seeds = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(params_path) as handle:
    fixed = json.load(handle)
runme.run_experiment(name, lambda: ({{}}, fixed, seeds),
                     simplicity_runner=serial, archive_experiment=False)
'''

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def run(experiment, params_template, seeds, horizon, population, perturb,
        consensus):
    """Build and run the experiment in one subprocess, through run_experiment."""
    import json
    import tempfile
    if params_template:
        with open(params_template) as handle:
            params = json.load(handle)
    else:
        import simplicity.settings_manager as sm
        params = dict(sm.read_standard_parameters_values())
    params['final_time'] = horizon
    params['population_size'] = population
    params['infected_individuals_at_start'] = max(1, population // 20)
    params['consensus'] = consensus
    params.pop('seed', None)        # the repeat supplies it, 0..seeds-1

    env = dict(os.environ)
    env['PYTHONHASHSEED'] = '0'
    env.pop('SIMPLICITY_PROFILE_DIR', None)

    handle, path = tempfile.mkstemp(suffix='.json')
    with os.fdopen(handle, 'w') as fh:
        json.dump(params, fh, indent=1)
    try:
        code = subprocess.call(
            [sys.executable, '-c', RUNNER.format(perturb=perturb),
             experiment, path, str(seeds)],
            cwd=REPO, env=env,
            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    finally:
        os.unlink(path)
    return code == 0


def outputs(experiment):
    root = os.path.join(REPO, 'Data', experiment, '04_Output')
    found = {}
    for dirpath, _, files in os.walk(root):
        if 'final_time.csv' in files:
            found[os.path.basename(dirpath)] = dirpath
    return found


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--params', default=DEFAULT_PARAMS)
    parser.add_argument('--seeds', type=int, default=3)
    parser.add_argument('--horizon', type=float, default=400.0)
    parser.add_argument('--population', type=int, default=1000)
    parser.add_argument('--consensus', default='distribution',
                        choices=['argmax', 'distribution'])
    args = parser.parse_args()

    if args.params and not os.path.isfile(args.params):
        sys.exit(f'no template parameters file at {args.params}')

    plain = 'zz_fitness_invariance_plain_#99'
    nudged = 'zz_fitness_invariance_nudged_#99'
    try:
        for experiment, perturb in ((plain, 'no'), (nudged, 'yes')):
            shutil.rmtree(os.path.join(REPO, 'Data', experiment),
                          ignore_errors=True)
            label = 'fitness nudged by one ulp' if perturb == 'yes' else 'baseline'
            print(f'  running {args.seeds} seeds, {label} ...',
                  end='', flush=True)
            ok = run(experiment, args.params, args.seeds, args.horizon,
                     args.population, perturb, args.consensus)
            print(' ok' if ok else ' FAILED')
            if not ok:
                sys.exit('a simulation failed; nothing to compare')

        before, after = outputs(plain), outputs(nudged)
        keys = sorted(set(before) & set(after))
        print(f'\nfitness perturbed by one ulp, {len(keys)} seeds, '
              f'{args.consensus} consensus')
        check('both sides produced the same seeds', len(keys), args.seeds)

        epidemic_differs = genealogy_differs = identical = 0
        for key in keys:
            names = sorted(set(os.listdir(before[key]))
                           | set(os.listdir(after[key])))
            for name in names:
                if not name.endswith('.csv'):
                    continue
                a = os.path.join(before[key], name)
                b = os.path.join(after[key], name)
                if not (os.path.isfile(a) and os.path.isfile(b)):
                    continue
                with open(a, 'rb') as f1, open(b, 'rb') as f2:
                    same = f1.read() == f2.read()
                if same:
                    identical += 1
                elif name in EPIDEMIC_FILES:
                    epidemic_differs += 1
                    print(f'      EPIDEMIC differs: {key}/{name}')
                else:
                    genealogy_differs += 1

        print(f'    {identical} identical, {genealogy_differs} genealogy '
              f'differing, {epidemic_differs} epidemic differing')
        check('the epidemic is immune to a one-ulp fitness change',
              epidemic_differs, 0)
        if genealogy_differs:
            print('    [expected] the genealogy moved: fitness steers the '
                  'parent draw, which is what a genealogy is')
        else:
            print('    [note] the genealogy did not move either -- this run was '
                  'too short or too small for the perturbation to reach a '
                  'draw boundary, so it does not exercise much. Try '
                  '--horizon 800 --population 2000.')
    finally:
        for experiment in (plain, nudged):
            shutil.rmtree(os.path.join(REPO, 'Data', experiment),
                          ignore_errors=True)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
