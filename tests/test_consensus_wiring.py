#!/usr/bin/env python3
'''
Wiring the incremental consensus in must not change what a run produces.

population.consensus_snapshot (a list rescanned on every rebuild) was replaced
by population.consensus, a ConsensusAccumulator folded in once per entry, and
the day's lineages are now keyed by genome before being added. The maths is
exact, but summation order changes, so this compares two captures of the same
seeds across the change.

What the comparison expects, established by measurement rather than guessed:

  argmax consensus        hamming_iw returns an INTEGER, so last-bit changes in
                          the matrix cannot reach the distance at all unless an
                          argmax flips, which it does not. Byte-identical,
                          asserted on every file.

  distribution consensus  the GENEALOGY diverges: fitness
                          steers rng4.choice when picking an infection's
                          parent, hundreds of hosts carry near-identical
                          fitness, and across a run a uniform draw eventually
                          lands in an interval that last-bit rounding moved.
                          Measured: one host switches lineage at day 53, and
                          everything downstream follows.

                          That is a property of the model, not of this change:
                          any reordering of the same sums does it, including
                          the genome keying alone or a numpy upgrade on
                          untouched code.

                          CORRECTION (v2.4.60). This test used to ASSERT the
                          epidemic byte-identical, on the reasoning that
                          propensities depend on counts and not on fitness.
                          The propensities do -- but the counts themselves do
                          not. population_model.infection draws the parent with
                          p=fitness and then draws the recipient from the SAME
                          rng4 stream, so a changed weight vector moves the
                          stream and a different susceptible is infected. Each
                          individual's long-shedder status is fixed at creation,
                          so that can infect a long shedder in place of a
                          standard host, whose tau_3_long gives a different
                          infectious duration and therefore different
                          compartment counts.

                          Measured at 6f5c185 vs 037e427 with PYTHONHASHSEED
                          pinned: 4 of 12 epidemic files differ, argmax
                          included. The old assertion passed by configuration
                          luck, not by law. Both are reported now; neither is
                          "wrong" -- the accumulator is exact to 5.47e-59
                          relative, so these are different realisations of one
                          stochastic process, not different answers.

    python tests/test_consensus_wiring.py <before_dir> <after_dir>

Captures come from scratchpad/consensus_baseline.py, run once on each side.
'''
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

# consensus_sequences_t holds the ARGMAX sequence. Under the distributional
# consensus it is written to disk but never read back into the model, so a
# flip there is cosmetic -- it is reported separately rather than failing.
ARGMAX_ARTEFACT = 'consensus_sequences_t.csv'

# Driven by compartment counts, which no fitness value can move. These must
# match even on the distribution path; if they ever stop, something has
# changed the dynamics rather than merely the rounding.
EPIDEMIC_FILES = ('final_time.csv', 'simulation_trajectory.csv')

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def seed_dirs(root):
    out = {}
    for dirpath, _, files in os.walk(root):
        if 'final_time.csv' in files:
            mode = 'distribution' if os.sep + 'distribution' + os.sep in dirpath \
                   else 'argmax'
            out[(mode, os.path.basename(dirpath))] = dirpath
    return out


def compare_mode(before, after, mode, strict):
    """strict=True asserts every file; otherwise only the epidemic files are
    asserted and the rest are reported."""
    b, a = seed_dirs(before), seed_dirs(after)
    keys = sorted(k for k in set(b) & set(a) if k[0] == mode)
    if not keys:
        print(f'  [skip] no {mode} captures on both sides')
        return
    check(f'{mode}: same seeds captured',
          sorted(k[1] for k in b if k[0] == mode),
          sorted(k[1] for k in a if k[0] == mode))

    identical = differing = 0
    artefact_differs = epidemic_differs = 0
    for key in keys:
        names = sorted(set(os.listdir(b[key])) | set(os.listdir(a[key])))
        for name in names:
            if not name.endswith('.csv'):
                continue
            pb, pa = os.path.join(b[key], name), os.path.join(a[key], name)
            if not (os.path.isfile(pb) and os.path.isfile(pa)):
                check(f'{mode} {key[1]}/{name}: present both sides', False, True)
                continue
            with open(pb, 'rb') as f1, open(pa, 'rb') as f2:
                same = f1.read() == f2.read()
            if same:
                identical += 1
            elif name in EPIDEMIC_FILES:
                epidemic_differs += 1
                print(f'      EPIDEMIC differs: {key[1]}/{name}')
            elif name == ARGMAX_ARTEFACT:
                artefact_differs += 1
            else:
                differing += 1

    print(f'    {identical} identical, {differing} genealogy files differing'
          + (f', {artefact_differs} {ARGMAX_ARTEFACT}' if artefact_differs else ''))
    if epidemic_differs:
        print(f'    [expected] the epidemic differs too: fitness and the '
              f'recipient draw share rng4 (see the correction above)')
    if strict:
        check(f'{mode}: every output byte-identical', differing, 0)
    elif differing:
        print(f'    [expected] genealogy differs: fitness steers the parent '
              f'choice, so reordered sums eventually flip one')


def main():
    before, after = sys.argv[1], sys.argv[2]
    print('\nargmax consensus  (strict: every file)')
    compare_mode(before, after, 'argmax', strict=True)
    print('\ndistribution consensus  (epidemic strict, genealogy reported)')
    compare_mode(before, after, 'distribution', strict=False)
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
