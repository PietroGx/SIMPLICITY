#!/usr/bin/env bash
# ============================================================================
# run_long_nsr_local_fit_test.sh -- the closed-loop test of the two-stage
# long-shedder NSR calibration, for one calib_long cohort (default SOT) or for
# all three (--scenario all).
#
# THE CLAIM UNDER TEST
#
# cal_1 fits OSR = A * NSR_long**B over the whole swept range and inverts it
# at target_osr_long. On run #3 that range spans 100x, the power law is only
# locally true, and the inversion lands 19% above what the measured means say
# at the crossing (0.000733807 against 0.000615). The proposal is to keep the
# coarse grid for locating the root and add a STAGE 2: a narrow grid within
# +/-1.5x of it, where A*x**B is locally valid, fitted and inverted exactly as
# the current code does.
#
# THE TEST IS THE CONFIRMATION RUN, NOT THE FIT
#
# A narrow fit with a healthy R^2 proves nothing -- the whole point of this
# finding is that R^2 and the fitted exponent stay put while the inverted
# LEVEL drifts by tens of percent (tests/analyze_long_nsr_local_fit.py prints
# that table). So this simulates both candidate rates and measures what each
# actually produces:
#
#   PASS 1  at NSR_cal (the stage-2 answer) the measured mean OSR brackets
#           the target 0.00205.
#   PASS 2  at the cohort's stored global-fit value the measured mean OSR is
#           displaced from the target IN THE DIRECTION THAT COHORT'S OWN LOCAL
#           EXPONENT PREDICTS, by a statistically resolvable amount. That is
#           the falsifiable half.
#
# [2] IS DIRECTIONAL, NOT ONE-SIDED. It used to ask only whether the old value
# read ABOVE the target, which is SOT's direction and not every cohort's: SOT's
# stored value sits +19.3% above its local root, HIV's -2.3% below it and
# edge_case's -2.0% below. For the two below, the theory predicts a measured
# OSR BELOW the target, so an "above" test was false however the simulation
# came out and reported FAIL for cohorts behaving exactly as predicted. The
# direction is now the sign of each cohort's own stored-vs-root displacement,
# read off the displacement block, never a fixed side and never keyed on a
# cohort's name. SOT's direction is still "above", so SOT passes under exactly
# the conditions it did before.
#
# THE THREE OUTCOMES OF [2], each reported in its own words rather than as one
# boolean, with the predicted gap printed beside the measured one so a reader
# sees whether the model of the error was right and not merely whether a flag
# flipped:
#
#   displaced as predicted    PASS -- the old value is genuinely wrong.
#   lands on the target       FAIL -- the INFORMATIVE NULL. This cohort had no
#                             error worth fixing, which is the result that
#                             kills the premise for it. Reported, not
#                             suppressed.
#   displaced the OTHER way   FAIL -- the third outcome, and the one that
#                             matters most: the local-exponent model of the
#                             error is itself wrong, so no cohort's number can
#                             be trusted until the sign is explained.
#
# A NULL IS ONLY EVIDENCE IF THE TEST COULD HAVE SEEN THE EFFECT. The predicted
# gap in units of the source run's SEM near the target is 14.6 for SOT, 9.0 for
# edge_case and 2.8 for HIV. A bootstrap 95% interval is +/-1.96 SEM, so a true
# displacement of d SEM is missed with probability about Phi(1.96-d): 0.1% at
# d=5 but 20% at d=2.8. Below MIN_SEPARATION_SEM (5.0, in
# tests/analyze_long_nsr_local_fit.py) the run says a null there is
# INCONCLUSIVE rather than letting it read as a clean FAIL.
#
# Each cohort's report opens with the displacement block: the stored value
# against the root, the OSR gap that implies at the local exponent, that gap in
# SEM, and the direction [2] will test. Read it before the verdict.
#
# WHAT IT RUNS
#
#   phase 1  calib_long__<cohort>          STEPS x SEEDS  narrow grid
#   phase 2  calib_long__<cohort>_confirm  2     x SEEDS  NSR_cal and the old
#                                                         global-fit value
#
# both phases per cohort, and with --scenario all that is every calib_long
# cohort the config enumerates -- SOT, HIV_low+HIV_high, edge_case. The list is
# not written here: it comes from impact_long_shedders_config's own
# unique_long_tau_r_long_pairs and group_name_for_pair, so a cohort added there
# is picked up and none can be named differently.
#
# Both phases are built by impact_long_shedders_config.build_cal1_settings --
# the cohort's own (tau_3_long, R_long) group, its own window, its own fixed
# params. No parameter dict is hand-built here and no scenario name is
# invented; the confirmation phase is the cohort's group name with a
# "_confirm" suffix, because a second phase needs a group of its own.
#
# The experiment is test_long_nsr_local_fit_#N, its own namespace: it can
# never collide with an impact_long_shedders_#N. YOU CHOOSE N. This script
# will not pick one -- a local Data/ is not the cluster's Data/, so nothing
# here can know which numbers are free. It refuses to run against an
# experiment number that already exists.
#
# ONE EXPERIMENT PER COHORT, so --scenario all wants three numbers. Phase 1
# dispatches with append=False, which rewrites the experiment from simulation
# id 0; cohorts sharing a number would make the first wipe the others, and
# re-running one cohort would then have to append into an experiment that
# already carries state, which jobs.write_repeats refuses. Separate numbers
# keep each cohort independently re-runnable and the free-number gate
# meaningful per cohort. Under --runner slurm each cohort gets its own
# orchestrator job.
#
# Usage, from the repo root:
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1 --check-only
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1 --yes
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1 --seeds 50 --steps 7
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1 --scenario edge_case
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1,2,3 --scenario all
#   bash tests/run_long_nsr_local_fit_test.sh --stamp 1 --workers 2
#
# --check-only runs every gate AND rehearses both phases of every cohort with
# the dispatch stubbed out, printing the exact grids and fixed params that
# would be submitted. It submits nothing.
#
# THE TWO STORED NUMBERS ARE READ PER COHORT, NOT PASSED IN.
#
#   the provisional root  log-log interpolation of --source-exp's own measured
#                         means (default impact_long_shedders_#3), through
#                         tests/analyze_long_nsr_local_fit.py's
#                         local_interpolation -- not through its global fit;
#   the old global-fit    the frozen <exp>_calibrated_long_nsr.csv row for this
#   NSR_long              cohort, through that tool's stored_global_fit, which
#                         reads it with long_nsr_calibration_plot's own
#                         read_calibrated_long_nsr, keyed on (tau_3_long,
#                         R_long) the way cal_2 reads it. Never refitted.
#
# Both are keyed on the SAME resolved cohort, so --scenario moves both at once.
# They used to be --old-nsr and --root, defaulting to SOT's values: changing
# --scenario alone simulated one cohort and judged it against another cohort's
# calibrated value, and the fit would still have looked healthy. Those flags
# remain as overrides, and the driver REFUSES one that contradicts the cohort
# it is given with, naming both and naming the cohort the value does belong to.
# They are refused outright with more than one cohort, where at least one
# pairing would have to go unchecked.
#
# The cohort's (tau_3_long, R_long) is also compared against the source run's
# own group key before anything is read, so a cohort cannot be simulated with
# another cohort's stored numbers even if the group names agree.
#
# --workers caps how many seeded simulations the multiprocessing runner holds
# at once. The default is per cohort, from a memory budget: the extended-window
# cohort (edge_case, 1095 d) has been measured near 1.57 GB per seed, so four
# want ~6.5 GB where a development box has ~4 free.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

STAMPS=""
SEEDS=30
STEPS=6
SPAN=1.5
TARGET=0.00205
# empty: read per cohort from the source run's own stored data. Only a value
# the user typed ends up in these, and the driver refuses one that contradicts
# the cohort it is given with.
OLD_NSR=""
PROV_ROOT=""
WORKERS=""
SCENARIO="SOT"
SOURCE_EXP='impact_long_shedders_#3'
RUNNER="slurm"
BOOT=2000
YES=0
CHECK_ONLY=0
JOB_TIME="${SIMPLICITY_LOCALFIT_TIME:-24:00:00}"
JOB_MEM="${SIMPLICITY_LOCALFIT_MEM:-8G}"
SIM_MEM="${SIMPLICITY_SLURM_MEM:-}"
SIM_TIME="${SIMPLICITY_SLURM_TIME:-}"

while [ $# -gt 0 ]; do
    case "$1" in
        --stamp)     STAMPS="$2"; shift ;;
        --seeds)       SEEDS="$2"; shift ;;
        --steps)       STEPS="$2"; shift ;;
        --span)        SPAN="$2"; shift ;;
        --target-osr)  TARGET="$2"; shift ;;
        --old-nsr)     OLD_NSR="$2"; shift ;;
        --scenario)    SCENARIO="$2"; shift ;;
        --source-exp)  SOURCE_EXP="$2"; shift ;;
        --root)        PROV_ROOT="$2"; shift ;;
        --workers)     WORKERS="$2"; shift ;;
        --runner)      RUNNER="$2"; shift ;;
        --boot)        BOOT="$2"; shift ;;
        --slurm-mem)   SIM_MEM="$2"; shift ;;
        --slurm-time)  SIM_TIME="$2"; shift ;;
        --yes|-y)      YES=1 ;;
        --check-only)  CHECK_ONLY=1 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

# Every driver invocation carries the same optional flags; built once so the
# checks, the rehearsal, the foreground run and the submitted job cannot pass
# a different set.
OPT_ARGS=(--seeds "$SEEDS" --steps "$STEPS" --span "$SPAN"
          --target-osr "$TARGET" --source-exp "$SOURCE_EXP"
          --runner "$RUNNER" --boot "$BOOT")
if [ -n "$SIM_MEM" ];   then OPT_ARGS+=(--slurm-mem "$SIM_MEM"); fi
if [ -n "$SIM_TIME" ];  then OPT_ARGS+=(--slurm-time "$SIM_TIME"); fi
if [ -n "$WORKERS" ];   then OPT_ARGS+=(--workers "$WORKERS"); fi
if [ -n "$OLD_NSR" ];   then OPT_ARGS+=(--old-nsr "$OLD_NSR"); fi
if [ -n "$PROV_ROOT" ]; then OPT_ARGS+=(--root "$PROV_ROOT"); fi

. tests/_suite.sh
resolve_python

# the one threshold, read from the module that applies it -- never a second copy
MIN_SEP=$("$PY" -c "import sys
sys.path[:0] = ['.', 'tests', 'scripts/experiments_pipeline']
import analyze_long_nsr_local_fit as la
print(f'{la.MIN_SEPARATION_SEM:g}')")

# The driver. Called twice: once with --rehearse (dispatch stubbed, nothing
# simulated) and once for real, from inside the submitted job.
driver () {
    "$PY" -u - "$@" <<'PYEOF'
import os
import sys
import argparse

import numpy as np
import pandas as pd

sys.path[:0] = ['.', 'tests', 'scripts/experiments_pipeline']

import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm
from simplicity.runme import run_experiment_phase
import simplicity.dir_manager as dm
import simplicity.tuning.evolutionary_rate as er
import simplicity.settings_manager as sm

import analyze_long_nsr_local_fit as la
from impact_long_shedders_config import (
    SCENARIOS, TAU_ROUND, CAL1_ISOLATED_FIXED_PARAMS, CAL1_FINAL_TIME,
    DEFAULT_CONSENSUS, build_cal1_settings, cal1_final_time, calib_long_group,
    group_name_for_pair, print_fixed_params, set_calibration_output_env,
    set_slurm_resource_env, unique_long_tau_r_long_pairs, USER_FIXED_PARAMS)

EXP_NAME = 'test_long_nsr_local_fit'
CONFIRM_SUFFIX = '_confirm'
RUNNERS = {'serial': simplicity.runners.serial,
           'multiprocessing': simplicity.runners.multiprocessing,
           'slurm': simplicity.runners.slurm}

# An explicit --old-nsr must BE the cohort's stored value: it is read from a
# frozen table, so anything else is a different cohort's row or a typo. A
# --root is provisional by construction and a reader may well type a rounded
# form of it, so it is held to a band wide enough for rounding and far too
# narrow to admit another cohort (the three roots differ by ~2x).
OLD_NSR_RTOL = 1e-6
ROOT_RTOL = 0.02

# Fallback band for --root when the source run has no stored per-seed frame for
# this cohort, so the derived root cannot be compared against: the frozen table
# is still readable, and a provisional root IS the cohort's stored value give
# or take the bias under test -- 19.3% for SOT and 2.3% for HIV on run #3. The
# nearest cross-cohort pairing of that run's roots and stored values is 32%
# out, so this separates them without refusing a real one.
ROOT_VS_STORED_MAX = 0.25

# Worker count for the multiprocessing runner. Peak RSS per seeded simulation
# is set by how many intra-host lineages a cohort accumulates, which tracks its
# simulated window: the cohorts on the standard CAL1_FINAL_TIME window stay
# small, while an extended-window cohort (edge_case, 1095 d) has been measured
# at ~1.57 GB. 4 of those want ~6.5 GB and a development box has ~4 free.
WORKER_CAP = 4
LOCAL_MEM_BUDGET_GB = float(os.environ.get('SIMPLICITY_LOCALFIT_BUDGET_GB', 4.0))
PEAK_GB_EXTENDED_WINDOW = 1.57
PEAK_GB_STANDARD_WINDOW = 0.5
WORKERS_ENV = 'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_MULTIPROCESS'


def peak_gb(final_time):
    return (PEAK_GB_EXTENDED_WINDOW if final_time > CAL1_FINAL_TIME
            else PEAK_GB_STANDARD_WINDOW)


def cohort_workers(final_time, budget=None):
    """How many seeded simulations of THIS cohort fit in the memory budget."""
    budget = LOCAL_MEM_BUDGET_GB if budget is None else budget
    return max(1, min(WORKER_CAP, int(budget // peak_gb(final_time))))


# =============================================================================
# COHORTS -- enumerated by the config, never listed here
# =============================================================================
def cohorts(sp):
    """Every cal_1 cohort, through the config's own pair enumeration and its
    own naming (unique_long_tau_r_long_pairs / group_name_for_pair /
    calib_long_group), so this test cannot invent a cohort the pipeline does
    not fit, name one differently, or miss one that is added."""
    out = []
    for pair in unique_long_tau_r_long_pairs(sp):
        name = group_name_for_pair(pair, sp)
        out.append({'cohort': name,
                    'group': calib_long_group(name),
                    'pair': pair,
                    'tau_3_long': pair[0],
                    'R_long': pair[1],
                    'final_time': cal1_final_time(pair[0], sp),
                    'scenarios': [s['name'] for s in SCENARIOS
                                  if s['name'] in name.split('+')]})
    return out


def resolve_cohort(selector, sp):
    """One cohort from a selector: its own name ('HIV_low+HIV_high'), any
    scenario name inside it ('HIV_low'), or its full group name. Replaces the
    old substring match of the scenario against the source run's group names,
    which would have accepted a selector naming no cohort of its own."""
    found = [c for c in cohorts(sp)
             if selector in (c['cohort'], c['group']) or
             selector in c['scenarios']]
    if len(found) != 1:
        raise SystemExit(
            f"--scenario '{selector}' does not name exactly one cal_1 cohort "
            f"(matched {[c['cohort'] for c in found]}). Cohorts: " +
            '; '.join(f"{c['cohort']} [scenarios {'+'.join(c['scenarios'])}]"
                      for c in cohorts(sp)) + "; or 'all'.")
    return found[0]


def source_group(source_exp, cohort):
    """The source run's group for this cohort, with the cohort KEY checked
    against the config's own.

    This is the repo's signature failure guarded: the old NSR_long and the
    provisional root are both read under this group, and if the stored group
    carries a different (tau_3_long, R_long) than the cohort about to be
    simulated, the test would judge one cohort against another's numbers and
    still report a healthy fit."""
    present = la.long_groups(source_exp)
    if cohort['group'] not in present:
        raise SystemExit(
            f"{source_exp} has no group '{cohort['group']}' for cohort "
            f"'{cohort['cohort']}'. calib_long groups present: {present}")
    stored = la.cohort_of(source_exp, cohort['group'])
    if stored != cohort['pair']:
        raise SystemExit(
            f"COHORT MISMATCH: this pipeline's '{cohort['cohort']}' is "
            f"(tau_3_long, R_long)={cohort['pair']}, but {source_exp} group "
            f"'{cohort['group']}' is {stored}. The stored NSR_long and the "
            f"provisional root would be read for a different cohort than the "
            f"one simulated. Refusing.")
    return cohort['group']


# =============================================================================
# THE TWO STORED NUMBERS, PER COHORT
# =============================================================================
def stored_old_nsr(source_exp, group, target, model):
    """(stored NSR_long, {A, B}) for one cohort, through
    analyze_long_nsr_local_fit.stored_global_fit -- which reads the frozen
    <exp>_calibrated_long_nsr.csv through c1.read_calibrated_long_nsr, keyed
    on (tau_3_long, R_long) the way cal_2 reads it, plus the run's own A/B.
    Never refitted here: the baseline is the number the run used."""
    return la.stored_global_fit(source_exp, group, model, target)


def derived_root(source_exp, group, target):
    """(root, how, grid, loc): the provisional root by log-log interpolation
    of the source run's own measured means, through
    analyze_long_nsr_local_fit's local_interpolation -- not its global fit.

    Never extracts: this runs on a login node under --check-only, and
    re-extracting a cohort costs minutes to half an hour."""
    df_src, src = la.stored_or_cached_frame(source_exp, group)
    if df_src is None:
        return None, None, None, None
    grid = la.grid_summary(df_src)
    loc = la.local_interpolation(grid[la.NSR_COL].values.astype(float),
                                 grid['mean_osr'].values.astype(float),
                                 target)
    how = (f'log-log local interpolation of {source_exp} {group} ({src}), '
           f'local exponent B={loc["exponent"]:.4f}')
    return loc['nsr_loglog'], how, grid, loc


def no_frame_message(source_exp, group):
    return (f"{source_exp} {group} has no stored per-seed OSR frame. Run\n"
            f"  python tests/analyze_long_nsr_local_fit.py "
            f"--experiment {source_exp} --group {group}\n"
            f"first (it caches one), or pass --root.")


def refuse_mismatch(flag, given, derived, cohort, rtol, owners, noun):
    """An explicit override that contradicts the cohort's stored value is a
    mismatched pair -- the exact failure --old-nsr/--root as free-standing
    flags invited. Refuse, naming both, and name the cohort the given value
    does belong to when one does."""
    if given is None or derived is None:
        return
    if abs(given / derived - 1.0) <= rtol:
        return
    whose = [name for name, val in owners.items()
             if val is not None and name != cohort
             and abs(given / val - 1.0) <= rtol]
    raise SystemExit(
        f"MISMATCHED PAIR: {flag} {given!r} with --scenario '{cohort}'.\n"
        f"  cohort '{cohort}' has {noun} {float(derived)!r} "
        f"({100 * (given / derived - 1):+.1f}% off, tolerance "
        f"{100 * rtol:g}%)." +
        (f"\n  {given!r} is cohort '{whose[0]}'s {noun}."
         if whose else '') +
        f"\n  Drop {flag} and it is read from the stored data, or pass the "
        f"value this cohort was calibrated to. Refusing: simulating one "
        f"cohort against another's number is a confident PASS/FAIL about a "
        f"pairing that was never checked.")


def phase_settings(nsr_values, seeds, want, sp, suffix=''):
    """cal_1's own builder, subset to one cohort's group.

    The swept values go in through `ranges`, which the builder turns into
    np.geomspace(min, max, steps) -- so a 2-step range IS exactly its two
    endpoints, which is how the confirmation phase asks for two specific
    rates without hand-building anything. Checked below rather than assumed.
    """
    vals = sorted(float(v) for v in nsr_values)
    ranges = {'min': vals[0], 'max': vals[-1], 'steps': len(vals)}
    varying, fixed, n_seeds = build_cal1_settings(
        seeds, ranges, R=CAL1_ISOLATED_FIXED_PARAMS['R'],
        ih_virus_emergence_rate=USER_FIXED_PARAMS['IH_virus_emergence_rate'],
        consensus=DEFAULT_CONSENSUS)()

    kept = [g for g in varying['_scenario_groups']
            if (round(g['tau_3_long'], TAU_ROUND),
                round(g['R_long'], TAU_ROUND)) == want]
    if len(kept) != 1:
        raise SystemExit(f"no single cal_1 group matches cohort {want}: "
                         f"{kept}")
    group = dict(kept[0])
    built = [float(v) for v in group['nucleotide_substitution_rate_long']]
    if not np.allclose(built, vals, rtol=1e-12):
        raise SystemExit(f"builder produced {built}, asked for {vals}")
    if suffix:
        group['name'] = group['name'] + suffix
    return (lambda: ({'_scenario_groups': [group]}, fixed, n_seeds)), group, fixed


def mean_ci(values, n_boot, seed=0):
    """Mean of the measured OSR with a bootstrap 95% interval over seeds."""
    v = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    means = np.array([v[rng.integers(0, len(v), len(v))].mean()
                      for _ in range(n_boot)]) if n_boot > 0 else None
    out = {'mean': float(v.mean()), 'n': len(v),
           'sd': float(v.std(ddof=1)) if len(v) > 1 else float('nan')}
    out['sem'] = out['sd'] / np.sqrt(out['n'])
    if means is not None:
        out['lo'] = float(np.percentile(means, 2.5))
        out['hi'] = float(np.percentile(means, 97.5))
    else:
        out['lo'] = out['mean'] - 1.96 * out['sem']
        out['hi'] = out['mean'] + 1.96 * out['sem']
    return out


def at_value(df, nsr):
    """Non-outlier OSR rows at the grid point nearest `nsr`."""
    present = sorted(df[la.NSR_COL].unique())
    hit = min(present, key=lambda v: abs(v / nsr - 1.0))
    if abs(hit / nsr - 1.0) > 1e-6:
        raise SystemExit(f"no grid point at NSR_long {nsr}; present: {present}")
    sub = df[(df[la.NSR_COL] == hit) & (df['is_outlier'] == 0)]
    return sub[la.OSR_COL].values


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stamp', type=int, default=None)
    p.add_argument('--seeds', type=int, default=30)
    p.add_argument('--steps', type=int, default=6)
    p.add_argument('--span', type=float, default=1.5)
    p.add_argument('--target-osr', type=float, default=0.00205)
    p.add_argument('--old-nsr', type=float, default=None)
    p.add_argument('--scenario', default='SOT')
    p.add_argument('--source-exp', default='impact_long_shedders_#3')
    p.add_argument('--root', type=float, default=None)
    p.add_argument('--runner', default='slurm', choices=sorted(RUNNERS))
    p.add_argument('--boot', type=int, default=2000)
    p.add_argument('--model', default='exp')
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--slurm-mem', default=os.environ['SIMPLICITY_SLURM_MEM'])
    p.add_argument('--slurm-time', default=os.environ['SIMPLICITY_SLURM_TIME'])
    p.add_argument('--rehearse', action='store_true')
    p.add_argument('--list-cohorts', action='store_true',
                   help="print the cal_1 cohort names --scenario all covers")
    args = p.parse_args()

    sp = sm.read_standard_parameters_values()
    if args.list_cohorts:
        for c in cohorts(sp):
            print(c['cohort'])
        return 0
    if args.stamp is None:
        raise SystemExit("--stamp is required")

    # calibration reads four files; do not write the other three
    set_calibration_output_env()
    # per-SIMULATION resources, like every stage script's two flags; the
    # orchestrator job's own --time/--mem are set in the shell around this.
    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    numbered = f'{EXP_NAME}_#{args.stamp}'

    cohort = resolve_cohort(args.scenario, sp)
    tau_3_long, r_long = cohort['pair']
    src_group = source_group(args.source_exp, cohort)

    # ---- the two stored numbers, both read for THIS cohort ---------------
    # Keyed on the same resolved cohort, from the same source run, in one
    # place: a --scenario change moves both, so they cannot drift apart.
    stored_nsr, g_params = stored_old_nsr(args.source_exp, src_group,
                                          args.target_osr, args.model)
    root_derived, how_derived, grid, loc = derived_root(
        args.source_exp, src_group, args.target_osr)

    # The other cohorts' stored values, only when there is an override to name
    # against them. Cheap: one frozen table, and frames that are never
    # extracted.
    others_nsr, others_root = {}, {}
    if args.old_nsr is not None or args.root is not None:
        present = la.long_groups(args.source_exp)
        for c in cohorts(sp):
            if c['group'] not in present:
                continue
            try:
                others_nsr[c['cohort']] = la.stored_global_fit(
                    args.source_exp, c['group'], args.model, args.target_osr)[0]
            except SystemExit:
                others_nsr[c['cohort']] = None
            others_root[c['cohort']] = derived_root(
                args.source_exp, c['group'], args.target_osr)[0]

    refuse_mismatch('--old-nsr', args.old_nsr, stored_nsr, cohort['cohort'],
                    OLD_NSR_RTOL, others_nsr, 'stored NSR_long')
    refuse_mismatch('--root', args.root, root_derived, cohort['cohort'],
                    ROOT_RTOL, others_root, 'provisional root')

    old_nsr = stored_nsr if args.old_nsr is None else args.old_nsr
    how_old = (f'frozen {os.path.basename(la.c1._calibrated_long_nsr_path(args.source_exp))}'
               f' row for cohort {cohort["pair"]}' if args.old_nsr is None
               else f'given with --old-nsr (matches the stored {stored_nsr!r})')

    if args.root is not None and root_derived is None:
        # no frame to contradict it, so the frozen table is the only check left
        off = abs(args.root / stored_nsr - 1.0)
        if off > ROOT_VS_STORED_MAX:
            whose = [n for n, v in others_root.items()
                     if v is not None and n != cohort['cohort']
                     and abs(args.root / v - 1.0) <= ROOT_RTOL]
            raise SystemExit(
                f"MISMATCHED PAIR: --root {args.root!r} with --scenario "
                f"'{cohort['cohort']}'.\n"
                f"  that cohort's stored NSR_long is {float(stored_nsr)!r}, and "
                f"{args.root!r} is {100 * (args.root / stored_nsr - 1):+.1f}% "
                f"from it -- outside the {100 * ROOT_VS_STORED_MAX:g}% a "
                f"provisional root for this cohort can be." +
                (f"\n  {args.root!r} is cohort '{whose[0]}'s provisional "
                 f"root." if whose else '') + "\n  " +
                no_frame_message(args.source_exp, src_group).replace('\n', '\n  '))

    if args.root is not None:
        root, how = args.root, 'given with --root'
        if root_derived is None:
            how += (f' (no stored frame to derive one; checked only against '
                    f'the frozen {stored_nsr!r})')
    elif root_derived is None:
        raise SystemExit(no_frame_message(args.source_exp, src_group))
    else:
        root, how = root_derived, how_derived

    workers = (cohort_workers(cohort['final_time']) if args.workers is None
               else args.workers)
    os.environ[WORKERS_ENV] = str(workers)

    nsr_grid = la.stage2_grid(root, args.span, args.steps)
    floor = la.resolution_floor(tau_3_long)

    print()
    print('=' * 78)
    print(f' long-shedder NSR local-fit test -- {numbered}')
    print('=' * 78)
    print(f'cohort             {cohort["cohort"]}: tau_3_long={tau_3_long} d, '
          f'R_long={r_long}')
    print(f'cohort group       {src_group}  '
          f'(key {cohort["pair"]} matches {args.source_exp}\'s)')
    print(f'target OSR         {args.target_osr}')
    print(f'provisional root   {root:.9f}   [{how}]')
    print(f'stage-2 grid       {args.steps} points, +/-{args.span:g}x: '
          f'{nsr_grid[0]:.9f} .. {nsr_grid[-1]:.9f}')
    print('                   [' + ', '.join(f'{v:.9f}' for v in nsr_grid) + ']')
    print(f'old global-fit NSR {old_nsr:.9f}   '
          f'({100*(old_nsr/root-1):+.1f}% vs the provisional root)')
    print(f'                   [{how_old}]')
    print(f'resolution floor   {floor:.8f} subs/site/y  '
          f'(target is {args.target_osr/floor:.2f}x it)')
    print(f'seeds per point    {args.seeds}')
    print(f'simulations        {args.steps} x {args.seeds} (grid) + '
          f'2 x {args.seeds} (confirmation) = '
          f'{(args.steps + 2) * args.seeds} seeded simulations')
    print(f'runner             {args.runner}')
    how_workers = 'given with --workers' if args.workers is not None else (
        f'{LOCAL_MEM_BUDGET_GB:g} GB budget / ~{peak_gb(cohort["final_time"]):g} '
        f'GB peak at a {cohort["final_time"]} d window, capped at {WORKER_CAP}')
    print(f'workers            {workers} ({how_workers}) '
          f'-- read by the multiprocessing runner only')
    print(f'per-task slurm     mem={args.slurm_mem}, time={args.slurm_time}')
    if args.runner == 'slurm':
        want_gb = peak_gb(cohort['final_time'])
        have = args.slurm_mem.upper()
        have_gb = (float(have[:-1]) if have.endswith('G')
                   else float(have[:-1]) / 1024 if have.endswith('M')
                   else float('nan'))
        if np.isfinite(have_gb) and have_gb < want_gb:
            print(f'  WARNING          --slurm-mem {args.slurm_mem} is below '
                  f'this cohort\'s measured ~{want_gb:g} GB peak; Slurm will '
                  f'kill the task. Pass --slurm-mem {int(np.ceil(want_gb))+1}G.')
    print(f'writes             {numbered} only. {args.source_exp} is read '
          f'for its per-seed OSR frame and its frozen calibrated-NSR table '
          f'and nothing else -- no fit, nothing written back to it.')

    # ---- what the confirmation run can resolve --------------------------
    # Criterion [2] is DIRECTIONAL: the direction it predicts is the sign of
    # this cohort's own stored-value displacement from its local root, and
    # nothing else. Its magnitude is that displacement through the local
    # exponent, and that only separates the two rates if it clears the noise.
    print()
    print('--- displacement: does the new method move this cohort? ---')
    disp = 100 * (old_nsr / root - 1)
    pred_gap, sep_pred = float('nan'), float('nan')
    print(f'  stored NSR_long vs provisional root   {disp:+.2f}%')
    if loc is not None:
        b_loc = loc['exponent']
        pred_gap = 100 * ((old_nsr / root) ** b_loc - 1.0)
        near = grid.loc[(grid[la.NSR_COL] >= loc['lo']) &
                        (grid[la.NSR_COL] <= loc['hi'])]
        sem_pct = (100 * near['sem'] / near['mean_osr']).max() if not near.empty \
            else float('nan')
        print(f'  local exponent B                      {b_loc:.4f}  '
              f'(stored global fit B={g_params["B"]:.4f})')
        print(f'  implied gap in measured OSR           {pred_gap:+.2f}%  '
              f'= ({old_nsr/root:.4f})**B')
        n_near = int(near['n'].max()) if not near.empty else 0
        print(f'  SEM near the target, {args.source_exp}   {sem_pct:.2f}% '
              f'of the mean at {n_near} seeds')
        if np.isfinite(sem_pct) and sem_pct:
            sep_pred = abs(pred_gap) / sem_pct
        print(f'  separation                            '
              f'{f"{sep_pred:.1f} SEM" if np.isfinite(sep_pred) else "unknown"}')
    # the predicted direction, from the sign of the displacement above
    side = la.SIDE[0 if disp == 0 else (1 if disp > 0 else -1)]
    print(f'  [2] predicts the measured OSR lands   {side.upper()} the target')
    if np.isfinite(sep_pred) and sep_pred < la.MIN_SEPARATION_SEM:
        print(f'  NOTE  that is only {sep_pred:.1f} SEM, under the '
              f'{la.MIN_SEPARATION_SEM:.1f} SEM needed to call a null a')
        print('        refutation. If [2] comes back null here it is '
              'INCONCLUSIVE -- a')
        print('        displacement of the predicted size would often be '
              'missed at this noise.')

    # CLAUDE.md's stub-the-dispatch rehearsal: everything but the simulating.
    dispatch = ((lambda *a, **k: None) if args.rehearse
                else run_experiment_phase)

    # ---- phase 1: the narrow grid ----------------------------------------
    settings_1, group_1, fixed_1 = phase_settings(
        nsr_grid, args.seeds, cohort['pair'], sp)
    print()
    print(f'--- phase 1: narrow grid, group {group_1["name"]} ---')
    print_fixed_params(fixed_1)
    print(f'  final_time (group override): {group_1["final_time"]}')
    print(f'  NSR_long: {[round(v, 9) for v in group_1["nucleotide_substitution_rate_long"]]}')
    if args.rehearse:
        print('  [rehearsal] run_experiment_phase stubbed; nothing simulated.')

    dispatch(numbered, settings_1, simplicity_runner=RUNNERS[args.runner],
             append=False)

    if not args.rehearse:
        df1 = la.extract_per_seed(numbered, group_1['name'])
        g1 = la.grid_summary(df1)
        print()
        print(f'{"NSR_long":>12} {"mean OSR":>11} {"SD":>10} {"n":>4} {"x floor":>8}')
        for _, r in g1.iterrows():
            print(f'{r[la.NSR_COL]:12.9f} {r["mean_osr"]:11.6f} '
                  f'{r["sd_osr"]:10.6f} {int(r["n"]):4d} '
                  f'{r["mean_osr"]/floor:8.2f}')
        ok1 = df1[df1['is_outlier'] == 0]
        fit = er.fit_observed_substitution_rate_regressor(
            ok1, args.model, parameter_name=la.NSR_COL)
        nsr_cal = float(er.compute_calibrated_parameter(
            args.model, fit, args.target_osr))
        b_narrow = float(fit.params['B'].value)
        print()
        print(f'narrow fit         A={fit.params["A"].value:.8g}  B={b_narrow:.6f}')
        print(f'NSR_cal            {nsr_cal:.9f}  '
              f'({100*(nsr_cal/root-1):+.1f}% vs the provisional root, '
              f'{100*(nsr_cal/old_nsr-1):+.1f}% vs the old global fit)')
    else:
        nsr_cal, b_narrow = root, float('nan')
        print(f'  [rehearsal] NSR_cal unknown; using the provisional root '
              f'{root:.9f} to shape phase 2.')

    # ---- phase 2: the confirmation runs ----------------------------------
    settings_2, group_2, fixed_2 = phase_settings(
        [nsr_cal, old_nsr], args.seeds, cohort['pair'], sp,
        suffix=CONFIRM_SUFFIX)
    print()
    print(f'--- phase 2: confirmation, group {group_2["name"]} ---')
    print(f'  NSR_long: {[round(v, 9) for v in group_2["nucleotide_substitution_rate_long"]]}')
    if args.rehearse:
        print('  [rehearsal] run_experiment_phase stubbed; nothing simulated.')
        print()
        print('Rehearsal complete. Nothing was simulated and nothing was '
              'submitted.')
        # one line the shell's per-cohort rehearsal table reads back
        sep_txt = (f'{sep_pred:.1f}' + ('*' if sep_pred < la.MIN_SEPARATION_SEM
                                        else '')
                   if np.isfinite(sep_pred) else '?')
        print(f'SUMMARY|{cohort["cohort"]}|{numbered}|{float(old_nsr)!r}|'
              f'{float(root)!r}|{disp:+.2f}|{workers}|{side}|'
              f'{pred_gap:+.2f}|{sep_txt}')
        return 0

    dispatch(numbered, settings_2, simplicity_runner=RUNNERS[args.runner],
             append=True)

    df2 = la.extract_per_seed(numbered, group_2['name'])
    cal = mean_ci(at_value(df2, nsr_cal), args.boot, seed=1)
    old = mean_ci(at_value(df2, old_nsr), args.boot, seed=2)

    print()
    print('=' * 78)
    print(' confirmation')
    print('=' * 78)
    print(f'{"NSR_long":>14} {"mean OSR":>11} {"SD":>10} {"SEM":>10} '
          f'{"boot 95% CI":>26} {"n":>4}')
    for label, nsr, s in (('NSR_cal', nsr_cal, cal),
                          ('old global fit', old_nsr, old)):
        print(f'{nsr:14.9f} {s["mean"]:11.6f} {s["sd"]:10.6f} {s["sem"]:10.6f} '
              f'  {s["lo"]:.6f} .. {s["hi"]:.6f} {s["n"]:4d}   {label}')
    print(f'{"target":>14} {args.target_osr:11.6f}')

    pass1 = cal['lo'] <= args.target_osr <= cal['hi']
    predicted = (args.target_osr * (old_nsr / nsr_cal) ** b_narrow
                 if np.isfinite(b_narrow) else float('nan'))
    # The predicted gap comes from the displacement block, so its sign and the
    # direction [2] tests are the same number. Only with no source frame to
    # derive a local exponent from does this run's own narrow fit stand in.
    pred_pct = (pred_gap if np.isfinite(pred_gap) else
                100 * (predicted / args.target_osr - 1)
                if np.isfinite(predicted) else float('nan'))
    v = la.directional_verdict(args.target_osr, old['mean'], old['lo'],
                               old['hi'], old['sem'], disp, pred_pct, sep_pred)
    pass2 = v['pass2']

    print()
    print(f'[1] {"PASS" if pass1 else "FAIL"}  mean OSR at NSR_cal brackets '
          f'the target: {cal["lo"]:.6f} .. {cal["hi"]:.6f} '
          f'{"contains" if pass1 else "does NOT contain"} '
          f'{args.target_osr}  '
          f'(mean is {100*(cal["mean"]/args.target_osr-1):+.2f}% off target)')
    for line in v['lines']:
        print(line)
    if np.isfinite(predicted) and np.isfinite(pred_gap):
        print(f'      cross-check  this run\'s own narrow exponent '
              f'B={b_narrow:.3f}, pivoted at NSR_cal, predicts '
              f'{predicted:.6f} ({100*(predicted/args.target_osr-1):+.2f}%)')

    cal_note = ('' if pass1 else ' NSR_cal does not land on the target either, '
                'so something upstream of the inversion is off as well.')
    print()
    if v['opposite']:
        print(f'RESULT: THE MODEL OF THE ERROR IS WRONG. The old value is '
              f'displaced from the target, but the OPPOSITE way to the '
              f'prediction: {pred_pct:+.2f}% predicted, '
              f'{v["meas_pct"]:+.2f}% measured. This is the outcome that '
              f'matters most of the three -- it is the local-exponent account '
              f'of the error that has failed, not one cohort\'s number, so '
              f'the stage-2 grid is not justified by it anywhere. Do not '
              f'recalibrate; explain the sign first.' + cal_note)
        rc = 1
    elif pass1 and pass2:
        print('RESULT: the two-stage calibration is confirmed. The narrow-grid '
              'inversion lands on the target and the global-fit value does '
              'not -- it is displaced by about the amount, and in the '
              'direction, the local exponent predicts.')
        rc = 0
    elif pass2:
        print('RESULT: INCONCLUSIVE. The old value is displaced as predicted, '
              'but the narrow-grid answer does not land on the target either. '
              'The stage-2 window or the target itself is wrong; check the '
              'resolution floor line above before trusting either rate.')
        rc = 1
    elif v['null'] and v['low_power']:
        weakest = min(x for x in (sep_pred, v['sep_run']) if np.isfinite(x))
        print(f'RESULT: INCONCLUSIVE, UNDERPOWERED. The old value lands on the '
              f'target, but the predicted displacement was only '
              f'{weakest:.1f} '
              f'SEM, under the {la.MIN_SEPARATION_SEM:.1f} SEM this test '
              f'needs. So this run cannot tell "this cohort had no error" '
              f'from "it had one this many seeds could not see". Raise '
              f'--seeds and re-run before reading it either way.' + cal_note)
        rc = 1
    elif v['null'] and pass1:
        print(f'RESULT: PREMISE NOT SUPPORTED. Both rates produce an OSR on '
              f'the target, so the {disp:+.1f}% difference between them does '
              f'not show up in what is measured. This is the informative '
              f'null: this cohort had no error worth fixing, and the case for '
              f'a stage-2 grid does not rest on it. Do not recalibrate on '
              f'this result.')
        rc = 1
    elif v['null']:
        print('RESULT: INCONCLUSIVE. Neither rate lands on the target. '
              'Something upstream of the inversion is off; this test cannot '
              'settle the question as run.')
        rc = 1
    else:
        print('RESULT: [2] NOT TESTED. The stored value sits on this cohort\'s '
              'local root, so the local exponent predicts no displacement and '
              'there is nothing for [2] to confirm or refute.' + cal_note)
        rc = 1
    print()
    print(f'experiment: {dm.get_experiment_dir(numbered)}')
    tag = ('PASS' if v['pass2'] else 'OPPOSITE' if v['opposite']
           else ('NULL*' if v['low_power'] else 'NULL') if v['null'] else 'n/a')
    # one line the shell's per-cohort summary reads back
    print(f'VERDICT|{cohort["cohort"]}|{numbered}|'
          f'{"PASS" if pass1 else "FAIL"}|{tag}|'
          f'{disp:+.2f}|{100*(cal["mean"]/args.target_osr-1):+.2f}|'
          f'{100*(old["mean"]/args.target_osr-1):+.2f}|{pred_pct:+.2f}|'
          f'{v["outcome"]}|{nsr_cal!r}')
    return rc


if __name__ == '__main__':
    sys.exit(main())
PYEOF
}

# ============================== the submitted job ===========================
# One cohort per job: the re-exec below and the sbatch wrap both name a single
# cohort and a single --stamp, so this branch never sees 'all'.
if [ -n "${SIMPLICITY_LOCALFIT_JOB:-}" ]; then
    if [ "$SCENARIO" = "all" ]; then
        echo "internal error: a job was given --scenario all" >&2; exit 2
    fi
    echo "node        : $(hostname)"
    echo "slurm job   : ${SLURM_JOB_ID:-none}"
    echo "started     : $(date -Is)"
    # 'rc=$?' on its own line never runs under set -e: a failing driver would
    # take the script with it and the finish line would not print.
    if driver --stamp "$STAMPS" --scenario "$SCENARIO" "${OPT_ARGS[@]}"
    then
        rc=0
    else
        rc=$?
    fi
    echo "finished    : $(date -Is)"
    exit $rc
fi

# ============================ the cohorts to run ============================
if [ "$SCENARIO" = "all" ]; then
    COHORT_LIST=$(driver --list-cohorts) || {
        echo "could not enumerate the cal_1 cohorts" >&2; exit 1; }
else
    COHORT_LIST="$SCENARIO"
fi
mapfile -t COHORTS < <(printf '%s\n' "$COHORT_LIST")
IFS=',' read -r -a STAMP_ARR <<< "$STAMPS"

# ================================== checks ==================================
echo
echo "================================================================"
echo " long-shedder NSR local-fit test -- checks"
echo "================================================================"

fail=0
note() { printf "  %-40s %s\n" "$1" "$2"; }
bad()  { printf "  %-40s %s\n" "$1" "$2"; fail=1; }

note "interpreter" "$PY"
note "repo root" "$ROOT"

note "cohorts" "${#COHORTS[@]}: ${COHORTS[*]}"

# ONE --stamp PER COHORT. Each cohort is its own experiment because phase 1
# dispatches with append=False, which rewrites the experiment from simulation
# id 0: cohorts sharing a number would make the first one wipe the others, and
# re-running a single cohort would have to append into an experiment that
# already carries state -- which write_repeats refuses. Separate numbers keep
# each cohort independently re-runnable and keep the gate below meaningful per
# cohort. The script still picks no number: a local Data/ is not the cluster's.
if [ -z "$STAMPS" ]; then
    bad "--stamp" "MISSING -- you choose it; this script will not"
    echo "      A local Data/ is not the cluster's Data/, so nothing here can"
    echo "      know which experiment numbers are free. Pass --stamp N"
    echo "      (one per cohort, comma-separated: --stamp 1,2,3)."
elif [ "${#STAMP_ARR[@]}" -ne "${#COHORTS[@]}" ]; then
    bad "--stamp" "${#STAMP_ARR[@]} given, ${#COHORTS[@]} cohorts to run"
    echo "      One experiment per cohort. Pass --stamp with"
    echo "      ${#COHORTS[@]} comma-separated numbers, one for each of:"
    echo "      ${COHORTS[*]}"
elif [ "$(printf '%s\n' "${STAMP_ARR[@]}" | sort -u | wc -l)" \
       -ne "${#STAMP_ARR[@]}" ]; then
    bad "--stamp" "repeated: $STAMPS"
    echo "      One experiment per cohort means one NUMBER per cohort;"
    echo "      a repeat would make one cohort overwrite another."
elif printf '%s\n' "${STAMP_ARR[@]}" | grep -qv '^[0-9][0-9]*$'; then
    bad "--stamp" "not all integers: $STAMPS"
else
    for i in "${!COHORTS[@]}"; do
        note "experiment" "test_long_nsr_local_fit_#${STAMP_ARR[$i]}  (${COHORTS[$i]})"
    done
fi

case "$RUNNER" in serial|multiprocessing|slurm) note "runner" "$RUNNER" ;;
    *) bad "runner" "unknown: $RUNNER" ;;
esac

# Nothing may be appended to an experiment that already exists: phase 1 writes
# it from scratch (append=False) and would renumber a group that already
# carries state.
if [ -n "$STAMPS" ] && [ "${#STAMP_ARR[@]}" -eq "${#COHORTS[@]}" ]; then
    DATA_DIR=$("$PY" -I -c "
import sys; sys.path.insert(0, '.')
import simplicity.dir_manager as dm
print(dm.get_data_dir())" 2>/dev/null) || DATA_DIR=""
    if [ -z "$DATA_DIR" ]; then
        bad "experiment path" "could not be resolved"
    else
        # Not dm.get_experiment_dir: that raises unless the directory is
        # already there, which is the case this is checking for.
        for n in "${STAMP_ARR[@]}"; do
            EXP_DIR="$DATA_DIR/test_long_nsr_local_fit_#$n"
            if [ -d "$EXP_DIR" ]; then
                bad "experiment number free" "ALREADY EXISTS: $EXP_DIR"
                echo "      Pick another --stamp, or move that directory aside."
            else
                note "experiment number free" "yes (no $EXP_DIR)"
            fi
        done
    fi
fi

# A --root or --old-nsr can only ever belong to ONE cohort, so with more than
# one to run at least one pairing would go unchecked -- the exact mismatch
# these flags used to invite. The per-cohort values are read from the source
# run instead; the driver refuses a single-cohort override that contradicts it.
if [ "${#COHORTS[@]}" -gt 1 ]; then
    for flag_val in "--root:$PROV_ROOT" "--old-nsr:$OLD_NSR"; do
        if [ -n "${flag_val#*:}" ]; then
            bad "${flag_val%%:*}" "REFUSED with ${#COHORTS[@]} cohorts: ${flag_val#*:}"
            echo "      That value belongs to one cohort. Running"
            echo "      ${COHORTS[*]}"
            echo "      against it would judge at least two of them by another"
            echo "      cohort's number. Drop it, or name one --scenario."
        fi
    done
fi

if [ -n "$PROV_ROOT" ]; then
    note "provisional root" "$PROV_ROOT (given with --root)"
else
    note "provisional root" "per cohort from $SOURCE_EXP (see the rehearsals)"
fi
if [ -n "$OLD_NSR" ]; then
    note "old global-fit NSR" "$OLD_NSR (given with --old-nsr)"
else
    note "old global-fit NSR" "per cohort from ${SOURCE_EXP}'s frozen calibrated_long_nsr.csv"
fi

SIMS=$(( (STEPS + 2) * SEEDS ))
TOTAL_SIMS=$(( SIMS * ${#COHORTS[@]} ))
note "grid" "$STEPS points, +/-${SPAN}x, $SEEDS seeds each"
note "confirmation" "2 points (NSR_cal, the old global fit), $SEEDS seeds each"
note "simulation cost" "$SIMS per cohort x ${#COHORTS[@]} = $TOTAL_SIMS seeded simulations"
note "orchestrator job" "--time=$JOB_TIME --mem=$JOB_MEM per cohort (it blocks on both phases)"

if [ "$RUNNER" = "slurm" ] && [ "$CHECK_ONLY" -ne 1 ]; then
    if command -v sbatch >/dev/null 2>&1; then
        note "sbatch" "$(command -v sbatch)"
    else
        bad "sbatch" "NOT FOUND -- this needs a real Slurm controller"
        echo "      Pass --runner serial to run it on this machine instead."
    fi
fi

# The rehearsal, ONCE PER COHORT: the cohort resolved against the source run's
# own group key, its two stored numbers read, and both phases built through the
# config's own builder with run_experiment_phase stubbed. It runs before EVERY
# submission, not only under --check-only: it is what proves the source run is
# readable and the settings functions work, while that still costs nothing. It
# reads a stored per-seed frame and refuses to extract one, so it stays cheap
# enough for a login node.
REHEARSALS=""
if [ "$fail" -eq 0 ]; then
    for i in "${!COHORTS[@]}"; do
        echo
        echo "------- rehearsal (dispatch stubbed): ${COHORTS[$i]} -------"
        if out=$(driver --stamp "${STAMP_ARR[$i]}" \
                        --scenario "${COHORTS[$i]}" "${OPT_ARGS[@]}" --rehearse \
                        2>&1); then
            echo "$out"
            REHEARSALS+="$(printf '%s\n' "$out" | grep '^SUMMARY|')"$'\n'
        else
            echo "$out"
            fail=1
        fi
    done
fi

echo
if [ "$fail" -ne 0 ]; then
    echo "  One or more checks FAILED. Nothing submitted."
    exit 1
fi

if [ "${#COHORTS[@]}" -gt 1 ]; then
    echo "---------------- rehearsed, per cohort ----------------"
    printf "  %-22s %-12s %14s %14s %9s %4s %9s %9s %10s\n" \
        cohort experiment "stored NSR" "prov. root" "stored-root" "wrk" \
        "[2] pred" "pred gap" "sep (SEM)"
    printf '%s' "$REHEARSALS" | while IFS='|' read -r _ coh exp old root disp wrk side gap sep; do
        [ -n "$coh" ] || continue
        printf "  %-22s %-12s %14.9f %14.9f %8s%% %4s %9s %8s%% %10s\n" \
            "$coh" "${exp#test_long_nsr_local_fit}" "$old" "$root" "$disp" \
            "$wrk" "$side" "$gap" "$sep"
    done
    echo
    echo "  [2] pred   the direction criterion [2] tests for this cohort: the"
    echo "             sign of its own stored-vs-root displacement, carried"
    echo "             into OSR by its own local exponent. Not a fixed side."
    echo "  sep (SEM)  that gap in units of the source run's SEM near the"
    echo "             target. '*' marks a cohort under the $MIN_SEP SEM this"
    echo "             test needs, where a null is inconclusive rather than a"
    echo "             refutation of the premise."
    echo
fi
echo "  All checks passed."

[ "$CHECK_ONLY" -eq 1 ] && { echo "  --check-only: submitting nothing."; exit 0; }

if [ "$YES" -ne 1 ]; then
    echo
    read -r -p "Submit the test? $TOTAL_SIMS seeded simulations, ${#COHORTS[@]} cohort(s) x 2 phases. [y/N] " reply
    case "$reply" in
        y|Y|yes|YES) ;;
        *) echo "  Nothing submitted."; exit 0 ;;
    esac
fi

# One cohort per invocation, whichever way it runs: the cohorts are independent
# and each one's phase 2 depends only on its own phase 1.
mkdir -p Data/logs/tests
LOG_TIME="$(date +%Y%m%d_%H%M%S)"

if [ "$RUNNER" != "slurm" ]; then
    echo
    echo "  runner=$RUNNER: running here, in the foreground, cohort by cohort."
    rc=0
    for i in "${!COHORTS[@]}"; do
        LOGFILE="$ROOT/Data/logs/tests/long_nsr_local_fit_${LOG_TIME}_${STAMP_ARR[$i]}.log"
        echo
        echo "=================== ${COHORTS[$i]} ==================="
        echo "  log $LOGFILE"
        if ! SIMPLICITY_LOCALFIT_JOB=1 bash "$0" \
                --stamp "${STAMP_ARR[$i]}" --scenario "${COHORTS[$i]}" \
                --seeds "$SEEDS" --steps "$STEPS" --span "$SPAN" \
                --target-osr "$TARGET" --source-exp "$SOURCE_EXP" \
                --runner "$RUNNER" --boot "$BOOT" \
                ${SIM_MEM:+--slurm-mem "$SIM_MEM"} \
                ${SIM_TIME:+--slurm-time "$SIM_TIME"} \
                ${WORKERS:+--workers "$WORKERS"} \
                ${OLD_NSR:+--old-nsr "$OLD_NSR"} \
                ${PROV_ROOT:+--root "$PROV_ROOT"} 2>&1 | tee "$LOGFILE"; then
            rc=1
        fi
        # the driver's rc is lost through the pipe; the VERDICT line is not
        grep -q '^VERDICT|' "$LOGFILE" || rc=1
    done

    echo
    echo "================================================================"
    echo " per cohort"
    echo "================================================================"
    printf "  %-22s %-12s %9s %6s %9s %10s %10s %10s\n" \
        cohort experiment "stored-root" "[1]" "[2]" "cal vs tgt" \
        "old pred" "old meas"
    for i in "${!COHORTS[@]}"; do
        LOGFILE="$ROOT/Data/logs/tests/long_nsr_local_fit_${LOG_TIME}_${STAMP_ARR[$i]}.log"
        line=$(grep '^VERDICT|' "$LOGFILE" 2>/dev/null | tail -1)
        if [ -z "$line" ]; then
            printf "  %-22s %-12s %s\n" "${COHORTS[$i]}" \
                "#${STAMP_ARR[$i]}" "no VERDICT line -- it did not finish"
            rc=1
            continue
        fi
        IFS='|' read -r _ coh exp p1 p2 disp calpct oldpct predpct _ <<< "$line"
        printf "  %-22s %-12s %8s%% %6s %9s %9s%% %9s%% %9s%%\n" \
            "$coh" "${exp#test_long_nsr_local_fit}" "$disp" "$p1" "$p2" \
            "$calpct" "$predpct" "$oldpct"
        case "$p2" in
            OPPOSITE) N_OPP=$((${N_OPP:-0} + 1)) ;;
            NULL\*)   N_WEAK=$((${N_WEAK:-0} + 1)) ;;
            NULL)     N_NULL=$((${N_NULL:-0} + 1)) ;;
        esac
    done
    echo
    echo "  [1] mean OSR at NSR_cal brackets the target."
    echo "  [2] mean OSR at the old global-fit value is displaced from the"
    echo "      target in the direction that cohort's OWN local exponent"
    echo "      predicts, by a resolvable amount. Compare 'old pred' with"
    echo "      'old meas': the measurement is the second, the model of the"
    echo "      error it tests is the first."
    echo "        PASS      displaced as predicted -- the old value is wrong."
    echo "        NULL      lands on the target. The informative null: that"
    echo "                  cohort had no error worth fixing, which kills the"
    echo "                  premise for it. Do not recalibrate it."
    echo "        NULL*     the same, but the predicted separation was under"
    echo "                  $MIN_SEP SEM -- inconclusive, not a refutation."
    echo "        OPPOSITE  displaced the OTHER way. This refutes the model of"
    echo "                  the error itself and matters more than either of"
    echo "                  the other two outcomes. Read that log first."
    echo
    if [ "${N_OPP:-0}" -gt 0 ]; then
        echo "  OVERALL: ${N_OPP} cohort(s) came out displaced the OPPOSITE way"
        echo "           to the prediction. Read those logs before anything"
        echo "           else: what has failed is the account of the error,"
        echo "           not one cohort's number."
    elif [ "${N_NULL:-0}" -gt 0 ]; then
        echo "  OVERALL: ${N_NULL} cohort(s) landed ON the target at the old"
        echo "           value -- the informative null. The old value was"
        echo "           already right there; do not recalibrate those."
    elif [ "${N_WEAK:-0}" -gt 0 ]; then
        echo "  OVERALL: ${N_WEAK} cohort(s) are INCONCLUSIVE -- their null is"
        echo "           at a predicted separation under $MIN_SEP SEM, which"
        echo "           this many seeds cannot resolve either way."
    elif [ "$rc" -eq 0 ]; then
        echo "  OVERALL: every cohort met both criteria."
    else
        echo "  OVERALL: not every cohort met both criteria -- see the table."
    fi
    exit $rc
fi

# One orchestrator job PER COHORT: the slurm runner inside it submits each
# phase's array and polls it to completion, so phase 2 cannot start before
# phase 1's answer exists. Nothing heavy runs on the login node. Separate jobs
# rather than one doing all three in sequence: they are independent, so three
# share the queue instead of stacking into one wall clock.
echo
for i in "${!COHORTS[@]}"; do
    COH="${COHORTS[$i]}"
    N="${STAMP_ARR[$i]}"
    OUT="$ROOT/Data/logs/tests/long_nsr_local_fit_${LOG_TIME}_${N}_%j.log"
    JOB=$(sbatch --parsable \
        --job-name="long_nsr_local_fit_$N" \
        --time="$JOB_TIME" --mem="$JOB_MEM" --cpus-per-task=1 \
        --output="$OUT" --error="$OUT" \
        --wrap "cd '$ROOT' && SIMPLICITY_PYTHON='$PY' SIMPLICITY_LOCALFIT_JOB=1 \
                bash tests/run_long_nsr_local_fit_test.sh --stamp $N \
                --scenario '$COH' \
                --seeds $SEEDS --steps $STEPS --span $SPAN --target-osr $TARGET \
                --source-exp '$SOURCE_EXP' --runner $RUNNER --boot $BOOT \
                ${SIM_MEM:+--slurm-mem $SIM_MEM} ${SIM_TIME:+--slurm-time $SIM_TIME} \
                ${WORKERS:+--workers $WORKERS} \
                ${OLD_NSR:+--old-nsr $OLD_NSR} \
                ${PROV_ROOT:+--root $PROV_ROOT}") || JOB=""
    if [ -z "$JOB" ]; then
        echo "  sbatch refused the job for $COH. Nothing more is submitted." >&2
        exit 1
    fi
    LOGFILE="${OUT/\%j/$JOB}"
    echo "  $COH"
    echo "    experiment  test_long_nsr_local_fit_#$N"
    echo "    job         $JOB"
    echo "    log         $LOGFILE"
    JOBS="${JOBS:+$JOBS,}$JOB"
done
echo
echo "  watch      tail -f $ROOT/Data/logs/tests/long_nsr_local_fit_${LOG_TIME}_*"
echo "  cancel     scancel $JOBS"
echo "  result     per cohort, the [1] and [2] lines and the VERDICT line at"
echo "             the end of its log. [2] is directional: the displacement"
echo "             block above it names the direction that cohort's own local"
echo "             exponent predicts, and whether the predicted gap is big"
echo "             enough ($MIN_SEP SEM) for a null to mean anything."
