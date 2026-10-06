"""Per-cell overrides for the unbound pipeline, applied without editing scripts/.

The pipeline hardcodes the configuration the grid is supposed to be choosing:

    impact_long_shedders_config.py:150,203   population_size 1000, I0 10 / 50
    impact_long_shedders_unbound_config.py   SCENARIOS[*]["R_long"] = 1.1
    impact_long_shedders_unbound_config.py   UNBOUND_CAL2_FINAL_TIME = 365
    impact_long_shedders_unbound_config.py   NSR_RANGES[*]["steps"]

run_impact_long_shedders_unbound_pipeline runs each stage as its own SUBPROCESS,
so patching the config from the driver would reach nothing. Python imports
sitecustomize at startup in every one of them, which is the one hook that does.

Set by tests/test_calibrated_grid.py, read here:

    SIMPLICITY_CELL_POPULATION    population_size for cal_1, cal_2 and production
    SIMPLICITY_CELL_I0            infected_individuals_at_start
    SIMPLICITY_CELL_RLONG         R_long for every long-shedder scenario
    SIMPLICITY_CELL_CAL2_TIME     UNBOUND_CAL2_FINAL_TIME
    SIMPLICITY_CELL_NSR_STEPS     points in both NSR sweeps

With none of them set this module does nothing at all.

The config is patched the instant its module finishes executing, by wrapping the
loader -- not by polling or by importing it here first, either of which would
race the pipeline's own import and silently apply to nothing.
"""
import os
import sys

_KEYS = ('SIMPLICITY_CELL_POPULATION', 'SIMPLICITY_CELL_I0',
         'SIMPLICITY_CELL_RLONG', 'SIMPLICITY_CELL_CAL2_TIME',
         'SIMPLICITY_CELL_NSR_STEPS')
_TARGET = 'impact_long_shedders_unbound_config'
# Population and I0 live in the BOUND config and are imported into the unbound
# one by reference, so patching the dicts there covers both.
_BOUND = 'impact_long_shedders_config'


def _apply_population(module):
    population = os.environ.get('SIMPLICITY_CELL_POPULATION')
    i0 = os.environ.get('SIMPLICITY_CELL_I0')
    if population is None and i0 is None:
        return
    for name in ('USER_FIXED_PARAMS', 'CAL1_ISOLATED_FIXED_PARAMS',
                 'CAL2_FIXED_PARAMS'):
        params = getattr(module, name, None)
        if not isinstance(params, dict):
            continue
        if population is not None and 'population_size' in params:
            params['population_size'] = int(population)
        if i0 is not None and 'infected_individuals_at_start' in params:
            params['infected_individuals_at_start'] = int(i0)


def _apply_unbound(module):
    r_long = os.environ.get('SIMPLICITY_CELL_RLONG')
    if r_long is not None:
        for scenario in getattr(module, 'SCENARIOS', []):
            # control carries R_long None and must keep it: the builder reads
            # that as "this cohort does not shed long" rather than as a rate.
            if scenario.get('R_long') is not None:
                scenario['R_long'] = float(r_long)

    cal2_time = os.environ.get('SIMPLICITY_CELL_CAL2_TIME')
    if cal2_time is not None and hasattr(module, 'UNBOUND_CAL2_FINAL_TIME'):
        # build_cal2_settings reads this as a module global at call time, so
        # rebinding the attribute is enough.
        module.UNBOUND_CAL2_FINAL_TIME = int(float(cal2_time))

    steps = os.environ.get('SIMPLICITY_CELL_NSR_STEPS')
    ranges = getattr(module, 'NSR_RANGES', None)
    if steps is not None and isinstance(ranges, dict):
        for entry in ranges.values():
            if isinstance(entry, dict) and 'steps' in entry:
                entry['steps'] = int(steps)

    _apply_population(module)


def _install():
    import importlib.abc

    handlers = {_TARGET: _apply_unbound, _BOUND: _apply_population}

    class Finder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            handler = handlers.get(fullname)
            if handler is None:
                return None
            for finder in list(sys.meta_path):
                if finder is self:
                    continue
                found = finder.find_spec(fullname, path, target)
                if found is not None:
                    break
            else:
                return None
            loader = getattr(found, 'loader', None)
            if loader is None:
                return found
            original = loader.exec_module

            def exec_module(module):
                original(module)
                try:
                    handler(module)
                except Exception:
                    # never let an override take a pipeline stage down; a cell
                    # that ran unpatched shows up in the report as the wrong
                    # population, which is loud enough
                    pass

            try:
                loader.exec_module = exec_module
            except (AttributeError, TypeError):
                pass
            return found

    sys.meta_path.insert(0, Finder())


if any(os.environ.get(key) for key in _KEYS):
    try:
        _install()
    except Exception:
        pass
