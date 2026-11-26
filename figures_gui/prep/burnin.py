"""The burn-in day and the run-length gates, shared by the clade build, the
transmission table and the diagnostics in tests/. Moved here from
scripts/nature_plots/fig2_preprocess_data.py when that folder was retired (v2.5.20)."""
import numpy as np

import simplicity.output_manager as om


# The window's LEFT EDGE. "phi" reads it from each simulation; a number is a
# fixed cutoff in days, which is what this was until v2.4.112.
#
# phi is not an approximation of the constant -- it is the quantity the
# constant was averaging. phenotype/update.py:39 scores a lineage as
#
#     infected_fraction * distance_from_weighted_consensus
#   + non_infected_fraction / active_lineages_n
#
# with infected_fraction = min((diagnosed + recovered)/size, 1) -- the SAME
# expression. Below 1 the selection term is damped and every lineage scores on
# 1/active_lineages_n. The 200 came from averaging the day that stops being
# true over 250 simulations of ONE configuration; across the calibrated grid
# (#101-#136) the same day ranges from 158 to 605, so the average describes no
# particular run. At the chosen configuration (N=3000, I0=90, R=1.05) its
# median is 284 -- 84 days of neutral evolution inside a 200-day window.
BURNIN_WINDOW = "phi"

# The fixed-cutoff alternative, and what every figure before v2.4.112 used.
# Kept so both can be produced and compared: the cache key carries the window,
# so they never collide.
BURNIN_CUTOFF_DAYS = 200

# Simulations shorter than MIN_FINAL_TIME are early-extinction runs and are
# dropped entirely (14 of 243 in unbound #1). MIN_WINDOW_DAYS then drops what
# is left with too little record ABOVE its own cutoff -- under "phi" a run
# reaching phi=1 on day 340 and stopping on day 400 has 60 usable days, and
# nothing in its length says so. Applied in BOTH modes, so the only difference
# between them is where the left edge sits.
MIN_FINAL_TIME = 365
MIN_WINDOW_DAYS = 365


def window_tag(window=None):
    """Cache-key fragment naming the window, so two modes never share a file."""
    window = BURNIN_WINDOW if window is None else window
    return "phi" if window == "phi" else f"{float(window):g}"


def phi_day(ssod):
    """The day this simulation completes burn-in: the first time at which
    (diagnosed + recovered)/N reaches 1, read from its own trajectory.

    The same rule as test_config_grid.measure(), which is where the grid's
    "phi=1 day" column comes from -- so the figure's window and the
    configuration grid cannot disagree about what burn-in means.

    N comes from the trajectory's first row, where nobody is diagnosed or
    recovered yet and susceptibles + infected is therefore the whole
    population. A wrong N moves the cutoff with no sign that it has, so it is
    not inferred any other way.

    None if the run never gets there: that run is neutral from start to finish
    and does not belong in a measurement of selection.
    """
    df = om.read_simulation_trajectory(ssod)
    if df.empty:
        return None
    size = float(df["susceptibles"].iloc[0] + df["infected"].iloc[0])
    if size <= 0:
        return None
    phi = ((df["diagnosed"].astype(float) + df["recovered"].astype(float))
           / size).to_numpy()
    hit = np.flatnonzero(phi >= 1.0)
    return float(df["time"].iloc[hit[0]]) if hit.size else None


def window_cutoff(ssod, window=None):
    """One simulation's left edge under the configured window, or None."""
    window = BURNIN_WINDOW if window is None else window
    return phi_day(ssod) if window == "phi" else float(window)
