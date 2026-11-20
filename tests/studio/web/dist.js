/* Per-run values, for the distribution overlay and the sensitivity lines.

   ----------------------------------------------------------------- why

   Every number panels.py sends is already a mean with a bootstrap CI over runs
   (panels.py:146-156). The per-run values exist on the way there -- panel 1
   builds per[cls]["frac_burden"] as a list over runs at panels.py:207 and then
   collapses it -- but nothing in the payload carries them, so no chart can show
   a distribution and a sensitivity line cannot carry a band.

   They CANNOT be added by editing panels.py. The panel memo key is
   [tag, parts, compute.code_hash(), panels.code_hash()] (panels.py:43-55), and
   panels.code_hash() is the sha256 of panels.py's own bytes (:27-32) -- so one
   byte changed there discards 1.93 GB of memoised answers and costs about an
   hour to warm back. tests/studio/dist.py therefore IMPORTS panels and reuses
   its functions without touching the file, memoising under its own tags through
   panels._memo. panels.code_hash() never moves; every existing memo stays good.

   ------------------------------------------------------- the two endpoints

   GET /api/dist?n=<1|2|3>&<the same build and read parameters as /api/panel>

     {kind: "dist", n, params, cached,
      seeds:   {scenario: [seed, ...]},        // the run behind each index
      series:  {"<metric key>": {scenario: [v per run, ...]}},
      ci:      {"<derived key>": {scenario: {n, mean, lo, hi}}}}

   `ci` carries only the DERIVED keys -- the sums and ratios formed inside a run
   (frac_burden_long, frac_burden_majority_long, intensity_ratio), built with
   panels._boot over the same per-run list so it is the identical estimator. The
   keys that mirror a panel field already have their interval in the panel
   payload. A chart must NOT form its own interval off `series`: a spread over
   50 runs shown where a 2000-resample bootstrap CI of the mean is shown
   everywhere else puts two different quantities under one whisker.

   Arrays are index-aligned with `seeds`, and a run with no value for a metric
   holds null at its index rather than being dropped -- so a point can always be
   traced back to its seed, which is the first question asked of an outlier.

   Metric keys are flat dotted paths naming the aggregate they sit under, so a
   tab looks up the key for the CI it is already drawing:

     n=1   frac_burden.<class>   frac_clades.<class>   per_clade.<class>
           -- every share is a FRACTION, 0..1. panels.py formed percentages
           until this pass; the studio converts at the display layer only, so
           a value and its label can no longer disagree about the unit.
           total_burden   total_clades
           frac_burden_long             sum over classes j >= 1
           frac_burden_majority_long    sum over classes j > k/2

     n=2   who_transmits   where.<ll|ls|sl|ss>
           intensity.<long|standard>   intensity_ratio

     n=3   rise_days   latent_days
           by_class_rise.<class>   by_class_latent.<class>

   The two summed panel-1 keys are not decoration: a class CI cannot be added to
   another class CI, so the long-born aggregate has to be formed inside each run
   and bootstrapped from there. Same for intensity_ratio, which is a per-run
   ratio and not one mean divided by another.

   GET /api/sens?axis=<k|mt|win|prom>&<build and read parameters>

     {kind: "sens", axis, pinned: {...}, levels: [...], operating: <level>,
      scenarios: {scenario: {metric: [{n, mean, lo, hi} per level, ...]}}}

   The panel-4 sensitivity line, regrouped and with intervals. panel4's own
   by_axis cannot serve it: P4_AXES is ["k","win","prom"] (panels.py:442), m_T
   is a build parameter that is not in that list, so _p4_by_axis (:503-523)
   leaves it unpinned and by_axis["k"] comes back with 21 levels for three
   values of k -- seven m_T rows per k, ordered lexicographically by directory
   name so m0 lands after m0.1. /api/sens pins every axis but the one asked for.

   `operating` is resolved from compute.DEFAULTS across every build parameter,
   m_T included, and NOT from panels._p4_operating_key (panels.py:458-465) --
   which matches on P4_AXES only, so it ignores m_T and answers
   k7_m0.001_... , m_T = 0.1% rather than the operating point's 1%. That key
   feeds _p4_summary's at_operating_point (panels.py:432-437), so panel 4's
   summary numbers are at the wrong m_T. Same root cause, one layer up.

   ---------------------------------------------------------------- caching

   Both go through apiGetSlot, so a repeat is free and a superseded request is
   aborted. A tab asks for the distribution only when its own box is ticked:
   nothing fetches per-run values to then not draw them.
*/

/* The per-run series for one panel, or null when the server has none -- an
   older store, or a build dist.py has not warmed. A tab must render without it:
   the overlay is an addition to a chart that already stands up alone. */
/* An endpoint this store does not serve is remembered. api.js caches only
   successes (api.js:96), and a tab whose box is ticked re-renders on every
   slider move -- so a 404 would be refetched each time. A supersede is NOT an
   absence: it means a newer render is already on its way. Cleared by a page
   reload, which is when dist.py could have appeared. */
const DIST_ABSENT = new Set();

async function distFor(n) {
  if (DIST_ABSENT.has("dist")) return null;
  try {
    const d = await apiGetSlot(`/api/dist?n=${n}&` + buildQuery());
    return d && d.series ? d : null;
  } catch (e) {
    if (!isSuperseded(e)) DIST_ABSENT.add("dist");
    return null;                 /* no endpoint, no build: draw the CIs alone */
  }
}

/* The runs for one metric key, laid out per category to match a chart's
   `vals`. `cats` are the category values the chart draws; `key` maps a
   category to its metric key, or is a plain string for a single-category
   chart. Returns undefined when there is nothing, which every chart in
   charts.js reads as "no points". */
function distRuns(dist, scenario, cats, key) {
  if (!dist || !dist.series) return undefined;
  const of = typeof key === "function" ? key : (() => key);
  const out = cats.map(c => {
    const s = dist.series[of(c)];
    return (s && s[scenario]) || null;
  });
  return out.some(r => r && r.length) ? out : undefined;
}

/* The same, for a chart whose CATEGORIES are the scenarios -- a forest plot or
   a scenario scatter -- where distRuns cannot help because it holds the
   scenario fixed. `key` is the single metric key for all of them. */
function distRunsByScenario(dist, cats, key) {
  if (!dist || !dist.series) return undefined;
  const s = dist.series[key];
  if (!s) return undefined;
  const out = cats.map(sc => s[sc] || null);
  return out.some(r => r && r.length) ? out : undefined;
}

/* The sensitivity table for one axis of the build grid, or null. */
async function sensFor(axis) {
  if (DIST_ABSENT.has("sens")) return null;
  try {
    const d = await apiGetSlot(`/api/sens?axis=${encodeURIComponent(axis)}&`
                               + buildQuery());
    return d && d.levels ? d : null;
  } catch (e) {
    if (!isSuperseded(e)) DIST_ABSENT.add("sens");
    return null;
  }
}

/* Whether the distribution box may be ticked for a given plot form. A stack of
   means has no per-run form that stays a stack, and a count is not an estimate,
   so those forms disable the box instead of ignoring it. */
const DIST_CAPABLE = {
  stacked: false, denominators: false, days: false, heatmap: false,
};

function distAllowed(plot) {
  return DIST_CAPABLE[plot] !== false;
}
