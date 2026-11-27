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
   hour to warm back. figures_gui/dist.py therefore IMPORTS panels and reuses
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

   ---------------------------------------------------------------- caching

   It goes through apiGetSlot, so a repeat is free and a superseded request is
   aborted. A tab asks for the distribution only when its own box is ticked:
   nothing fetches per-run values to then not draw them.
*/
