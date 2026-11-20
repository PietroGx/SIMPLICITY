/* Episode detection in the browser.

   This is the one place in the studio that deliberately duplicates logic
   already written in python: `smooth`, `episodesOf` and `riseOf` here mirror
   `smooth`, `episodes_of` and `rise_of` in tests/studio/compute.py.

   The duplication is the price of the sliders. The clade build holds episodes
   at ONE set of episode rules, so moving the smoothing window or the
   prominence bar would otherwise need a twelve-minute rebuild before anything
   could be seen. Computing them here from the run's own trajectories makes the
   rules immediate -- which is what the tab is for.

   The risk is the obvious one: two implementations of one definition can
   drift, and a drift would show the panels one number and this tab another.
   tests/studio/check_episodes.js pins them together by running both over the
   same run and comparing episode counts and rise days.

   The functions below were taken from tests/artifacts/fig2_studio.html
   unchanged; `findPeaks` is its own reimplementation of
   scipy.signal.find_peaks's prominence test, which is what compute.py calls.
*/

function smooth(y, times, winDays){
  const n = y.length;
  if (n < 2) return Float64Array.from(y);
  const step = (times[n-1] - times[0]) / (n - 1);
  let w = Math.max(1, Math.round(winDays / Math.max(step, 1e-9)));
  if (w % 2 === 0) w += 1;
  const off = (w - 1) >> 1;
  const pre = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) pre[i+1] = pre[i] + y[i];
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++){
    const a = Math.max(0, i - off), b = Math.min(n - 1, i + off);
    out[i] = (pre[b+1] - pre[a]) / (b - a + 1);
  }
  return out;
}

function findPeaks(ys, minProm){
  const n = ys.length, out = [];
  for (let i = 1; i < n - 1; i++){
    if (!(ys[i] > ys[i-1])) continue;
    let j = i;
    while (j + 1 < n && ys[j+1] === ys[i]) j++;
    if (j + 1 >= n || ys[j+1] > ys[i]){ i = j; continue; }
    const mid = (i + j) >> 1;
    let lmin = ys[mid];
    for (let k = mid; k >= 0 && ys[k] <= ys[mid]; k--) if (ys[k] < lmin) lmin = ys[k];
    let rmin = ys[mid];
    for (let k = mid; k < n && ys[k] <= ys[mid]; k++) if (ys[k] < rmin) rmin = ys[k];
    if (ys[mid] - Math.max(lmin, rmin) >= minProm) out.push(mid);
    i = j;
  }
  return out;
}

/* The censored peak, mirroring clade_preprocess.episodes_of: the clade's
   maximum, when it rose to it by the prominence bar and never falls back by
   that bar before the record ends. Placed at the first sample within
   CENSOR_TOL of the maximum, so ~1e-16 noise on a plateau at 1.0 cannot move
   it. -1 when there is none, or when the maximum is the last sample. */
const TROUGH_PROM = 1e-4, CENSOR_TOL = 1e-9;
function censoredPeak(ys, h, minProm){
  const n = ys.length;
  let m = 0;
  while (m < n && !(ys[m] >= h * (1 - CENSOR_TOL))) m++;
  if (m >= n - 1) return -1;
  let lmin = Infinity, rmin = Infinity;
  for (let k = 0; k <= m; k++) if (ys[k] < lmin) lmin = ys[k];
  for (let k = m; k < n; k++) if (ys[k] < rmin) rmin = ys[k];
  return h - lmin >= minProm && h - rmin < minProm ? m : -1;
}

function episodesOf(ys, prom){
  let h = 0;
  for (let i = 0; i < ys.length; i++) if (ys[i] > h) h = ys[i];
  if (!(h > 0)) return {h:0, eps:[]};
  const pk = findPeaks(ys, prom * h).map(i => [i, false]);
  const cm = censoredPeak(ys, h, prom * h);
  if (cm >= 0) pk.push([cm, true]);
  pk.sort((a, b) => a[0] - b[0]);
  if (!pk.length) return {h, eps:[]};
  const neg = new Float64Array(ys.length);
  for (let i = 0; i < ys.length; i++) neg[i] = -ys[i];
  const mins = [0, ...findPeaks(neg, TROUGH_PROM * h)];
  const eps = pk.map(([i, cens]) => {
    let t = 0;
    for (const m of mins){ if (m < i) t = m; else break; }
    return {a:t, b:i, cens};
  });
  return {h, eps};
}

function riseOf(ys, times, e, lo, hi){
  const ylo = ys[e.a], yhi = ys[e.b];
  if (!(yhi > ylo)) return null;
  const A = ylo + lo * (yhi - ylo), B = ylo + hi * (yhi - ylo);
  let ia = e.a, ib = e.b;
  for (let i = e.a; i <= e.b; i++) if (ys[i] >= A){ ia = i; break; }
  for (let i = e.a; i <= e.b; i++) if (ys[i] >= B){ ib = i; break; }
  return {...e, ia, ib, days: times[ib] - times[ia], amp: yhi - ylo, top: yhi};
}

function nearIdx(v, list){
  let b = 0;
  for (let i = 1; i < list.length; i++) if (Math.abs(list[i] - v) < Math.abs(list[b] - v)) b = i;
  return b;
}

/* The clade's chosen episode. The artifact made this a four-way control
   (largest / weighted / first / fastest); "first" was definition 2, dropped at
   tests/artifacts/STUDIO_PLAN.md:211 for a measured confound with lifetime of
   +0.41 to +0.69 against +0.14 to +0.36 for this one. The decision stands, so
   this is a constant rather than a knob. */
const EPISODE_PICK = "largest";   // superseded below: see the pick in
                                 // prepareEpisodes, which mirrors compute.py

/* Per-clade trajectories, smoothed, with every episode measured.

   Measured on the window, drawn over the whole run: the cells share the
   stacked plot's x axis, so a peak in a cell is the same day as a peak above.
   Ported from the artifact's `compute()`; the only changes are that the
   settings come from the `settings` object and the pick rule is the constant
   above.
*/
function prepareEpisodes(view, daily) {
  const {D, T, C} = view;
  const {ordered, colMap} = colourMapFor(T, C);
  /* fig2's pool: non-root clades carrying at least one defining mutation.

     Whose count of defining mutations, though. The browser derives `ndef` from
     the export's per-lineage mutation counts; the build derives it from the
     clade metadata, over the segment between the parent clade's root and the
     clade's own root. Under rule B those disagree, because a refused lineage
     folds into its ancestor and brings its mutations with it -- on control
     seed 0 at m_T=0.012% the build's pool is 32 clades and the browser's own
     filter left 9. So when a build is in hand, the pool is the BUILD's. */
  const dailyMap = daily ? new Map(daily.clades.map(c => [c.root, c])) : null;
  const ids = [...C.alive].filter(r => {
    if (dailyMap) {
      const rec = dailyMap.get(r);
      return !!rec && !rec.founding && rec.ndef > 0;
    }
    const i = C.info.get(r);
    return !i.isRoot && i.ndef > 0;
  });
  const nt = D.times.length;
  const F = cladeSeries(D, C, ids);
  const phi = settings.burnin ? burninDayOf(D) : null;

  /* The measurement axis. runs/ samples 120 time points about 9.2 days apart,
     so a 30-day rise is three samples and every crossing day is quantised;
     `daily` is the same record at 1-day spacing, served by /api/run_daily,
     which is the axis compute.py measures on. A clade is matched by the index
     of its root lineage, which both sides take from ph row order.

     The Muller plot above is still drawn from runs/ -- that is a shape, not a
     number, and export_clade_explorer.py says its frequencies "feed a stacked
     area plot, nothing numerical". Both span the same days, so the cells and
     the stacked plot still line up in day space. */
  const dailyOf = dailyMap;
  const ft = dailyOf ? daily.t : Array.from(D.times);
  const ntm = ft.length;
  let j0 = 0;
  if (dailyOf) j0 = daily.j0 || 0;      // the server's own slice index
  else if (phi !== null) while (j0 < ntm - 1 && ft[j0] < phi) j0++;
  const times = ft.slice(j0);
  const rows = [];
  let unserved = 0;
  for (let i = 0; i < ids.length; i++) {
    const rec = dailyOf ? dailyOf.get(ids[i]) : null;
    if (dailyOf && !rec) { unserved++; continue; }
    const yAll = rec ? rec.y : Array.from(F.subarray(i * nt, i * nt + nt));
    const y = yAll.slice(j0);
    let peak = 0;
    for (const v of y) if (v > peak) peak = v;
    if (!(peak > 0)) continue;
    /* The smoothed line is drawn here; the EPISODES are measured by
       compute.episodes_of and compute.rise_of on the server, because this file
       used to reimplement scipy's prominence test and chose a different trough
       -- on HIV_high seed 0 Clade_0.366.5.11, index 496 against compute.py's
       600, a 71-day difference in the same clade's rise. One definition now
       has one implementation; `episodesOf` below stays only for the fallback
       axis, where no build is available to measure on. */
    const ys = smooth(y, times, settings.win);
    let h, rs, nFound;
    if (rec) {
      h = rec.h;
      rs = rec.eps || [];
      nFound = rec.n_found ?? rs.length;
    } else {
      const got = episodesOf(ys, settings.prom);
      h = got.h;
      nFound = got.eps.length;
      rs = got.eps.map(e => riseOf(ys, times, e, settings.lo, settings.hi))
             .filter(Boolean);
    }
    let main = null, value = null;
    if (rs.length) {
      main = rec ? (rs[rec.main] || rs.find(e => e.emg) || rs[0])
                 : rs.reduce((a, b) => b.amp > a.amp ? b : a);
      value = main.days;
    }
    const ampShare = rs.length ? main.amp / rs.reduce((t, e) => t + e.amp, 0) : 0;
    let peakAll = 0;
    for (const v of yAll) if (v > peakAll) peakAll = v;
    rows.push({id: ids[i], inf: C.info.get(ids[i]), col: colMap.get(ids[i]),
               yAll, y, ys, times, h, peak, peakAll, eps: rs, nFound, main, value, ampShare});
  }
  rows.sort((a, b) => b.peak - a.peak);
  let top = 0;
  for (const r of rows) top = Math.max(top, r.h, settings.burnin ? r.peak : r.peakAll);
  return {rows, ids, ordered, colMap, nt: ntm, phi,
          daily: !!dailyOf, unserved,
          view: {ft, nt: ntm, j0, top: top || 1,
                 phiIdx: phi !== null ? dayToIndex(ft, phi) : null, phi}};
}
