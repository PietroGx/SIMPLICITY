/* The painters: the artifact's charts, unchanged in what they draw.

   Each function was ported from tests/artifacts/fig2_studio.html with three
   seams rewired and nothing else touched:

     the canvas   `sizeCanvas` in five copies became draw.js's `fit`
     the colours  every literal became a css() token, so dark mode is free
     the state    `S.rule` / `S.col` / `S.ax` became the `settings` object

   `paintMuller` replaces the artifact's TWO Muller painters (its lines 1834 and
   2789). They shared the band algorithm exactly and differed only in which
   canvas they wrote to, their height, and whether the pre-burn-in stretch is
   washed out -- so those three are arguments now.

   Colours: a clade's colour is either its index on a golden-angle hue wheel, or
   its origin. Which, and how origin is decided, are the `colourBy` and
   `originRule` settings -- the Ratio / Majority / Root choice restored from the
   artifact. The three rules disagree, so the control stays visible.
*/

const MAX_BANDS = 70;        // clades with their own band; the rest collect in grey

function labelOf(inf){
  if (inf.isRoot) return "founder";
  return settings.originRule === "root" ? inf.root : inf.vote;
}

function colourFor(r, inf, ordIdx){
  return settings.colourBy === "origin" ? originColour(inf) : cladeColour(r, ordIdx);
}

function cladeColour(k, i){
  const h = (i * 137.508) % 360, s = 52 + ((i * 37) % 26), l = 48 + ((i * 17) % 20);
  return `hsl(${h.toFixed(1)} ${s}% ${l}%)`;
}

function originColour(inf){
  if (inf.isRoot) return css("--founder");
  if (settings.originRule === "frac"){
    const p = inf.p, a = css("--std"), b = css("--long");
    return `color-mix(in srgb, ${b} ${(p*100).toFixed(0)}%, ${a})`;
  }
  const l = labelOf(inf);
  return l === "long" ? css("--long") : l === "mixed" ? css("--mixed")
       : l === "founder" ? css("--founder") : css("--std");
}

function colourMapFor(T, C){
  const ordered = [...C.alive].sort((a, b) => T.rank[a] - T.rank[b]);
  const colMap = new Map();
  ordered.forEach((r, i) => colMap.set(r, colourFor(r, C.info.get(r), i)));
  return {ordered, colMap};
}

function paintTree(g, w, h, D, T, C, colMap, lw){
  g.clearRect(0, 0, w, h);
  const n = D.n, par = D.parent;
  const xv = settings.treeAxis === "time" ? D.tem : D.depth;
  let xmax = 0;
  for (let i = 0; i < n; i++) if (xv[i] > xmax) xmax = xv[i];
  xmax = xmax || 1;
  const padL = 2, padR = 2, padT = 4, padB = 4;
  const sx = v => padL + (v / xmax) * (w - padL - padR);
  const sy = v => padT + (v / Math.max(1, T.nTips - 1)) * (h - padT - padB);

  // batch by colour: a few thousand strokes at most
  const paths = new Map();
  for (let i = 0; i < n; i++){
    const p = par[i];
    if (p < 0) continue;
    const c = colMap.get(C.clade[i]) || css("--ink-3");
    let pth = paths.get(c);
    if (!pth){ pth = new Path2D(); paths.set(c, pth); }
    const x0 = sx(xv[p]), x1 = sx(xv[i]), y0 = sy(T.y[p]), y1 = sy(T.y[i]);
    pth.moveTo(x0, y0); pth.lineTo(x0, y1); pth.lineTo(x1, y1);
  }
  g.lineWidth = lw != null ? lw
              : T.nTips > 4000 ? 0.5 : T.nTips > 1200 ? 0.7 : 1.1;
  g.lineJoin = "round";
  for (const [c, pth] of paths){ g.strokeStyle = c; g.stroke(pth); }
  return {xmax, sx};
}

function paintOne(cv, r, hCss, big, V){
  const {g, w, h} = fit(cv, hCss);
  g.clearRect(0, 0, w, h);
  const nt = V.nt, j0 = V.j0;
  const padL = 2, padR = 2, padT = 5, padB = 5;
  /* the stacked plot's x, sample for sample, and one ceiling for every cell
     unless the reader asks for each clade's own */
  const top = (settings.yScale === "own" ? Math.max(r.h, settings.burnin ? r.peak : r.peakAll) : V.top) || 1;
  const sx = j => padL + (nt < 2 ? 0 : j / (nt - 1)) * (w - padL - padR);
  const sy = v => h - padB - (v / top) * (h - padT - padB);

  if (settings.showRaw){
    g.strokeStyle = css("--ink-3"); g.globalAlpha = .35; g.lineWidth = 1;
    g.beginPath();
    for (let j = 0; j < nt; j++) j ? g.lineTo(sx(j), sy(r.yAll[j])) : g.moveTo(sx(j), sy(r.yAll[j]));
    g.stroke(); g.globalAlpha = 1;
  }
  const col = r.col || css("--ink-3");
  g.strokeStyle = col; g.lineWidth = big ? 1.6 : 1.3;
  g.beginPath();
  for (let j = 0; j < r.ys.length; j++)
    j ? g.lineTo(sx(j0 + j), sy(r.ys[j])) : g.moveTo(sx(j0 + j), sy(r.ys[j]));
  g.stroke();

  const rise = css("--rise");
  for (const e of r.eps){
    // the measured stretch: lo-crossing to hi-crossing of THIS episode
    g.strokeStyle = rise; g.lineWidth = big ? 4 : 3; g.lineCap = "round";
    g.globalAlpha = (r.main && e === r.main) ? 1 : .45;
    g.beginPath();
    for (let j = e.ia; j <= e.ib; j++)
      j > e.ia ? g.lineTo(sx(j0 + j), sy(r.ys[j])) : g.moveTo(sx(j0 + j), sy(r.ys[j]));
    if (e.ib === e.ia){ g.lineTo(sx(j0 + e.ia) + .01, sy(r.ys[e.ia])); }
    g.stroke(); g.lineCap = "butt";
    // the episode's own ends: trough it rose from, peak it reached
    g.fillStyle = rise;
    g.beginPath(); g.arc(sx(j0 + e.a), sy(r.ys[e.a]), big ? 3 : 2.1, 0, 7); g.fill();
    g.beginPath(); g.arc(sx(j0 + e.b), sy(r.ys[e.b]), big ? 3.4 : 2.4, 0, 7); g.fill();
    g.globalAlpha = 1;
  }
  /* what the metric never looked at, washed out exactly as in the plot above */
  if (settings.burnin && V.phiIdx != null){
    g.save(); g.globalAlpha = .62; g.fillStyle = css("--bg");
    g.fillRect(0, 0, sx(V.phiIdx), h); g.restore();
    drawBurnin(g, sx(V.phiIdx), w, h, V.phi, true);
  }
  return {sx, sy};
}

/* Clade frequency as stacked bands over the run. One painter for both the
   clustering tab and the sweep-time tab.

     washBurnin  shade the stretch the metric never looks at, rather than crop
                 it: the shape to the left is the reason the window is there

   Returns what the caller needs for its own subtitle, so this function touches
   no DOM but the canvas it was given. */
/* `excluded` is a Set of clade ids that m_F has dropped: they keep their band,
   because the clade is still there, and lose their colour.

   `selected` is one clade id, which this draws at full strength over a washed
   plot. The band geometry it needs for that is the geometry the stack already
   computes, so it is kept rather than thrown away, and `hit()` in the returned
   object answers which clade is under a point. */
function paintMuller(cv, {D, C, ordered, colMap, height, washBurnin, excluded,
                          selected}) {
  const {g, w, h} = fit(cv, height);
  const times = D.times, nt = times.length;
  const padL = 2, padR = 2, padT = 4, padB = 4;
  const sx = j => padL + (nt < 2 ? 0 : j / (nt - 1)) * (w - padL - padR);

  const F = cladeSeries(D, C, ordered);
  const nc = ordered.length;
  const peak = ordered.map((_, i) => {
    let m = 0;
    for (let j = 0; j < nt; j++) if (F[i * nt + j] > m) m = F[i * nt + j];
    return m;
  });
  const rank = ordered.map((_, i) => i).sort((a, b) => peak[b] - peak[a]);
  const own = new Set(rank.slice(0, MAX_BANDS));

  const base = new Float32Array(nt);
  const grey = css("--ink-3");
  const sy = v => padT + v * (h - padT - padB);
  /* Every band that was drawn, bottom and top per sample. The overflow band
     carries id null: it is many clades at once, so it can be hit but not
     named. */
  const bands = [];
  const bandPath = b => {
    g.beginPath();
    g.moveTo(sx(0), sy(b.lo[0]));
    for (let j = 1; j < nt; j++) g.lineTo(sx(j), sy(b.lo[j]));
    for (let j = nt - 1; j >= 0; j--) g.lineTo(sx(j), sy(b.hi[j]));
    g.closePath();
  };
  const drawBand = (vals, colour, id) => {
    const lo = Float32Array.from(base);
    for (let j = 0; j < nt; j++) base[j] = Math.min(1, base[j] + vals[j]);
    const b = {id: id === undefined ? null : id, colour,
               lo, hi: Float32Array.from(base)};
    bands.push(b);
    bandPath(b); g.fillStyle = colour; g.fill();
  };

  const other = new Float32Array(nt);
  for (let i = 0; i < nc; i++) {
    const vals = F.subarray(i * nt, i * nt + nt);
    if (own.has(i)) drawBand(vals, (excluded && excluded.has(ordered[i]))
                                   ? grey : (colMap.get(ordered[i]) || grey),
                             ordered[i]);
    else for (let j = 0; j < nt; j++) other[j] += vals[j];
  }
  let anyOther = false;
  for (let j = 0; j < nt; j++) if (other[j] > 1e-6) { anyOther = true; break; }
  if (anyOther) drawBand(other, grey);

  /* The pop-out: wash the whole plot back, then redraw the one band over it.
     Washing rather than recolouring every other band keeps their boundaries
     visible, so the selected clade is read in its place in the stack. */
  const sel = selected === null || selected === undefined ? null
            : bands.find(b => b.id === selected);
  if (sel) {
    g.save(); g.globalAlpha = 0.8; g.fillStyle = css("--bg");
    g.fillRect(0, 0, w, h); g.restore();
    bandPath(sel); g.fillStyle = sel.colour; g.fill();
    g.strokeStyle = css("--ink"); g.lineWidth = 1; g.stroke();
  }

  const day = burninDayOf(D);
  if (day !== null) {
    const ji = dayToIndex(times, day);
    if (ji !== null) {
      if (washBurnin) {
        g.save(); g.globalAlpha = 0.62; g.fillStyle = css("--bg");
        g.fillRect(0, 0, sx(ji), h); g.restore();
      }
      drawBurnin(g, sx(ji), w, h, day, false);
    }
  }
  /* Which clade's band holds a point, in CSS pixels -- which is what fit()
     leaves the context in, so a getBoundingClientRect offset drops straight
     in. null above a clipped stack, below it, or on the overflow band. */
  const hit = (px, py) => {
    if (nt < 2) return null;
    const j = Math.max(0, Math.min(nt - 1, Math.round(
      ((px - padL) / Math.max(1e-9, w - padL - padR)) * (nt - 1))));
    const v = (py - padT) / Math.max(1e-9, h - padT - padB);
    for (const b of bands)
      if (b.id !== null && v >= b.lo[j] && v < b.hi[j]) return b.id;
    return null;
  };

  return {nClades: nc, shown: Math.min(nc, MAX_BANDS), anyOther, bands, hit,
          selected: sel ? sel.id : null, sx, sy,
          day0: times[0], day1: times[nt - 1]};
}
