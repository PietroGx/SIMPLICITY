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

function colourFor(r, inf, ordIdx, n){
  return settings.colourBy === "origin" ? originColour(inf) : cladeColour(r, ordIdx, n);
}

/* settings.cladePal: [value, label, cycle or ramp]. A cycle repeats, lighter then
   darker on each pass; a ramp runs along the tree order, top to bottom. */
const CLADE_PALETTES = [
  /* colour-blind safe only: the hue wheel, Tableau 10 and Spectral are gone (v2.5.21+) */
  ["okabe", "Okabe–Ito", {cycle: ["#e69f00", "#56b4e9", "#009e73",
    "#f0e442", "#0072b2", "#d55e00", "#cc79a7", "#999999"]}],
  ["viridis", "Viridis, by tree order", {ramp: ["#440154", "#482878", "#3e4989", "#31688e",
    "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"]}],
];
/* a stored palette that is no longer offered falls back to Okabe–Ito */
if (!CLADE_PALETTES.some(([v]) => v === settings.cladePal)) settings.cladePal = "okabe";

function cladeColour(k, i, n){
  const pal = (CLADE_PALETTES.find(([v]) => v === settings.cladePal) || [])[2];
  if (pal && pal.cycle){
    const c = pal.cycle, base = c[i % c.length], pass = Math.floor(i / c.length) % 3;
    return pass === 0 ? base : mixHex(base, pass === 1 ? "#ffffff" : "#000000", pass === 1 ? 0.4 : 0.3);
  }
  if (pal && pal.ramp){
    const r = pal.ramp, t = (n > 1 ? i / (n - 1) : 0.5) * (r.length - 1), j = Math.min(r.length - 2, Math.floor(t));
    return mixHex(r[j], r[j + 1], t - j);
  }
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
  ordered.forEach((r, i) => colMap.set(r, colourFor(r, C.info.get(r), i, ordered.length)));
  return {ordered, colMap};
}

/* `selected`, a clade id: its lineages drawn thicker on top, the rest washed */
function paintTree(g, w, h, D, T, C, colMap, lw, selected){
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
  const sel = selected === null || selected === undefined || !colMap.has(selected) ? null : selected;
  if (sel !== null) g.globalAlpha = 0.14;
  for (const [c, pth] of paths){ g.strokeStyle = c; g.stroke(pth); }
  g.globalAlpha = 1;
  let nSel = 0;
  if (sel !== null){
    const pth = new Path2D();
    for (let i = 0; i < n; i++){
      const p = par[i];
      if (p < 0 || C.clade[i] !== sel) continue;
      nSel++;
      const x0 = sx(xv[p]), x1 = sx(xv[i]), y0 = sy(T.y[p]), y1 = sy(T.y[i]);
      pth.moveTo(x0, y0); pth.lineTo(x0, y1); pth.lineTo(x1, y1);
    }
    g.lineWidth = Math.max(1.6, 2.4 * g.lineWidth);
    g.strokeStyle = colMap.get(sel); g.stroke(pth);
  }
  return {xmax, sx, nSel};
}

/* What each rail card moves, by the layers the painters draw: a focused card
   draws its layers emphasised and fades the rest. */
const SW_FOCUS = {
  curve: ["curve", "raw"], peak: ["eps", "peaks"], start: ["eps", "starts", "launch"],
  breakout: ["beps", "bpeak", "F"], hidden: ["band", "anchor", "launch"],
  cryptic: ["band", "launch", "ceil"],
};
const focusOn = (focus, k) => !focus || (SW_FOCUS[focus] || []).includes(k);
const FADED = .15;

function paintOne(cv, r, hCss, big, V, showEps = true, focus = null){
  const {g, w, h} = fit(cv, hCss);
  g.clearRect(0, 0, w, h);
  const nt = V.nt, j0 = V.j0;
  const padL = 2, padR = 2, padT = 5, padB = 5;
  /* the stacked plot's x, sample for sample, and one ceiling for every cell
     unless the reader asks for each clade's own */
  const top = (settings.yScale === "own" ? Math.max(r.h, settings.burnin ? r.peak : r.peakAll) : V.top) || 1;
  const x0 = V.x0 || 0;                     // the first sample on screen
  const sx = j => padL + (nt - x0 < 2 ? 0 : (j - x0) / (nt - 1 - x0)) * (w - padL - padR);
  const sy = v => h - padB - (v / top) * (h - padT - padB);

  const jE = r.times ? emergedFrom(r.times, r.tEm) : 0;   // into ys; j0 + jE into yAll
  const on = k => focusOn(focus, k);
  /* the raw daily record is an area, the kernel's output a line over it */
  if (settings.showRaw || focus === "curve"){
    const yb = sy(0);
    g.fillStyle = css("--ink-3"); g.globalAlpha = on("raw") ? (focus ? .32 : .2) : .06;
    g.beginPath(); g.moveTo(sx(j0 + jE), yb);
    for (let j = j0 + jE; j < nt; j++) g.lineTo(sx(j), sy(r.yAll[j]));
    g.lineTo(sx(nt - 1), yb); g.closePath(); g.fill(); g.globalAlpha = 1;
  }
  const col = r.col || css("--ink-3");
  /* the outbreak's own episode, matched by its peak day */
  const isBrk = e => !!(r.brk && isNum(r.brk.peakDay) && r.times
    && Math.abs(r.times[e.b] - r.brk.peakDay) < 1e-6);
  const epsOn = showEps || focus === "peak" || focus === "start";
  const eps = epsOn ? r.eps : focus === "breakout" ? r.eps.filter(isBrk) : [];
  /* with episodes on, the curve recedes and each episode, start to peak, is
     drawn over it in the same colour, bold */
  g.strokeStyle = col;
  g.lineWidth = focus === "curve" ? (big ? 2.4 : 2) : eps.length ? 1 : big ? 1.6 : 1.3;
  g.globalAlpha = focus ? (on("curve") ? 1 : FADED) : eps.length ? .35 : 1;
  g.beginPath();
  for (let j = jE; j < r.ys.length; j++)
    j > jE ? g.lineTo(sx(j0 + j), sy(r.ys[j])) : g.moveTo(sx(j0 + j), sy(r.ys[j]));
  g.stroke(); g.globalAlpha = 1;

  if (focus !== "curve") {
    g.lineWidth = big ? 3.2 : 2.6; g.lineCap = "round"; g.lineJoin = "round";
    for (const e of eps){
      g.globalAlpha = on("eps") || (isBrk(e) && on("beps")) ? 1 : FADED;
      g.beginPath();
      for (let j = e.a; j <= e.b; j++)
        j > e.a ? g.lineTo(sx(j0 + j), sy(r.ys[j])) : g.moveTo(sx(j0 + j), sy(r.ys[j]));
      if (e.b === e.a) g.lineTo(sx(j0 + e.a) + .01, sy(r.ys[e.a]));
      g.stroke();
    }
    g.lineCap = "butt"; g.globalAlpha = 1;
  }
  /* a focused card's points: every episode's peak, or every episode's start */
  const pts = (key, idx, fill) => { if (!focus || !(SW_FOCUS[focus] || []).includes(key)) return;
    for (const e of r.eps) { g.fillStyle = fill(); g.beginPath();
      g.arc(sx(j0 + e[idx]), sy(r.ys[e[idx]]), big ? 3.6 : 3, 0, 6.2832); g.fill();
      g.strokeStyle = css("--panel"); g.lineWidth = 1; g.stroke(); } };
  pts("peaks", "b", () => col);
  pts("starts", "a", () => css("--long"));
  /* Episodes card: each peak's width at half its prominence, the segment Min
     width is tested on; a peak it dropped keeps a dashed bar and a cross */
  if (focus === "peak") for (const [p, li, ri, yh, days, kept] of r.wid || []) {
    const xl = sx(j0 + li), xr = sx(j0 + ri), y = sy(yh), ink = css(kept ? "--ink" : "--ink-3");
    g.strokeStyle = ink; g.lineWidth = kept ? 1.6 : 1.2; g.setLineDash(kept ? [] : [3, 2]);
    g.beginPath(); g.moveTo(xl, y); g.lineTo(xr, y); g.stroke(); g.setLineDash([]);
    g.beginPath(); g.moveTo(xl, y - 3); g.lineTo(xl, y + 3); g.moveTo(xr, y - 3); g.lineTo(xr, y + 3); g.stroke();
    if (!kept) { const xp = sx(j0 + p), yp = sy(r.ys[p]), s = 3;
      g.lineWidth = 1.6; g.beginPath();
      g.moveTo(xp - s, yp - s); g.lineTo(xp + s, yp + s); g.moveTo(xp - s, yp + s); g.lineTo(xp + s, yp - s);
      g.stroke(); }
    const lab = `${Math.round(days)} d`;
    g.font = "9px " + css("--f-num");
    const tw = g.measureText(lab).width;
    if (xr - xl > tw + 4) { g.fillStyle = ink; g.textAlign = "center";
      g.fillText(lab, (xl + xr) / 2, Math.max(9, y - 3)); g.textAlign = "left"; }
  }
  return {sx, sy, w, h};
}

/* An outbreak's marks over a clade curve already drawn on X (day) and Y
   (frequency) inside box {l, r, t, b}: hidden-time band, first appearance,
   sweep start, launch and peak dots, F and f rules. One painter for the Sweep
   time cells and Figure 5's popup, so the two cannot drift. */
function drawBreakout(g, X, Y, box, b, {labels = false, r = 3.5, hidden = true, startLine = true,
                                        focus = null, showF = true} = {}) {
  /* a focused layer is drawn stronger, an unfocused one faded */
  const lay = k => { const on = focusOn(focus, k);
    g.globalAlpha = on ? 1 : FADED; return focus && on; };
  const rule = (v, dash, n, col) => { const y = Y(v); if (!(y >= box.t && y <= box.b)) return;
    const hi = lay(n);
    g.strokeStyle = css(col || (hi ? "--ink" : "--ink-3")); g.lineWidth = hi ? 1.8 : col ? 1.3 : 1; g.setLineDash(dash);
    g.beginPath(); g.moveTo(box.l, y); g.lineTo(box.r, y); g.stroke(); g.setLineDash([]);
    if (labels) { g.fillStyle = css(col || "--ink-3"); g.font = "9px " + css("--f-num"); g.textAlign = "right";
      g.fillText(`${n} ${brkPct(v)}`, box.r - 2, y - 3); g.textAlign = "left"; } };
  if (showF) rule(+settings.brkF, [4, 3], "F");
  g.globalAlpha = 1;
  if (!b) return;
  const xa = X(b.tem), xs = X(b.start), xh = X(isNum(b.hEnd) ? b.hEnd : b.start);
  if (hidden) { const hi = lay("band");
    g.fillStyle = css("--long"); g.globalAlpha = hi ? .3 : focus ? .05 : .14;
    g.fillRect(xa, box.t, Math.max(1, xh - xa), box.b - box.t); g.globalAlpha = 1; }
  const vline = (x, dash, k) => { const hi = lay(k);
    g.strokeStyle = css(hi ? "--ink" : "--ink-2"); g.lineWidth = hi ? 1.8 : 1; g.setLineDash(dash);
    g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.b); g.stroke(); g.setLineDash([]); };
  if (hidden) vline(xa, [], "anchor");
  /* the ceiling across the hidden time: f, which the curve stays under while hidden */
  const ceil = b.ceil !== undefined ? b.ceil : +settings.brkLo;
  if (hidden && isNum(ceil) && xh > xa) { const hi = lay("ceil"), y = Y(ceil);
    g.strokeStyle = css("--long"); g.lineWidth = hi ? 2 : 1.3; g.setLineDash([2, 3]);
    g.beginPath(); g.moveTo(xa, y); g.lineTo(xh, y); g.stroke(); g.setLineDash([]); }
  if (startLine) vline(xs, [3, 2], "launch");
  const dot = (x, y, col, k) => { const hi = lay(k), rr = hi ? r * 1.5 : r;
    g.fillStyle = col; g.beginPath(); g.arc(x, y, rr, 0, 6.2832); g.fill();
    g.strokeStyle = css(hi ? "--ink" : "--panel"); g.lineWidth = 1; g.stroke(); };
  dot(xs, Y(b.launch), css("--long"), "launch");
  if (isNum(b.peakDay)) dot(X(b.peakDay), Y(b.peak), css("--ink"), "bpeak");
  g.globalAlpha = 1;
}

/* drawBreakout's legend: one small SVG per mark, drawn as the canvas draws it */
function breakoutKey(curve = "clade frequency (smoothed)", anchor = "first appearance", tx = false) {
  const ink = css("--ink"), ink2 = css("--ink-2"), ink3 = css("--ink-3"), lng = css("--long");
  const sv = body => `<svg width="22" height="12" viewBox="0 0 22 12">${body}</svg>`;
  const hl = dash => sv(`<line x1="1" y1="6" x2="21" y2="6" stroke="${ink3}" stroke-width="1.2" stroke-dasharray="${dash}"/>`);
  const vl = dash => sv(`<line x1="11" y1="0" x2="11" y2="12" stroke="${ink2}" stroke-width="1.2" stroke-dasharray="${dash}"/>`);
  const dot = col => sv(`<circle cx="11" cy="6" r="3.5" fill="${col}" stroke="${css("--panel")}"/>`);
  const items = [
    ...(curve ? [[sv(`<line x1="1" y1="6" x2="21" y2="6" stroke="${ink}" stroke-width="1.6"/>`), curve]] : []),
    [sv(`<rect x="1" y="0" width="20" height="12" fill="${lng}" opacity=".18"/>`), "hidden time"],
    [vl("none"), anchor],
    ...(tx ? [[sv(`<line x1="11" y1="0" x2="11" y2="12" stroke="${ink2}" stroke-width="1.2" stroke-dasharray="1 2"/>`),
     "first transmission: intra-host hidden time before it, between-host hidden time after"],
    [sv(`<line x1="11" y1="0" x2="11" y2="12" stroke="${ink2}" stroke-width="1" stroke-dasharray="6 2 1 2"/>`),
     "first day above f: between-host hidden time ends"]] : []),
    [vl("3 2"), "outbreak sweep starts"],
    [dot(lng), "launch frequency"],
    [dot(ink), "outbreak peak"],
    [hl("4 3"), "F: big-sweep threshold"],
    [sv(`<line x1="1" y1="6" x2="21" y2="6" stroke="${lng}" stroke-width="1.4" stroke-dasharray="2 3"/>`),
     cryCeilText()],
  ];
  return `<div class="bkey">${items.map(([i, t]) => `<span>${i}${t}</span>`).join("")}</div>`;
}

/* Clade frequency as stacked bands over the run. One painter for both the
   clustering tab and the sweep-time tab.

     washBurnin  shade the stretch the metric never looks at, rather than crop
                 it: the shape to the left is the reason the window is there
     burninLine  draw the burn-in day's dashed line and label

   Returns what the caller needs for its own subtitle, so this function touches
   no DOM but the canvas it was given. */
/* `excluded` is a Set of clade ids that m_F has dropped: they keep their band,
   because the clade is still there, and lose their colour.

   `selected` is one clade id, which this draws at full strength over a washed
   plot. The band geometry it needs for that is the geometry the stack already
   computes, so it is kept rather than thrown away, and `hit()` in the returned
   object answers which clade is under a point. */
/* `hideExcluded`: excluded clades and the overflow are left out of the stack
   entirely instead of drawn grey, so it holds only the clades that are kept. */
function paintMuller(cv, {D, C, ordered, colMap, height, washBurnin, excluded,
                          selected, hideExcluded, burninLine = true, burninLabel = true}) {
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
  const out = i => !!(hideExcluded && excluded && excluded.has(ordered[i]));
  const rank = ordered.map((_, i) => i).filter(i => !out(i)).sort((a, b) => peak[b] - peak[a]);
  const own = new Set(rank.slice(0, MAX_BANDS));

  const base = new Float32Array(nt);
  const grey = css("--ink-3");
  /* a partial stack grows up from the floor, so its gap is at the top */
  const sy = v => hideExcluded ? h - padB - v * (h - padT - padB) : padT + v * (h - padT - padB);
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
    if (out(i)) continue;
    const vals = F.subarray(i * nt, i * nt + nt);
    if (own.has(i)) drawBand(vals, (excluded && excluded.has(ordered[i]))
                                   ? grey : (colMap.get(ordered[i]) || grey),
                             ordered[i]);
    else for (let j = 0; j < nt; j++) other[j] += vals[j];
  }
  let anyOther = false;
  for (let j = 0; j < nt; j++) if (other[j] > 1e-6) { anyOther = true; break; }
  if (hideExcluded) anyOther = false;
  if (anyOther) drawBand(other, grey);

  /* The pop-out: wash the whole plot back, then redraw the one band over it.
     Washing rather than recolouring every other band keeps their boundaries
     visible, so the selected clade is read in its place in the stack. */
  /* a selected clade the stack left out (excluded, or past MAX_BANDS) is drawn
     from the floor on its own, so a click always shows it */
  const lone = id => { const i = ordered.indexOf(id); if (i < 0) return null;
    return {id, colour: colMap.get(id) || grey, lo: new Float32Array(nt),
            hi: Float32Array.from(F.subarray(i * nt, i * nt + nt))}; };
  const sel = selected === null || selected === undefined ? null
            : bands.find(b => b.id === selected) || lone(selected);
  /* on a floor-up stack the selected clade is lifted down to the floor, so its
     outline is its own frequency rather than its place between the others */
  const pop = sel && hideExcluded ? {lo: new Float32Array(nt),
    hi: sel.hi.map((v, j) => v - sel.lo[j])} : null;
  if (sel) {
    g.save(); g.globalAlpha = 0.8; g.fillStyle = css("--bg");
    g.fillRect(0, 0, w, h); g.restore();
    bandPath(pop || sel); g.fillStyle = sel.colour; g.fill();
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
      if (burninLine) drawBurnin(g, sx(ji), w, h, burninLabel ? day : null, false);
    }
  }
  /* Which clade's band holds a point, in CSS pixels -- which is what fit()
     leaves the context in, so a getBoundingClientRect offset drops straight
     in. null above a clipped stack, below it, or on the overflow band. */
  const hit = (px, py) => {
    if (nt < 2) return null;
    const t = Math.max(0, Math.min(nt - 1,
      ((px - padL) / Math.max(1e-9, w - padL - padR)) * (nt - 1)));
    const j = Math.min(nt - 2, Math.floor(t)), f = t - j;
    const at = a => a[j] + (a[j + 1] - a[j]) * f;
    const v = (hideExcluded ? h - padB - py : py - padT) / Math.max(1e-9, h - padT - padB);
    if (pop && v >= 0 && v < at(pop.hi)) return sel.id;
    for (const b of bands)
      if (b.id !== null && v >= at(b.lo) && v < at(b.hi)) return b.id;
    return null;
  };
  /* a click, from the canvas as laid out NOW: the page can narrow after the
     paint (a scrollbar arriving with the overview grid) and CSS stretches the
     drawing to fit, so the offset is scaled back to the drawn width */
  const hitEvent = ev => {
    const r = cv.getBoundingClientRect();
    return hit((ev.clientX - r.left) * w / Math.max(1, r.width),
               (ev.clientY - r.top) * h / Math.max(1, r.height));
  };

  return {nClades: nc, shown: Math.min(nc, MAX_BANDS), drawn: bands.filter(b => b.id !== null).length, anyOther, bands, hit, hitEvent,
          selected: sel ? sel.id : null, sx, sy,
          day0: times[0], day1: times[nt - 1]};
}
