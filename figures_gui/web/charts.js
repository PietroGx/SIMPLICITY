/* The chart forms a panel tab chooses between.

   draw.js keeps the primitives every tab has always used -- grouped bars, day
   lines, the burn-in rule. These are the forms the panels now PICK between, and
   they all take the same series shape, so a tab maps its payload once:

     {colour, label, vals: [CI|null, ...], runs: [[v, ...]|null, ...]}

   `vals[i]` is the bootstrap CI for category i, exactly as panels.py sends it:
   {n, mean, lo, hi}, with lo/hi null under three runs. `runs[i]` is the per-run
   values behind it, from /api/dist -- optional, and drawn only when the panel's
   distribution box is ticked. A chart handed no runs draws the interval alone;
   none of them ever fabricates a spread.

   The interval is a bootstrap CI of the MEAN (panels.py:146-156), so it is
   narrower than the run-to-run scatter by construction. That is the whole
   reason the distribution overlay exists: the whiskers answer "where is the
   mean", the points answer "how much do runs differ", and a figure that shows
   only the first invites the second to be read off it.
*/

const CHART_H = 250;

/* Deterministic jitter. A redraw must not move a point: these canvases repaint
   on every slider tick, and Math.random() would make the cloud shimmer. */
function jitter(i, span) {
  const t = ((Math.sin((i + 1) * 12.9898) * 43758.5453) % 1 + 1) % 1;
  return (t - 0.5) * span;
}

const isNum = v => v !== null && v !== undefined && isFinite(v);

/* The y window, over the intervals AND over the points when they are shown --
   a cloud clipped by the frame reads as no cloud at all. */
function yWindow(series, {points, zero, pct}) {
  /* a share is read against the whole: 0-100%. pct "fit" (panel 2) is 0 to
     10% above the largest value drawn, for shares too small to read on 0-100 */
  if (pct && pct !== "fit") return [0, 100];
  let lo = Infinity, hi = -Infinity;
  const see = v => { if (isNum(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } };
  series.forEach(s => {
    (s.vals || []).forEach(b => { if (b) { see(b.mean); see(b.lo); see(b.hi); } });
    if (points) (s.runs || []).forEach(r => (r || []).forEach(see));
  });
  if (!isFinite(lo)) return [0, 1];
  if (zero) lo = Math.min(0, lo);
  if (pct === "fit") return [0, hi > 0 ? hi * 1.1 : 1];
  const m = (hi - lo) * 0.08 || Math.abs(hi) * 0.1 || 1;
  return [zero ? lo : lo - m, hi + m];
}

/* Category labels along the x axis, thinned so they never overlap. */
function xLabels(g, w, h, pad, cats, slotW) {
  g.fillStyle = css("--ink-3");
  g.font = "10px " + css("--f-num");
  const every = Math.max(1, Math.ceil(44 / Math.max(slotW, 1)));
  cats.forEach((c, i) => {
    if (i % every) return;
    const t = String(c);
    const x = pad.l + slotW * (i + 0.5) - g.measureText(t).width / 2;
    g.fillText(t, x, h - pad.b + 12);
  });
}

function chartKey(series) {
  return keyRow(series.filter(s => s.label).map(s => [s.colour, s.label]));
}

/* ------------------------------------------------------------------ dots */

/* Dot and interval: the forest plot. x = category, one offset per series so
   two intervals in the same category do not sit on top of each other. */
function dots(cv, {cats, series, ylab, fmtY, height, points, zero, pct}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 10, t: 10, b: 30};
  const [y0, y1] = yWindow(series, {points, zero, pct});
  axes(g, w, h, pad, "", ylab, y0, y1, fmtY);
  const slot = (w - pad.l - pad.r) / Math.max(cats.length, 1);
  const Y = v => h - pad.b - (h - pad.b - pad.t) * ((v - y0) / (y1 - y0 || 1));
  const step = slot / (series.length + 1);
  cats.forEach((c, ci2) => series.forEach((s, si) => {
    const x = pad.l + slot * ci2 + step * (si + 1);
    const b = (s.vals || [])[ci2];
    if (points) {
      g.fillStyle = s.colour; g.globalAlpha = 0.28;
      ((s.runs || [])[ci2] || []).forEach((v, ri) => {
        if (!isNum(v)) return;
        g.beginPath();
        g.arc(x + jitter(ri, step * 0.7), Y(v), 1.6, 0, 6.2832);
        g.fill();
      });
      g.globalAlpha = 1;
    }
    if (!b || !isNum(b.mean)) return;
    if (isNum(b.lo) && isNum(b.hi)) {
      g.strokeStyle = s.colour; g.lineWidth = 1.5;
      g.beginPath(); g.moveTo(x, Y(b.lo)); g.lineTo(x, Y(b.hi));
      g.moveTo(x - 3, Y(b.lo)); g.lineTo(x + 3, Y(b.lo));
      g.moveTo(x - 3, Y(b.hi)); g.lineTo(x + 3, Y(b.hi)); g.stroke();
    }
    g.fillStyle = s.colour;
    g.beginPath(); g.arc(x, Y(b.mean), 3.2, 0, 6.2832); g.fill();
    g.strokeStyle = css("--panel"); g.lineWidth = 1; g.stroke();
  }));
  xLabels(g, w, h, pad, cats, slot);
  return {g, w, h};
}

/* --------------------------------------------------------------- scatter */

/* x against y, one marked point per group, arms = each axis's own interval.
   `points` draws that group's runs as a cloud behind it; `diagonal` rules y=x,
   which is what makes "above the line" a readable claim. `lines` rules any
   segment, [{x0, y0, x1, y1, colour, label}], or polyline ({pts: [[x, y], ...]}),
   dashed, clipped to the plot, label at its top end. `xticks`: x values labelled
   between the two end labels.

   groups: [{colour, label, x: CI, y: CI, xruns: [], yruns: []}] -- xruns and
   yruns are index-aligned by run, so a point is one run on both axes.
   Optional per group, each defaulting to the scenario-coloured dot: `shape`
   (draw.js SCENARIO_SHAPE values), `fill`, `stroke` (outline), `arm` (whisker
   colour), `alpha` (faded groups are drawn first), `tip` ({text, math} on hover). */
function scatterCI(cv, {groups, xlab, ylab, fmtX, fmtY, height, points,
                        diagonal, zero, pctX, pctY, xlim, ylim, vlines, hlines, lines, xticks}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 12, t: 12, b: 34};
  const win = (pick, runs) => {
    let lo = Infinity, hi = -Infinity;
    const see = v => { if (isNum(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } };
    groups.forEach(gr => {
      const b = gr[pick];
      if (b) { see(b.mean); see(b.lo); see(b.hi); }
      if (points) (gr[runs] || []).forEach(see);
    });
    if (!isFinite(lo)) return [0, 1];
    if (zero) lo = Math.min(0, lo);
    const m = (hi - lo) * 0.1 || Math.abs(hi) * 0.1 || 1;
    return [zero ? lo : lo - m, hi + m];
  };
  let [x0, x1] = win("x", "xruns");
  let [y0, y1] = win("y", "yruns");
  if (pctX) [x0, x1] = [0, 100];
  if (pctY) [y0, y1] = [0, 100];
  if (xlim) [x0, x1] = xlim;
  if (ylim) [y0, y1] = ylim;
  if (diagonal) {
    const lo = Math.min(x0, y0), hi = Math.max(x1, y1);
    x0 = y0 = lo; x1 = y1 = hi;
  }
  axes(g, w, h, pad, xlab, ylab, y0, y1, fmtY);
  const X = v => pad.l + (w - pad.l - pad.r) * ((v - x0) / (x1 - x0 || 1));
  const Y = v => h - pad.b - (h - pad.b - pad.t) * ((v - y0) / (y1 - y0 || 1));
  if (diagonal) {
    g.strokeStyle = css("--line"); g.lineWidth = 1; g.setLineDash([4, 4]);
    g.beginPath(); g.moveTo(X(x0), Y(y0)); g.lineTo(X(x1), Y(y1)); g.stroke();
    g.setLineDash([]);
  }
  /* reference x values: a dashed rule, an optional band, a label at the top */
  (vlines || []).forEach(v => {
    if (!isNum(v.x)) return;
    g.strokeStyle = g.fillStyle = v.colour || css("--ink-3");
    if (v.band && isNum(v.band[0]) && isNum(v.band[1])) {
      g.globalAlpha = 0.12;
      g.fillRect(X(v.band[0]), pad.t, Math.max(1, X(v.band[1]) - X(v.band[0])), h - pad.b - pad.t);
      g.globalAlpha = 1;
    }
    g.lineWidth = 1.2; g.setLineDash(v.dash || [5, 3]);
    g.beginPath(); g.moveTo(X(v.x), pad.t); g.lineTo(X(v.x), h - pad.b); g.stroke();
    g.setLineDash([]);
    if (v.label) {
      g.font = "600 9px " + css("--f-ui");
      const tw = g.measureText(v.label).width || 0;
      const right = X(v.x) + 3 + tw < w - pad.r;
      g.textAlign = right ? "left" : "right";
      g.fillText(v.label, X(v.x) + (right ? 3 : -3), pad.t + 9 + (v.row || 0) * 11);
      g.textAlign = "left";
    }
  });
  (hlines || []).forEach(v => {
    if (!isNum(v.y)) return;
    g.strokeStyle = g.fillStyle = v.colour || css("--ink-3");
    g.lineWidth = 1.2; g.setLineDash(v.dash || [5, 3]);
    g.beginPath(); g.moveTo(pad.l, Y(v.y)); g.lineTo(w - pad.r, Y(v.y)); g.stroke();
    g.setLineDash([]);
    if (v.label) { g.font = "600 9px " + css("--f-ui"); g.textAlign = "right";
      g.fillText(v.label, w - pad.r - 2, Y(v.y) - 3); g.textAlign = "left"; }
  });
  (lines || []).forEach(v => {
    g.strokeStyle = g.fillStyle = v.colour || css("--ink-3");
    g.lineWidth = 1.2; g.setLineDash(v.dash || [5, 3]);
    g.save(); g.beginPath(); g.rect(pad.l, pad.t, w - pad.l - pad.r, h - pad.t - pad.b); g.clip();
    const pts = v.pts || [[v.x0, v.y0], [v.x1, v.y1]];
    g.beginPath(); pts.forEach(([x, y], i) => i ? g.lineTo(X(x), Y(y)) : g.moveTo(X(x), Y(y)));
    g.stroke(); g.restore();
    g.setLineDash([]);
    if (v.label) { const top = pts.reduce((a, p) => p[1] > a[1] ? p : a);
      g.font = "600 9px " + css("--f-ui");
      g.fillText(v.label, X(top[0]) + 4, Math.max(pad.t + 9, Y(top[1]) + 9)); }
  });
  const hits = [];
  const al = gr => isNum(gr.alpha) ? gr.alpha : 1;
  [...groups.filter(gr => al(gr) < 1), ...groups.filter(gr => al(gr) >= 1)].forEach(gr => {
    const a = al(gr), fill = gr.fill || gr.colour;
    if (points) {
      const xs = gr.xruns || [], ys = gr.yruns || [];
      g.fillStyle = fill; g.globalAlpha = (gr.fill ? 0.35 : 0.25) * a;
      for (let i = 0; i < Math.max(xs.length, ys.length); i++) {
        const xv = xs.length ? xs[i] : (gr.x || {}).mean;
        const yv = ys.length ? ys[i] : (gr.y || {}).mean;
        if (!isNum(xv) || !isNum(yv)) continue;
        g.beginPath(); g.arc(X(xv), Y(yv), 1.8, 0, 6.2832); g.fill();
      }
      g.globalAlpha = 1;
    }
    const bx = gr.x, by = gr.y;
    if (!bx || !by || !isNum(bx.mean) || !isNum(by.mean)) return;
    g.globalAlpha = a;
    g.strokeStyle = gr.arm || gr.colour; g.lineWidth = 1.3;
    if (isNum(bx.lo)) { g.beginPath(); g.moveTo(X(bx.lo), Y(by.mean));
                        g.lineTo(X(bx.hi), Y(by.mean)); g.stroke(); }
    if (isNum(by.lo)) { g.beginPath(); g.moveTo(X(bx.mean), Y(by.lo));
                        g.lineTo(X(bx.mean), Y(by.hi)); g.stroke(); }
    g.fillStyle = fill;
    if (gr.shape) traceShape(g, gr.shape, X(bx.mean), Y(by.mean), gr.r || 4.5);
    else { g.beginPath(); g.arc(X(bx.mean), Y(by.mean), gr.r || 4, 0, 6.2832); }
    g.fill();
    g.strokeStyle = gr.stroke || css("--panel"); g.lineWidth = gr.r ? 0.5 : gr.stroke ? 1 : 1.2; g.stroke();
    g.globalAlpha = 1;
    if (gr.tip) hits.push({x: X(bx.mean), y: Y(by.mean), faded: a < 1, text: gr.tip});
  });
  cv._ptHit = hits;
  g.fillStyle = css("--ink-3"); g.font = "10px " + css("--f-num");
  const f = fmtX || (v => fmt(v, 1));
  g.fillText(f(x0), pad.l, h - 16);
  const lt = f(x1);
  g.fillText(lt, w - pad.r - g.measureText(lt).width, h - 16);
  g.textAlign = "center"; (xticks || []).forEach(v => g.fillText(f(v), X(v), h - 16)); g.textAlign = "left";
  return {g, w, h};
}

/* ----------------------------------------------------------- line across */

/* One line per series over an ordinal x -- composition class j, or one axis of
   the build grid. `ribbon` fills each series' interval; `points` draws one
   faint line per run behind it, so a smooth mean over noisy runs is visible as
   such. `mark` rules one level, for the operating point. */
function linesOverLevels(cv, {levels, series, xlab, ylab, fmtY, height,
                              ribbon, points, mark, zero, pct}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 12, t: 10, b: 30};
  const [y0, y1] = yWindow(series, {points, zero, pct});
  axes(g, w, h, pad, xlab, ylab, y0, y1, fmtY);
  const n = Math.max(levels.length, 1);
  const slot = (w - pad.l - pad.r) / n;
  const X = i => pad.l + slot * (i + 0.5);
  const Y = v => h - pad.b - (h - pad.b - pad.t) * ((v - y0) / (y1 - y0 || 1));
  if (mark !== undefined && mark !== null) {
    const i = levels.indexOf(mark);
    if (i >= 0) {
      g.strokeStyle = css("--warn"); g.lineWidth = 1; g.setLineDash([3, 3]);
      g.beginPath(); g.moveTo(X(i), pad.t); g.lineTo(X(i), h - pad.b); g.stroke();
      g.setLineDash([]);
    }
  }
  const path = (get, close) => {
    g.beginPath();
    let on = false;
    for (let i = 0; i < n; i++) {
      const v = get(i);
      if (!isNum(v)) { on = false; continue; }
      if (!on) { g.moveTo(X(i), Y(v)); on = true; } else g.lineTo(X(i), Y(v));
    }
    if (close) g.closePath();
  };
  /* per-run lines first, then the ribbon, then the mean on top */
  if (points) series.forEach(s => {
    const runs = s.runs || [];
    const nr = Math.max(0, ...runs.map(r => (r || []).length));
    g.strokeStyle = s.colour; g.lineWidth = 0.8; g.globalAlpha = 0.16;
    for (let r = 0; r < nr; r++) { path(i => (runs[i] || [])[r]); g.stroke(); }
    g.globalAlpha = 1;
  });
  if (ribbon) series.forEach(s => {
    const v = s.vals || [];
    if (!v.some(b => b && isNum(b.lo))) return;
    g.fillStyle = s.colour; g.globalAlpha = 0.13;
    g.beginPath();
    let on = false;
    for (let i = 0; i < n; i++) {
      const b = v[i];
      if (!b || !isNum(b.lo)) { on = false; continue; }
      if (!on) { g.moveTo(X(i), Y(b.lo)); on = true; } else g.lineTo(X(i), Y(b.lo));
    }
    for (let i = n - 1; i >= 0; i--) {
      const b = v[i];
      if (b && isNum(b.hi)) g.lineTo(X(i), Y(b.hi));
    }
    g.closePath(); g.fill(); g.globalAlpha = 1;
  });
  series.forEach(s => {
    g.strokeStyle = s.colour; g.lineWidth = 1.6;
    path(i => ((s.vals || [])[i] || {}).mean);
    g.stroke();
    g.fillStyle = s.colour;
    for (let i = 0; i < n; i++) {
      const b = (s.vals || [])[i];
      if (b && isNum(b.mean)) {
        g.beginPath(); g.arc(X(i), Y(b.mean), 2.4, 0, 6.2832); g.fill();
      }
    }
  });
  xLabels(g, w, h, pad, levels, slot);
  return {g, w, h};
}

/* --------------------------------------------------------------- stacked */

/* 100% stacked bars: one bar per category, one band per part, each part's mean
   over runs. No interval can be drawn on a stack and no per-run form of a
   stack stays a stack, which is why the distribution box is DISABLED wherever
   this form is selected rather than quietly ignored. Bands that reach 6% of
   the bar carry their label. */
function stacked(cv, {cats, parts, height, pct, fmtV}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 10, t: 10, b: 30};
  const tot = cats.map((c, i) =>
    parts.reduce((a, p) => a + (isNum(p.vals[i]) ? p.vals[i] : 0), 0));
  const hi = pct ? 100 : Math.max(1, ...tot) * 1.05;
  axes(g, w, h, pad, "", pct ? "% of the simulation" : "", 0, hi, v => fmt(v, 0));
  const slot = (w - pad.l - pad.r) / Math.max(cats.length, 1);
  const bw = Math.min(slot - 10, 54);
  const Y = v => h - pad.b - (h - pad.b - pad.t) * (v / (hi || 1));
  const f = fmtV || (v => fmt(v, 0));
  cats.forEach((c, i) => {
    const x = pad.l + slot * i + (slot - bw) / 2;
    const scale = pct && tot[i] > 0 ? 100 / tot[i] : 1;
    let base = 0;
    parts.forEach(p => {
      const v = isNum(p.vals[i]) ? p.vals[i] * scale : 0;
      if (v <= 0) return;
      g.fillStyle = p.colour;
      g.fillRect(x, Y(base + v), bw, Y(base) - Y(base + v));
      if (v / (hi || 1) > 0.06) {
        g.fillStyle = css("--panel");
        g.font = "9px " + css("--f-num");
        const t = f(v);
        g.fillText(t, x + bw / 2 - g.measureText(t).width / 2,
                   (Y(base) + Y(base + v)) / 2 + 3);
      }
      base += v;
    });
    g.strokeStyle = css("--line"); g.lineWidth = 1;
    g.strokeRect(x, Y(base), bw, Y(0) - Y(base));
  });
  xLabels(g, w, h, pad, cats, slot);
  return {g, w, h};
}

/* --------------------------------------------------------------- heatmap */

/* A small matrix, rows x cols, colour by value against `max`. Every cell
   carries its number, because a reader cannot take a value off a ramp. */
function heatcells(cv, {rows, cols, grid, height, fmtV, max, hue}) {
  const {g, w, h} = fit(cv, height || 150);
  const pad = {l: 56, r: 10, t: 18, b: 22};
  const vals = grid.flat().filter(isNum);
  const top = max || (vals.length ? Math.max(...vals) : 1) || 1;
  const cw = (w - pad.l - pad.r) / Math.max(cols.length, 1);
  const ch = (h - pad.t - pad.b) / Math.max(rows.length, 1);
  const base = hue || css("--accent");
  const f = fmtV || (v => fmt(v, 1));
  g.font = "10px " + css("--f-num");
  cols.forEach((c, j) => {
    g.fillStyle = css("--ink-3");
    const t = String(c);
    g.fillText(t, pad.l + cw * (j + 0.5) - g.measureText(t).width / 2, pad.t - 5);
  });
  rows.forEach((r, i) => {
    g.fillStyle = css("--ink-3");
    g.fillText(String(r), 4, pad.t + ch * (i + 0.5) + 3);
    cols.forEach((c, j) => {
      const v = (grid[i] || [])[j];
      const x = pad.l + cw * j, y = pad.t + ch * i;
      g.globalAlpha = isNum(v) ? 0.12 + 0.82 * Math.min(1, v / top) : 0.05;
      g.fillStyle = base;
      g.fillRect(x + 1, y + 1, cw - 2, ch - 2);
      g.globalAlpha = 1;
      g.strokeStyle = css("--line-2"); g.lineWidth = 1;
      g.strokeRect(x + 1, y + 1, cw - 2, ch - 2);
      if (!isNum(v)) return;
      const t = f(v);
      g.fillStyle = v / top > 0.55 ? css("--panel") : css("--ink");
      g.fillText(t, x + cw / 2 - g.measureText(t).width / 2, y + ch / 2 + 3);
    });
  });
  return {g, w, h};
}

/* ------------------------------------------------------------ run counts */

/* Counts, not estimates: denominator bars carry no interval because there is
   nothing to bootstrap -- n_with_rise is a count of runs, not a mean over
   them. One group per category, one bar per count. */
function countBars(cv, {cats, series, ylab, height}) {
  return bars(document.getElementById(cv.id) || cv, {
    cats, ylab, height: height || CHART_H,
    series: series.map(s => ({colour: s.colour,
      vals: s.vals.map(v => isNum(v) ? {mean: v, lo: null, hi: null} : null)}))});
}
