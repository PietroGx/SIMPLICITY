/* The chart forms the panel tabs draw with: a scatter with interval arms and a
   line over ordinal levels with an interval ribbon. Each interval is a
   bootstrap CI of the MEAN, {n, mean, lo, hi} exactly as panels.py sends it,
   with lo/hi null under three runs; no chart ever fabricates a spread. */

const CHART_H = 250;

/* Deterministic jitter. A redraw must not move a point: these canvases repaint
   on every slider tick, and Math.random() would make the cloud shimmer. */
function jitter(i, span) {
  const t = ((Math.sin((i + 1) * 12.9898) * 43758.5453) % 1 + 1) % 1;
  return (t - 0.5) * span;
}

const isNum = v => v !== null && v !== undefined && isFinite(v);

/* The y window, over the intervals --
   a cloud clipped by the frame reads as no cloud at all. */
function yWindow(series, {zero, pct}) {
  /* a share is read against the whole: 0-100%. pct "fit" (panel 2) is 0 to
     10% above the largest value drawn, for shares too small to read on 0-100 */
  if (pct && pct !== "fit") return [0, 100];
  let lo = Infinity, hi = -Infinity;
  const see = v => { if (isNum(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } };
  series.forEach(s => {
    (s.vals || []).forEach(b => { if (b) { see(b.mean); see(b.lo); see(b.hi); } });
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

/* --------------------------------------------------------------- scatter */

/* x against y, one marked point per group, arms = each axis's own interval.
   `diagonal` rules y=x,
   which is what makes "above the line" a readable claim. `lines` rules any
   segment, [{x0, y0, x1, y1, colour, label}], or polyline ({pts: [[x, y], ...]}),
   dashed, clipped to the plot, label at its top end. `xticks`: x values labelled
   between the two end labels.

   groups: [{colour, label, x: CI, y: CI}].
   Optional per group, each defaulting to the scenario-coloured dot: `shape`
   (a draw.js shapePoly name), `fill`, `stroke` (outline), `arm` (whisker
   colour), `alpha` (faded groups are drawn first), `tip` ({text, math} on hover). */
function scatterCI(cv, {groups, xlab, ylab, fmtX, fmtY, height,
                        diagonal, zero, pctX, pctY, xlim, ylim, vlines, hlines, lines, xticks}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 12, t: 12, b: 34};
  const win = pick => {
    let lo = Infinity, hi = -Infinity;
    const see = v => { if (isNum(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } };
    groups.forEach(gr => {
      const b = gr[pick];
      if (b) { see(b.mean); see(b.lo); see(b.hi); }
    });
    if (!isFinite(lo)) return [0, 1];
    if (zero) lo = Math.min(0, lo);
    const m = (hi - lo) * 0.1 || Math.abs(hi) * 0.1 || 1;
    return [zero ? lo : lo - m, hi + m];
  };
  let [x0, x1] = win("x");
  let [y0, y1] = win("y");
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
   the build grid. `ribbon` fills each series' interval; `mark` rules one
   level, for the operating point. */
function linesOverLevels(cv, {levels, series, xlab, ylab, fmtY, height,
                              ribbon, mark, zero, pct}) {
  const {g, w, h} = fit(cv, height || CHART_H);
  const pad = {l: 46, r: 12, t: 10, b: 30};
  const [y0, y1] = yWindow(series, {zero, pct});
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
  /* the ribbon first, then the mean on top */
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
