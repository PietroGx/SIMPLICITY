/* Drawing primitives shared by every tab: the canvas plumbing, the palette and
   the number formatting.

   One source of truth for three things that were previously repeated per
   chart:

     fit()        device-pixel ratio, size, clearing, and the transform, so a
                  painter never touches cv.width or devicePixelRatio
     css()        every colour comes from a CSS custom property, so the theme
                  lives in studio.css alone and dark mode needs no JS
     fmt()/ci()   one rendering of a number and of a bootstrap interval

   `bars` and `linesOverDays` are the generic charts. A phylogeny, a Muller plot
   and an episode grid are not expressible through them, so those live in
   paint.js as their own painters -- but they still take their canvas from
   fit() and their colours from css().
*/

const SCENARIOS = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"];
const SCENARIO_LABEL = {control: "control", SOT: "SOT", HIV_low: "HIV low",
                        HIV_high: "HIV high", edge_case: "edge case"};
/* The build grid's axes, filled by loadAxes() from /api/index, which reads
   them off sweep.py -- the script that builds the grid.

   Nothing here is hardcoded on purpose. A build parameter names a directory on
   disk, so a second copy of the list in the browser is a list that can drift
   from the grid, and the k slider already had: it offered [0,3,5,7,10,15,20]
   against a grid of [5,7,10]. Empty until the fetch lands, which is why
   index.html awaits it before the first render. */
const AXES = {};

async function loadAxes() {
  try {
    const ix = await apiGet("/api/index");
    if (ix && ix.axes) Object.assign(AXES, ix.axes);
  } catch (e) { /* the rail says so: a build slider with no axis is disabled */ }
}
/* m_T is stored as a fraction and shown as a percentage. One conversion,
   so a panel label and a slider cannot disagree about the unit. */
const mtPct = v => +(Number(v) * 100).toFixed(3);

/* Every share in a panel payload is a FRACTION, 0..1 -- panels.py formed
   percentages until v2.5.16 and the field names said `pct_`, so a value and
   its label could disagree about the unit. The server no longer converts; the
   studio does, here and only here.

   `asPct` scales one {n, mean, lo, hi} and keeps its nulls, which matter: lo
   and hi are null under three runs and a chart must not draw a whisker there.
   `asPctRuns` does the same for a per-run array from /api/dist. */
const asPct = b => !b ? b : {
  n: b.n,
  mean: b.mean === null || b.mean === undefined ? b.mean : b.mean * 100,
  lo: b.lo === null || b.lo === undefined ? b.lo : b.lo * 100,
  hi: b.hi === null || b.hi === undefined ? b.hi : b.hi * 100,
};

const asPctNum = v => v === null || v === undefined ? v : v * 100;

const asPctRuns = rs => !rs ? rs
  : rs.map(r => !r ? r : r.map(asPctNum));
/* m_F: the peak frequency an episode must reach to be counted. A READ
   parameter -- episodes are stored, so it filters a build it does not name. */
const MF_VALUES = [0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50];

function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
function scenarioColour(scenario) {
  return css("--sc" + Math.max(0, SCENARIOS.indexOf(scenario)));
}

/* Size a canvas to its container at the screen's pixel ratio and hand back a
   context already scaled to CSS pixels, so a painter works in CSS units. */
function fit(cv, cssHeight) {
  const w = cv.parentElement.clientWidth || 620;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(cssHeight * dpr);
  cv.style.height = cssHeight + "px";
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, cssHeight);
  return {g, w, h: cssHeight};
}

function fmt(v, dp) {
  if (v === null || v === undefined || !isFinite(v)) return "—";
  return Number(v).toFixed(dp === undefined ? 1 : dp);
}

/* A bootstrap interval as the panels return it: {mean, lo, hi} or nulls. */
function ci(b, dp, unit) {
  if (!b || b.mean === null || b.mean === undefined) return "—";
  const u = unit || "";
  if (b.lo === null || b.lo === undefined) return `${fmt(b.mean, dp)}${u}`;
  return `${fmt(b.mean, dp)}${u} (${fmt(b.lo, dp)}–${fmt(b.hi, dp)})`;
}

function errorBox(message) {
  const d = document.createElement("div");
  d.className = "err";
  d.textContent = message;
  return d;
}

function keyRow(items) {
  return `<div class="key">` + items.map(([colour, text]) =>
    `<span><i style="background:${colour}"></i>${text}</span>`).join("") + `</div>`;
}

function axes(g, w, h, pad, xlab, ylab, y0, y1, fmtY) {
  g.strokeStyle = css("--line"); g.lineWidth = 1;
  g.beginPath(); g.moveTo(pad.l, pad.t); g.lineTo(pad.l, h - pad.b);
  g.lineTo(w - pad.r, h - pad.b); g.stroke();
  g.fillStyle = css("--ink-3");
  g.font = "10px " + css("--f-num");
  const n = 4;
  for (let i = 0; i <= n; i++) {
    const v = y0 + (y1 - y0) * i / n;
    const y = h - pad.b - (h - pad.b - pad.t) * i / n;
    g.fillText((fmtY ? fmtY(v) : v.toFixed(0)), 3, y + 3);
    if (i) { g.strokeStyle = css("--line-2"); g.beginPath();
             g.moveTo(pad.l, y); g.lineTo(w - pad.r, y); g.stroke(); }
  }
  g.font = "10px " + css("--f-ui");
  if (xlab) g.fillText(xlab, pad.l, h - 3);
  if (ylab) { g.save(); g.translate(10, pad.t + 4); g.rotate(-Math.PI / 2);
              g.fillText(ylab, -(h - pad.b - pad.t), 0); g.restore(); }
}

/* Grouped bars with bootstrap whiskers: x = category, one series per scenario. */
function bars(cv, {cats, series, ylab, fmtY, height, pct}) {
  const {g, w, h} = fit(cv, height || 240);
  const pad = {l: 44, r: 8, t: 10, b: 30};
  let hi = 0;
  series.forEach(s => s.vals.forEach(b => {
    if (b && b.hi !== null && b.hi !== undefined) hi = Math.max(hi, b.hi);
    if (b && b.mean) hi = Math.max(hi, b.mean);
  }));
  if (pct) hi = Math.max(hi, 5);
  hi = hi > 0 ? hi * 1.12 : 1;
  axes(g, w, h, pad, "", ylab, 0, hi, fmtY);
  const iw = (w - pad.l - pad.r) / Math.max(cats.length, 1);
  const bw = Math.max(2, (iw - 6) / Math.max(series.length, 1));
  const Y = v => h - pad.b - (h - pad.b - pad.t) * (v / hi);
  cats.forEach((c, ci2) => {
    const x0 = pad.l + iw * ci2 + 3;
    series.forEach((s, si) => {
      const b = s.vals[ci2];
      if (!b || b.mean === null || b.mean === undefined) return;
      const x = x0 + bw * si;
      g.fillStyle = s.colour;
      g.fillRect(x, Y(b.mean), bw - 1, h - pad.b - Y(b.mean));
      if (b.lo !== null && b.lo !== undefined) {
        g.strokeStyle = css("--ink-2"); g.lineWidth = 1;
        const xm = x + (bw - 1) / 2;
        g.beginPath(); g.moveTo(xm, Y(b.lo)); g.lineTo(xm, Y(b.hi));
        g.moveTo(xm - 2, Y(b.lo)); g.lineTo(xm + 2, Y(b.lo));
        g.moveTo(xm - 2, Y(b.hi)); g.lineTo(xm + 2, Y(b.hi)); g.stroke();
      }
    });
    g.fillStyle = css("--ink-3");
    g.font = "10px " + css("--f-num");
    const t = String(c);
    g.fillText(t, x0 + (iw - 6) / 2 - g.measureText(t).width / 2, h - pad.b + 12);
  });
  return {g, w, h};
}

/* One line per series over a day axis. */
function linesOverDays(cv, {series, ylab, fmtY, height, xmax, marks}) {
  const {g, w, h} = fit(cv, height || 240);
  const pad = {l: 44, r: 8, t: 10, b: 26};
  let hi = 0, n = 0;
  series.forEach(s => { n = Math.max(n, s.y.length);
                        s.y.forEach(v => { if (isFinite(v)) hi = Math.max(hi, v); }); });
  hi = hi > 0 ? hi * 1.1 : 1;
  const XM = xmax || n;
  axes(g, w, h, pad, "day", ylab, 0, hi, fmtY);
  const X = d => pad.l + (w - pad.l - pad.r) * (d / Math.max(XM - 1, 1));
  const Y = v => h - pad.b - (h - pad.b - pad.t) * (v / hi);
  (marks || []).forEach(m => {
    g.strokeStyle = css("--warn"); g.setLineDash([3, 3]);
    g.beginPath(); g.moveTo(X(m), pad.t); g.lineTo(X(m), h - pad.b); g.stroke();
    g.setLineDash([]);
  });
  series.forEach(s => {
    g.strokeStyle = s.colour; g.lineWidth = 1.4; g.beginPath();
    let on = false;
    s.y.forEach((v, i) => {
      if (!isFinite(v)) { on = false; return; }
      const x = X(i), y = Y(v);
      if (!on) { g.moveTo(x, y); on = true; } else g.lineTo(x, y);
    });
    g.stroke();
  });
  g.fillStyle = css("--ink-3"); g.font = "10px " + css("--f-num");
  g.fillText("0", pad.l, h - 8);
  const lt = String(XM - 1);
  g.fillText(lt, w - pad.r - g.measureText(lt).width, h - 8);
  return {g, w, h};
}

/* The burn-in day: a haloed dashed rule with a label, `mini` for a grid cell
   where there is no room for one. Taken verbatim from
   tests/artifacts/fig2_studio.html, which carried two identical copies of it
   (lines 1707 and 2763); this is the only one. */
function drawBurnin(g, x, w, h, day, mini) {
  g.save();
  /* halo first, so the line reads over dense colour */
  g.globalAlpha = .8; g.strokeStyle = css("--panel"); g.lineWidth = 4;
  g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke();
  g.globalAlpha = 1; g.strokeStyle = css("--ink-2"); g.lineWidth = 1.25;
  g.setLineDash(mini ? [2, 2] : [5, 3]);
  g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke();
  g.setLineDash([]);
  if (mini) { g.restore(); return; }
  const t = `burn-in \u00b7 day ${Math.round(day)}`;
  g.font = "500 10px " + css("--f-num");
  const tw = g.measureText(t).width;
  const bx = (x + 6 + tw + 4 > w) ? x - 6 - tw - 3 : x + 6;
  g.globalAlpha = .88; g.fillStyle = css("--panel");
  g.fillRect(bx - 3, 1, tw + 6, 14);
  g.globalAlpha = 1; g.fillStyle = css("--ink-2");
  g.fillText(t, bx, 11);
  g.restore();
}
