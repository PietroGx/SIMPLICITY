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

/* The studio says "simulation" where the code and the server say "run"
   (studio only, for now). For text the server sends, such as a panel's note. */
function simWords(t) {
  return String(t == null ? "" : t).replace(/\bRUNS\b/g, "SIMULATIONS").replace(/\bRuns\b/g, "Simulations")
    .replace(/\bruns\b/g, "simulations").replace(/\bRun\b/g, "Simulation")
    .replace(/\brun\b/g, "simulation");
}

function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
/* settings.focus is a comma list of scenarios. When it names any, every other
   scenario fades to ~22% so the chosen ones stand out; the key fades with
   them, since it is drawn from here too. */
function focusSet() {
  const f = typeof settings !== "undefined" ? settings.focus || "" : "";
  return new Set(String(f).split(",").filter(Boolean));
}
/* scenarios with long shedders, one subpanel each; their colour unfaded, since a
   subpanel holds one scenario and has nothing to fade against */
const SUB_SCENARIOS = SCENARIOS.filter(sc => sc !== "control");
const scenarioRaw = sc => css("--sc" + Math.max(0, SCENARIOS.indexOf(sc)));
const scenarioHead = (sc, note = "") => `<h4>${shapeSvg(SCENARIO_SHAPE[sc], scenarioRaw(sc))}${SCENARIO_LABEL[sc]}`
  + `${note ? ` <small style="color:var(--ink-3);font-weight:400">${note}</small>` : ""}</h4>`;
/* figure panel heads: a letter, on a disc in the scenario's colour when `disc`,
   and the scenario's name, capitalised, in its colour when `col` */
const FIG_LETTERS = ["A", "B", "C", "D", "E", "F"];
const cap1 = t => t.charAt(0).toUpperCase() + t.slice(1);
/* white on a dark disc (control's near-black), the ink colour on the rest */
const discInk = c => { const m = /^#?([0-9a-f]{6})$/i.exec(String(c).trim()); if (!m) return "";
  const n = parseInt(m[1], 16), lum = (0.2126 * (n >> 16) + 0.7152 * (n >> 8 & 255) + 0.0722 * (n & 255)) / 255;
  return lum < 0.35 ? ";color:#fff" : ""; };
const figLetter = (sc, letter, disc) => `<b class="c4let${disc ? " on" : ""}"${disc
  ? ` style="background:${scenarioRaw(sc)}${discInk(scenarioRaw(sc))}"` : ""}>${letter}</b>`;
const figName = (sc, t, col) => col ? `<span style="color:${scenarioRaw(sc)}">${t}</span>` : t;
/* a dashed segment as spans, its dash lengths the canvas's, so html2canvas draws it */
const dashSwatch = (on, off, wpx, c = "var(--ink)") => {
  let out = "", x = 0;
  while (x < 20) { const l = Math.min(on, 20 - x);
    out += `<i style="position:absolute;left:${x}px;top:4px;width:${l}px;border-top:${wpx}px solid ${c}"></i>`; x += on + off; }
  return `<i class="c4dash">${out}</i>`;
};
function scenarioColour(scenario) {
  const c = css("--sc" + Math.max(0, SCENARIOS.indexOf(scenario)));
  const f = focusSet();
  return f.size && !f.has(scenario) && /^#[0-9a-f]{6}$/i.test(c) ? c + "38" : c;
}

/* A tab that fits the window: zero the given canvases, lay out, and return the
   pixels the window has left below everything else. The caller shares them out.
   Includes the .wrap's bottom padding; never negative. */
function fitFree(main, els) {
  els.forEach(e => { if (e) e.style.height = "0px"; });
  if (!main.getBoundingClientRect) return 0;
  const r = main.getBoundingClientRect();
  const top = r.top + (window.scrollY || 0);
  const wrap = main.closest ? main.closest(".wrap") : null;
  const pad = wrap ? parseFloat(getComputedStyle(wrap).paddingBottom) || 0 : 16;
  const last = main.lastElementChild;     // its bottom margin falls outside main's box
  const lm = last ? parseFloat(getComputedStyle(last).marginBottom) || 0 : 0;
  return Math.max(0, (window.innerHeight || 900) - top - (r.height || 0) - lm - pad - 2);
}
const clampH = (v, lo, hi) => Math.round(Math.max(lo, Math.min(hi, v)));

/* Size a canvas to its container at the screen's pixel ratio and hand back a
   context already scaled to CSS pixels, so a painter works in CSS units. */
/* set by tab_export.js while it redraws a figure at print resolution */
let EXPORT_DPR = 0;

function fit(cv, cssHeight) {
  const w = cv.parentElement.clientWidth || 620;
  const dpr = EXPORT_DPR || Math.min(window.devicePixelRatio || 1, 2);
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

/* Panel 1's scatter: scenario -> marker shape, composition class -> shade of
   j/k. One geometry for canvas and SVG key, so the two cannot drift. */
const SCENARIO_SHAPE = {control: "circle", SOT: "square", HIV_low: "triangle",
                        HIV_high: "diamond", edge_case: "cross"};
const FADE_ALPHA = 0.22;

function shapePoly(shape, x, y, r) {
  const R = r * 1.2, s = r * 0.85, d = r * 1.15, a = r * 1.1, b = r * 0.4;
  if (shape === "square") return [[x - s, y - s], [x + s, y - s], [x + s, y + s], [x - s, y + s]];
  if (shape === "triangle") return [[x, y - R], [x + R * 0.866, y + R * 0.5], [x - R * 0.866, y + R * 0.5]];
  if (shape === "diamond") return [[x, y - d], [x + d, y], [x, y + d], [x - d, y]];
  if (shape === "cross") return [[x - b, y - a], [x + b, y - a], [x + b, y - b], [x + a, y - b],
    [x + a, y + b], [x + b, y + b], [x + b, y + a], [x - b, y + a], [x - b, y + b],
    [x - a, y + b], [x - a, y - b], [x - b, y - b]];
  return null;
}
function traceShape(g, shape, x, y, r) {
  g.beginPath();
  const p = shapePoly(shape, x, y, r);
  if (!p) { g.arc(x, y, r * 0.95, 0, 6.2832); return; }
  p.forEach(([px, py], i) => i ? g.lineTo(px, py) : g.moveTo(px, py));
  g.closePath();
}
function shapeSvg(shape, fill) {
  const p = shapePoly(shape, 0, 0, 4.5);
  const body = p ? `<polygon points="${p.map(q => q.map(v => v.toFixed(2)).join(",")).join(" ")}"/>`
                 : `<circle r="4.27"/>`;
  return `<svg width="12" height="12" viewBox="-6 -6 12 12" style="fill:${fill};flex:none">${body}</svg>`;
}
/* scenario key by shape, in ink; fades with the focus like the colour key */
function shapeKey() {
  const f = focusSet();
  return `<div class="key">` + SCENARIOS.map(sc =>
    `<span${f.size && !f.has(sc) ? ` style="opacity:${FADE_ALPHA}"` : ""}>`
    + `${shapeSvg(SCENARIO_SHAPE[sc], "var(--ink-2)")}${SCENARIO_LABEL[sc]}</span>`).join("") + `</div>`;
}

/* The long-born share j/k of a class, 0..1; null for wt (0/0, undefined) and
   for other (n_defining is not k, so j/k does not apply). */
function classShare(c, k) {
  return /^\d+$/.test(c) && Number(k) > 0 ? Math.min(1, Number(c) / Number(k)) : null;
}
/* Green to red, as the Clustering tab's "Ratio" origin colouring
   (paint.js originColour): --std mixed into --long by the long-born share,
   in sRGB like its color-mix. The founding clade is --founder there too. */
function mixHex(a, b, t) {
  const p = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
  const A = p(a), B = p(b);
  return "#" + A.map((x, i) => Math.round(x + (B[i] - x) * t)
    .toString(16).padStart(2, "0")).join("");
}
const jShade = t => mixHex(css("--std"), css("--long"), Math.max(0, Math.min(1, t)));
function classShade(c, k) {
  const t = classShare(c, k);
  return t !== null ? jShade(t) : css(c === "wt" ? "--founder" : "--ndef");
}
function classLabel(c, k) {
  return c === "wt" ? "wild type" : c === "other" ? "ndef≠k" : `${c} long of ${k}`;
}
/* the shade's legend: the ramp with a notch at each class j/k, then wild type and ndef≠k */
function jShadeBar(k) {
  const kk = Number(k) || 0;
  const stops = [0, .25, .5, .75, 1].map(t => `${jShade(t)} ${t * 100}%`).join(",");
  const notches = kk > 0 ? Array.from({length: kk + 1}, (_, j) =>
    `<u style="left:${(100 * j / kk).toFixed(2)}%"></u>`).join("") : "";
  return `<div class="jbar"><div class="ramp">
      <b>long-born share of defining mutations (j/k)</b>
      <div class="bar" style="background:linear-gradient(90deg,${stops})">${notches}</div>
      <div class="ticks">${[0, 25, 50, 75, 100].map(v => `<span style="left:${v}%">${v}%</span>`).join("")}</div>
    </div>
    <span class="grey"><i style="background:${css("--founder")}"></i>wild type</span>
    <span class="grey"><i style="background:${css("--ndef")}"></i>ndef≠k</span></div>`;
}

/* How each axis is computed, shown on hovering its label. Keyed by the label
   text, so a label is specific enough to have one meaning. */
const AXIS_HELP = {
  "long-shedder prevalence (% of individuals)": "The scenario's long_shedders_ratio from config, as a percentage: the share of hosts who are long shedders. A design value, not a measurement, so it has no interval.",
  "long-shedder prevalence (%)": "The scenario's long_shedders_ratio from config, as a percentage: the share of hosts who are long shedders. A design value, not a measurement, so it has no interval.",
  "% of simulation frequency AUC, classes j ≥ 1": "Per simulation: the frequency AUC of every clade with at least one defining mutation that arose in a long shedder (classes j ≥ 1), divided by the simulation's total frequency AUC. Frequency AUC is the area under a clade's frequency curve, the trapezoid-rule integral of its frequency over the sampling times (stored as <code>burden</code>). Averaged over simulations; the burn-in toggle sets which days count.",
  "% of the simulation's frequency AUC": "Per simulation: the frequency AUC of this class's clades divided by the simulation's total frequency AUC (the founding clade included, as wild type, so the classes sum to 100%). Frequency AUC is the area under a clade's frequency curve, the trapezoid-rule integral of its frequency over the sampling times (stored as <code>burden</code>). Averaged over simulations.",
  "% of the simulation's clades": "Per simulation: the number of clades in this class divided by the number of clades in the simulation, founding clade included. Averaged over simulations.",
  "% of events with a long-shedder transmitter": "Per simulation: transmission events whose transmitter is a long shedder (long→long + long→standard) divided by all transmission events. Seeded infections are not events. With the burn-in toggle on, only events on or after the simulation's burn-in day count. Averaged over simulations.",
  "% of events in each cell": "Per simulation: the events in one transmitter→recipient cell (long→long, long→standard, standard→long, standard→standard) divided by all transmission events in the window. Averaged over simulations; the four cells sum to 100%.",
  "events / day": "Transmission events per day in each transmitter→recipient cell, averaged over the scenario's simulations day by day.",
  "events per host of that type": "Per simulation: events transmitted by long shedders divided by the number of long shedders infected in the simulation, and the same for standard hosts. Events are counted from the burn-in day when the toggle is on; the host counts are the whole simulation's. Averaged over simulations.",
  "events per long shedder ÷ events per standard host": "Per simulation: the per-host rate of long shedders divided by that of standard hosts. 1 is parity; above 1, each long shedder transmits more than each standard host.",
  "latent interval (days) · clades with a latency": "Per clade: days from the clade's first appearance to the start of its emergence episode (the episode holding its highest frequency). Averaged over clades in a simulation, then over simulations. A clade truncated at the burn-in day has no latency.",
  "rise time (days) · clades with a rise": "Per clade: days for the emergence episode to go from 10% to 90% of its own amplitude, on the smoothed daily record. Averaged over clades in a simulation, then over simulations.",
  "rise time (days)": "Per clade: days for the emergence episode to go from 10% to 90% of its own amplitude, on the smoothed daily record. Averaged over the class's clades in a simulation, then over simulations.",
  "latent interval (days)": "Per clade: days from the clade's first appearance to the start of its emergence episode. Averaged over the class's clades in a simulation, then over simulations.",
  "days · rise time and latent interval": "Rise time: days for the emergence episode to go from 10% to 90% of its own amplitude. Latent interval: days from first appearance to the start of that episode. Each averaged over clades in a simulation, then over simulations.",
  "clades (count)": "Counts of clades, summed over the scenario's simulations: clades with a rise, clades with a latency, and clades truncated at the burn-in day. Counts, not estimates, so no interval.",
  "% of simulation frequency AUC, majority-long classes": "Per simulation: the frequency AUC of clades whose defining mutations arose mostly in long shedders (classes j > k/2) divided by the simulation's total frequency AUC. Averaged over simulations, one value per build.",
  "mean rise time (days)": "Mean emergence rise time over simulations at this build: days from 10% to 90% of the emergence episode's own amplitude.",
  "mean latent interval (days)": "Mean latent interval over simulations at this build: days from a clade's first appearance to the start of its emergence episode.",
  "% of the simulation": "Each band is a share of the whole; a bar sums to 100%.",
  "day": "Days since the start of the simulation.",
  "k": "The shared-mutation threshold: how far a lineage must have diverged from its clade root to found a clade. Each value is its own clade build.",
  "m_T": "The minimum share of the simulation's transmissions a lineage's subtree must carry to found a clade. Each value is its own clade build.",
  "smoothing window": "The width of the centred moving average the daily frequency is read through before episodes are found.",
  "average cryptic outbreaks per simulation": "Per scenario: the mean, over its simulations, of how many cryptic outbreaks each one has. Each segment is one origin group's own mean, so the segments add up to the column.",
  "intra-host hidden time (days)": "Per outbreak: days from the clade's first appearance (its root lineage's emergence inside a host) to its first transmission, the earliest infection that inherited any of its lineages. How long the clade stayed inside the host where it arose.",
  "between-host hidden time (days)": "Per outbreak: days from the clade's first transmission to the first day its smoothed frequency rises above f. How long it spread from host to host before it became visible. Under cryptic rule 2 the two hidden times add up to the hidden time; negative when the clade rose above f before its first transmission.",
  "prominence": "How far a peak must stand above its neighbouring dip, as a fraction of the clade's own peak height, to count as an episode.",
};
/* The same axes as formulas. s is one simulation (S of them), c a clade, f_c(t)
   its frequency at sampling time t, B_c its frequency AUC (trapezoid rule), E_xy the events from
   host type x to y (l long, s standard), φ the burn-in day. */
const AXIS_MATH = {
  "long-shedder prevalence (% of individuals)": [String.raw`\text{prevalence} = \frac{N_{\text{long}}}{N}\times 100`, String.raw`\text{(config: long\_shedders\_ratio}\times 100\text{)}`],
  "long-shedder prevalence (%)": [String.raw`\text{prevalence} = \frac{N_{\text{long}}}{N}\times 100`, String.raw`\text{(config: long\_shedders\_ratio}\times 100\text{)}`],
  "% of simulation frequency AUC, classes j ≥ 1": [String.raw`y = \frac{1}{S}\sum_{s}\frac{\sum_{c:\,j_c\ge 1} B_c}{\sum_{c} B_c}\times 100`, String.raw`B_c = \sum_i \tfrac{1}{2}\big(f_c(t_i)+f_c(t_{i+1})\big)(t_{i+1}-t_i)`],
  "% of the simulation's frequency AUC": [String.raw`\text{share}_j = \frac{1}{S}\sum_{s}\frac{\sum_{c\in j} B_c}{\sum_{c} B_c}\times 100`, String.raw`B_c = \sum_i \tfrac{1}{2}\big(f_c(t_i)+f_c(t_{i+1})\big)(t_{i+1}-t_i)`],
  "% of the simulation's clades": [String.raw`\text{share}_j = \frac{1}{S}\sum_{s}\frac{n_j}{n}\times 100`, String.raw`n_j:\ \text{clades in class } j,\quad n:\ \text{all clades}`],
  "% of events with a long-shedder transmitter": [String.raw`y = \frac{1}{S}\sum_{s}\frac{E_{ll}+E_{ls}}{E_{ll}+E_{ls}+E_{sl}+E_{ss}}\times 100`, String.raw`E_{xy} = \sum_{t\ge\varphi} e_{xy}(t)`],
  "% of events in each cell": [String.raw`y_{xy} = \frac{1}{S}\sum_{s}\frac{E_{xy}}{E_{ll}+E_{ls}+E_{sl}+E_{ss}}\times 100`, String.raw`E_{xy} = \sum_{t\ge\varphi} e_{xy}(t)`],
  "events / day": [String.raw`y_{xy}(t) = \frac{1}{S}\sum_{s} e_{xy,s}(t)`],
  "events per host of that type": [String.raw`I_{\text{long}} = \frac{1}{S}\sum_{s}\frac{E_{ll}+E_{ls}}{N_{\text{long}}}`, String.raw`I_{\text{std}} = \frac{1}{S}\sum_{s}\frac{E_{sl}+E_{ss}}{N_{\text{std}}}`, String.raw`N:\ \text{hosts of that type infected in the simulation}`],
  "events per long shedder ÷ events per standard host": [String.raw`R = \frac{\overline{I_{\text{long}}}}{\overline{I_{\text{std}}}} = \frac{\frac{1}{S}\sum_s (E_{ll}+E_{ls})/N_{\text{long}}}{\frac{1}{S}\sum_s (E_{sl}+E_{ss})/N_{\text{std}}}`],
  "latent interval (days) · clades with a latency": [String.raw`L_c = t_{\text{start}}(\text{emergence episode}) - t_{\text{first}}(c)`, String.raw`y = \frac{1}{S}\sum_{s}\operatorname{mean}_c L_c`],
  "rise time (days) · clades with a rise": [String.raw`A = \tilde y_{\text{peak}} - \tilde y_{\text{start}}`, String.raw`T_c = t\big[\tilde y \ge \tilde y_{\text{start}}+0.9A\big] - t\big[\tilde y \ge \tilde y_{\text{start}}+0.1A\big]`, String.raw`y = \frac{1}{S}\sum_{s}\operatorname{mean}_c T_c`],
  "rise time (days)": [String.raw`A = \tilde y_{\text{peak}} - \tilde y_{\text{start}}`, String.raw`T_c = t\big[\tilde y \ge \tilde y_{\text{start}}+0.9A\big] - t\big[\tilde y \ge \tilde y_{\text{start}}+0.1A\big]`, String.raw`y_j = \frac{1}{S}\sum_{s}\operatorname{mean}_{c\in j} T_c`],
  "latent interval (days)": [String.raw`L_c = t_{\text{start}}(\text{emergence episode}) - t_{\text{first}}(c)`, String.raw`y_j = \frac{1}{S}\sum_{s}\operatorname{mean}_{c\in j} L_c`],
  "days · rise time and latent interval": [String.raw`T_c = t\big[\tilde y \ge \tilde y_{\text{start}}+0.9A\big] - t\big[\tilde y \ge \tilde y_{\text{start}}+0.1A\big]`, String.raw`L_c = t_{\text{start}}(\text{emergence episode}) - t_{\text{first}}(c)`, String.raw`\text{each: } \frac{1}{S}\sum_{s}\operatorname{mean}_c`],
  "clades (count)": [String.raw`n = \sum_{s} n_s`],
  "% of simulation frequency AUC, majority-long classes": [String.raw`y = \frac{1}{S}\sum_{s}\frac{\sum_{c:\,j_c > k/2} B_c}{\sum_{c} B_c}\times 100`, String.raw`B_c = \sum_i \tfrac{1}{2}\big(f_c(t_i)+f_c(t_{i+1})\big)(t_{i+1}-t_i)`],
  "mean rise time (days)": [String.raw`y = \frac{1}{S}\sum_{s}\operatorname{mean}_c T_c`],
  "mean latent interval (days)": [String.raw`y = \frac{1}{S}\sum_{s}\operatorname{mean}_c L_c`],
  "% of the simulation": [String.raw`\text{band} = \frac{\text{part}}{\sum \text{parts}}\times 100`],
  "k": [String.raw`\ell \text{ founds a clade if } d(\ell,\ \text{clade root}) \ge k`, String.raw`\text{and } \frac{T(\text{subtree } \ell)}{T(\text{simulation})} \ge m_T`],
  "m_T": [String.raw`\ell \text{ founds a clade if } \frac{T(\text{subtree } \ell)}{T(\text{simulation})} \ge m_T`, String.raw`\text{and } d(\ell,\ \text{clade root}) \ge k`],
  "smoothing window": [String.raw`\tilde y(t) = \operatorname{mean}\{\, f(\tau) : |\tau - t| \le w/2 \,\}`, String.raw`\text{(w in samples, forced odd; edges shrink)}`],
  "intra-host hidden time (days)": [String.raw`x_c = t_{\text{tx}}(c) - t_{\text{first}}(c)`],
  "between-host hidden time (days)": [String.raw`y_c = \min\{t \ge t_{\text{first}}(c) : \tilde y_c(t) > f\} - t_{\text{tx}}(c)`, String.raw`\text{rule 2: } x_c + y_c = \text{hidden time}`],
  "prominence": [String.raw`\text{peak } p \text{ counts if } \operatorname{prom}(p) \ge \text{prom}\times \max_t \tilde y(t)`],
};

function axisHelp(label) {
  if (AXIS_HELP[label]) return {text: AXIS_HELP[label], tex: AXIS_MATH[label] || null};
  if (/^j long of /.test(label)) return {text: "The composition class: how many of a clade's k defining mutations arose in a long-shedder host, 0 to k. ndef≠k holds clades whose defining-mutation count is not k; wild type is the founding clade.",
    tex: [String.raw`j_c = \#\{\text{defining mutations of } c \text{ that arose in a long shedder}\}`, String.raw`0 \le j_c \le k`]};
  return null;
}

/* One floating tip for every canvas; each canvas keeps its label boxes, and a
   scatter its points (scatterCI's `tip`), an axis label taking precedence. */
function axisTip(cv) {
  if (cv._axBound || !cv.addEventListener) return;
  cv._axBound = true;
  const tip = () => {
    let t = document.getElementById && document.getElementById("axtip");
    if (!t && document.createElement && document.body && document.body.append) {
      t = document.createElement("div"); t.id = "axtip"; t.className = "axtip";
      document.body.append(t);
    }
    return t;
  };
  cv.addEventListener("mousemove", ev => {
    const r = cv.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
    const ax = (cv._axHit || []).find(b => x >= b.x0 && x <= b.x1 && y >= b.y0 && y <= b.y1);
    let hit = ax, best = 8;
    if (!ax) (cv._ptHit || []).forEach(p => {
      const dd = Math.hypot(p.x - x, p.y - y) + (p.faded ? 4 : 0);
      if (dd < best) { best = dd; hit = p; }
    });
    const t = tip(); if (!t) return;
    cv.style.cursor = ax ? "help" : "";
    if (!hit) { t.style.display = "none"; return; }
    t.textContent = "";
    const p = document.createElement("div"); p.textContent = hit.text.text;
    t.append(p);
    if (hit.text.math) { const m = document.createElement("pre");
      m.className = "axmath"; m.textContent = hit.text.math; t.append(m); }
    /* LaTeX through KaTeX (index.html loads it); its source if KaTeX is absent */
    if (hit.text.tex) { const m = document.createElement("div");
      m.className = "axmath axtex"; t.append(m);
      hit.text.tex.forEach(src => { const line = document.createElement("div");
        m.append(line);
        try { if (window.katex) window.katex.render(src, line, {displayMode: true, throwOnError: false});
              else line.textContent = src; }
        catch (e) { line.textContent = src; } }); }
    t.style.display = "block";
    t.style.left = Math.min(ev.clientX + 12, window.innerWidth - 450) + "px";
    t.style.top = (ev.clientY + 14) + "px";
  });
  cv.addEventListener("mouseleave", () => { const t = tip(); if (t) t.style.display = "none"; });
}

/* The y ticks right-aligned against the axis and the y label rotated to their
   left, so neither overlaps; the margin widens to fit both. Mutates pad.l --
   every caller maps x only after this returns. */
function axes(g, w, h, pad, xlab, ylab, y0, y1, fmtY) {
  const n = 4;
  const tick = i => { const v = y0 + (y1 - y0) * i / n; return fmtY ? fmtY(v) : v.toFixed(0); };
  g.font = "10px " + css("--f-num");
  let tw = 0;
  for (let i = 0; i <= n; i++) tw = Math.max(tw, g.measureText(String(tick(i))).width || 0);
  const LAB = ylab ? 18 : 2;
  pad.l = Math.max(pad.l, Math.ceil(LAB + tw + 8));
  g.strokeStyle = css("--line"); g.lineWidth = 1;
  g.beginPath(); g.moveTo(pad.l, pad.t); g.lineTo(pad.l, h - pad.b);
  g.lineTo(w - pad.r, h - pad.b); g.stroke();
  g.fillStyle = css("--ink-3");
  g.textAlign = "right";
  for (let i = 0; i <= n; i++) {
    const y = h - pad.b - (h - pad.b - pad.t) * i / n;
    g.fillText(String(tick(i)), pad.l - 5, y + 3);
    if (i) { g.strokeStyle = css("--line-2"); g.beginPath();
             g.moveTo(pad.l, y); g.lineTo(w - pad.r, y); g.stroke(); }
  }
  g.textAlign = "left";
  g.font = "600 11px " + css("--f-ui");
  g.fillStyle = css("--ink-2");
  const hits = [];
  const cv = g.canvas;
  if (xlab) {
    const xw = g.measureText(xlab).width || 0;
    const cx = (pad.l + w - pad.r) / 2;
    g.textAlign = "center"; g.fillText(xlab, cx, h - 3); g.textAlign = "left";
    const help = axisHelp(xlab);
    if (help) hits.push({x0: cx - xw / 2, x1: cx + xw / 2, y0: h - 15, y1: h, text: help});
  }
  if (ylab) {
    /* fit the rotated label to the plot's height: smaller type first, then cut */
    const room = h - pad.t - pad.b;
    let text = ylab;
    for (const px of [11, 10, 9]) {
      g.font = `600 ${px}px ` + css("--f-ui");
      if ((g.measureText(text).width || 0) <= room) break;
    }
    while (text.length > 4 && (g.measureText(text).width || 0) > room)
      text = text.replace(/…$/, "").slice(0, -1) + "…";
    const yw = g.measureText(text).width || 0;
    const cy = (pad.t + h - pad.b) / 2;
    g.save(); g.translate(12, cy); g.rotate(-Math.PI / 2);
    g.textAlign = "center"; g.fillText(text, 0, 0); g.restore();
    const help = axisHelp(ylab);
    if (help) hits.push({x0: 0, x1: 16, y0: cy - yw / 2, y1: cy + yw / 2, text: help});
  }
  if (cv) { cv._axHit = hits; axisTip(cv); }
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
  hi = pct === "fit" ? (hi > 0 ? hi * 1.1 : 1) : pct ? 100 : hi > 0 ? hi * 1.12 : 1;
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
  if (mini || day == null) { g.restore(); return; }   /* day null: the line alone, its key elsewhere */
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
