/* Outbreak detection and Cryptic outbreaks: two tabs, one run each, one painter.

   Outbreak detection (id "sweep") finds each clade's episodes and picks its
   outbreak: smoothing, which peaks count, where an episode starts, F. Cryptic
   outbreaks (id "cryptic") takes those outbreaks as they are and asks which
   are cryptic: where hidden time starts, H, f. Both measure through
   /api/explore (figures_gui/explore.py) on the run's daily record, with no
   clade build behind them, so every rule answers immediately. Figure 5's
   default source measures every simulation with these same rules.

   The stacked plot is drawn from the per-seed export; the cells are measured
   and drawn on the daily record, on the same days.
*/

const SWEEP_MAX_CELLS = 60;
const DAILY_CAP = 200;          // compute.MAX_DAILY_CLADES

/* The clade each tab's grid was last scrolled to, so a redraw that did not
   change the selection does not scroll again. */
const sweepScrolledTo = {};

/* The last Muller painted on Outbreak detection, kept so check_app.js can pin
   the band hit test to the geometry that was actually drawn. */
let sweepMuller = null;

/* ------------------------------------------------- detection rules (sd*)

   These tabs' settings are their own (sd*), never the build parameters
   win/prom/lo: dragging the window to 45 days here must not point every other
   tab at a build that does not exist. SD_FIG5 is what Figure 5 measures with, not
   the build's rules: exploreQuery sends gate=start, dip=low; the build is trunc, any (explore.py FIG5). */
const SD_WIN = [1, 2, 3, 5, 7, 10, 15, 21, 30, 45, 60, 90, 120, 180];
const SD_PROM = [0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9];
const SD_TROUGH = [0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5];
const SD_FRAC = [0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5];
const BRK_LO = [0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2];
const SD_FIG5 = {sdKernel: "box", sdWidth: 0,
                 sdCens: 1, sdStart: "trough", sdTrough: 0.0001, sdPromRef: "peak",
                 sdPick: "first", sdGate: "start"};
const SD_LABEL = {sdKernel: "kernel", sdWin: "window",
  sdProm: "prominence", sdPromRef: "prominence reference", sdWidth: "width", sdCens: "censored peaks", sdStart: "episode start",
  sdTrough: "dip depth", sdPick: "outbreak pick",
  sdGate: "burn-in rule"};

/* Figure 5's rules: SD_FIG5 plus the build's own window and prominence */
const sweepFig5Patch = () => Object.assign({}, SD_FIG5,
  {sdWin: +settings.win, sdProm: +settings.prom});
const sweepOffRules = () => Object.entries(sweepFig5Patch())
  .filter(([k, v]) => settings[k] !== v).map(([k]) => k);

const sdEsc = v => String(v).replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"})[c]);
const sdMedian = a => {
  if (!a.length) return NaN;
  const s = a.slice().sort((x, y) => x - y), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};

/* Each card draws the mark it moves, as the cells draw it (paintOne, drawBreakout). */
const SW_PATH = "1,10 6,10 11,2.5 16,8 21,8";
const swIcon = body => `<svg class="mk-ic" width="22" height="12" viewBox="0 0 22 12">${body}</svg>`;
const swLine = (stroke = "var(--ink-2)") =>
  `<polyline points="${SW_PATH}" fill="none" style="stroke:${stroke}" stroke-width="1.3"/>`;
const swDot = (x, y, fill, r = 2.6) =>
  `<circle cx="${x}" cy="${y}" r="${r}" style="fill:${fill};stroke:var(--panel)"/>`;
const swEp = `<polyline points="6,10 11,2.5" fill="none" style="stroke:var(--std)" stroke-width="2.6" stroke-linecap="round"/>`;
const SW_MARK = {
  curve: swIcon(swLine("var(--std)")),
  peak: swIcon(swLine() + swEp),
  start: swIcon(swLine() + swEp + swDot(6, 10, "var(--long)", 2)),
  breakout: swIcon(`<line x1="1" y1="2.5" x2="21" y2="2.5" style="stroke:var(--ink-3)" stroke-dasharray="4 3"/>`
    + swLine() + swDot(11, 2.5, "var(--ink)", 2.8)),
  hidden: swIcon(`<rect x="3" y="0" width="13" height="12" style="fill:var(--long)" opacity=".18"/>`
    + `<line x1="3" y1="0" x2="3" y2="12" style="stroke:var(--ink-2)"/>`),
  cryptic: swIcon(`<rect x="1.5" y="1" width="19" height="10" rx="1.5" style="fill:var(--long)" opacity=".16"/>`),
};
const swKeyIcon = body => `<svg width="22" height="12" viewBox="0 0 22 12">${body}</svg>`;
/* The card a reader clicked, per tab: the cells emphasise what it moves */
const swFocus = {sweep: null, cryptic: null};
function swCard(mark, title, moves, children) {
  const el = document.createElement("div");
  el.className = "card sw-focusable";
  el.dataset.mark = mark;
  el.innerHTML = `<h2 data-hint="${sdEsc(moves + " Click the card to highlight it on the plots.")}">${
    SW_MARK[mark] || ""}${title}</h2>`;
  children.filter(Boolean).forEach(c => el.append(c));
  el.onclick = ev => {
    const t = ev.target;
    if (t.closest && t.closest("input, select, button, label, .seg, .val") && !t.closest("h2")) return;
    const rail = el.closest(".rail"), mode = rail && rail.id.replace(/_rail$/, "");
    if (!(mode in swFocus)) return;
    swFocus[mode] = swFocus[mode] === mark ? null : mark;
    swMarkFocus(rail, mode);
    renderTab(tabById(mode));
  };
  return el;
}
function swMarkFocus(rail, mode) {
  rail.querySelectorAll(".sw-focusable").forEach(c =>
    c.classList.toggle("sw-on", c.dataset.mark === swFocus[mode]));
}

/* The two tabs: their DOM prefix, and which view settings are their own */
const SW_MODES = {
  sweep:   {p: "rl", cryptic: false, show: "sdShowDet", eps: "sdEpsDet"},
  cryptic: {p: "cr", cryptic: true,  show: "sdShow",    eps: "sdEps"},
  crystats: {p: "cs", cryptic: true,  show: "sdShow",    eps: "sdEps"},
};

/* --------------------------------------------------------------- rails

   Every explanation is a hover popup (controls.js wireHintPopups), so the
   rails hold controls only and fit one screen. The view controls sit above
   the cells, in the main area (swTools). */

const SW_WIDTH = [0, 14, 21, 30, 45, 60, 90];

function swSimCard() {
  const seed = sliderControl({setting: "seed", label: "Seed", min: 0, max: 49, step: 1,
      format: v => String(v).padStart(4, "0"),
      hint: "One simulation, the independent unit of every analysis. Clades are "
          + "defined at the k and m<sub>T</sub> shown top right, set on Clustering. "
          + "Red ticks: simulations that never reach the burn-in day."});
  seed.classList.add("sw-seed");
  return card("Simulation", [
    selectControl({setting: "scenario", label: "Scenario",
      options: SCENARIOS.map(s => [s, SCENARIO_LABEL[s]]),
      hint: "Scenarios are never pooled: each is its own design point."}),
    seed,
    swBurninControl(),
  ]);
}

/* One menu for two settings: burnin (shared with every tab) and sdGate */
function swBurninControl() {
  const cur = !settings.burnin ? "off" : settings.sdGate;
  const opts = [["start", "Outbreak starts after it"], ["trunc", "Not circulating at it"],
                ["off", "Ignore it"]];
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `<label>Burn-in</label>
    <select>${opts.map(([v, t]) => `<option value="${v}"${v === cur ? " selected" : ""}>${t}</option>`).join("")}</select>
    <div class="hint">Which clades the burn-in day, the first day everyone has been
      infected once, lets through. <b>Outbreak starts after it</b>: measured on the
      whole record, a clade counts only if its outbreak begins after the burn-in day;
      one that broke out before it is left out. <b>Not circulating at it</b> (Figure 5):
      measured from the burn-in day, a clade already above 0.01% then is left out, even
      if its outbreak comes later. <b>Ignore it</b>: the whole record, every clade.</div>`;
  el.querySelector("select").addEventListener("change", ev => {
    const v = ev.target.value;
    setSettings(v === "off" ? {burnin: false} : {burnin: true, sdGate: v});
  });
  return el;
}

function swRailStart(rail) {
  railCompact(rail, {toggle: false});
  rail.classList.add("hintpop-on");
  wireHintPopups();
  rail.append(swSimCard());
}

function detectionRail(rail) {
  swRailStart(rail);
  rail.append(swCard("curve", "1 · Smoothing",
    "The coloured line. Every mark is measured on it, not on the raw record.", [
    selectControl({setting: "sdKernel", label: "Kernel",
      options: SD_KERNELS,
      hint: "How the daily record is smoothed. <b>Moving average</b>: equal weights. "
          + "<b>Gaussian</b>: weights the centre, keeps short dips. "
          + "<b>Median</b>: ignores one-day spikes."}),
    sliderControl({setting: "sdWin", label: "Window", values: SD_WIN,
      format: v => v + " d",
      hint: "How many days the kernel spans. Wider merges nearby peaks."}),
  ]));
  rail.append(swCard("peak", "2 · Episodes",
    "Which peaks count. Each episode is drawn bold from its start to its peak.", [
    sliderControl({setting: "sdProm", label: "Prominence",
      values: SD_PROM, format: v => Number(v).toFixed(2),
      hint: "How far a peak must rise above the higher of its two neighbouring lows. "
          + "Lower finds more, smaller peaks."}),
    segmentControl({setting: "sdPromRef", label: "Prominence relative to",
      options: [
        ["peak", "Each peak", "A share of the peak's own height: a small peak before a big one can count."],
        ["clade", "Clade maximum", "A share of the clade's highest frequency: a small peak before a big one "
          + "rarely counts. The rule before v2.5.25."]],
      hint: "What the prominence is a share of."}),
    sliderControl({setting: "sdWidth", label: "Min width",
      values: SW_WIDTH, format: v => v ? v + " d" : "any",
      hint: "Drops peaks narrower than this, measured at half their prominence."}),
    checkControl({setting: "sdCens", label: "Peak still high at the end",
      hint: "Counts a clade still at its maximum when the record ends. Its peak is a lower bound."}),
  ]));
  const start = settings.sdStart;
  rail.append(swCard("start", "3 · Episode start",
    "Where each episode begins: the red dot.", [
    segmentControl({setting: "sdStart", label: "Rule",
      options: [
        ["trough", "Last dip", "The low point just before the peak. With the launch rule it leaves "
          + "no cryptic outbreaks: the clade was above that low earlier."],
        ["frac", "Below %", "The last day the clade was under a share of the peak, set below."],
        ["base", "Base", "scipy's prominence base: the lowest point back to a taller peak. Often far left."]],
      hint: "Where each episode begins."}),
    start === "frac" ? sliderControl({setting: "sdFrac", label: "Below",
      values: SD_FRAC, format: v => fmt(100 * v, v < 0.01 ? 1 : 0) + "% of peak",
      hint: "The episode starts on the last day the clade was under this share of its peak."}) : null,
    start === "trough" ? sliderControl({setting: "sdTrough", label: "Dip depth",
      values: SD_TROUGH, format: v => v === 0 ? "any" : fmt(100 * v, v < 0.01 ? 2 : 0) + "% of peak",
      hint: "How deep a dip must be to start an episode, as a share of the clade's peak."}) : null,
  ]));
  rail.append(swCard("breakout", "4 · The outbreak",
    "Which episode is the outbreak. F is the dashed line.", [
    sliderControl({setting: "brkF", label: "Peak ≥ <i>F</i>", values: BRK_F,
      format: v => `${Math.round(100 * v)}%`,
      hint: "An episode can be the outbreak only if its peak reaches F, as a share of all virus that day."}),
    segmentControl({setting: "sdPick", label: "Which",
      options: [
        ["first", "First ≥ F", "The first episode to reach F. Later peaks do not matter."],
        ["largest", "Largest ≥ F", "The tallest episode to reach F. Earlier peaks fall in its hidden time."]],
      hint: "Which episode reaching F is the outbreak."}),
  ]));
}

function crypticRail(rail) {
  swRailStart(rail);
  crypticCards(rail);
}
function crypticCards(rail) {
  rail.append(swCard("hidden", "1 · Hidden time",
    "The red band: from the solid line to the outbreak's start, the red dot.", [
    (() => { const d = document.createElement("div");
      d.innerHTML = `<div class="hint" style="display:block">From <b>emergence in host</b>, when the
        clade's root lineage arose inside a host.</div>`; return d; })(),
  ]));
  const under = settings.brkRule === "under";
  rail.append(swCard("cryptic", "2 · Cryptic if",
    "Tinted cells: cryptic outbreaks under the rule chosen.", [
    segmentControl({setting: "brkRule", label: "Rule",
      options: [
        ["launch", "Rule 1", "Hidden until the outbreak starts. With the box below ticked, the curve "
          + "must also stay under its launch level while hidden."],
        ["under", "Rule 2", "Hidden while the curve stays under f: hidden time ends the first day "
          + "it rises above f."]],
      hint: "What hidden means."}),
    sliderControl({setting: "brkH", label: "Hidden ≥ <i>H</i>", values: BRK_H,
      format: v => v + " d",
      hint: "Cryptic only if the hidden time, the red band, lasts at least H."}),
    sliderControl({setting: "brkLo", label: under ? "Under <i>f</i>" : "<i>f</i> (Figure 5)", values: BRK_LO,
      format: v => brkPct(v),
      hint: (under ? "Rule 2's ceiling: the clade is hidden while its smoothed frequency stays at or under f. "
          : "Under rule 1, f decides no cryptic outbreak. ")
          + "Figure 5's between-host hidden time ends on the first day the smoothed frequency rises above f."}),
    under ? null : checkControl({setting: "brkQuiet", label: "No rise above the launch before it",
      hint: "Hidden time counts only if the curve never rose above the launch level, "
          + "the red dotted line, during it."}),
  ]));
}

/* The view controls, above the cells: rebuilt on every render, so they
   always show the current values */
function swTools(mode, box, brkOn) {
  const M = SW_MODES[mode], O = SW_OV[mode], ov = swOvOn(mode);
  const grid = !ov && settings.episodeView === "grid";
  const cellOpts = M.cryptic ? [["cryptic", "Cryptic"], ["breakout", "Outbreaks"], ["all", "All"]]
                             : [["breakout", "Outbreaks"], ["all", "All"]];
  [
    O ? checkControl({setting: O.on, label: "Whole scenario",
      hint: "Every outbreak in every simulation of this scenario, one cell each, "
          + "under the rules in force. Click a cell to open its simulation."}) : null,
    ov ? segmentControl({setting: "sdOvCols", label: "Columns",
      options: [4, 6, 8, 10].map(n => [n, String(n)])}) : null,
    ov ? segmentControl({setting: O.show, label: "Show", options: O.shows,
      hint: "Which cells to show."}) : null,
    ov && settings[O.show] === O.filtered ? segmentControl({setting: O.why, label: "Why",
      options: swOvWhys(O), hint: "Show only the cells left out for one reason."}) : null,
    brkOn && !ov ? segmentControl({setting: M.show, label: "Show", options: cellOpts,
      hint: M.cryptic ? "Cryptic outbreaks first, then other outbreaks, then the rest."
                      : "Clades with an outbreak first, then the rest."}) : null,
    checkControl({setting: M.eps, label: "Episodes",
      hint: "Every episode drawn bold in the clade's colour, from its start to its "
          + "peak; the rest of the curve fades."}),
    checkControl({setting: "showRaw", label: "Raw",
      hint: "The unsmoothed daily record in grey behind the smoothed line: how "
          + "much the kernel is doing."}),
    ov ? null : segmentControl({setting: "episodeView", label: "Layout",
      options: [["grid", "Cells"], ["one", "Overlay"]],
      hint: "<b>Cells</b>: one plot per clade. <b>Overlay</b>: every clade on one "
          + "axis; click a band above to mark one."}),
    grid ? segmentControl({setting: "episodeCols", label: "Columns",
      options: [1, 2, 3, 4, 5, 6].map(n => [n, String(n)])}) : null,
    grid || ov ? segmentControl({setting: "yScale", label: "Height",
      options: [["shared", "Shared"], ["own", "Own"]],
      hint: "<b>Shared</b>: one ceiling, so cells compare. <b>Own</b>: each "
          + "clade's own peak, to see a small clade's shape."}) : null,
    segmentControl({setting: "colourBy", label: "Colour",
      options: [["origin", "Origin"], ["clade", "Clade"]],
      hint: "<b>Origin</b>: long-shedder or standard origin. <b>Clade</b>: one "
          + "arbitrary hue per clade."}),
  ].filter(Boolean).forEach(c => box.append(c));
}

/* --------------------------------------------------------- main area */

function swChip(mode, offRules) {
  const m = swPresetMatch();
  const reset = `<button type="button" class="sw-reset" id="${SW_MODES[mode].p}_reset">Reset to defaults</button>`;
  const chip = m ? `<span class="sw-chip saved" data-hint="Saved settings">${sdEsc(m)}</span>`
    : offRules.length ? `<span class="sw-chip off" data-hint="${sdEsc("Changed from the defaults: "
        + offRules.map(k => SD_LABEL[k]).join(", ") + ".")}">changed</span>`
    : `<span class="sw-chip">defaults</span>`;
  return chip + (offRules.length ? reset : "") + swPresetBar(mode);
}
function swStartText() {
  return settings.sdStart === "frac"
    ? `the last day below ${fmt(100 * settings.sdFrac, settings.sdFrac < 0.01 ? 1 : 0)}% of the peak`
    : settings.sdStart === "trough"
    ? `the last dip ≥ ${fmt(100 * settings.sdTrough, settings.sdTrough < 0.01 ? 2 : 0)}% of peak deep`
    : "the prominence base";
}
function swRuleLine(mode, offRules) {
  const kern = {box: "moving average", gauss: "Gaussian", median: "median"}[settings.sdKernel];
  const brk = `the <b>${settings.sdPick}</b> episode peaking ≥ <b>${brkPct(+settings.brkF)}</b>`;
  if (!SW_MODES[mode].cryptic)
    return `<div class="sw-rule hintpop-on">${swChip(mode, offRules)} Episodes are peaks standing
      <b>${Number(settings.sdProm).toFixed(2)}</b> of ${promRefText()} above their
      dips${+settings.sdWidth ? `, at least <b>${settings.sdWidth} d</b> wide` : ""}, on a
      <b>${kern} ${settings.sdWin} d</b>. Each starts at <b>${swStartText()}</b>.
      The outbreak is ${brk}. <a href="#cryptic">Which are cryptic →</a></div>`;
  return `<div class="sw-rule hintpop-on">${swChip(mode, offRules)} <b>Rule ${cryUnder() ? 2 : 1}</b>: hidden from
    <b>emergence in host</b> to ${cryUnder() ? `the first day above <b>f = ${brkPct(+settings.brkLo)}</b>` : "the outbreak's start"}.
    Cryptic: hidden ≥ <b>${settings.brkH} d</b>${!cryUnder() && settings.brkQuiet
      ? ", <b>never above its launch level before it</b>" : ""}.
    <span class="sw-sub">Outbreaks from <a href="#sweep">← Outbreak detection</a>:
    ${brk.replace(/<\/?b>/g, "")}, starting at ${swStartText()}.</span></div>`;
}

function swStats(mode, prep) {
  const rows = prep.rows;
  const stat = ([k, v, s]) => `<div class="stat"><div class="k">${k}</div><div class="v">${v}`
    + (s ? ` <small>${s}</small>` : "") + `</div></div>`;
  const brks = rows.filter(r => r.brk);
  const withEp = rows.filter(r => r.main);
  const med = (a, f) => a.length ? f(sdMedian(a)) : "—";
  const cells = !SW_MODES[mode].cryptic || !prep.brkOn ? [
    ["clades", rows.length.toLocaleString(), ""],
    ["with an episode", withEp.length.toLocaleString(), ""],
    ["episodes / clade", rows.length
      ? (rows.reduce((t, r) => t + r.eps.length, 0) / rows.length).toFixed(2) : "—", ""],
    ...(prep.brkOn ? [
      ["outbreaks", brks.length.toLocaleString(), `peak ≥ ${brkPct(+settings.brkF)}`],
      ["median launch", med(brks.map(r => r.brk.launch), brkPct), ""],
      ["median outbreak rise", med(brks.map(r => r.brk.rise).filter(isNum), Math.round), "days"],
    ] : [["median clade rise", med(withEp.map(r => r.value), Math.round), "days"]]),
  ] : (() => {
    const cry = brks.filter(r => r.brk.cry), zero = brks.filter(r => r.brk.hidden <= 0);
    const gap = rows.filter(r => isNum(r.tEm) && isNum(r.tTx)).map(r => r.tTx - r.tEm);
    return [
      ["outbreaks", brks.length.toLocaleString(), `peak ≥ ${brkPct(+settings.brkF)}`],
      ["cryptic", cry.length.toLocaleString(), `of ${brks.length}`],
      ["median hidden", med(brks.map(r => r.brk.hidden), Math.round), "days"],
      ["median launch", med(brks.map(r => r.brk.launch), brkPct), ""],
      ["hidden 0", zero.length.toLocaleString(), `of ${brks.length}, start ≤ anchor`],
      ["emergence → 1st tx", med(gap, v => fmt(v, 1)), "days, median"],
    ];
  })();
  return `<div class="stats">${cells.map(stat).join("")}</div>`;
}

function swLegend(mode, showEps) {
  const ink2 = "var(--ink-2)", ink3 = "var(--ink-3)";
  const vl = dash => swKeyIcon(`<line x1="11" y1="0" x2="11" y2="12" style="stroke:${ink2}" stroke-dasharray="${dash}"/>`);
  const hl = dash => swKeyIcon(`<line x1="1" y1="6" x2="21" y2="6" style="stroke:${ink3}" stroke-dasharray="${dash}"/>`);
  const dot = col => swKeyIcon(swDot(11, 6, col, 3.5));
  const cry = SW_MODES[mode].cryptic;
  const items = [
    ...(showEps ? [[swKeyIcon(`<polyline points="${SW_PATH}" fill="none" style="stroke:var(--ink-2)" stroke-width="1" opacity=".4"/>`
                    + `<polyline points="6,10 11,2.5" fill="none" style="stroke:var(--ink-2)" stroke-width="2.6" stroke-linecap="round"/>`),
                    "episode, start to peak"]] : []),
    ...(cry ? [[swKeyIcon(`<rect x="1" y="0" width="20" height="12" style="fill:var(--long)" opacity=".18"/>`), "hidden time"],
               [vl("none"), "emergence in host"]] : []),
    [dot("var(--long)"), "outbreak starts, at its launch frequency"],
    [dot("var(--ink)"), "outbreak peak"],
    [swKeyIcon(swDot(6, 6, "var(--long)", 2.4) + swDot(16, 6, "var(--ink)", 2.4)), "start and peak of another episode reaching <i>F</i>"],
    ...(cry ? [] : [[hl("4 3"), "<i>F</i>"]]),
    ...(cry && cryCeilText() ? [[swKeyIcon(`<line x1="1" y1="6" x2="21" y2="6" style="stroke:var(--long)" stroke-width="1.4" stroke-dasharray="2 3"/>`),
               cryCeilText()]] : []),
    ...(cry ? [[`<i class="sw-crysw"></i>`, "cryptic"]] : []),
  ];
  return `<div class="bkey">${items.map(([i, t]) => `<span>${i}${t}</span>`).join("")}</div>`;
}

function sweepKey(anyOther) {
  if (settings.colourBy === "origin")
    return keyRow([[css("--long"), "long origin"], [css("--std"), "standard origin"],
                   [css("--mixed"), "mixed (tie)"]]
      .concat(anyOther ? [[css("--ink-3"), "smaller clades"]] : []));
  return anyOther ? `<div class="key"><span><i style="background:${css("--ink-3")}"></i>smaller clades</span></div>` : "";
}

/* drawBreakout on a cell's own axes: sx maps a daily sample index, so a day
   goes through dayToIndex, clamped to the record */
function sweepBreakout(g, {sx, sy, w, h}, r, V, rad, hidden, focus = null) {
  const X = d => { const i = dayToIndex(V.ft, d);
    return sx(Math.max(V.x0 || 0, i !== null ? i : d < V.ft[0] ? 0 : V.nt - 1)); };
  /* every other episode reaching F, its start and peak, smaller: the candidates
     the Which rule chose between */
  const F = +settings.brkF, j0 = V.j0 || 0, x0 = V.x0 || 0;
  const pk = r.brk && isNum(r.brk.peakDay) ? r.brk.peakDay : null;
  const dot = (j, col, rr) => { const y = r.ys[j];
    g.fillStyle = col; g.beginPath(); g.arc(sx(Math.max(x0, j0 + j)), sy(y), rr, 0, 6.2832); g.fill();
    g.strokeStyle = css("--panel"); g.lineWidth = 1; g.stroke(); };
  /* a clade left without an outbreak: the episode the Which rule would have
     picked, full size, so the reason it was dropped shows */
  const over = r.ys && r.times ? (r.eps || []).filter(e => r.ys[e.b] >= F) : [];
  const would = r.brk || !over.length ? null : settings.sdPick === "largest"
    ? over.reduce((m, e) => r.ys[e.b] > r.ys[m.b] ? e : m) : over[0];
  over.forEach(e => {
    if (pk !== null && Math.abs(r.times[e.b] - pk) < 1e-6) return;
    const rr = e === would ? rad : rad * 0.7;
    g.globalAlpha = focus ? 0.5 : 1;
    dot(e.a, css("--long"), rr); dot(e.b, css("--ink"), rr); g.globalAlpha = 1; });
  drawBreakout(g, X, sy, {l: 2, r: w - 2, t: 5, b: h - 5}, r.brk, {r: rad, hidden, startLine: false, focus, showF: !hidden});
}

/* Seeds that never reach the burn-in day (/api/runs phi null): red ticks on
   the seed slider, and a red tag beside the seed when it is one of them */
async function swSeedTags(mode) {
  const ctl = document.querySelector(`#${mode}_rail .sw-seed`);
  if (!ctl) return;
  let runs;
  try { runs = await apiGet("/api/runs"); } catch (e) { return; }
  const list = (runs.scenarios || {})[settings.scenario] || [];
  const bad = list.map((r, i) => r.phi == null ? i : -1).filter(i => i >= 0);
  const input = ctl.querySelector("input[type=range]");
  const n = input ? +input.max : 49;
  let ticks = ctl.querySelector(".sw-ticks");
  if (!ticks) { ticks = document.createElement("div"); ticks.className = "sw-ticks"; input.after(ticks); }
  ticks.innerHTML = bad.map(i => `<i style="left:calc(7px + (100% - 14px) * ${i / Math.max(1, n)})"></i>`).join("");
  ticks.title = bad.length ? `never reach the burn-in day: ${bad.map(i => String(i).padStart(4, "0")).join(", ")}` : "";
  const lab = ctl.querySelector(".lab");
  let tag = lab.querySelector(".sw-nophi");
  if (bad.includes(+settings.seed)) {
    if (!tag) { tag = document.createElement("span"); tag.className = "sw-nophi"; tag.textContent = "no burn-in"; lab.append(tag); }
  } else if (tag) tag.remove();
}

/* ---------------------------------------------------------- presets
   The shared data-settings sets (presets.js), shown in this tab's rule line */
const swPresetsLoad = presetsLoad;
const swPresetMatch = presetMatch;
function swApply(patch, mode) { presetApply(patch); }
const swPresetBar = mode => presetBarHtml(SW_MODES[mode].p);
const swPresetWire = mode => presetWire(SW_MODES[mode].p, () => renderTab(tabById(mode)));

/* drawBurnin's short dashed line, with a label where the cell has room */
function swBurninMark(g, x, w, h) {
  drawBurnin(g, x, w, h, null, true);
  g.save(); g.font = "8px " + css("--f-num"); g.fillStyle = css("--ink-3");
  const t = "burn-in", tw = g.measureText(t).width;
  if (x + 3 + tw < w - 2) g.fillText(t, x + 3, 9);
  else if (x - 3 - tw > 2) { g.textAlign = "right"; g.fillText(t, x - 3, 9); }
  g.restore();
}

/* Why a cell is greyed out: it has no outbreak under the rules in force */
function swWhyGrey(r, V) {
  const F = +settings.brkF;
  if (isNum(r.preBurn)) return `<b>Left out: its outbreak began before the burn-in day.</b> `
    + `It starts on day ${Math.round(r.preBurn)}; the burn-in day is day ${Math.round(V.phi)}, `
    + `and under <b>Outbreak starts after it</b> only later outbreaks count.`;
  if (r.trunc) return `<b>No outbreak: already circulating at the burn-in day</b>`
    + `${isNum(V.phi) ? ` (day ${Math.round(V.phi)})` : ""}. It emerged before the `
    + `measured window, so there is no emergence to time hidden time from, and Figure 5 leaves it out.`;
  if (r.eps.length) {
    const top = Math.max(...r.eps.map(e => r.ys[e.b]));
    return `<b>No outbreak: no episode peaks at F or above.</b> Its tallest `
      + `episode peaks at ${brkPct(top)}, under <i>F</i> = ${brkPct(F)}. Lower `
      + `<i>F</i> on <b>4 · The outbreak</b> to count it.`;
  }
  let iMax = 0;
  r.ys.forEach((v, j) => { if (v > r.ys[iMax]) iMax = j; });
  const atEnd = iMax >= r.ys.length - 2 && !settings.sdCens;
  return `<b>No outbreak: no episode.</b> No peak of its smoothed curve stands `
    + `${Number(settings.sdProm).toFixed(2)} of ${promRefText()} above its dips`
    + `${settings.sdWidth ? ` and is ${settings.sdWidth} d wide` : ""}.`
    + (atEnd ? " Its maximum is at the end of the record and <b>Peak still high at "
             + "the end</b> is unticked." : "")
    + ` The plot above stacks only clades with an episode, so it leaves this one `
    + `out; click the cell to draw it there.`;
}

/* Cells in rank order, each rank by peak; the Show filter cuts from the end.
   Detection ranks outbreak, rest; cryptic ranks cryptic, outbreak, rest. */
function sweepCells(mode, rows, brkOn) {
  const cry = SW_MODES[mode].cryptic, show = settings[SW_MODES[mode].show];
  const rank = r => r.brk ? (cry && r.brk.cry ? 0 : 1) : 2;
  const keep = !brkOn || show === "all" ? 2 : show === "breakout" ? 1 : 0;
  return rows.filter(r => rank(r) <= keep)
    .sort((a, b) => rank(a) - rank(b) || b.peak - a.peak);
}

/* Per simulation, per scenario: outbreaks and cryptic outbreaks under the tabs'
   rules, from /api/explore_panel with breakoutOf -- the pick the cells use. A
   simulation with no burn-in day is left out while the burn-in is in force. */
function swScenCounts(d) {
  const out = {};
  SCENARIOS.forEach(sc => {
    const runs = ((d.scenarios || {})[sc] || {}).runs || [];
    const brk = [], cry = [];
    let skipped = 0;
    runs.forEach(r => {
      if (settings.burnin && r.phi == null) { skipped++; return; }
      const gate = sdGateDay(settings.burnin, r.phi);
      let nb = 0, nc = 0;
      r.clades.forEach(([, , tem, eps, , , , tCross]) => {
        const b = breakoutOf(eps.map(([aDay, top, amp, days, pre]) => ({aDay, top, amp, days, pre})),
                             tem, settings.sdPick, gate, tCross);
        if (b) { nb++; if (b.cry) nc++; }
      });
      brk.push(nb); cry.push(nc);
    });
    out[sc] = {brk, cry, skipped};
  });
  return out;
}
const swMean = a => a.length ? a.reduce((s, x) => s + x, 0) / a.length : NaN;
const swSd = a => { if (a.length < 2) return NaN; const m = swMean(a);
  return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1)); };
const swP = p => !isNum(p) ? "—" : p < 0.001 ? "&lt; 0.001" : p.toFixed(3);

/* Run fn with `patch` in force over `settings`, without notifying anyone.
   fn must be synchronous: nothing else may see the swapped values. */
function swWithSettings(patch, fn) {
  if (!patch) return fn();
  const old = {};
  Object.keys(patch).forEach(k => { old[k] = settings[k]; settings[k] = patch[k]; });
  try { return fn(); } finally { Object.assign(settings, old); }
}

/* One measurement: counts per simulation and the tests, under the settings in
   force, or under `patch` over them (a saved set). Only the live one takes the
   slot, so the saved one never aborts it. */
async function swCryMeasure(patch) {
  const url = "/api/explore_panel?" + swWithSettings(patch, () => exploreQuery({mf: 0}));
  const d = patch ? await apiGet(url, {label: "measuring under the saved settings"})
                  : await apiGetSlot(url);
  const C = swWithSettings(patch, () => swScenCounts(d));
  const scs = SCENARIOS.filter(sc => C[sc].brk.length);
  const metrics = {cryptic: Object.fromEntries(scs.map(sc => [sc, C[sc].cry])),
                   breakouts: Object.fromEntries(scs.map(sc => [sc, C[sc].brk]))};
  let T = null, terr = null;
  try {
    const r = await fetch("/api/count_tests", {method: "POST", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({metrics})});
    T = await r.json();
    if (!r.ok) { terr = T.error || `HTTP ${r.status}`; T = null; }
  } catch (e) { terr = e.message; }
  return {C, scs, T, terr};
}

/* what a matrix cell says: not significant, or significant with the row higher/lower */
function swVerdict(x, mA, mB) {
  if (!x || !isNum(x.p_holm)) return null;
  if (x.p_holm >= 0.05) return "ns";
  return mA > mB ? "hi" : mA < mB ? "lo" : "eq";
}
const SW_ARROW = {hi: "▲", lo: "▼", eq: "", ns: ""};
const swSigned = (v, d) => (v > 0 ? "+" : v < 0 ? "−" : "±") + fmt(Math.abs(v), d);

/* How a pair's verdict moved from the saved set: its direction reversed, it
   became significant, or it stopped being. Direction is the sign of the mean
   difference, and counts only where the pair is significant on either side. */
function swChange(v, rv, d, rd) {
  if (!v || !rv) return null;
  const sig = v !== "ns", rsig = rv !== "ns";
  if ((sig || rsig) && d * rd < 0) return "flip";
  if (sig && !rsig) return "gain";
  if (!sig && rsig) return "loss";
  return null;
}
const SW_CHG = {flip: ["⇅ flip", "direction reversed"], gain: ["new", "became significant"],
                loss: ["lost", "no longer significant"]};
let swFadeSame = true;

let swCryTok = 0;
async function renderCryScen(main, mode = "crystats") {
  const tok = ++swCryTok;
  await swPresetsLoad();
  const head = swRuleLine(mode, sweepOffRules());
  const wire = () => { swPresetWire(mode);
    const reset = document.getElementById(`${SW_MODES[mode].p}_reset`);
    if (reset) reset.onclick = () => swApply(sweepFig5Patch(), mode); };
  main.innerHTML = head + `<div class="sw-sub">Measuring every simulation under these rules…</div>`;
  wire();
  /* off a saved set, the same measurement under it, to mark what moved */
  const refName = presetRef(), live = refName && refName !== presetMatch();
  const [cur, ref] = await Promise.all([swCryMeasure(null),
    live ? swCryMeasure(presetOf(refName)).catch(() => null) : null]);
  if (tok !== swCryTok) return;
  const {C, scs, T, terr} = cur;
  const R = ref && ref.T ? ref : null;
  const key2 = key => key === "cryptic" ? "cry" : "brk";
  const ms = a => `${fmt(swMean(a), 2)} <small>± ${fmt(swSd(a), 2)}</small>`;
  const dMean = (sc, k) => { if (!R || !R.C[sc]) return "";
    const v = swMean(C[sc][k]) - swMean(R.C[sc][k]);
    return Math.abs(v) < 0.005 ? "" : ` <span class="sw-d" title="change from the saved ${sdEsc(refName)}">${swSigned(v, 2)}</span>`; };
  const share = c => { const nb = c.brk.reduce((s, x) => s + x, 0), nc = c.cry.reduce((s, x) => s + x, 0);
    return {nb, nc, f: nb ? 100 * nc / nb : NaN}; };
  const rows = scs.map(sc => { const c = C[sc], sh = share(c);
    const rs = R && R.C[sc] ? share(R.C[sc]) : null;
    const dsh = rs && isNum(sh.f) && isNum(rs.f) && Math.abs(sh.f - rs.f) >= 0.05
      ? ` <span class="sw-d" title="change from the saved ${sdEsc(refName)}, percentage points">${swSigned(sh.f - rs.f, 1)} pp</span>` : "";
    const dn = R && R.C[sc] && R.C[sc].brk.length !== c.brk.length
      ? ` <span class="sw-d">was ${R.C[sc].brk.length}</span>` : "";
    return `<tr><th>${SCENARIO_LABEL[sc]}</th><td>${c.brk.length}${c.skipped
        ? ` <span class="sw-nophi" title="left out: no burn-in day">−${c.skipped}</span>` : ""}${dn}</td>
      <td>${ms(c.brk)}${dMean(sc, "brk")}</td><td>${ms(c.cry)}${dMean(sc, "cry")}</td>
      <td>${isNum(sh.f) ? fmt(sh.f, 1) + "%" : "—"}${dsh} <small>${sh.nc}/${sh.nb}</small></td></tr>`; }).join("");
  const nChg = {flip: 0, gain: 0, loss: 0};
  const matrix = (key, title) => {
    const t = T && T[key];
    if (!t) return "";
    const rt = R && R.T[key];
    const index = pairs => { const P = {};
      (pairs || []).forEach(x => { P[x.a + "|" + x.b] = P[x.b + "|" + x.a] = x; }); return P; };
    const P = index(t.pairs), RP = index(rt && rt.pairs);
    const mean = (M, sc) => M.C[sc] ? swMean(M.C[sc][key2(key)]) : NaN;
    const cell = (a, b) => { if (a === b) return `<td class="dg"></td>`;
      const x = P[a + "|" + b]; if (!x) return "<td>—</td>";
      const v = swVerdict(x, mean(cur, a), mean(cur, b));
      const sig = v !== "ns", cls = sig ? "sig " + (v === "hi" || v === "lo" ? v : "") : "";
      let tip = `U = ${x.U}, p = ${x.p.toPrecision(3)} unadjusted`, was = "", tag = "", ch = null;
      const rx = RP[a + "|" + b];
      if (rt && rx) {
        const rv = swVerdict(rx, mean(R, a), mean(R, b));
        ch = swChange(v, rv, mean(cur, a) - mean(cur, b), mean(R, a) - mean(R, b));
        if (ch && a < b) nChg[ch]++;
        tip += `; saved ${refName}: p = ${swP(rx.p_holm).replace("&lt;", "<")} Holm-adjusted`;
        if (ch) {
          tip += ` -- ${SW_CHG[ch][1]}`;
          tag = `<span class="sw-tag">${SW_CHG[ch][0]}</span>`;
          was = `<span class="sw-was">was ${swP(rx.p_holm)}${rv !== "ns" ? " " + SW_ARROW[rv] : ""}</span>`;
        }
      }
      const mark = ch ? ` chg-${ch}` : rt && swFadeSame ? " same" : "";
      return `<td class="${cls}${mark}" title="${sdEsc(tip)}">${tag}${swP(x.p_holm)} ${sig ? SW_ARROW[v] : ""}${was}</td>`; };
    const kw = t.kruskal, rkw = rt && rt.kruskal;
    const kwP = k => k.p < 0.001 ? "&lt; 0.001" : "= " + swP(k.p);
    return `<div class="sw-mx"><h4>${title}</h4>
      <div class="sw-sub">Kruskal–Wallis across ${kw ? kw.groups : 0} scenarios: H = ${kw ? fmt(kw.H, 2) : "—"},
        p ${kw ? kwP(kw) : "= —"}${rkw ? ` <span class="sw-d">(saved: H = ${fmt(rkw.H, 2)}, p ${kwP(rkw)})</span>` : ""}</div>
      <table class="sw-tab"><tr><th></th>${scs.map(sc => `<th>${SCENARIO_LABEL[sc]}</th>`).join("")}</tr>
      ${scs.map(a => `<tr><th>${SCENARIO_LABEL[a]}</th>${scs.map(b => cell(a, b)).join("")}</tr>`).join("")}</table></div>`;
  };
  const mxs = matrix("cryptic", "Cryptic outbreaks per simulation") + matrix("breakouts", "Outbreaks per simulation");
  const cmp = R ? `<div class="sw-cmp">Compared with the saved <b>${sdEsc(refName)}</b>
      (significant: Holm-adjusted p &lt; 0.05):
      ${Object.keys(SW_CHG).map(k => `<span class="sw-key chg-${k}"><span class="sw-tag">${SW_CHG[k][0]}</span></span>
        ${SW_CHG[k][1]} <b>${nChg[k]}</b>`).join(" · ")}
      <label class="sw-fade"><input type="checkbox" id="cs_fade"${swFadeSame ? " checked" : ""}> fade unchanged</label>
      <div class="sw-sub">Counts are pairs; each pair appears twice in its matrix. A marked cell shows its
        saved p under it. Other p shifts are on hover only.</div></div>`
    : live ? `<div class="sw-cmp">Could not measure under the saved <b>${sdEsc(refName)}</b>, so nothing is marked.</div>` : "";
  main.innerHTML = head + `
    <div class="plot"><div class="ph"><h3>Outbreaks by scenario</h3>
      <span class="sub">per simulation · k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></div>
      <table class="sw-tab sw-sum"><tr><th>Scenario</th><th>Simulations</th><th>Outbreaks / simulation</th>
        <th>Cryptic / simulation</th><th>Cryptic share</th></tr>${rows}</table>
      <div class="sw-sub">Mean ± standard deviation over simulations. Cryptic share pools every
        outbreak in the scenario.${settings.burnin ? " Simulations that never reach the burn-in day are left out (red count)." : ""}</div></div>
    <div class="plot"><div class="ph"><h3>Do scenarios differ?</h3>
      <span class="sub">two-sided Mann–Whitney U per pair, Holm-adjusted over ${scs.length * (scs.length - 1) / 2} pairs</span></div>
      ${terr ? `<div class="note err">The tests did not run: ${sdEsc(terr)}</div>` : ""}
      ${cmp}
      <div class="sw-mxs">${mxs}</div>
      <div class="sw-sub">Each cell: the Holm-adjusted p for the row scenario against the column
        scenario; where p &lt; 0.05 it is red with ▲ if the row's mean is higher, green with ▼ if lower.
        Hover a cell for U and the unadjusted p. The simulation is the unit: one count per run.</div></div>`;
  wire();
  const fade = document.getElementById("cs_fade");
  if (fade) fade.onchange = ev => { swFadeSame = ev.target.checked;
    main.querySelectorAll(".sw-mxs td.sig, .sw-mxs td").forEach(td => {
      if (td.matches(".dg") || /chg-/.test(td.className) || !td.title.includes("saved")) return;
      td.classList.toggle("same", swFadeSame); }); };
}

/* ------------------------------------------------------- overview
   Every outbreak of the scenario, one cell each, drawn as the one-simulation
   cells are, so a glance down the grid shows how the rules behave. Runs load
   one at a time and their cells are appended as each arrives. */
const SW_WHY = [
  ["all", "All", "Every clade without an outbreak."],
  ["below", "below F", "Has episodes, but none peaks at F or above. Black dot: its tallest peak."],
  ["pre", "before burn-in", "Its outbreak started on or before the burn-in day, so it is left out."],
  ["at", "at burn-in", "Already circulating at the burn-in day: no emergence to time from."],
  ["nopeak", "no peak", "No peak passes the episode rules."]];
const CR_WHY = [
  ["all", "All", "Every outbreak that is not cryptic."],
  ["short", "hidden < H", "Hidden for less than H days."],
  ["rose", "rose while hidden", "Rule 1 with the launch rule: hidden long enough, but the curve rose above its launch level while hidden."]];
/* the scenario view's settings, per tab */
const SW_OV = {
  sweep: {on: "sdOverview", show: "sdOvShow", why: "sdOvWhy", eps: "sdEpsDet", whys: SW_WHY,
          shows: [["breakout", "Outbreaks", "Every clade the rules picked an outbreak for."],
                  ["none", "No outbreak", "Every other clade, each saying why."]], filtered: "none"},
  cryptic: {on: "crOverview", show: "crOvShow", why: "crOvWhy", eps: "sdEps", whys: CR_WHY,
          shows: [["cryptic", "Cryptic", "Only the cryptic outbreaks."],
                  ["breakout", "Outbreaks", "Every outbreak, cryptic ones tinted."],
                  ["notcry", "Not cryptic", "Outbreaks that are not cryptic, each saying why."]], filtered: "notcry"},
};
const swOvOn = mode => SW_OV[mode] && !!settings[SW_OV[mode].on];
/* the reasons that can occur under the rules in force: rule 2 has no launch rule */
const swOvWhys = O => O.whys.filter(([v]) => !(v === "rose" && cryUnder()));

async function swOverview(mode, main) {
  const O = SW_OV[mode], P = SW_MODES[mode].p, cry = SW_MODES[mode].cryptic;
  const sc = settings.scenario, F = +settings.brkF, H = +settings.brkH;
  const show = settings[O.show], showEps = !!settings[O.eps];
  const filtered = show === O.filtered;
  const whyKey = cry
    ? r => r.brk.hidden < H ? "short" : "rose"
    : r => isNum(r.preBurn) ? "pre" : r.trunc ? "at" : r.eps.length ? "below" : "nopeak";
  const whyText = Object.fromEntries(O.whys.map(([v, t]) => [v, t]));
  const why = r => whyText[whyKey(r)];
  const only = filtered ? settings[O.why] : "all";
  const keep = cry
    ? r => r.brk && (show === "breakout" || (show === "cryptic") === r.brk.cry)
    : r => !r.brk === filtered;
  const cols = +settings.sdOvCols || 6;
  const cropped = !!settings.burnin && settings.sdGate === "trunc";
  const what = cry ? {cryptic: "Every cryptic outbreak", breakout: "Every outbreak", notcry: "Every outbreak that is not cryptic"}[show]
    : filtered ? "Every clade without an outbreak" : "Every outbreak";
  main.innerHTML = swRuleLine(mode, sweepOffRules()) + `<div id="${P}_ovStats"></div>
    <div class="plot"><div class="ph"><h3>${what} in ${SCENARIO_LABEL[sc]}</h3>
      <span class="sub" id="${P}_ovSub">k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></div>
      <div class="sw-tools hintpop-on" id="${P}_tools"></div>
      ${filtered || (cry && show === "breakout") ? `<div class="whykey">${
        cry && show === "breakout" ? `<span><i class="sw-crysw"></i>cryptic</span>` : ""}${
        swOvWhys(O).filter(([v]) => v !== "all").map(([v, t, h]) =>
        `<span title="${sdEsc(h)}"><i class="why-${v}"></i>${cry && !filtered ? "not cryptic: " + t : t}</span>`).join("")}</div>` : ""}
      <div class="rgrid hintpop-on" id="${P}_ovGrid"></div>
      ${swLegend(mode, showEps)}
      <div class="cap">One cell per ${what.replace(/^Every /, "")}, simulation by simulation${
        filtered && !cry ? ", simulations that never reach the burn-in day left out" : ""}: the clade's smoothed
        frequency from its emergence to the end of the record,
        ${settings.yScale === "own" ? "each to its own peak" : "0–100% in every cell"}. Header: simulation,
        the clade's size in lineages, ${filtered ? "why, also the outline's colour" : cry ? "its hidden time" : "the outbreak's rise"};
        hover for more. Click a cell to open its simulation with the clade selected.</div></div>`;
  swTools(mode, document.getElementById(`${P}_tools`), true);
  swPresetWire(mode);
  const reset = document.getElementById(`${P}_reset`);
  if (reset) reset.onclick = () => swApply(sweepFig5Patch(), mode);
  const runs = (await apiGet("/api/runs", {slot: "runs", label: "simulation index"})).scenarios[sc] || [];
  const host = document.getElementById(`${P}_ovGrid`);
  if (!host) return;
  host.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
  const cellH = cols >= 8 ? 56 : cols >= 6 ? 66 : 80;
  const picks = [];
  host.onclick = ev => { const c = ev.target.closest && ev.target.closest(".rcell");
    const x = c ? picks[+c.dataset.i] : null;
    if (!x) return;
    /* the rail's seed slider is built once: rebuild it to show the seed opened */
    DATA_TABS.forEach(id => { const t = tabById(id); if (t) t.railBuilt = false; });
    setSettings({seed: x.seed, clade: x.id, [O.on]: 0}); };
  const brks = [], acc = {done: 0, nophi: 0, clades: 0, failed: 0},
        whyN = Object.fromEntries(swOvWhys(O).filter(([v]) => v !== "all").map(([, t]) => [t, 0]));
  const stats = () => {
    const el = document.getElementById(`${P}_ovStats`);
    if (!el) return;
    const med = (a, f) => a.length ? f(sdMedian(a)) : "—";
    const used = acc.done - acc.nophi - acc.failed;
    const nCry = brks.filter(b => b.cry).length;
    el.innerHTML = `<div class="stats">${[
      ["simulations", `${acc.done}${acc.done < runs.length ? ` of ${runs.length}` : ""}`,
        [acc.nophi ? `${acc.nophi} never reach burn-in` : "", acc.failed ? `${acc.failed} unread` : ""].filter(Boolean).join(", ")],
      ["outbreaks", brks.length.toLocaleString(), `peak ≥ ${brkPct(F)}`],
      ...(cry ? [
      ["cryptic", nCry.toLocaleString(), brks.length ? `${Math.round(100 * nCry / brks.length)}% of outbreaks` : ""],
      ["cryptic / simulation", used ? (nCry / used).toFixed(2) : "—", ""],
      ["median hidden", med(brks.map(b => b.hidden), Math.round), "days"],
      ...Object.entries(whyN).map(([k, v]) => [k, v.toLocaleString(), "not cryptic"]),
    ] : filtered ? [
      ["no outbreak", Object.values(whyN).reduce((t, v) => t + v, 0).toLocaleString(), `of ${acc.clades} clades`],
      ...Object.entries(whyN).map(([k, v]) => [k, v.toLocaleString(), ""]),
    ] : [
      ["outbreaks / simulation", used ? (brks.length / used).toFixed(2) : "—", ""],
      ["outbreaks / clade", acc.clades ? (brks.length / acc.clades).toFixed(2) : "—", `of ${acc.clades} clades`],
      ["median launch", med(brks.map(b => b.launch), brkPct), ""],
      ["median outbreak rise", med(brks.map(b => b.rise).filter(isNum), Math.round), "days"],
    ])].map(([k, v, s]) => `<div class="stat"><div class="k">${k}</div><div class="v">${v}`
      + (s ? ` <small>${s}</small>` : "") + `</div></div>`).join("")}</div>`;
  };
  stats();
  for (let i = 0; i < runs.length; i++) {
    if (!document.body.contains(host)) return;   // re-rendered under us
    let prep;
    try {
      const [v, d] = await Promise.all([runView(sc, i, settings.k, settings.mt),
        apiGet("/api/explore?" + exploreQuery({scenario: sc, seed: i, mf: 0}))]);
      prep = prepareEpisodes(v, d);
    } catch (e) {
      if (isSuperseded(e)) return;
      acc.failed++; acc.done++; stats(); continue;
    }
    if (!document.body.contains(host)) return;
    const nophi = !prep.brkOn || prep.fig5Off === "nophi";
    acc.nophi += nophi; acc.clades += nophi ? 0 : prep.rows.length; acc.done++;
    const rows = nophi ? [] : prep.rows.filter(r => keep(r) && (only === "all" || whyKey(r) === only))
      .sort((a, b) => b.peak - a.peak);
    if (!nophi) prep.rows.forEach(r => {
      if (cry ? r.brk && !r.brk.cry : !r.brk) whyN[why(r)]++; });
    brks.push(...prep.rows.filter(r => r.brk && !nophi).map(r => r.brk));
    const V0 = prep.view;
    const V = Object.assign({}, V0, {top: 1, x0: cropped && V0.phiIdx != null ? V0.phiIdx : 0});
    const frag = document.createElement("div");
    frag.innerHTML = rows.map(r => {
      const k = picks.push({seed: i, id: r.id}) - 1, b = r.brk, lab = labelOf(r.inf);
      const out = cry ? !b.cry : !b;
      const hint = `${sdEsc(runs[i].seed)} · ${lab}` + (!cry && !b ? `<br>${swWhyGrey(r, prep.view)}`
        : cry ? `<br>hidden ${Math.round(b.hidden)} d: ${b.cry ? "cryptic" : "not cryptic, " + why(r)}` : "");
      const corner = cry ? `hidden ${Math.round(b.hidden)} d`
        : b ? `rise ${isNum(b.rise) ? Math.round(b.rise) : "—"} d` : why(r);
      return `<div class="rcell${out ? " why-" + whyKey(r) : cry ? " cry" : ""}" data-i="${k}" data-hint="${sdEsc(hint)}">
        <div class="rh"><b>${String(i).padStart(4, "0")} · ${r.inf.size}L</b>
          <span class="d">${corner}</span></div>
        <canvas></canvas>
        <div class="rf"><span>${b ? `launch ${brkPct(b.launch)} · peak ${(100 * b.peak).toFixed(0)}%`
          : `peak ${(100 * r.ys.reduce((m, v) => v > m ? v : m, 0)).toFixed(0)}%`}</span><span>${r.eps.length} ep${r.eps.length === 1 ? "" : "s"}</span></div></div>`;
    }).join("");
    const cells = [...frag.children];
    cells.forEach(c => host.append(c));
    rows.forEach((r, j) => {
      const cv = cells[j].querySelector("canvas");
      const Vr = Object.assign({}, V, {x0: Math.max(V.x0, V.j0 + emergedFrom(r.times, r.tEm))});
      const p = paintOne(cv, r, cellH, false, Vr, showEps, null);
      if (V.phiIdx != null && Vr.x0 < V.phiIdx) swBurninMark(cv.getContext("2d"), p.sx(V.phiIdx), p.w, p.h);
      const g = cv.getContext("2d");
      sweepBreakout(g, p, r, Vr, 2.8, cry);
      /* below F: its tallest peak, to see how far under F it stayed */
      if (!cry && !r.brk && whyKey(r) === "below" && r.eps.length) {
        const e = r.eps.reduce((m, x) => r.ys[x.b] > r.ys[m.b] ? x : m);
        g.fillStyle = css("--ink"); g.beginPath();
        g.arc(p.sx(Math.max(Vr.x0, Vr.j0 + e.b)), p.sy(r.ys[e.b]), 2.8, 0, 6.2832); g.fill();
        g.strokeStyle = css("--panel"); g.lineWidth = 1; g.stroke();
      }
    });
    stats();
  }
}

async function renderSweepMode(mode, main) {
  const M = SW_MODES[mode], P = M.p, cry = M.cryptic;
  main.classList.add("fit");
  let dailyErr = null;
  const q = exploreQuery({scenario: settings.scenario, seed: settings.seed, mf: 0});
  await swPresetsLoad();
  if (swOvOn(mode)) return swOverview(mode, main);
  const [runv, daily] = await Promise.all([
    runView(settings.scenario, settings.seed, settings.k, settings.mt),
    apiGetSlot("/api/explore?" + q)
      .then(d => d, e => { if (isSuperseded(e)) throw e; dailyErr = e.message; return null; }),
  ]);
  const prep = prepareEpisodes(runv, daily);
  /* measured from the burn-in day, the cells start there */
  /* Figure 5's burn-in rule measures from the burn-in day, so the cells start
     there; the "start" rule measures the whole record */
  const cropped = !!settings.burnin && settings.sdGate === "trunc" && prep.view.phiIdx != null;
  const V = cropped ? Object.assign({}, prep.view, {x0: prep.view.phiIdx}) : prep.view;
  const offRules = sweepOffRules();
  const rows = prep.rows;
  const cells = sweepCells(mode, rows, prep.brkOn);
  const shown = cells.slice(0, SWEEP_MAX_CELLS);
  const gridMode = settings.episodeView === "grid";
  const showEps = !!settings[M.eps];
  const focus = swFocus[mode];
  const rail = document.getElementById(`${mode}_rail`);
  if (rail) swMarkFocus(rail, mode);
  const show = settings[M.show];

  const offNote = prep.fig5Off === "nophi"
    ? `<div class="note nophi"><span class="sw-nophi">no burn-in day</span> <strong>No outbreaks
       in this simulation.</strong> (diagnosed + recovered) / N never reaches 1, so it has no
       burn-in day and Figure 5 takes nothing from it. Set <em>Burn-in</em> to <em>Ignore it</em>
       to see its outbreaks.</div>`
    : prep.fig5Off ? `<div class="note"><strong>Not Figure 5's numbers.</strong>
       Measured ${settings.burnin ? "from the burn-in day" : "over the whole record"},
       every clade counted. Set <em>Burn-in</em> in Simulation.</div>` : "";
  const errNote = prep.daily ? "" : `<div class="note err"><strong>No daily record.</strong>
       The detector could not read this run (${sdEsc(dailyErr || "no answer")}); the
       cells fall back to the 9.2-day export and no outbreak is drawn.</div>`;
  const unserved = prep.daily && prep.unserved ? `<div class="sw-sub">${prep.unserved}
       clade${prep.unserved === 1 ? "" : "s"} carry no daily series and are not shown.</div>` : "";
  const what = !prep.brkOn || show === "all" ? `${rows.length} clades`
    : show === "cryptic" ? `${cells.length} cryptic of ${rows.length} clades`
    : `${cells.length} outbreaks of ${rows.length} clades`;
  const empty = show === "cryptic" ? "cryptic outbreak" : "outbreak";

  main.innerHTML = errNote + offNote + swRuleLine(mode, offRules) + swStats(mode, prep) + `
    <div class="plot"><div class="ph"><h3>Lineage frequency by clade</h3>
      <span class="sub" id="${P}_mullerSub"></span></div>
      <canvas id="${P}_muller"></canvas>
      <div class="axis"><span id="${P}_mX0"></span><span>${gridMode
        ? "each cell below starts at its clade's emergence"
        : cropped ? `the plot below starts at the burn-in, day ${Math.round(V.phi)}`
        : "the plot below shares this axis"}</span>
        <span id="${P}_mX1"></span></div>
      <div id="${P}_mullerKey"></div></div>
    <div class="plot"><div class="ph">
      <h3>${gridMode ? (cry ? "Each clade, its hidden time marked" : "Each clade, its outbreak marked")
                     : "Every clade on one axis"}</h3>
      <span class="sub">${shown.length < cells.length ? `${shown.length} shown · ` : ""}${what}
        · k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></div>
      <div class="sw-tools hintpop-on" id="${P}_tools"></div>
      ${unserved}
      ${gridMode ? (shown.length ? `<div class="rgrid hintpop-on${focus ? " sw-f-" + focus : ""}" id="${P}_grid"></div>`
                    : `<div class="sw-empty" id="${P}_grid">No ${empty} in this simulation under these rules.</div>`)
                 : `<canvas id="${P}_one"></canvas>`}
      ${prep.brkOn ? swLegend(mode, showEps) : ""}
    </div>`;

  swTools(mode, document.getElementById(`${P}_tools`), prep.brkOn);
  swPresetWire(mode);
  swSeedTags(mode);
  const reset = document.getElementById(`${P}_reset`);
  if (reset) reset.onclick = () => swApply(sweepFig5Patch(), mode);

  if (gridMode && shown.length) {
    const grid = document.getElementById(`${P}_grid`);
    grid.style.gridTemplateColumns = `repeat(${settings.episodeCols}, minmax(0, 1fr))`;
    grid.innerHTML = shown.map((r, i) => {
      const lab = labelOf(r.inf);
      const b = r.brk;
      const none = isNum(r.preBurn) ? "before burn-in" : r.trunc ? "at burn-in"
        : r.eps.length ? "no outbreak" : "no peak";
      const corner = !prep.brkOn ? (r.main ? `rise ${Math.round(r.value)} d` : "no peak")
        : !b ? none
        : cry ? `hidden ${Math.round(b.hidden)} d${b.cry ? " · <b>cryptic</b>" : ""}`
        : `outbreak · rise ${isNum(b.rise) ? Math.round(b.rise) : "—"} d`;
      const tip = !b ? "" : cry ? "from emergence to the outbreak's start"
        : "the outbreak episode's rise days";
      const foot = b
        ? `launch ${brkPct(b.launch)} · peak ${(100 * b.peak).toFixed(1)}%`
          + (cry && showEps && isNum(b.rise) ? ` · rise ${Math.round(b.rise)} d` : "")
        : `peak ${(100 * r.peak).toFixed(1)}%` + (r.trunc ? " · circulating at the burn-in day"
          : !r.eps.length && r.nFound ? ` · under m_F = ${mtPct(settings.mf)}%` : "");
      const from = isNum(r.tEm) ? `from day ${Math.round(Math.max(r.tEm, cropped ? V.phi : -Infinity))}` : "";
      const tail = [from, showEps || !cry ? `${r.eps.length} ep${r.eps.length === 1 ? "" : "s"}` : ""]
        .filter(Boolean).join(" · ");
      const why = b || !prep.brkOn ? "" : swWhyGrey(r, prep.view);
      return `<div class="rcell${b || !prep.brkOn ? "" : " nobrk"}${cry && b && b.cry ? " cry" : ""}${
        r.id === settings.clade ? " sel" : ""}" data-i="${i}"${why ? ` data-hint="${sdEsc(why)}"` : ""}>
        <div class="rh"><b title="${lab}">${lab} · ${r.inf.size}L</b>
          <span class="d" title="${tip}">${corner}</span></div>
        <canvas></canvas>
        <div class="rf"><span>${foot}</span><span>${tail}</span></div></div>`;
    }).join("");
  }
  const mcv = document.getElementById(`${P}_muller`);
  document.getElementById(`${P}_mullerKey`).innerHTML = sweepKey(false);
  const cellCvs = gridMode && shown.length
    ? [...document.getElementById(`${P}_grid`).querySelectorAll(".rcell canvas")] : [];
  const oneCv = gridMode ? null : document.getElementById(`${P}_one`);
  const free = fitFree(main, [mcv, oneCv, ...cellCvs]);
  const cols = settings.episodeCols, nRows = Math.ceil(shown.length / cols);
  const mulH = clampH(free * 0.28, 70, 150);
  const cellH = clampH((free - mulH) / Math.max(1, nRows),
    cols >= 3 ? 50 : 80, cols >= 5 ? 110 : cols >= 3 ? 150 : 220);
  const oneH = clampH(free - mulH, 160, 520);
  const m = paintMuller(mcv,
    {D: runv.D, C: runv.C, ordered: prep.ordered, colMap: prep.colMap,
     height: mulH, washBurnin: !!settings.burnin,
     selected: settings.clade, hideExcluded: true,
     excluded: new Set(prep.ordered.filter(id => !rows.some(r => r.id === id && r.eps.length)))});
  /* Click a band to pop it out of the stack; click it again, or outside, to clear */
  if (!cry) sweepMuller = m;
  mcv.style.cursor = "pointer";
  mcv.onclick = ev => {
    const id = m.hitEvent(ev);
    setSetting("clade", id !== null && id === settings.clade ? null : id);
  };
  const selRow = rows.find(r => r.id === settings.clade);
  const selInf = settings.clade === null ? null : runv.C.info.get(settings.clade);
  const selNote =
    settings.clade === null ? ""
    : selRow ? ` · showing ${labelOf(selRow.inf)} · ${selRow.inf.size}L`
      + (shown.includes(selRow) ? "" : " (hidden by the Cells filter)")
    : !selInf ? " · the selected clade is not in this partition"
    : selInf.isRoot ? " · the founding clade, which carries no metric"
    : !selInf.ndef ? " · no defining mutations, so it carries no metric"
    : " · no daily series for this clade";
  document.getElementById(`${P}_mullerSub`).textContent = "click a band to find its cell" + selNote;
  document.getElementById(`${P}_mullerKey`).innerHTML = sweepKey(m.anyOther);
  document.getElementById(`${P}_mX0`).textContent = `day ${m.day0}`;
  document.getElementById(`${P}_mX1`).textContent = `day ${m.day1}`;

  if (gridMode && shown.length) {
    const grid = document.getElementById(`${P}_grid`);
    shown.forEach((r, i) => { if (!cellCvs[i]) return;
      /* each cell from its clade's emergence, never before the burn-in it is measured from */
      const Vr = Object.assign({}, V, {x0: Math.max(V.x0 || 0, V.j0 + emergedFrom(r.times, r.tEm))});
      const p = paintOne(cellCvs[i], r, cellH, false, Vr, showEps, focus);
      /* the burn-in day, on a cell that starts before it: a clade born earlier */
      if (V.phiIdx != null && Vr.x0 < V.phiIdx) swBurninMark(cellCvs[i].getContext("2d"), p.sx(V.phiIdx), p.w, p.h);
      if (prep.brkOn) sweepBreakout(cellCvs[i].getContext("2d"), p, r, Vr, 2.6, cry, focus); });
    grid.onclick = ev => {
      const c = ev.target.closest && ev.target.closest(".rcell");
      const r = c ? shown[+c.dataset.i] : null;
      if (r) setSetting("clade", r.id === settings.clade ? null : r.id);
    };
    /* Only when the selection CHANGED: scrolling on every redraw would fight
       with the slider the redraw came from. */
    if (settings.clade === null) sweepScrolledTo[mode] = null;
    else if (settings.clade !== sweepScrolledTo[mode]) {
      const el = grid.querySelector(".rcell.sel");
      if (el && el.scrollIntoView) el.scrollIntoView({block: "nearest"});
      sweepScrolledTo[mode] = settings.clade;
    }
  } else if (!gridMode) {
    const cv = document.getElementById(`${P}_one`);
    const {g, w, h} = fit(cv, oneH);
    const x0 = V.x0 || 0;
    const sx = j => 2 + (V.nt - x0 < 2 ? 0 : (j - x0) / (V.nt - 1 - x0)) * (w - 4);
    const sy = v => h - 5 - (v / V.top) * (h - 10);
    shown.forEach(r => {
      const dim = settings.clade !== null && r.id !== settings.clade;
      g.strokeStyle = (!r.brk && prep.brkOn) || !r.eps.length ? css("--ink-3")
        : cry && r.brk && r.brk.cry ? css("--long") : (r.col || css("--ink-3"));
      g.globalAlpha = Math.min(1, (dim ? 0.12 : (r.brk || !prep.brkOn ? (showEps ? 0.35 : 0.6) : 0.22))
        * (focus === "curve" ? 1.6 : focus ? 0.5 : 1));
      g.lineWidth = r.id === settings.clade ? 2 : focus === "curve" ? 1.8 : 1;
      g.beginPath();
      const jE = emergedFrom(r.times, r.tEm);
      for (let j = jE; j < r.ys.length; j++)
        j > jE ? g.lineTo(sx(V.j0 + j), sy(r.ys[j])) : g.moveTo(sx(V.j0 + j), sy(r.ys[j]));
      g.stroke();
    });
    g.globalAlpha = 1;
    if (showEps || focus === "peak" || focus === "start") shown.forEach(r => r.eps.forEach(e => {
      const dim = settings.clade !== null && r.id !== settings.clade;
      g.strokeStyle = (!r.brk && prep.brkOn) ? css("--ink-3")
        : cry && r.brk && r.brk.cry ? css("--long") : (r.col || css("--ink-3"));
      g.lineWidth = 2.4; g.lineCap = "round"; g.lineJoin = "round";
      g.globalAlpha = dim ? 0.1 : (r.brk || !prep.brkOn ? 1 : 0.4);
      g.beginPath();
      for (let j = e.a; j <= e.b; j++)
        j > e.a ? g.lineTo(sx(V.j0 + j), sy(r.ys[j])) : g.moveTo(sx(V.j0 + j), sy(r.ys[j]));
      g.stroke(); g.lineCap = "butt";
    }));
    g.globalAlpha = 1;
    const sel = prep.brkOn ? shown.find(r => r.id === settings.clade) : null;
    if (V.phiIdx != null && x0 < V.phiIdx) drawBurnin(g, sx(V.phiIdx), w, h, V.phi, false);
    if (sel) sweepBreakout(g, {sx, sy, w, h}, sel, V, 3.5, cry, focus);
  }
}

registerTab({
  id: "sweep", group: "explore", title: "Outbreak detection", viewId: "viewRise",
  lede: "Find each clade's episodes and pick its outbreak, on one simulation. "
      + "Figure 5 measures every simulation with these rules.",
  rail: detectionRail,
  render: main => renderSweepMode("sweep", main),
});
registerTab({
  id: "cryptic", group: "explore", title: "Cryptic outbreaks", viewId: "viewCryptic",
  lede: "Of the outbreaks found on Outbreak detection, which stayed hidden: where "
      + "hidden time starts, and how long and how low makes one cryptic.",
  rail: crypticRail,
  render: main => renderSweepMode("cryptic", main),
});
/* Every simulation by scenario: no scenario or seed to pick, so the rail keeps
   only the burn-in rule and the cryptic rules */
function crystatsRail(rail) {
  railCompact(rail, {toggle: false});
  rail.classList.add("hintpop-on");
  wireHintPopups();
  rail.append(card("Simulations", [swBurninControl()]));
  crypticCards(rail);
}
registerTab({
  id: "crystats", group: "explore", title: "Cryptic by scenario", viewId: "viewCryStats",
  lede: "Outbreaks and cryptic outbreaks per simulation in every scenario, under the "
      + "rules set on Outbreak detection and Cryptic outbreaks, and whether the "
      + "scenarios differ.",
  rail: crystatsRail,
  render: main => renderCryScen(main, "crystats"),
});
