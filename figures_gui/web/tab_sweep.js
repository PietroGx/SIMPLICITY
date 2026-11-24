/* Breakout detection and Cryptic breakouts: two tabs, one run each, one painter.

   Breakout detection (id "sweep") finds each clade's episodes and picks its
   breakout: smoothing, which peaks count, where an episode starts, F. Cryptic
   breakouts (id "cryptic") takes those breakouts as they are and asks which
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

/* The last Muller painted on Breakout detection, kept so check_app.js can pin
   the band hit test to the geometry that was actually drawn. */
let sweepMuller = null;

/* ------------------------------------------------- detection rules (sd*)

   These tabs' settings are their own (sd*), never the build parameters
   win/prom/lo: dragging the window to 45 days here must not point every other
   tab at a build that does not exist. SD_FIG5 is what Figure 5 measured with;
   at those values explore.py IS the build's measurement. */
const SD_WIN = [1, 2, 3, 5, 7, 10, 15, 21, 30, 45, 60, 90, 120, 180];
const SD_PROM = [0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9];
const SD_TROUGH = [0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5];
const SD_FRAC = [0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5];
const SD_FIG5 = {sdKernel: "box", sdSpace: "lin", sdWidth: 0,
                 sdCens: 1, sdStart: "trough", sdTrough: 0.0001,
                 sdAnchor: "emergence", sdPick: "first", sdGate: "start"};
const SD_LABEL = {sdKernel: "kernel", sdSpace: "log smoothing", sdWin: "window",
  sdProm: "prominence", sdWidth: "width", sdCens: "censored peaks", sdStart: "episode start",
  sdTrough: "dip depth", sdAnchor: "hidden-time start", sdPick: "breakout pick",
  sdGate: "burn-in rule"};

/* Figure 5's rules: SD_FIG5 plus the build's own window and prominence */
const sweepFig5Patch = () => Object.assign({}, SD_FIG5,
  {sdWin: +settings.win, sdProm: +settings.prom});
const sweepOffRules = () => Object.entries(sweepFig5Patch())
  .filter(([k, v]) => settings[k] !== v).map(([k]) => k);

const sdPct = v => v === 0 ? "0" : brkPct(v);
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
  const opts = [["start", "Breakout starts after it"], ["trunc", "Not circulating at it"],
                ["off", "Ignore it"]];
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `<label>Burn-in</label>
    <select>${opts.map(([v, t]) => `<option value="${v}"${v === cur ? " selected" : ""}>${t}</option>`).join("")}</select>
    <div class="hint">Which clades the burn-in day, the first day everyone has been
      infected once, lets through. <b>Breakout starts after it</b>: measured on the
      whole record, a clade counts only if its breakout begins after the burn-in day;
      one that broke out before it is left out. <b>Not circulating at it</b> (Figure 5):
      measured from the burn-in day, a clade already above 0.01% then is left out, even
      if its breakout comes later. <b>Ignore it</b>: the whole record, every clade.</div>`;
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
    "The coloured line. Every mark is measured on it, never on the raw daily record.", [
    selectControl({setting: "sdKernel", label: "Kernel",
      options: SD_KERNELS,
      hint: "<b>Moving average</b> (Figure 5): every day in the window weighs the "
          + "same. <b>Gaussian</b> weights the centre, so a short dip survives. "
          + "<b>Median</b> ignores single-day spikes."}),
    sliderControl({setting: "sdWin", label: "Window", values: SD_WIN,
      format: v => v + " d",
      hint: "How many days the kernel averages over (the Gaussian's full width "
          + "at half maximum). Wider merges neighbouring episodes and spreads a "
          + "rise up to half a window before the clade emerged. Figure 5: 15 d."}),
    segmentControl({setting: "sdSpace", label: "Smooth in",
      options: [["lin", "Frequency"], ["log", "log₁₀"]],
      hint: "Smoothing in log₁₀ space weighs a climb from 0.01% to 1% as much as "
          + "the last doubling, so early, small episodes survive. Moved 6–7 of "
          + "235–425 breakouts on 50 runs."}),
  ]));
  rail.append(swCard("peak", "2 · Episodes",
    "Which peaks count as episodes: each one is a bold stretch of the curve, from its start to its peak. "
    + "Selected, each peak shows the width Min width is tested on, at half its prominence; "
    + "a peak it drops keeps a dashed bar and a cross.", [
    sliderControl({setting: "sdProm", label: "Prominence",
      values: SD_PROM, format: v => Number(v).toFixed(2),
      hint: "How far a peak must stand above the lower of the dips either side, "
          + "as a fraction of the clade's own highest frequency. Lower finds more, "
          + "smaller episodes. Figure 5: 0.10."}),
    sliderControl({setting: "sdWidth", label: "Min width",
      values: SW_WIDTH, format: v => v ? v + " d" : "any",
      hint: "Drops a peak narrower than this at half its prominence. Below 14 d "
          + "the 15-day window already does this; at 30 d it moved 39 of 425 "
          + "breakouts on 50 runs. Click this card to see each peak's width."}),
    checkControl({setting: "sdCens", label: "Peak still high at the end",
      hint: "A clade holding its maximum to the end of the record never comes "
          + "down, so no peak finder returns it. Ticked (Figure 5), it counts, "
          + "flagged censored: its rise is real, its peak a lower bound. So does one "
          + "still rising on the last day, its peak taken on that day."}),
  ]));
  const start = settings.sdStart;
  rail.append(swCard("start", "3 · Episode start",
    "Where each episode begins: the left end of its bold stretch. On the breakout "
    + "this is the red dot, where launch frequency is read and hidden time ends.", [
    segmentControl({setting: "sdStart", label: "Rule",
      options: [["trough", "Last dip"], ["frac", "Below %"], ["base", "Base"]],
      hint: "<b>Last dip</b> (Figure 5): the last local minimum before the peak "
          + "that is at least the depth below deep, and lies at least the "
          + "prominence below the peak, so a wobble at the top of a plateau never "
          + "starts an episode. <b>Below %</b>: the last day "
          + "the clade was under that share of the peak. <b>Base</b>: scipy's "
          + "prominence base, the lowest point before a taller peak; it reaches "
          + "far left, which is why Figure 5 does not use it."}),
    start === "frac" ? sliderControl({setting: "sdFrac", label: "Below",
      values: SD_FRAC, format: v => fmt(100 * v, v < 0.01 ? 1 : 0) + "% of peak",
      hint: "The start is the last day the clade was under this share of the "
          + "peak. Under about 1% that is within a day of emergence, when one "
          + "host carries it, so hidden time collapses to 0; 2–5% catches the "
          + "low stretch before the rise."}) : null,
    start === "trough" ? sliderControl({setting: "sdTrough", label: "Dip depth",
      values: SD_TROUGH, format: v => v === 0 ? "any" : fmt(100 * v, v < 0.01 ? 2 : 0) + "% of peak",
      hint: "How deep a dip must be to count as a start, as a share of the "
          + "clade's own peak. Deeper skips shallow wiggles and moves the start "
          + "left. Figure 5: 0.01%, the smallest that stops noise choosing it "
          + "(clade_preprocess.py:84)."}) : null,
  ]));
  rail.append(swCard("breakout", "4 · The breakout",
    "The episode Figure 5 times: the black dot at its peak. F is the dashed line.", [
    sliderControl({setting: "brkF", label: "Peak ≥ <i>F</i>", values: BRK_F,
      format: v => brkPct(v),
      hint: "An episode is a breakout only if it peaks at or above this share of "
          + "all virus that day. Shared with Figure 5."}),
    segmentControl({setting: "sdPick", label: "Which",
      options: [["first", "First ≥ F"], ["largest", "Largest ≥ F"]],
      hint: "<b>First</b> (Figure 5): the first episode to reach F; earlier, "
          + "smaller ones count as hidden. <b>Largest</b>: the tallest episode "
          + "reaching F; every episode before it counts as hidden."}),
  ]));
}

function crypticRail(rail) {
  swRailStart(rail);
  crypticCards(rail);
}
function crypticCards(rail) {
  rail.append(swCard("hidden", "1 · Hidden time",
    "The red band: from the solid line to the breakout's start (the red dot). "
    + "The breakout itself is set on Breakout detection.", [
    segmentControl({setting: "sdAnchor", label: "From",
      options: [["emergence", "Emergence"], ["transmission", "1st transmission"]],
      hint: "<b>Emergence</b> (Figure 5): when the clade's root lineage arose "
          + "inside a host, its Time_emergence. <b>1st transmission</b>: the "
          + "earliest infection that inherited any of the clade's lineages, when "
          + "it first left a host. The gap between them is time spent in one host."}),
  ]));
  rail.append(swCard("cryptic", "2 · Cryptic if",
    "Tinted cells: the red band at least H long, and, with the box ticked, the curve "
    + "never above the red dot inside it.", [
    sliderControl({setting: "brkH", label: "Hidden ≥ <i>H</i>", values: BRK_H,
      format: v => v + " d",
      hint: "A breakout is cryptic only if its hidden time, the red band, is at "
          + "least this long. Shared with Figure 5."}),
    checkControl({setting: "brkQuiet", label: "No rise above the launch before it",
      hint: "Ticked: hidden time only counts if the clade never rose above its "
          + "launch level (the red dot, drawn as the red dotted line across the band) "
          + "between the anchor and the breakout's start. A clade with an evident "
          + "earlier rise was not hidden, so it is not cryptic. Unticked: hidden "
          + "time alone decides. With the start rule <b>Last dip</b> the curve always "
          + "comes down into the red dot from above, so ticked excludes nearly every "
          + "breakout; <b>Below %</b> does not have that problem. Shared with Figure 5."}),
  ]));
}

/* The view controls, above the cells: rebuilt on every render, so they
   always show the current values */
function swTools(mode, box, brkOn) {
  const M = SW_MODES[mode], grid = settings.episodeView === "grid";
  const cellOpts = M.cryptic ? [["cryptic", "Cryptic"], ["breakout", "Breakouts"], ["all", "All"]]
                             : [["breakout", "Breakouts"], ["all", "All"]];
  [
    brkOn ? segmentControl({setting: M.show, label: "Show", options: cellOpts,
      hint: M.cryptic ? "Cryptic breakouts first, then other breakouts, then the rest."
                      : "Clades with a breakout first, then the rest."}) : null,
    checkControl({setting: M.eps, label: "Episodes",
      hint: "Every episode drawn bold in the clade's colour, from its start to its "
          + "peak; the rest of the curve fades."}),
    checkControl({setting: "showRaw", label: "Raw",
      hint: "The unsmoothed daily record in grey behind the smoothed line: how "
          + "much the kernel is doing."}),
    segmentControl({setting: "episodeView", label: "Layout",
      options: [["grid", "Cells"], ["one", "Overlay"]],
      hint: "<b>Cells</b>: one plot per clade. <b>Overlay</b>: every clade on one "
          + "axis; click a band above to mark one."}),
    grid ? segmentControl({setting: "episodeCols", label: "Columns",
      options: [1, 2, 3, 4, 5, 6].map(n => [n, String(n)])}) : null,
    grid ? segmentControl({setting: "yScale", label: "Height",
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
      <b>${Number(settings.sdProm).toFixed(2)}</b> of the clade's height above their
      dips${+settings.sdWidth ? `, at least <b>${settings.sdWidth} d</b> wide` : ""}, on a
      <b>${kern} ${settings.sdWin} d</b>${
      settings.sdSpace === "log" ? " in log₁₀" : ""}. Each starts at <b>${swStartText()}</b>.
      The breakout is ${brk}. <a href="#cryptic">Which are cryptic →</a></div>`;
  const anchor = settings.sdAnchor === "emergence" ? "emergence in host" : "first transmission";
  return `<div class="sw-rule hintpop-on">${swChip(mode, offRules)} Hidden from <b>${anchor}</b> to the
    breakout's start. Cryptic: hidden ≥ <b>${settings.brkH} d</b>${settings.brkQuiet
      ? ", <b>never above its launch level before it</b>" : ""}.
    <span class="sw-sub">Breakouts from <a href="#sweep">← Breakout detection</a>:
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
      ["breakouts", brks.length.toLocaleString(), `peak ≥ ${brkPct(+settings.brkF)}`],
      ["median launch", med(brks.map(r => r.brk.launch), brkPct), ""],
      ["median breakout rise", med(brks.map(r => r.brk.rise).filter(isNum), Math.round), "days"],
    ] : [["median clade rise", med(withEp.map(r => r.value), Math.round), "days"]]),
  ] : (() => {
    const cry = brks.filter(r => r.brk.cry), zero = brks.filter(r => r.brk.hidden <= 0);
    const gap = rows.filter(r => isNum(r.tEm) && isNum(r.tTx)).map(r => r.tTx - r.tEm);
    return [
      ["breakouts", brks.length.toLocaleString(), `peak ≥ ${brkPct(+settings.brkF)}`],
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
               [vl("none"), settings.sdAnchor === "emergence" ? "emergence in host" : "first transmission"]] : []),
    [dot("var(--long)"), "breakout starts, at its launch frequency"],
    [dot("var(--ink)"), "breakout peak"],
    ...(cry ? [] : [[hl("4 3"), "<i>F</i>"]]),
    ...(cry && settings.brkQuiet ? [[swKeyIcon(`<line x1="1" y1="6" x2="21" y2="6" style="stroke:var(--long)" stroke-width="1.4" stroke-dasharray="2 3"/>`),
               "launch level: the curve must stay under it"]] : []),
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

/* Why a cell is greyed out: it has no breakout under the rules in force */
function swWhyGrey(r, V) {
  const F = +settings.brkF;
  if (isNum(r.preBurn)) return `<b>Left out: its breakout began before the burn-in day.</b> `
    + `It starts on day ${Math.round(r.preBurn)}; the burn-in day is day ${Math.round(V.phi)}, `
    + `and under <b>Breakout starts after it</b> only later breakouts count.`;
  if (r.trunc) return `<b>No breakout: already circulating at the burn-in day</b>`
    + `${isNum(V.phi) ? ` (day ${Math.round(V.phi)})` : ""}. It emerged before the `
    + `measured window, so there is no emergence to time hidden time from, and Figure 5 leaves it out.`;
  if (r.eps.length) {
    const top = Math.max(...r.eps.map(e => r.ys[e.b]));
    return `<b>No breakout: no episode peaks at F or above.</b> Its tallest `
      + `episode peaks at ${brkPct(top)}, under <i>F</i> = ${brkPct(F)}. Lower `
      + `<i>F</i> on <b>4 · The breakout</b> to count it.`;
  }
  let iMax = 0;
  r.ys.forEach((v, j) => { if (v > r.ys[iMax]) iMax = j; });
  const atEnd = iMax >= r.ys.length - 2 && !settings.sdCens;
  return `<b>No breakout: no episode.</b> No peak of its smoothed curve stands `
    + `${Number(settings.sdProm).toFixed(2)} of its height above its dips`
    + `${settings.sdWidth ? ` and is ${settings.sdWidth} d wide` : ""}.`
    + (atEnd ? " Its maximum is at the end of the record and <b>Peak still high at "
             + "the end</b> is unticked." : "")
    + ` The plot above stacks only clades with an episode, so it leaves this one `
    + `out; click the cell to draw it there.`;
}

/* Cells in rank order, each rank by peak; the Show filter cuts from the end.
   Detection ranks breakout, rest; cryptic ranks cryptic, breakout, rest. */
function sweepCells(mode, rows, brkOn) {
  const cry = SW_MODES[mode].cryptic, show = settings[SW_MODES[mode].show];
  const rank = r => r.brk ? (cry && r.brk.cry ? 0 : 1) : 2;
  const keep = !brkOn || show === "all" ? 2 : show === "breakout" ? 1 : 0;
  return rows.filter(r => rank(r) <= keep)
    .sort((a, b) => rank(a) - rank(b) || b.peak - a.peak);
}

/* Per simulation, per scenario: breakouts and cryptic breakouts under the tabs'
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
      r.clades.forEach(([, , tem, eps]) => {
        const b = breakoutOf(eps.map(([aDay, top, amp, days, pre]) => ({aDay, top, amp, days, pre})),
                             tem, settings.sdPick, gate);
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

async function renderCryScen(main, mode = "crystats") {
  await swPresetsLoad();
  const head = swRuleLine(mode, sweepOffRules());
  const wire = () => { swPresetWire(mode);
    const reset = document.getElementById(`${SW_MODES[mode].p}_reset`);
    if (reset) reset.onclick = () => swApply(sweepFig5Patch(), mode); };
  main.innerHTML = head + `<div class="sw-sub">Measuring every simulation under these rules…</div>`;
  wire();
  const d = await apiGetSlot("/api/explore_panel?" + exploreQuery({mf: 0}));
  const C = swScenCounts(d);
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
  const ms = a => `${fmt(swMean(a), 2)} <small>± ${fmt(swSd(a), 2)}</small>`;
  const rows = scs.map(sc => { const c = C[sc], nb = c.brk.reduce((s, x) => s + x, 0),
      nc = c.cry.reduce((s, x) => s + x, 0);
    return `<tr><th>${SCENARIO_LABEL[sc]}</th><td>${c.brk.length}${c.skipped
        ? ` <span class="sw-nophi" title="left out: no burn-in day">−${c.skipped}</span>` : ""}</td>
      <td>${ms(c.brk)}</td><td>${ms(c.cry)}</td>
      <td>${nb ? fmt(100 * nc / nb, 1) + "%" : "—"} <small>${nc}/${nb}</small></td></tr>`; }).join("");
  const matrix = (key, title) => {
    const t = T && T[key];
    if (!t) return "";
    const P = {}, mean = sc => swMean(C[sc][key === "cryptic" ? "cry" : "brk"]);
    t.pairs.forEach(x => { P[x.a + "|" + x.b] = P[x.b + "|" + x.a] = x; });
    const cell = (a, b) => { if (a === b) return `<td class="dg"></td>`;
      const x = P[a + "|" + b]; if (!x) return "<td>—</td>";
      const sig = x.p_holm < 0.05, hi = mean(a) > mean(b) ? "▲" : mean(a) < mean(b) ? "▼" : "";
      return `<td class="${sig ? "sig " + (hi === "▲" ? "hi" : hi === "▼" ? "lo" : "") : ""}" title="U = ${x.U}, p = ${x.p.toPrecision(3)} unadjusted">${swP(x.p_holm)} ${sig ? hi : ""}</td>`; };
    const kw = t.kruskal;
    return `<div class="sw-mx"><h4>${title}</h4>
      <div class="sw-sub">Kruskal–Wallis across ${kw ? kw.groups : 0} scenarios: H = ${kw ? fmt(kw.H, 2) : "—"},
        p ${kw ? (kw.p < 0.001 ? "&lt; 0.001" : "= " + swP(kw.p)) : "= —"}</div>
      <table class="sw-tab"><tr><th></th>${scs.map(sc => `<th>${SCENARIO_LABEL[sc]}</th>`).join("")}</tr>
      ${scs.map(a => `<tr><th>${SCENARIO_LABEL[a]}</th>${scs.map(b => cell(a, b)).join("")}</tr>`).join("")}</table></div>`;
  };
  main.innerHTML = head + `
    <div class="plot"><div class="ph"><h3>Breakouts by scenario</h3>
      <span class="sub">per simulation · k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></div>
      <table class="sw-tab sw-sum"><tr><th>Scenario</th><th>Simulations</th><th>Breakouts / simulation</th>
        <th>Cryptic / simulation</th><th>Cryptic share</th></tr>${rows}</table>
      <div class="sw-sub">Mean ± standard deviation over simulations. Cryptic share pools every
        breakout in the scenario.${settings.burnin ? " Simulations that never reach the burn-in day are left out (red count)." : ""}</div></div>
    <div class="plot"><div class="ph"><h3>Do scenarios differ?</h3>
      <span class="sub">two-sided Mann–Whitney U per pair, Holm-adjusted over ${scs.length * (scs.length - 1) / 2} pairs</span></div>
      ${terr ? `<div class="note err">The tests did not run: ${sdEsc(terr)}</div>` : ""}
      <div class="sw-mxs">${matrix("cryptic", "Cryptic breakouts per simulation")}
        ${matrix("breakouts", "Breakouts per simulation")}</div>
      <div class="sw-sub">Each cell: the Holm-adjusted p for the row scenario against the column
        scenario; where p &lt; 0.05 it is red with ▲ if the row's mean is higher, green with ▼ if lower.
        Hover a cell for U and the unadjusted p. The simulation is the unit: one count per run.</div></div>`;
  wire();
}

async function renderSweepMode(mode, main) {
  const M = SW_MODES[mode], P = M.p, cry = M.cryptic;
  main.classList.add("fit");
  let dailyErr = null;
  /* Detection times nothing from an anchor, so it asks with Figure 5's: which
     clades carry a breakout must not move with the cryptic tab's choice */
  const q = exploreQuery(Object.assign({scenario: settings.scenario, seed: settings.seed, mf: 0},
                                       cry ? {} : {anchor: "emergence"}));
  await swPresetsLoad();
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
  const anchorWord = settings.sdAnchor === "emergence" ? "emergence" : "1st transmission";

  const offNote = prep.fig5Off === "nophi"
    ? `<div class="note nophi"><span class="sw-nophi">no burn-in day</span> <strong>No breakouts
       in this simulation.</strong> (diagnosed + recovered) / N never reaches 1, so it has no
       burn-in day and Figure 5 takes nothing from it. Set <em>Burn-in</em> to <em>Ignore it</em>
       to see its breakouts.</div>`
    : prep.fig5Off ? `<div class="note"><strong>Not Figure 5's numbers.</strong>
       Measured ${settings.burnin ? "from the burn-in day" : "over the whole record"},
       every clade counted. Set <em>Burn-in</em> in Simulation.</div>` : "";
  const errNote = prep.daily ? "" : `<div class="note err"><strong>No daily record.</strong>
       The detector could not read this run (${sdEsc(dailyErr || "no answer")}); the
       cells fall back to the 9.2-day export and no breakout is drawn.</div>`;
  const unserved = prep.daily && prep.unserved ? `<div class="sw-sub">${prep.unserved}
       clade${prep.unserved === 1 ? "" : "s"} carry no daily series and are not shown.</div>` : "";
  const what = !prep.brkOn || show === "all" ? `${rows.length} clades`
    : show === "cryptic" ? `${cells.length} cryptic of ${rows.length} clades`
    : `${cells.length} breakouts of ${rows.length} clades`;
  const empty = show === "cryptic" ? "cryptic breakout" : "breakout";

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
      <h3>${gridMode ? (cry ? "Each clade, its hidden time marked" : "Each clade, its breakout marked")
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
        : r.eps.length ? "no breakout" : "no peak";
      const corner = !prep.brkOn ? (r.main ? `rise ${Math.round(r.value)} d` : "no peak")
        : !b ? none
        : cry ? `hidden ${Math.round(b.hidden)} d${b.cry ? " · <b>cryptic</b>" : ""}`
        : `breakout · rise ${isNum(b.rise) ? Math.round(b.rise) : "—"} d`;
      const tip = !b ? "" : cry ? `from ${anchorWord} to the breakout's start`
        : "the breakout episode's rise days";
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
  id: "sweep", group: "explore", title: "Breakout detection", viewId: "viewRise",
  lede: "Find each clade's episodes and pick its breakout, on one simulation. "
      + "Figure 5 measures every simulation with these rules.",
  rail: detectionRail,
  render: main => renderSweepMode("sweep", main),
});
registerTab({
  id: "cryptic", group: "explore", title: "Cryptic breakouts", viewId: "viewCryptic",
  lede: "Of the breakouts found on Breakout detection, which stayed hidden: where "
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
  lede: "Breakouts and cryptic breakouts per simulation in every scenario, under the "
      + "rules set on Breakout detection and Cryptic breakouts, and whether the "
      + "scenarios differ.",
  rail: crystatsRail,
  render: main => renderCryScen(main, "crystats"),
});
