/* Panel 2: ported from panel2.html.

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_panel2 when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.

   One metric and one plot form at a time, both from the rail. `metric2` gates
   the plot list -- `who_transmits` has no 2x2 to put in a heat map, and only
   `where` carries a day curve -- so a pair that cannot be drawn is not
   selectable rather than silently blank.

   Canvas ids are p2_* : shell.js:138 builds all ten tab views into ONE
   document, so a bare `c1` here and a bare `c1` in tab_panel3.js were the same
   element, and whichever rendered second painted over the other.
*/

const P2_CELLS = ["ll", "ls", "sl", "ss"];
const P2_CELL = {ll: "long → long", ls: "long → standard",
                 sl: "standard → long", ss: "standard → standard"};
const p2CellColour = () => ({ll: css("--long"), ls: css("--sc3"),
                             sl: css("--sc2"), ss: css("--std")});

const P2_METRICS = [
  ["who_transmits", "Who transmits",
   "The percentage of all transmission events whose "
   + "<strong>transmitter</strong> is a long shedder: "
   + "(<code>ll</code>+<code>ls</code>) over all four cells "
   + "(panels.py:263-264). The denominator is <em>events</em>, so it rises with "
   + "long-shedder prevalence on its own."],
  ["where", "Where it goes",
   "The 2×2 split of events by the host type at each end, every cell as a "
   + "percentage of all events (panels.py:265-266). <code>ll</code> is "
   + "long→long, <code>ls</code> long→standard, <code>sl</code> standard→long, "
   + "<code>ss</code> standard→standard."],
  ["intensity", "Intensity",
   "Transmission events <strong>per host of that type</strong>: "
   + "(<code>ll</code>+<code>ls</code>) over <code>n_long</code> and "
   + "(<code>sl</code>+<code>ss</code>) over <code>n_standard</code> "
   + "(panels.py:267-270). The denominator is <em>hosts</em>, not events, so "
   + "prevalence cannot inflate it. <strong>This is the figure "
   + "candidate.</strong>"],
];

const P2_PLOTS = {
  who_transmits: [
    ["dots", "Dot and interval",
     "One point per scenario with its bootstrap CI of the mean over runs."],
    ["bars", "Grouped bars",
     "The same five numbers as bars from zero — the form the panel shipped in."],
    ["prev", "Against prevalence",
     "The percentage against each scenario's long-shedder prevalence. The "
     + "climb here is the confound, not the finding: read Intensity instead."],
  ],
  where: [
    ["heatmap", "2×2 heat map",
     "One matrix per scenario, rows the transmitter and columns the recipient, "
     + "on a shared 0–100% ramp so cells compare across scenarios."],
    ["stacked", "100% stacked",
     "One bar per scenario, one band per cell: the composition of the events "
     + "and nothing about their number."],
    ["bars", "Grouped bars",
     "The four cells side by side per scenario, each with its own interval — "
     + "the only form of this metric that carries one."],
    ["days", "Over days",
     "The four cells as daily counts for one scenario, averaged over runs "
     + "inside the panel (panels.py:284-289)."],
  ],
  intensity: [
    ["ratio", "Long ÷ standard",
     "Events per long shedder divided by events per standard host. One is "
     + "parity; above one means each long shedder transmits more."],
    ["paired", "Long and standard",
     "Both per-capita rates side by side, each with its own interval, so the "
     + "ratio can be read off its two parts."],
    ["prev", "Ratio against prevalence",
     "The per-capita ratio against long-shedder prevalence. Flat is the claim: "
     + "a per-host rate that does not track how many hosts there are."],
  ],
};

/* The ratio is formed from two MEANS, so it gets no interval: dividing CI
   bounds is not a CI of a ratio. /api/dist's `intensity_ratio` is the per-run
   ratio and is the only honest source of one; it is drawn as points. */
function p2Ratio(s) {
  const L = ((s || {}).intensity || {}).long;
  const S = ((s || {}).intensity || {}).standard;
  const miss = v => !v || v.mean === null || v.mean === undefined;
  if (miss(L) || miss(S) || !S.mean) return null;
  return {n: Math.min(L.n || 0, S.n || 0), mean: L.mean / S.mean,
          lo: null, hi: null};
}

/* Per-run values laid out over the SCENARIO axis: dist.js's distRuns maps
   categories within one scenario, and `paired` has the scenarios as its
   categories instead. */
function p2RunsOverScenarios(dist, key) {
  const s = dist && dist.series && dist.series[key];
  if (!s) return undefined;
  const out = SCENARIOS.map(sc => s[sc] || null);
  return out.some(r => r && r.length) ? out : undefined;
}

/* settings.js:47 lists metric2 and plot2 as RAIL_SETTINGS, so shell.js:94-97
   clears and rebuilds this rail when either moves -- and only then, so a drag
   elsewhere in the rail does not lose its focus. The metric pick runs first, so
   its correction of plot2 is in place before the plot pick and the
   distribution box read it. */
function rail_panel2(rail) {
  if (!rail) return;
  rail.append(card("Measure", [
    pickControl({setting: "metric2", label: "What is measured",
      options: P2_METRICS,
      hint: "All three come from one transmission table, counted by the host "
          + "type at each end. Independent of <i>k</i> and of the clustering."}),
    pickControl({setting: "plot2", label: "Plot form",
      options: P2_PLOTS[settings.metric2] || P2_PLOTS.intensity}),
  ]));
  rail.append(card("Runs", [railDistribution("dist2", settings.plot2)]));
  rail.append(card("Window", [
    checkControl({setting: "burnin",
      label: "Count only events on or after the burn-in day",
      hint: "With per-day counts the burn-in cut is a sum range, so it is exact "
          + "and any window is available. 5 of 250 runs never reach the burn-in "
          + "and drop out when this is on."}),
  ]));
}

const p2Plot = (h, sub, body, key, cap) => `
  <div class="plot"><div class="ph"><h3>${h}</h3>
    <span class="sub">${sub}</span></div>
    ${body}${key || ""}${cap ? `<div class="cap">${cap}</div>` : ""}</div>`;

const P2_CONFOUND =
  "<code>who_transmits</code> counts events, so it rises mechanically with "
  + "prevalence and cannot separate &ldquo;long shedders matter&rdquo; from "
  + "&ldquo;there are more of them&rdquo; — the confound the scenario design "
  + "exists to break. SOT and HIV low hold prevalence at 1.0% and move duration "
  + "(63 d against 109 d); HIV low and HIV high hold duration at 109 d and move "
  + "prevalence (1.0% against 12.0%). Scenarios are never pooled.";

const P2_RATIO_NOTE =
  "The point is one mean divided by another, so it carries <strong>no "
  + "interval</strong>: dividing two CIs is not a CI of a ratio. "
  + "<code>intensity_ratio</code> in /api/dist is the per-run ratio and is the "
  + "only honest source of one; until it lands the point stands alone, and the "
  + "cloud it then draws need not centre on the point. control has no long "
  + "shedders at all (<code>n=0</code>), so its ratio is absent, not zero.";

async function render_panel2(main) {

  const metric = settings.metric2, plot = settings.plot2;
  const points = !!settings.dist2 && distAllowed(plot);

  let d; try { d = await apiGetSlot("/api/panel?n=2&" + buildQuery()); }
  catch (e) { main.innerHTML = ""; main.append(errorBox(e.message)); return; }
  const dist = points ? await distFor(2) : null;

  /* panel 2 carries no prevalence; it is a scenario design constant and
     /api/panel?n=1 already reports it, so this is one cached request and not a
     recompute. It needs a clade build, which panel 2 itself does not, so a
     missing one leaves the x axis empty rather than invented. */
  let prev = null;
  if (plot === "prev") {
    try { prev = (await apiGetSlot("/api/panel?n=1&" + buildQuery())).scenarios; }
    catch (e) { prev = null; }
  }
  const prevOf = sc => {
    const p = ((prev || {})[sc] || {}).prevalence;
    return p === null || p === undefined ? null : {mean: asPctNum(p), lo: null, hi: null};
  };
  const noPrev = prev ? "" : "<strong>Prevalence could not be read</strong> "
    + "from panel 1 at this build, so the x axis is empty rather than invented. ";

  const CC = p2CellColour();
  const scKey = keyRow(SCENARIOS.map(sc => [scenarioColour(sc), SCENARIO_LABEL[sc]]));
  const cellKey = keyRow(P2_CELLS.map(c => [CC[c], P2_CELL[c]]));
  const sOf = sc => d.scenarios[sc] || {};
  /* who_transmits and the where cells arrive as fractions; intensity is events
     per host and is not a share */
  const whoOf = sc => asPct(sOf(sc).who_transmits);
  const cv1 = () => document.getElementById("p2_c1");

  let head = "", sub = "", body = `<canvas id="p2_c1"></canvas>`;
  let key = "", cap = "", paint = () => {};

  if (metric === "who_transmits") {
    head = "Who transmits";
    sub = "% of all events whose transmitter is a long shedder";
    const vals = sc => [whoOf(sc) || null];
    if (plot === "prev") {
      key = scKey;
      cap = `Each scenario's percentage against its own long-shedder prevalence.
        ${noPrev}SOT and HIV low sit at the same 1.0% and differ only in
        duration, so they overlap on x by design. ${P2_CONFOUND}`;
      paint = () => scatterCI(cv1(), {
        xlab: "long-shedder prevalence (%)", ylab: "% of events", height: 280,
        points, zero: 1, fmtY: v => fmt(v, 1), fmtX: v => fmt(v, 1),
        groups: SCENARIOS.map(sc => ({colour: scenarioColour(sc),
          label: SCENARIO_LABEL[sc], x: prevOf(sc), y: whoOf(sc),
          yruns: ((distRuns(dist, sc, ["w"], "who_transmits") || [])[0] || []).map(asPctNum)}))});
    } else if (plot === "bars") {
      key = scKey;
      cap = `Bars from zero, the form the panel shipped in. This is the
        host-level counterpart to panel 1's mutation-level claim.
        ${P2_CONFOUND}`;
      paint = () => bars(cv1(), {
        cats: ["% of events"], ylab: "%", pct: 1, height: 280,
        series: SCENARIOS.map(sc => ({colour: scenarioColour(sc),
          vals: vals(sc)}))});
    } else {
      key = scKey;
      cap = `Dot and bootstrap interval, the five scenarios in one category.
        Compare against the long-shedder prevalence of each: 0%, 1.0%, 1.0%,
        12.0%, 0.1%. ${P2_CONFOUND}`;
      paint = () => dots(cv1(), {
        cats: ["% of events"], ylab: "% of events", height: 280,
        points, zero: 1, pct: 1, fmtY: v => fmt(v, 1),
        series: SCENARIOS.map(sc => ({colour: scenarioColour(sc),
          label: SCENARIO_LABEL[sc], vals: vals(sc),
          runs: asPctRuns(distRuns(dist, sc, ["w"], "who_transmits"))}))});
    }
  }

  if (metric === "where") {
    head = "Where the virus goes";
    sub = "the four cells as % of all events";
    const cellOf = (sc, c) => asPct((sOf(sc).where || {})[c]) || null;
    const cellMean = (sc, c) => {
      const b = cellOf(sc, c);
      return b && b.mean !== null && b.mean !== undefined ? b.mean : null;
    };
    if (plot === "heatmap") {
      sub = "rows = transmitter, columns = recipient · % of all events";
      body = `<div style="display:grid;gap:10px;
        grid-template-columns:repeat(auto-fit,minmax(190px,1fr))">`
        + SCENARIOS.map(sc => `<div><div style="font:11px var(--f-ui);
            color:var(--ink-3);padding:0 0 2px 4px">${SCENARIO_LABEL[sc]}</div>
          <canvas id="p2_c1_${sc}"></canvas></div>`).join("") + `</div>`;
      cap = `One matrix per scenario on a shared 0–100% ramp, so a cell can be
        compared across scenarios; every cell also carries its number, because
        a value cannot be read off a ramp. The diagonal says whether long
        shedders mainly infect each other — a different question from how much
        they transmit. A cell of a matrix has room for neither an interval nor
        per-run points, so the distribution box is disabled here.`;
      paint = () => SCENARIOS.forEach(sc => heatcells(
        document.getElementById("p2_c1_" + sc), {
          rows: ["long", "standard"], cols: ["long", "standard"], height: 130,
          max: 100, fmtV: v => fmt(v, 2),
          grid: [[cellMean(sc, "ll"), cellMean(sc, "ls")],
                 [cellMean(sc, "sl"), cellMean(sc, "ss")]]}));
    } else if (plot === "stacked") {
      sub = "composition of the events, one bar per scenario";
      key = cellKey;
      cap = `The composition of the events and nothing about their number: every
        bar is the same height whatever the scenario's transmission count. The
        four cells already sum to 100% inside each run, so normalising changes
        no value — it makes the reading explicit. No interval can be drawn on a
        stack and no per-run form of a stack stays a stack, so the distribution
        box is disabled rather than ignored.`;
      paint = () => stacked(cv1(), {
        cats: SCENARIOS.map(sc => SCENARIO_LABEL[sc]), pct: 1, height: 280,
        fmtV: v => fmt(v, 1),
        parts: P2_CELLS.map(c => ({colour: CC[c],
          vals: SCENARIOS.map(sc => cellMean(sc, c))}))});
    } else if (plot === "days") {
      sub = "mean events per day across runs, one scenario";
      body = `<select id="p2_scsel" style="max-width:220px">${SCENARIOS.map(sc =>
          `<option value="${sc}">${SCENARIO_LABEL[sc]}</option>`).join("")}</select>
        <canvas id="p2_c1"></canvas>`;
      key = cellKey;
      cap = `Resolved over days, one scenario at a time. The day curves are
        averaged across runs inside the panel (panels.py:284-289), so there are
        no per-run curves to overlay — the distribution box is disabled here,
        and per-run curves would be about a million floats besides.
        <strong>The burn-in toggle reaches these curves differently</strong>
        from every other form on this tab: the curve is accumulated from the
        whole record and is never cut at the burn-in day (panels.py:272-277),
        but a run that never reaches the burn-in is dropped from the average
        along with everything else (panels.py:255-257). So the toggle moves the
        denominator here and the day range everywhere else — the x axis below
        still starts at day 0 with the toggle on.`;
      paint = () => {
        const sel = document.getElementById("p2_scsel");
        const curve = () => {
          const pick = sel && d.series[sel.value] ? sel.value : SCENARIOS[0];
          const s = d.series[pick];
          if (!s) return;
          linesOverDays(cv1(), {
            ylab: "events / day", xmax: s.days, height: 280,
            series: P2_CELLS.map(c => ({colour: CC[c], y: s.mean[c] || []}))});
        };
        if (sel) sel.addEventListener("change", curve);
        curve();
      };
    } else {
      key = cellKey;
      cap = `The only form of this metric that carries an interval: four cells
        per scenario, each a bootstrap CI of the mean over runs. Whether long
        shedders seed the general population or mainly each other — a different
        question from how much they transmit.`;
      paint = () => bars(cv1(), {
        cats: SCENARIOS.map(sc => SCENARIO_LABEL[sc]),
        ylab: "% of events", pct: 1, height: 280,
        series: P2_CELLS.map(c => ({colour: CC[c],
          vals: SCENARIOS.map(sc => cellOf(sc, c))}))});
    }
  }

  if (metric === "intensity") {
    head = "Intensity per host";
    sub = "transmission events per host of that type";
    if (plot === "paired") {
      sub = "events per long shedder and per standard host";
      key = keyRow([[css("--long"), "per long shedder"],
                    [css("--std"), "per standard host"]]);
      cap = `Both per-capita rates with their own intervals, so the ratio can be
        read off its two parts rather than taken on trust. Separates
        &ldquo;each long shedder transmits more&rdquo; from &ldquo;there are
        more of them&rdquo;: the denominator is hosts, not events, which is why
        this metric is the figure candidate. control has no long shedders, so it
        shows the standard rate alone.`;
      paint = () => dots(cv1(), {
        cats: SCENARIOS.map(sc => SCENARIO_LABEL[sc]),
        ylab: "events per host", height: 280, points, zero: 1,
        fmtY: v => fmt(v, 2),
        series: [["long", css("--long"), "per long shedder"],
                 ["standard", css("--std"), "per standard host"]]
          .map(([h, colour, label]) => ({colour, label,
            vals: SCENARIOS.map(sc => (sOf(sc).intensity || {})[h] || null),
            runs: p2RunsOverScenarios(dist, "intensity." + h)}))});
    } else if (plot === "prev") {
      sub = "the long ÷ standard ratio against long-shedder prevalence";
      key = scKey;
      cap = `The per-capita ratio against prevalence. ${noPrev}A flat line is
        the claim: a per-host rate that does not track how many hosts there are.
        Contrast the same form under <em>Who transmits</em>, which climbs with
        prevalence by construction. ${P2_RATIO_NOTE}`;
      paint = () => scatterCI(cv1(), {
        xlab: "long-shedder prevalence (%)", ylab: "long ÷ standard",
        height: 280, points, fmtY: v => fmt(v, 2), fmtX: v => fmt(v, 1),
        groups: SCENARIOS.map(sc => ({colour: scenarioColour(sc),
          label: SCENARIO_LABEL[sc], x: prevOf(sc), y: p2Ratio(sOf(sc)),
          yruns: (distRuns(dist, sc, ["r"], "intensity_ratio") || [])[0] || []}))});
    } else {
      sub = "events per long shedder ÷ events per standard host";
      key = keyRow(SCENARIOS.map(sc => [scenarioColour(sc), SCENARIO_LABEL[sc]])
        .concat([[css("--ink-3"), "parity = 1"]]));
      cap = `One is parity; above one means each long shedder transmits more
        than each standard host. The grey point is parity itself, plotted so the
        y window always contains it — all five scenarios sit below it here.
        ${P2_RATIO_NOTE} This is the panel's figure candidate: ${P2_CONFOUND}
        <code>intensity</code> is per-capita, so it breaks that confound.`;
      paint = () => dots(cv1(), {
        cats: ["long ÷ standard"], ylab: "ratio", height: 280,
        points, fmtY: v => fmt(v, 2),
        series: SCENARIOS.map(sc => ({colour: scenarioColour(sc),
          label: SCENARIO_LABEL[sc], vals: [p2Ratio(sOf(sc))],
          runs: distRuns(dist, sc, ["r"], "intensity_ratio")}))
          .concat([{colour: css("--ink-3"), label: "parity = 1",
                    vals: [{n: 0, mean: 1, lo: null, hi: null}]}])});
    }
  }

  const rows = SCENARIOS.map(sc => { const s = sOf(sc); return `<tr>
    <td>${SCENARIO_LABEL[sc]}</td>
    <td class="n">${s.runs ?? 0}</td>
    <td class="n">${ci(asPct(s.who_transmits), 2, "%")}</td>
    <td class="n">${ci((s.intensity || {}).long, 2)}</td>
    <td class="n">${ci((s.intensity || {}).standard, 2)}</td>
    <td class="n">${ci(p2Ratio(s), 2)}</td></tr>`; }).join("");

  main.innerHTML = `
  <div class="note">${d.note}.</div>
  ${p2Plot(head, sub, body, key, cap)}
  <div class="plot"><div class="ph"><h3>Summary</h3>
    <span class="sub">all three metrics, whichever one is plotted</span></div>
    <table class="t"><thead><tr><th>scenario</th><th>runs</th>
      <th>caused by long shedders</th><th>events per long shedder</th>
      <th>events per standard host</th><th>long ÷ standard</th>
      </tr></thead><tbody>${rows}</tbody></table>
    <div class="cap">The last column has no interval, for the reason the ratio
      plot gives. control's per-long-shedder rate is absent because the scenario
      has no long shedders to divide by.</div></div>`;
  paint();

}


registerTab({
  id: "panel2", num: 5, title: "Panel 2", viewId: "viewP2",
  lede: "A <strong>transmission event</strong> is one infection: a transmitter, a recipient, a day and the lineage passed on. All three metrics come from one <strong>transmission table</strong> — per-day 2×2 counts of events by the host type on each side — counted by the host type at each end. Independent of <i>k</i> and of the clustering entirely.",
  rail: rail_panel2,
  render: render_panel2,
});
