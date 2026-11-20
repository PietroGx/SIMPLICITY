/* Panel 3: ported from panel3.html.

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_panel3 when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.

   The four fixed bar charts became one picked form, settings.plot3.

   Canvas ids carry a p3_ prefix. The shell builds every tab view in ONE
   document (shell.js, `tabs.forEach(viewElement)`), so a bare "c1" resolved to
   panel 2's canvas of the same name and this panel drew into it: panel 1's
   uniquely named canvases took 50/50/36 fillRect calls and panel 3 took zero.

   railDistribution greys its box from settings.plot3, so the rail has to be
   rebuilt when the form changes -- which settings.js's RAIL_SETTINGS and
   shell.js do. Rebuilding it from render() instead would rebuild it on every
   build-slider drag and take the handle's focus with it.
*/

const P3_PLOTS = [
  ["forest", "Dot and interval",
   "Both metrics as a mean with its bootstrap interval, the five scenarios on "
   + "one day axis. The standard form for a mean, and it does not imply a zero "
   + "the way a bar does."],
  ["rise_latency", "Rise against latency",
   "One point per scenario: the latent interval on x against the rise time on "
   + "y, each arm that axis's own interval. Puts the waiting time and the "
   + "sweep speed on one pair of axes."],
  ["across", "Line across j",
   "Both metrics over the composition classes, one line per scenario with its "
   + "interval as a ribbon. <i>j</i> is how many of a clade's <i>k</i> defining "
   + "mutations arose in a long shedder."],
  ["denominators", "The three counts",
   "Clades with a rise, clades with a latency, and clades truncated at the "
   + "burn-in day. Counts of clades, not estimates, so they carry no interval."],
];

function rail_panel3(rail) {
  rail.innerHTML = "";
  /* innerHTML alone empties it in a browser. The headless harness stores
     innerHTML without parsing (tests/artifacts/fakedom.js:1), so appended
     cards have to leave the child list too, or a rebuilt rail shows two of
     every control to check_app.js. */
  if (Array.isArray(rail.children)) rail.children.length = 0;
  rail.append(card("Plot", [
    pickControl({setting: "plot3", label: "Plot", options: P3_PLOTS}),
    railDistribution("dist3", settings.plot3),
  ]));
  rail.append(railClade());
  const r = document.createElement("div");
  r.className = "card";
  r.innerHTML = `<h2>Episode rules</h2><div class="hint">Set on the
  <a href="frequency.html">Frequency</a> page. Currently
  <code>win=${settings.win}d prom=${settings.prom} ${settings.lo}–${settings.hi} ${settings.window}</code>.</div>`;
  rail.append(r);
}

/* distRuns holds the scenario fixed, so a chart whose CATEGORIES are the
   scenarios asks once per scenario and takes the one column back. */
const p3_runs = (dist, sc, key) => (distRuns(dist, sc, [key], key) || [null])[0];

function p3_runsBySc(dist, key) {
  const out = SCENARIOS.map(sc => p3_runs(dist, sc, key));
  return out.some(r => r && r.length) ? out : undefined;
}

const p3_classLab = c => c === "wt" ? "wild type" : c === "other" ? "ndef≠k" : c;

async function render_panel3(main, rail){

  let d; try { d = await apiGetSlot("/api/panel?n=3&" + buildQuery()); }
  catch(e){ main.innerHTML=""; main.append(errorBox(e.message +
    " — build these rules on the Export page.")); return; }

  const plot = settings.plot3;
  const points = !!settings.dist3 && distAllowed(plot);
  const dist = points ? await distFor(3) : null;
  const S = sc => d.scenarios[sc] || {};
  const scLabs = SCENARIOS.map(sc => SCENARIO_LABEL[sc]);
  const scKey = keyRow(SCENARIOS.map(sc=>[scenarioColour(sc),SCENARIO_LABEL[sc]]));
  const classes = [...new Set(SCENARIOS.flatMap(sc =>
    Object.keys(S(sc).by_class_rise || {})))]
    .sort((a,b)=> (a==="wt")-(b==="wt") || (a==="other")-(b==="other") || (+a)-(+b));

  const tbl = SCENARIOS.map(sc => { const s = S(sc);
    return `<tr><td>${SCENARIO_LABEL[sc]} <small style="color:var(--ink-3)">${fmt(asPctNum(s.prevalence),1)}% · ${s.duration ?? "—"}d</small></td>
      <td class="n">${s.runs ?? 0}</td>
      <td class="n">${ci(s.rise_days,1," d")}</td>
      <td class="n">${s.n_with_rise ?? 0}</td>
      <td class="n">${ci(s.latent_days,1," d")}</td>
      <td class="n">${s.n_with_latency ?? 0}</td>
      <td class="n">${s.n_truncated ?? 0}</td>
      <td class="n">${s.n_censored ?? "—"}</td></tr>`; }).join("");

  const head = `<div class="note">${d.note}</div>`;
  const table = `
  <div class="plot"><div class="ph"><h3>Summary</h3></div>
    <table class="t"><thead><tr><th>scenario</th><th>runs</th><th>rise</th>
      <th>clades with a rise</th><th>latency</th><th>clades with a latency</th>
      <th>truncated at burn-in</th><th>still high at the end</th></tr></thead><tbody>${tbl}</tbody></table>
    <div class="cap">The truncation columns are why the two metrics cannot share
      a denominator: a clade already rising at the burn-in day keeps its rise
      time and loses its latency. <strong>Still high at the end</strong> counts
      emergence episodes the record cuts off before the clade comes back down —
      a sweep to fixation among them. Their rise is fully observed and is
      counted; only their peak is a lower bound. A dash means the build predates
      the flag.</div></div>`;

  if (plot === "rise_latency") {
    const groups = SCENARIOS.map(sc => ({colour: scenarioColour(sc),
      label: SCENARIO_LABEL[sc], x: S(sc).latent_days, y: S(sc).rise_days,
      xruns: p3_runs(dist, sc, "latent_days") || [],
      yruns: p3_runs(dist, sc, "rise_days") || []}));
    main.innerHTML = head + `
  <div class="plot"><div class="ph"><h3>Rise time against latent interval</h3>
    <span class="sub">one point per scenario · arms are each axis's own bootstrap CI of the mean</span></div>
    <canvas id="p3_c1"></canvas>
    ${scKey}
    <div class="cap"><strong>Measured:</strong> two means per scenario — days
      from first appearance to the start of the sweep on x, days from 10% to 90%
      of that episode's own amplitude on y — each with a bootstrap interval over
      runs.<br>
      <strong>Read:</strong> down and to the right would be a clade that waited
      longer <em>and then</em> swept faster, which is what chronic shedding
      predicts: more time inside one host accumulating divergence, then a
      lineage that arrives already fit. That is the reading, not the
      measurement. Five scenarios are five design points, not a trend, and the
      two axes do not share a denominator — x is over clades with a latency,
      y over clades with a rise (Summary below).</div></div>` + table;
    scatterCI(document.getElementById("p3_c1"), {groups, points, height: 300,
      xlab: "latent interval (days) · clades with a latency",
      ylab: "rise time (days) · clades with a rise"});

  } else if (plot === "across") {
    const mk = key => SCENARIOS.map(sc => ({colour: scenarioColour(sc),
      label: SCENARIO_LABEL[sc],
      vals: classes.map(c => (S(sc)[key] || {})[c]),
      runs: distRuns(dist, sc, classes, c => `${key}.${c}`)}));
    const levels = classes.map(p3_classLab);
    main.innerHTML = head + `
  <div class="plot"><div class="ph"><h3>Rise time by composition class</h3>
    <span class="sub">days from 10% to 90% of the episode's own amplitude</span></div>
    <canvas id="p3_c1"></canvas>
    ${scKey}
    <div class="cap">Classes run 0…${settings.k}, then <strong>ndef≠k</strong>
      (clades whose defining-mutation count is not <i>k</i>) and
      <strong>wild type</strong> (the founding clade), the order
      <code>panels.classes_for</code> defines. A class a scenario never produced
      has no point and the line breaks there rather than being interpolated
      across it.</div></div>
  <div class="plot"><div class="ph"><h3>Latent interval by composition class</h3>
    <span class="sub">days from first appearance to the start of the sweep</span></div>
    <canvas id="p3_c2"></canvas>
    <div class="cap">The waiting time is what chronic shedding would plausibly
      lengthen: more time inside one host accumulating divergence before the
      lineage ever competes in the population. Drawn over clades with a latency,
      a smaller set than the chart above — a truncated clade keeps its rise time
      and loses its latency.</div></div>` + table;
    linesOverLevels(document.getElementById("p3_c1"), {levels, points,
      series: mk("by_class_rise"), ribbon: true, height: 260,
      xlab: `j long of ${settings.k}`, ylab: "days"});
    linesOverLevels(document.getElementById("p3_c2"), {levels, points,
      series: mk("by_class_latent"), ribbon: true, height: 260,
      xlab: `j long of ${settings.k}`, ylab: "days"});

  } else if (plot === "denominators") {
    const series = [
      {colour: css("--accent"), label: "clades with a rise",
       vals: SCENARIOS.map(sc => S(sc).n_with_rise)},
      {colour: css("--wt"), label: "clades with a latency",
       vals: SCENARIOS.map(sc => S(sc).n_with_latency)},
      {colour: css("--warn"), label: "truncated at the burn-in day",
       vals: SCENARIOS.map(sc => S(sc).n_truncated)},
      {colour: css("--ink-3"), label: "still high at the end",
       vals: SCENARIOS.map(sc => S(sc).n_censored ?? null)},
    ];
    main.innerHTML = head + `
  <div class="plot"><div class="ph"><h3>The four counts</h3>
    <span class="sub">clades, pooled over the runs of each scenario</span></div>
    <canvas id="p3_c1"></canvas>
    ${chartKey(series)}
    <div class="cap"><strong>Clades with a rise</strong>: clades whose emergence
      episode yielded a 10–90% rise.
      <strong>Clades with a latency</strong>: those that also have a measurable
      wait before that sweep.
      <strong>Truncated at the burn-in day</strong>: designated clades already
      under way at the burn-in day, so their first appearance is not in the
      window — the latency has no start to measure from, and the clade is
      dropped from the latency only. It keeps its rise time, which is measured
      inside the episode and needs no origin. That is why the first two counts
      differ, and why this panel has to state both rather than let one stand in
      for the other.<br>
      These are <strong>counts of clades</strong>, not estimates over runs:
      there is nothing to bootstrap, so no bar carries an interval and the
      distribution overlay is greyed out. The <i>n</i> behind the two intervals
      elsewhere on this tab is a count of <strong>runs</strong> (the Summary's
      runs column), and it differs between the two metrics as well.</div></div>`
      + table;
    countBars(document.getElementById("p3_c1"), {cats: scLabs, series,
      ylab: "clades", height: 260});

  } else {
    const series = [
      {colour: css("--accent"), label: "rise time · clades with a rise",
       vals: SCENARIOS.map(sc => S(sc).rise_days),
       runs: p3_runsBySc(dist, "rise_days")},
      {colour: css("--wt"), label: "latent interval · clades with a latency",
       vals: SCENARIOS.map(sc => S(sc).latent_days),
       runs: p3_runsBySc(dist, "latent_days")},
    ];
    main.innerHTML = head + `
  <div class="plot"><div class="ph"><h3>Rise time and latent interval</h3>
    <span class="sub">days · mean over runs, whiskers a bootstrap CI of the mean</span></div>
    <canvas id="p3_c1"></canvas>
    ${chartKey(series)}
    <div class="cap">A day count has no meaningful zero — neither quantity is a
      share of anything — so the axis does not start at one, and a bar, whose
      length is read from zero, would imply it did. A dot with its interval puts
      the difference between two scenarios in the gap between two dots instead
      of in the difference between two bar lengths.<br>
      <strong>Rise time</strong> is days from 10% to 90% of the emergence
      episode's own amplitude. The <strong>latent interval</strong> is days from
      the clade's first appearance to the start of that sweep — the waiting time
      is what chronic shedding would plausibly lengthen: more time inside one
      host accumulating divergence before the lineage ever competes in the
      population. The two series do not share a denominator: rise is over clades
      with a rise, latency over clades with a latency, and the Summary below
      gives both per scenario.</div></div>` + table;
    dots(document.getElementById("p3_c1"), {cats: scLabs, series, points,
      ylab: "days", height: 260});
  }

}


registerTab({
  id: "panel3", num: 6, title: "Panel 3", viewId: "viewP3",
  lede: "How long a clade lurked, then how fast it went. The <strong>emergence episode</strong> is the one containing the clade's highest frequency; its rise is measured from the dip immediately before the sweep, which is what makes it the speed of a growth event rather than a disguised lifetime.",
  rail: rail_panel3,
  render: render_panel3,
});
