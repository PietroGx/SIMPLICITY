/* Panel 4: ported from panel4.html.

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_panel4 when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.

   The sensitivity lines read /api/sens, never this payload's `by_axis`:
   P4_AXES is ["k","win","prom"] (panels.py:442), so _p4_by_axis (:503-523)
   leaves m_T unpinned and by_axis["k"] comes back with 21 levels for three
   values of k. See p4SensNote().
*/

/* Which build axis the sensitivity lines walk. */
const P4_AXIS_OPTS = [
  ["k", "shared-mutation threshold k",
   "How far a lineage must have diverged from its clade root to found a clade. "
   + "It moves the partition AND the class labels, so it is the one axis that "
   + "changes what the majority-long classes are."],
  ["mt", "minimum transmissions m_T",
   "The share of the simulation's transmissions a lineage's subtree must carry to "
   + "found a clade. A refused lineage folds into its ancestor, so this merges "
   + "clades rather than relabelling them."],
  ["win", "smoothing window",
   "The centred moving average the frequency track is read through. Wider "
   + "merges neighbouring episodes, which lengthens a rise by joining two."],
  ["prom", "peak prominence",
   "How far a peak must stand above its neighbouring dip, as a fraction of "
   + "the clade's own height. It decides what counts as an episode at all."],
];
const P4_AXIS_NAME = {k: "k", mt: "m_T", win: "smoothing window",
                      prom: "prominence"};
const p4Level = (ax, v) =>
  ax === "mt" ? mtPct(v) + "%" : ax === "win" ? fmt(v, 0) + " d"
  : ax === "prom" ? fmt(v, 2) : String(v);
/* the one share among the metrics: it arrives as a fraction and is shown as a
   percentage, so the server's "(fraction)" label is restated for display */
const P4_SHARE = "frac_burden_majority_long";
const p4Val = (m, v) => m === P4_SHARE ? asPctNum(v) : v;
const p4Label = (M, m) => String(M[m] || m).replace("(fraction)", "(%)");
/* m_T off a build's directory name, which is the build's own identity --
   by_axis reports a level per row without saying which m_T that row is. */
const p4KeyMt = key => (String(key).match(/_m([0-9.]+)_/) || [])[1];

function rail_panel4(rail) {
rail.append(railFocus());
rail.append(railWindow());
rail.append(card("Sensitivity", [
  pickControl({setting: "axis4", label: "Axis to vary",
    options: P4_AXIS_OPTS,
    hint: "The other three build parameters are pinned at the operating point "
        + "(k=7, m_T=1%, window 15 d, prominence 0.10) and named on the chart. "
        + "One at a time, not a marginal over the grid."}),
]));
}

/* Why the lines are missing, with the numbers. Shown in place of the chart so
   nothing silently draws by_axis against a duplicated axis. */
function p4SensNote(ax) {
  return `The three sensitivity lines need <code>/api/sens</code>, which is
    not answering yet — <code>sensFor("${ax}")</code> in
    <code>web/dist.js</code> returned null. This panel's own
    <code>by_axis</code> cannot stand in: <code>P4_AXES = ["k","win","prom"]</code>
    (<code>panels.py:442</code>) omits m_T, so <code>_p4_by_axis</code>
    (<code>panels.py:503-523</code>) leaves it unpinned and
    <code>by_axis["k"]</code> carries 21 levels for three values of k — seven
    m_T rows per k, ordered lexicographically by directory name, so m_T=0 lands
    after m_T=0.1. At k=5, HIV high, mean emergence rise, those seven rows read
    20.9, 22.6, 23.8, 24.3, 25.9, 28.2, 10.7 days: rise climbing with m_T and
    then the 10.7 of m_T=0, not instability in k. The tables below are that
    grouping, each column named with the m_T it is actually at.
    <code>/api/sens</code> pins every build axis but the one asked for, sorts
    the levels numerically, and carries the interval per level, which
    <code>_panel4_row</code> (<code>panels.py:388-393</code>) drops.`;
}

async function render_panel4(main, rail, {idp = "p4", compact = false} = {}){

  let d; try { d = await apiGetSlot("/api/panel?n=4&" + buildQuery()); }
  catch(e){ main.innerHTML=""; main.append(errorBox(e.message)); return; }
  if (!d.rows.length){ main.innerHTML=""; main.append(errorBox(
    "no clade builds cached. Build at least two settings on the Export page.")); return; }
  const ax = settings.axis4;
  const sens = await sensFor(ax);
  const haveSens = !!(sens && sens.axis === ax && (sens.levels || []).length);
  const head = `<tr><th>build</th>${SCENARIOS.map(s=>`<th>${SCENARIO_LABEL[s]}</th>`).join("")}</tr>`;
  const body = m => d.rows.map(r => {
    if (r.error) return `<tr><td class="n">${r.key}</td><td colspan="5">${r.error}</td></tr>`;
    const p = r.params;
    return `<tr><td class="n">k=${p.k} m_T=${mtPct(p.mt)}% win=${p.win} prom=${p.prom} ${p.window}</td>`
      + SCENARIOS.map(s => `<td class="n">${fmt(p4Val(m, (r.scenarios[s]||{})[m]), 1)}</td>`).join("")
      + `</tr>`;
  }).join("");
  /* The surface, before the per-build table: the span over the whole grid is
     the honest width of a claim, and the one-at-a-time rows say which axis
     moves it. A number quoted at one setting means little without both. */
  const M = d.metrics || {};
  const sumRows = Object.keys(M).map(m => {
    const cells = SCENARIOS.map(s => {
      const v = ((d.summary || {})[s] || {})[m];
      if (!v || v.n === 0) return `<td class="n">—</td>`;
      const at = v.at_operating_point;
      return `<td class="n">${at === null || at === undefined ? "—" : fmt(p4Val(m, at), 1)}`
        + `<br><small>${fmt(p4Val(m, v.min), 1)}–${fmt(p4Val(m, v.max), 1)} over ${v.n}</small></td>`;
    }).join("");
    return `<tr><td class="n">${p4Label(M, m)}</td>${cells}</tr>`;
  }).join("");

  /* ---------------------------------------------- the sensitivity facets */

  const mKeys = Object.keys(M);
  const sensBlock = () => {
    const nm = P4_AXIS_NAME[ax] || ax;
    if (!haveSens) return `<div class="plot"><div class="ph">
      <h3>Sensitivity over ${nm}</h3>
      <span class="sub">not drawn</span></div>
      <div class="cap">${p4SensNote(ax)}</div></div>`;
    const pin = Object.entries(sens.pinned || {})
      .map(([k2, v]) => `${P4_AXIS_NAME[k2] || k2}=${p4Level(k2, v)}`).join(" · ");
    const facets = mKeys.map((m, i) => `<div class="ph"><h3>${p4Label(M, m)}</h3>
      <span class="sub">${i === mKeys.length - 1 ? nm + " along the bottom" : ""}</span></div>
      <canvas id="${idp}s${i + 1}"></canvas>`).join("");
    return `<div class="plot"><div class="ph">
      <h3>Sensitivity over ${nm}</h3>
      <span class="sub">${pin ? pin + " pinned" : "nothing pinned"}
        · ${sens.levels.length} level${sens.levels.length === 1 ? "" : "s"}</span></div>
      ${facets}
      ${keyRow(SCENARIOS.map(s => [scenarioColour(s), SCENARIO_LABEL[s]]))}
      <div class="cap">Three facets over one x axis, not fifteen lines on one:
        burden is a percentage and the other two are days, so a shared y would
        be meaningless. The dashed rule is the operating point. The band is the
        bootstrap CI of the mean over simulations — without it a two-day wiggle
        invites being read as a real shift, which is the whole reason the
        intervals matter here. The other build parameters are pinned at the
        operating point (${pin || "—"}), so this is a one-at-a-time
        sensitivity, not a marginal over the grid. The two day facets are
        scaled to their own range; only the percentage facet starts at zero.
      </div></div>`;
  };

  const drawSens = () => {
    if (!haveSens) return;
    const levels = sens.levels.map(v => p4Level(ax, v));
    /* The operating point is compute.DEFAULTS; SETTING_DEFAULTS is the
       browser's one copy of it (settings.js:30). Compared as formatted
       levels, so 15 and 15.0 are the same stop. */
    const mark = p4Level(ax, SETTING_DEFAULTS[ax]);
    mKeys.forEach((m, i) => {
      const cv = document.getElementById(idp + "s" + (i + 1));
      if (!cv) return;
      const pct = m === P4_SHARE;
      linesOverLevels(cv, {
        levels, mark, ribbon: true, height: 190, zero: pct, pct,
        xlab: i === mKeys.length - 1 ? (P4_AXIS_NAME[ax] || ax) : "",
        ylab: pct ? "% of simulation burden, majority-long classes"
          : m === "rise_days" ? "mean rise time (days)" : "mean latent interval (days)",
        fmtY: v => fmt(v, 0),
        series: SCENARIOS.map(s => ({
          colour: scenarioColour(s), label: SCENARIO_LABEL[s],
          vals: (((sens.scenarios || {})[s] || {})[m] || []).map(b => pct ? asPct(b) : b),
        })),
      });
    });
  };

  /* --------------------------------------------- the by_axis tables, kept */

  const axisBlock = ax2 => {
    const a = (d.by_axis || {})[ax2];
    if (!a || !a.levels.length) return "";
    const pinned = Object.entries(a.pinned || {})
      .map(([k2, v]) => `${k2}=${v}`).join(", ");
    const col = i => {
      const mt = p4KeyMt((a.keys || [])[i]);
      return `${ax2}=${a.levels[i]}${mt === undefined ? ""
        : `<br><small>m_T=${mtPct(mt)}%</small>`}`;
    };
    const rows = Object.keys(M).map(m => {
      const cells = a.levels.map((lv, i) => {
        const per = SCENARIOS.map(s =>
          fmt(p4Val(m, ((a.scenarios[s] || {})[m] || [])[i]), 1)).join(" / ");
        return `<td class="n">${per}</td>`;
      }).join("");
      return `<tr><td class="n">${p4Label(M, m)}</td>${cells}</tr>`;
    }).join("");
    return `<div class="plot"><div class="ph">
      <h3>Moving ${ax2} alone</h3><span class="sub">${pinned} held at the
      operating point · each cell is ${SCENARIOS.map(s =>
        SCENARIO_LABEL[s]).join(" / ")}</span></div>
      <table class="t"><thead><tr><th>metric</th>${a.levels.map((lv, i) =>
        `<th>${col(i)}</th>`).join("")}</tr></thead>
      <tbody>${rows}</tbody></table>
      <div class="cap">m_T is not one of the axes this grouping pins
        (<code>panels.py:442</code>), so ${ax2}=${a.levels[0]} repeats once per
        m_T value on disk and the columns are ordered by directory name, not
        numerically. Read the second line of each header before comparing
        columns. The chart above takes its levels from
        <code>/api/sens</code> instead.</div></div>`;
  };

  const opMt = p4KeyMt(d.operating_point);

  if (compact) { main.innerHTML = sensBlock().replace(/<div class="cap">[\s\S]*?<\/div>/g, ""); drawSens(); return; }
  main.innerHTML = `
  <div class="note">${simWords(d.note)}</div>
  ${sensBlock()}
  <div class="plot"><div class="ph"><h3>At the operating point, and across the grid</h3>
    <span class="sub">${d.operating_point
      ? "operating point k=7, m_T=1%, win=15, prom=0.1"
      : "the operating-point build is NOT on disk — only the range is shown"}
      · grid k ${(d.axes.k || []).join("/")} · win ${(d.axes.win || []).join("/")}
      · prom ${(d.axes.prom || []).join("/")}</span></div>
    <table class="t"><thead><tr><th>metric</th>${SCENARIOS.map(s =>
      `<th>${SCENARIO_LABEL[s]}</th>`).join("")}</tr></thead>
    <tbody>${sumRows}</tbody></table>
    <div class="cap">The big number is the operating point; beneath it the
      range over every build on disk, and how many builds that is. A range
      wider than the difference between scenarios means the grid, not the
      biology, is doing the talking.${opMt === undefined ? "" :
      ` <strong>The big number is at m_T=${mtPct(opMt)}%</strong>, not
      necessarily the 1% of the operating point:
      <code>_p4_operating_key</code> (<code>panels.py:458-465</code>) matches on
      k, win and prom only, so it returns whichever m_T sorts first —
      <code>${d.operating_point}</code>.`}</div></div>
  ${["k", "win", "prom"].map(axisBlock).join("")}
  <div class="plot"><div class="ph"><h3>Figure 4 — burden in the majority-long classes</h3>
    <span class="sub">% of the simulation's burden</span></div>
    <table class="t"><thead>${head}</thead><tbody>${body(P4_SHARE)}</tbody></table>
    <div class="cap">Majority-long means more than half of a clade's defining
      mutations arose in a long shedder, so the classes summed depend on k.</div></div>
  <div class="plot"><div class="ph"><h3>Figure 5 — emergence rise time</h3>
    <span class="sub">days</span></div>
    <table class="t"><thead>${head}</thead><tbody>${body("rise_days")}</tbody></table></div>
  <div class="plot"><div class="ph"><h3>Figure 5 — latent interval</h3>
    <span class="sub">days</span></div>
    <table class="t"><thead>${head}</thead><tbody>${body("latent_days")}</tbody></table></div>
  <div class="plot"><div class="ph"><h3>Simulations contributing</h3></div>
    <table class="t"><thead>${head}</thead><tbody>${body("runs")}</tbody></table>
    <div class="cap">A row with fewer simulations is not a weaker effect, it is a
      different sample — check this before comparing rows.</div></div>`;
  drawSens();

}


registerTab({
  id: "panel4", group: "supp", num: "S1", title: "Clustering sensitivity", viewId: "viewP4",
  lede: "Does anything on Figures 4 and 5 survive changing the clade definition or the episode rules? One row per cached clade build. This is not an independent choice — it is those figures, recomputed per build.",
  rail: rail_panel4,
  render: render_panel4,
});
