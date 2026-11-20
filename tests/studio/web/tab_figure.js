/* Figure: ported from figure.html.

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_figure when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.
*/

function rail_figure(rail) {
rail.append(railClade());
const note = document.createElement("div");
note.className = "card";
note.innerHTML = `<h2>Read-only</h2><div class="hint">This page composes; it
  decides nothing. Change <i>k</i> or <i>m</i><sub>T</sub> on
  <a href="clustering.html">Clustering</a>, the episode rules on
  <a href="frequency.html">Frequency</a>.</div>`;
rail.append(note);
}

async function render_figure(main){

  let p1, p2, p3;
  try {
    [p1, p2, p3] = await Promise.all([
      apiGetSlot("/api/panel?n=1&" + buildQuery()), apiGetSlot("/api/panel?n=2&" + buildQuery()),
      apiGetSlot("/api/panel?n=3&" + buildQuery())]);
  } catch(e){ main.innerHTML=""; main.append(errorBox(e.message +
    " — build the current settings on the Export page.")); return; }
  main.innerHTML = `
  <div class="grid2">
    <div class="plot"><div class="ph"><h3>A — composition against burden</h3>
      <span class="sub">% of run burden</span></div><canvas id="a"></canvas>
      ${keyRow(SCENARIOS.map(sc=>[scenarioColour(sc),SCENARIO_LABEL[sc]]))}</div>
    <div class="plot"><div class="ph"><h3>B — who transmits</h3>
      <span class="sub">% of events caused by a long shedder</span></div>
      <canvas id="b"></canvas></div>
    <div class="plot"><div class="ph"><h3>C — emergence speed and latency</h3>
      <span class="sub">days</span></div><canvas id="c"></canvas>
      <div class="cap">Left bar of each pair is the rise, right is the latent
        interval. They have different denominators.</div></div>
    <div class="plot"><div class="ph"><h3>D — sensitivity</h3>
      <span class="sub">builds cached</span></div>
      <div id="d" class="cap">…</div></div>
  </div>
  <div class="note">k=${settings.k} · m_T=${mtPct(settings.mt)}% · win=${settings.win}d · prom=${settings.prom} ·
    ${settings.window === "post" ? "post burn-in" : "whole record"} ·
    ${settings.burnin ? "burn-in excluded from counts" : "whole run counted"}</div>`;
  bars(document.getElementById("a"), {cats:p1.classes.map(c=>c==="wt"?"wt":c==="other"?"≠k":c),
    ylab:"% burden", pct:1, height:210,
    series:SCENARIOS.map(sc=>({colour:scenarioColour(sc),
      vals:p1.classes.map(c=>asPct((((p1.scenarios[sc]||{}).classes||{})[c]||{}).frac_burden))}))});
  bars(document.getElementById("b"), {cats:SCENARIOS.map(s=>SCENARIO_LABEL[s]),
    ylab:"%", pct:1, height:210,
    series:[{colour:css("--long"),
      vals:SCENARIOS.map(sc=>asPct((p2.scenarios[sc]||{}).who_transmits))}]});
  bars(document.getElementById("c"), {cats:SCENARIOS.map(s=>SCENARIO_LABEL[s]),
    ylab:"days", height:210,
    series:[{colour:css("--accent"),vals:SCENARIOS.map(sc=>(p3.scenarios[sc]||{}).rise_days)},
            {colour:css("--wt"),    vals:SCENARIOS.map(sc=>(p3.scenarios[sc]||{}).latent_days)}]});
  const ix = await apiGetSlot("/api/index");
  document.getElementById("d").innerHTML = (ix.builds||[]).length > 1
    ? `${(ix.builds||[]).length} builds cached — see <a href="panel4.html">Panel 4</a>.`
    : `Only one clade build is cached, so there is nothing to compare yet.
       Build a second setting on the <a href="export.html">Export</a> page.`;

}


registerTab({
  id: "figure", num: 8, title: "Figure", viewId: "viewFig",
  lede: "The four panels as they would be assembled, read-only. Everything here is set on the pages upstream.",
  rail: rail_figure,
  render: render_figure,
});
