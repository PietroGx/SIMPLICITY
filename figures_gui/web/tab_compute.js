/* Compute: ported from export.html (the old Export tab; figure export is tab_export.js now).

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_export when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.
*/

function rail_compute(rail) {
  rail.append(railClade());
}

/* While a job runs the tab re-reads /api/jobs so the log grows on screen. */
let exportPollTimer = null;

async function render_compute(main){

  let d, j;
  try { d = await apiGetSlot("/api/index"); j = await apiGetSlot("/api/jobs"); }
  catch(e){ main.innerHTML=""; main.append(errorBox(e.message)); return; }
  const running = !!j.running;
  const builds = (d.builds||[]).map(b => {
    const p = b.params || {};
    const cur = p.k==settings.k && p.win==settings.win && Math.abs(p.prom-settings.prom)<1e-9
             && p.window==settings.window;
    return `<tr${cur?' class="tot"':''}><td class="n">${b.key}</td>
      <td class="n">${Object.values(b.counts||{}).reduce((a,b)=>a+b,0)}</td>
      <td class="n">${(b.failed||[]).length}</td>
      <td>${cur?"current settings":""}</td></tr>`;
  }).join("") || `<tr><td colspan="4">nothing built yet</td></tr>`;
  const tx = d.transmissions;
  main.innerHTML = `
  <div class="stats">
    <div class="stat"><div class="k">clade builds</div><div class="v">${(d.builds||[]).length}</div></div>
    <div class="stat"><div class="k">transmissions</div><div class="v">${tx?Object.values(tx.counts||{}).reduce((a,b)=>a+b,0):0}<small> simulations</small></div></div>
    <div class="stat"><div class="k">code hash</div><div class="v" style="font-size:.8rem">${d.code}</div></div>
    <div class="stat"><div class="k">job</div><div class="v" style="font-size:.8rem">${running?"running":"idle"}</div></div>
  </div>
  <div class="note">A cache key carries the parameters <strong>and</strong> a hash
    of the computing code, so a fix to the maths cannot be served from a stale
    build.</div>
  <div class="plot"><div class="ph"><h3>Run a pass</h3>
    <span class="sub">${running ? j.running.name : "idle"}</span></div>
    <div class="grid2">
      <button class="wide prim" id="b3" ${running?"disabled":""}>
        Build clades at k=${settings.k}</button>
      <button class="wide" id="btx" ${running?"disabled":""}>
        Rebuild transmissions — Figure 3, ~10 min</button></div>
    <div class="cap">The clade pass reads the daily record and clusters with the
      pipeline's own function; expect roughly ten minutes for 250 simulations at two
      processes. A second request is refused, not queued.</div></div>
  <div class="plot"><div class="ph"><h3>Job log</h3></div>
    <pre class="log" id="log">${(j.tail||[]).join("\n") || "no output"}</pre></div>
  <div class="plot"><div class="ph"><h3>Cached clade builds</h3></div>
    <table class="t"><thead><tr><th>key</th><th>simulations</th><th>failed</th><th></th>
      </tr></thead><tbody>${builds}</tbody></table></div>`;
  document.getElementById("b3").onclick = async () => {
    /* every build parameter, taken from BUILD_SETTINGS rather than listed
       here: m_T joined that list in v2.5.5 and this call did not notice, so
       the button rebuilt a different build than the one the page was reading. */
    const body = {procs: 2};
    for (const n of BUILD_SETTINGS) body[n] = settings[n];
    const r = await post("/api/rebuild/clades", body);
    if (r.error || r.started === false) alert(r.error || r.msg);
    apiInvalidate();
    renderTab(tabById("compute"));
  };
  document.getElementById("btx").onclick = async () => {
    const r = await post("/api/rebuild/transmissions", {procs:2});
    if (r.error || r.started === false) alert(r.error || r.msg);
    apiInvalidate();
    renderTab(tabById("compute"));
  };
  clearTimeout(exportPollTimer);
  if (running) exportPollTimer = setTimeout(() => {
    apiInvalidate("/api/jobs");
    renderTab(tabById("compute"));
  }, 3000);

}


registerTab({
  id: "compute", group: "tools", title: "Compute", viewId: "viewExp",
  lede: "Runs a compute pass with the current settings and caches the result. These endpoints are loopback-only: a pass takes two of four cores for minutes, and this box also serves production vhosts.",
  rail: rail_compute,
  render: render_compute,
});
