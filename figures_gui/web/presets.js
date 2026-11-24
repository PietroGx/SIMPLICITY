/* Saved data settings: one named set of every rule the exploration tabs decide
   (clade definition, burn-in, breakout detection), stored by server.py in
   Data/figures_studio/detection_presets.json so it outlives the browser. Every
   exploration tab loads and saves the same sets; the figure tabs only read them. */

const DATA_PRESET_KEYS = ["k", "mt", "originRule", "burnin",
  "sdKernel", "sdSpace", "sdWin", "sdProm", "sdWidth", "sdCens", "sdStart",
  "sdTrough", "sdFrac", "sdAnchor", "sdPick", "sdGate", "brkF", "brkH", "brkQuiet"];
const DATA_TABS = ["clustering", "across", "sweep", "cryptic", "crystats"];
const presetEsc = v => String(v).replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"})[c]);

let dataPresets = null;
async function presetsLoad() {
  if (dataPresets) return dataPresets;
  try { const r = await fetch("/api/presets", {cache: "no-store"}); dataPresets = r.ok ? await r.json() : {}; }
  catch (e) { dataPresets = {}; }
  return dataPresets;
}
async function presetsPost(name, values) {
  const r = await fetch("/api/presets", {method: "POST", headers: {"Content-Type": "application/json"},
    body: JSON.stringify({name, settings: values})});
  const d = await r.json();
  if (!r.ok) throw new Error(d.error || `HTTP ${r.status}`);
  dataPresets = d;
}
const presetNow = () => Object.fromEntries(DATA_PRESET_KEYS.map(k => [k, settings[k]]));
/* a set saved before it held every key matches on the keys it has */
const presetMatch = () => Object.keys(dataPresets || {}).find(n => {
  const p = dataPresets[n], ks = DATA_PRESET_KEYS.filter(k => k in p);
  return ks.length && ks.every(k => String(p[k]) === String(settings[k]));
});
/* every exploration tab's rail shows these values, so all of them rebuild */
function presetApply(patch) {
  setSettings(patch);
  DATA_TABS.forEach(id => { const t = tabById(id); if (t) t.railBuilt = false; });
  const t = tabById(activeTabId); if (t) renderTab(t);
}
function presetBarHtml(P) {
  const names = Object.keys(dataPresets || {}).sort(), m = presetMatch();
  return `<span class="sw-presets"><select id="${P}_pload" title="Load a saved set of data settings">
      <option value="">Load saved…</option>${names.map(n =>
        `<option value="${presetEsc(n)}"${n === m ? " selected" : ""}>${presetEsc(n)}</option>`).join("")}</select>
    <button type="button" class="sw-reset" id="${P}_psave">Save…</button>${m
      ? `<button type="button" class="sw-reset" id="${P}_pdel">Delete</button>` : ""}</span>`;
}
function presetWire(P, rerender) {
  const el = id => document.getElementById(`${P}_${id}`);
  const fail = e => alert(`Could not save: ${e.message}`);
  if (el("pload")) el("pload").onchange = ev => {
    const n = ev.target.value;
    if (n && dataPresets[n]) presetApply(dataPresets[n]);
  };
  if (el("psave")) el("psave").onclick = () => {
    const n = (prompt("Save these data settings as:", presetMatch() || "") || "").trim();
    if (!n || (dataPresets[n] && !confirm(`Replace the saved "${n}"?`))) return;
    presetsPost(n, presetNow()).then(rerender, fail);
  };
  if (el("pdel")) el("pdel").onclick = () => {
    const n = presetMatch();
    if (n && confirm(`Delete the saved "${n}"?`)) presetsPost(n, null).then(rerender, fail);
  };
}
/* the rail card for the tabs whose rule line has no room for the bar */
function presetCard(P) {
  const el = document.createElement("div");
  el.className = "card";
  const paint = () => {
    const m = presetMatch();
    el.innerHTML = `<h2>Saved settings</h2>${presetBarHtml(P)}
      <div class="hint" style="display:block">${m ? `Matches <b>${presetEsc(m)}</b>. ` : ""}A saved set holds
        the clade definition, the burn-in and the breakout rules: every exploration tab
        loads and saves the same sets, and the figures read whatever is in force.</div>`;
    presetWire(P, () => { paint(); });
  };
  presetsLoad().then(paint);
  el._paint = paint;
  return el;
}
