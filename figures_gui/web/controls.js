/* Rail controls, built from three factories instead of repeated markup.

   Every control in the studio is a slider, a segmented choice or a checkbox,
   and each one does the same four things: show a label, show its current
   value, write the setting, and carry a hint explaining what it does. The
   artifact wrote that markup out once per control -- about fifteen times. Here
   the factories own it, so a change to how a control looks or how it writes a
   setting happens in one place.

   Every factory writes through `setSetting`, so no control keeps private state
   and the shell's render loop hears every change the same way.
*/

/* Put the current value onto `values`, the grid's own axis for this setting.

   A saved setting can sit between the stops: the rail used to let prom land
   anywhere on 0.01..0.50 against a grid of [0.10, 0.20], and the handle then
   showed index 0 while the readout said 0.50. Nearest wins, and the setting is
   written, so what the rail shows is what the server will be asked for. */
function snapToAxis(setting, values) {
  const cur = settings[setting];
  if (!values || !values.length || values.includes(cur)) return cur;
  const v = values.every(x => typeof x === "number") && typeof cur === "number"
    ? values.reduce((a, b) => Math.abs(b - cur) < Math.abs(a - cur) ? b : a)
    : values[0];
  setSetting(setting, v);
  return v;
}

/* Every build parameter onto the grid, before the first render.

   Snapping when a control is BUILT is too late: the rails are built lazily, one
   per tab, so a saved prom=0.5 would survive until the tab owning that slider
   was opened -- and until then every panel asks for a build that is not there. */
function snapBuildSettings() {
  for (const n of BUILD_SETTINGS) snapToAxis(n, AXES[n]);
}

/* A labelled range slider. `values` makes it a slider over a fixed list -- for
   every build parameter, the grid's own axis, so the slider cannot stop on a
   value that has no build. An axis of one value is shown and disabled; an axis
   the server did not report leaves the control dead rather than inventing a
   list, which is how the k slider came to offer [0,3,5,7,10,15,20] against a
   grid of [5,7,10]. */
function sliderControl({setting, label, hint, min, max, step, values, format}) {
  const show = format || (v => String(v));
  const axis = values && values.length ? values : null;
  const dead = !!values && !axis;
  const current = axis ? snapToAxis(setting, axis) : settings[setting];
  const idx = axis ? axis.indexOf(current) : current;
  const off = dead || (axis && axis.length < 2);
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    <div class="lab"><label>${label}</label>
      <span class="val">${show(current)}</span></div>
    <input type="range" min="${axis ? 0 : min}" max="${axis ? axis.length - 1 : max}"
           step="${axis ? 1 : (step || 1)}" value="${idx}"${off ? " disabled" : ""}>
    ${dead ? `<div class="hint err">The server did not report this axis of the
       build grid, so there is nothing to choose from.</div>` : ""}
    ${hint ? `<div class="hint">${hint}</div>` : ""}`;
  const input = el.querySelector("input");
  const readout = el.querySelector(".val");
  input.addEventListener("input", () => {
    const v = axis ? axis[+input.value] : +input.value;
    readout.textContent = show(v);
    setSetting(setting, v);
  });
  return el;
}

/* A segmented button group. `options` is [[value, label], ...]. `axis`, where
   the setting names a build, keeps only the options the grid was built over. */
function segmentControl({setting, label, hint, options, axis}) {
  if (axis && axis.length) {
    options = options.filter(([v]) => axis.includes(v));
    snapToAxis(setting, axis);
  }
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    ${label ? `<div class="lab"><span>${label}</span></div>` : ""}
    <div class="seg">${options.map(([v, t, h]) =>
      `<button type="button" data-v="${v}" aria-pressed="${settings[setting] === v}"${
        h ? ` data-hint="${String(h).replace(/"/g, "&quot;")}"` : ""}>${t}</button>`
      ).join("")}</div>
    ${hint ? `<div class="hint">${hint}</div>` : ""}`;
  el.querySelectorAll("button").forEach(b => b.addEventListener("click", () => {
    const raw = b.dataset.v;
    const v = isNaN(Number(raw)) || raw === "" ? raw : Number(raw);
    setSetting(setting, v);
    el.querySelectorAll("button").forEach(o =>
      o.setAttribute("aria-pressed", String(o === b)));
  }));
  return el;
}

/* `disabled` is for a box that CANNOT apply to what is selected -- the
   distribution overlay under a stacked bar chart, which has no per-run form.
   Shown and greyed, with `why` saying so, rather than quietly doing nothing. */
function checkControl({setting, label, hint, disabled, why}) {
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    <label class="chk"><input type="checkbox" ${settings[setting] ? "checked" : ""}${disabled ? " disabled" : ""}>
      <span>${label}</span></label>
    ${disabled && why ? `<div class="hint">${why}</div>` : ""}
    ${hint ? `<div class="hint">${hint}</div>` : ""}`;
  el.querySelector("input").addEventListener("change", ev =>
    setSetting(setting, ev.target.checked ? 1 : 0));
  return el;
}

/* A choice whose options need explaining: the metric a panel measures, or the
   form it is plotted in. `options` is [[value, label, explain], ...], and the
   explanation of whatever is selected is shown under the control and swapped
   on change -- the alternative is a hint listing all five at once, which is
   how a rail stops being read.

   An option list filtered down to one is shown and disabled, so a metric that
   admits a single form says so instead of looking like a free choice. */
function pickControl({setting, label, options, hint}) {
  const opts = options.filter(Boolean);
  const vals = opts.map(o => o[0]);
  if (vals.length && !vals.includes(settings[setting])) setSetting(setting, vals[0]);
  const explain = v => (opts.find(o => o[0] === v) || [])[2] || "";
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    <label>${label}</label>
    <select${opts.length < 2 ? " disabled" : ""}>${opts.map(([v, t]) =>
      `<option value="${v}"${settings[setting] === v ? " selected" : ""}>${t}</option>`
      ).join("")}</select>
    <div class="hint" data-role="explain">${explain(settings[setting])}</div>
    ${hint ? `<div class="hint">${hint}</div>` : ""}`;
  const box = el.querySelector("[data-role=explain]");
  el.querySelector("select").addEventListener("change", ev => {
    box.innerHTML = explain(ev.target.value);
    setSetting(setting, ev.target.value);
  });
  return el;
}

function selectControl({setting, label, hint, options}) {
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    <label>${label}</label>
    <select>${options.map(([v, t]) =>
      `<option value="${v}"${settings[setting] === v ? " selected" : ""}>${t}</option>`
      ).join("")}</select>
    ${hint ? `<div class="hint" data-role="hint">${hint}</div>` : ""}`;
  el.querySelector("select").addEventListener("change", ev =>
    setSetting(setting, ev.target.value));
  return el;
}

function buttonControl({label, hint, onClick, primary}) {
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `<button type="button" class="wide${primary ? " prim" : ""}">${label}</button>
    ${hint ? `<div class="hint">${hint}</div>` : ""}`;
  el.querySelector("button").addEventListener("click", onClick);
  return el;
}

function card(title, children) {
  const el = document.createElement("div");
  el.className = "card";
  if (title) el.innerHTML = `<h2>${title}</h2>`;
  children.filter(Boolean).forEach(c => el.append(c));
  return el;
}

/* The clade palette, with a strip of its first twelve colours. */
function paletteControl() {
  const el = document.createElement("div");
  el.className = "ctl";
  const strip = () => Array.from({length: 12}, (_, i) =>
    `<i style="background:${cladeColour(null, i, 12)}"></i>`).join("");
  el.innerHTML = `<label>Clade palette</label>
    <select>${CLADE_PALETTES.map(([v, t]) =>
      `<option value="${v}"${settings.cladePal === v ? " selected" : ""}>${t}</option>`).join("")}</select>
    <div class="palstrip" data-role="strip">${strip()}</div>
    <div class="hint">For colour by clade, here and in Figure 2. A cycled palette
      repeats lighter, then darker; a ramp runs down the tree, so nearby clades
      get nearby colours.</div>`;
  const box = el.querySelector("[data-role=strip]");
  el.querySelector("select").addEventListener("change", ev => {
    setSetting("cladePal", ev.target.value);
    box.innerHTML = strip();
  });
  return el;
}

/* A compact rail: hints hidden behind one toggle, so every control fits on
   screen. On a tab whose main carries .fit, the captions go with them. Per
   page view, not stored. */
let RAIL_EXPLAIN = false;
function railCompact(rail, {toggle = true} = {}) {
  rail.classList.add("compact");
  if (!toggle) return;
  const b = document.createElement("button");
  b.type = "button";
  b.className = "hinttog";
  const lab = () => RAIL_EXPLAIN ? "Hide explanations" : "Show explanations";
  b.textContent = lab();
  b.addEventListener("click", () => {
    RAIL_EXPLAIN = !RAIL_EXPLAIN;
    document.body.classList.toggle("explain", RAIL_EXPLAIN);
    document.querySelectorAll(".hinttog").forEach(x => { x.textContent = lab(); });
    const tab = tabById(activeTabId);
    if (tab) renderTab(tab);
  });
  rail.append(b);
}

/* Hover popups: inside a `.hintpop-on` container, a control's own hint -- or
   any element's data-hint -- appears beside it on hover, or on a tap on touch.
   One popup for the page, wired once. */
let HINT_POP = null;
function wireHintPopups() {
  if (HINT_POP || typeof document.addEventListener !== "function") return;
  HINT_POP = document.createElement("div");
  HINT_POP.className = "hintpop";
  HINT_POP.hidden = true;
  document.body.append(HINT_POP);
  let cur = null;
  const source = t => {
    const box = t && t.closest ? t.closest(".hintpop-on") : null;
    if (!box) return null;
    const d = t.closest("[data-hint]");
    if (d && box.contains(d)) return {el: d, html: d.dataset.hint};
    const c = t.closest(".ctl");
    const h = c ? [...c.children].find(x => x.classList.contains("hint") && !x.classList.contains("err")) : null;
    return h ? {el: c, html: h.innerHTML} : null;
  };
  const hide = () => { HINT_POP.hidden = true; cur = null; };
  const show = s => {
    cur = s.el;
    HINT_POP.innerHTML = s.html;
    HINT_POP.hidden = false;
    const r = s.el.getBoundingClientRect(), p = HINT_POP.getBoundingClientRect();
    const W = innerWidth, H = innerHeight;
    let left, top;
    if (r.left > p.width + 16) { left = r.left - p.width - 10; top = r.top; }
    else { left = Math.min(r.left, W - p.width - 8); top = r.bottom + 6; }
    HINT_POP.style.left = Math.max(8, left) + "px";
    HINT_POP.style.top = Math.max(8, Math.min(top, H - p.height - 8)) + "px";
  };
  document.addEventListener("mouseover", ev => {
    const s = source(ev.target);
    if (!s) { if (cur) hide(); return; }
    if (s.el !== cur) show(s);
  });
  document.addEventListener("touchstart", ev => {
    const s = source(ev.target);
    s && s.el !== cur ? show(s) : hide();
  }, {passive: true});
  document.addEventListener("scroll", hide, true);
}

/* ------------------------------------------------- composed rail cards */

/* Scenario and seed. Used by the tabs that show one run at a time. The seed
   maximum is the run count the index reports, not a hardcoded 49. */
function railSimulation(maxSeed) {
  return card("Simulation", [
    selectControl({setting: "scenario", label: "Scenario",
      options: SCENARIOS.map(s => [s, SCENARIO_LABEL[s]]),
      hint: "Scenarios are never pooled: each one is its own design point."}),
    sliderControl({setting: "seed", label: "Seed", min: 0,
      max: Math.max(0, (maxSeed || 50) - 1), step: 1,
      format: v => String(v).padStart(4, "0"),
      hint: "One simulation. The simulation is the independent unit everywhere in the studio."}),
  ]);
}

/* The clade definition: both halves of rule B, plus the burn-in toggle. */
function railClade({origin = false} = {}) {
  return card("Clade definition", [
    sliderControl({setting: "k", label: "Shared-mutation threshold <i>k</i>",
      values: AXES.k,
      hint: "A lineage founds a clade when it has diverged at least <i>k</i> "
          + "mutations from its clade root.<br>"
          + "<strong>This names a clade build</strong>, so the stops are the "
          + "grid's own <i>k</i> axis."}),
    sliderControl({setting: "mt", label: "Minimum transmissions <i>m<sub>T</sub></i>",
      values: AXES.mt, format: v => mtPct(v) + "%",
      hint: "…<strong>and</strong> it and its descendants carry at least this "
          + "share of the simulation's transmission events. A share, not a "
          + "count: the scenarios run 11,846 to 27,510 transmissions, so a "
          + "fixed count would be a different filter in each. Stored as a "
          + "fraction, shown here as a percentage.<br>"
          + "<strong>This names a clade build.</strong> Under rule B a refused "
          + "lineage is folded into its ancestor, so m<sub>T</sub> changes the "
          + "partition and cannot be applied over an existing build — it is "
          + "one of the values the grid was built over, and a value with no build "
          + "shows as not computed."}),
    checkControl({setting: "burnin",
      label: "Exclude the burn-in from frequency and infection counts",
      hint: "Clades and their defining mutations are never windowed — a clade's "
          + "defining mutations were laid down once."}),
    origin && originControl("Origin rule"),
  ]);
}

/* A figure tab's rail never sets data: this card says where each rule is set. */
function dataCard(html) {
  const el = document.createElement("div");
  el.className = "card";
  el.innerHTML = `<h2>Data</h2><div class="hint" style="display:block">${html}
    Save and load whole sets of them on any exploration tab.</div>`;
  return el;
}
const DATA_CLADE = `The clade definition (<i>k</i>, <i>m</i><sub>T</sub>, burn-in) and the origin
  rule are set on <a href="#clustering">Clustering</a>.`;

/* Scenario toggles in each scenario's own colour: the picked ones are drawn in
   full and the rest faded (scenarioColour). None picked is all in full. */
function railFocus({hide = false} = {}) {
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `
    <div class="sctog">${SCENARIOS.map(sc => `<button type="button" data-v="${sc}"
      aria-pressed="${focusSet().has(sc)}" style="--c:${css("--sc" + SCENARIOS.indexOf(sc))}"
      ><i></i>${SCENARIO_LABEL[sc]}</button>`).join("")}</div>
    <div class="hint">${hide ? "Pick one or several to show only those."
      : "Pick one or several to make them stand out; the others are faded, not removed."}
      None picked shows every scenario.</div>`;
  const btns = el.querySelectorAll("button");
  const paint = () => { const f = focusSet();
    btns.forEach(b => b.setAttribute("aria-pressed", String(f.has(b.dataset.v)))); };
  btns.forEach(b => b.addEventListener("click", () => {
    const f = focusSet();
    f.has(b.dataset.v) ? f.delete(b.dataset.v) : f.add(b.dataset.v);
    setSetting("focus", SCENARIOS.filter(sc => f.has(sc)).join(","));
    paint();
  }));
  return card("Scenarios", [el]);
}

function originControl(label) {
  return segmentControl({setting: "originRule", label, options: [
      ["frac", "Ratio"], ["vote", "Majority"], ["root", "Root"]],
      hint: "<strong>Ratio</strong> is the long-born share of the clade's "
          + "defining mutations. <strong>Majority</strong> is the pipeline's own "
          + "vote over those mutations, which cannot tie on an odd <i>k</i>. "
          + "<strong>Root</strong> is the host type of the clade's root lineage "
          + "alone, ignoring the rest."});
}

/* A clade id names one lineage of ONE run under ONE partition, so a selection
   cannot outlive a change to either. Registered before the shell's own
   listener, so the stale id is gone before anything redraws. */
function wireCladeSelection() {
  const resets = ["scenario", "seed", "k", "mt"];
  onSettingChange(name => { if (resets.includes(name)) setSetting("clade", null); });
}

/* The rise band writes two settings, so it needs its own wiring: hi = 1 - lo. */
function wireRiseBand() {
  onSettingChange(name => { if (name === "lo") setSetting("hi", +(1 - settings.lo).toFixed(2)); });
}
