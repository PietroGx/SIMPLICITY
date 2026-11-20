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
    <div class="seg">${options.map(([v, t]) =>
      `<button type="button" data-v="${v}" aria-pressed="${settings[setting] === v}">${t}</button>`
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

/* The distribution overlay box, the same in every panel that offers it. `plot`
   is the panel's selected form; dist.js decides which forms can carry points. */
function railDistribution(setting, plot) {
  const ok = distAllowed(plot);
  return checkControl({
    setting, disabled: !ok,
    label: "Show the individual runs",
    why: "This plot form cannot carry them: a stack of means has no per-run "
       + "form that stays a stack, and a count of runs is not an estimate.",
    hint: ok ? "The whiskers are a bootstrap CI of the <strong>mean</strong> "
             + "(2000 resamples over runs), so they are narrower than the runs "
             + "themselves by construction. This draws the runs behind them."
             : "",
  });
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
      hint: "One run. The run is the independent unit everywhere in the studio."}),
  ]);
}

/* The clade definition: both halves of rule B, plus the burn-in toggle. */
function railClade() {
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
          + "share of the run's transmission events. A share, not a "
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
  ]);
}

/* How a clade's origin is decided. Restored from the artifact: the three rules
   disagree, and which is right is still open (BACKLOG, k-parity). */
function railOrigin() {
  return card("Origin", [
    segmentControl({setting: "originRule", options: [
      ["frac", "Ratio"], ["vote", "Majority"], ["root", "Root"]],
      hint: "<strong>Ratio</strong> is the long-born share of the clade's "
          + "defining mutations. <strong>Majority</strong> is the pipeline's own "
          + "vote over those mutations, which cannot tie on an odd <i>k</i>. "
          + "<strong>Root</strong> is the host type of the clade's root lineage "
          + "alone, ignoring the rest."}),
  ]);
}

/* Episode detection: the rules that decide what counts as a rise. */
function railEpisodes() {
  return card("Episode rules", [
    sliderControl({setting: "win", label: "Smoothing window",
      values: AXES.win, format: v => v + " d",
      hint: "Centred moving average. Wider merges neighbouring episodes — and "
          + "at 1 d there is no smoothing, so a rise collapses to the sampling "
          + "interval.<br>"
          + "<strong>This names a clade build</strong>: the stops are the "
          + "grid's own window axis."}),
    sliderControl({setting: "prom", label: "Peak prominence",
      values: AXES.prom, format: v => Number(v).toFixed(2),
      hint: "As a fraction of the clade's own peak height.<br>"
          + "<strong>This names a clade build</strong>, and the grid holds two "
          + "values. Off them the daily record has no build, the tab falls "
          + "back to the 9.2-day export and m<sub>F</sub> stops applying."}),
    sliderControl({setting: "mf", label: "Minimum peak <i>m<sub>F</sub></i>",
      values: MF_VALUES, format: v => mtPct(v) + "%",
      hint: "An episode is counted only if its peak reaches this share of the "
          + "whole viral population. Prominence is relative to the clade's own "
          + "height, so a clade that never leaves 1% still gets peaks; this is "
          + "the absolute floor.<br>"
          + "The clade is kept either way — only the episode stops counting. "
          + "At 0 every episode counts, which is what the build recorded."}),
    sliderControl({setting: "lo", label: "Rise band",
      values: AXES.lo,
      format: v => `${Number(v).toFixed(2)}–${(1 - v).toFixed(2)}`,
      hint: "Both ends are fractions of <strong>that episode's own "
          + "amplitude</strong>, peak minus its start. Kept at 0.10–0.90 so the "
          + "measure stays comparable to the final-approach definition."}),
    segmentControl({setting: "window", label: "Window", axis: AXES.window,
      options: [["post", "Post burn-in"], ["all", "Whole record"]],
      hint: "Selection strength is proportional to the fraction of the "
          + "population already infected, reaching full weight at the burn-in "
          + "day, so a rise before it grew partly under drift."}),
  ]);
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
