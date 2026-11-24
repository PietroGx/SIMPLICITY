"use strict";
// Every value is validated and derived by the server, through the pipeline's
// own code (experiment_settings.py, impact_long_shedders_config.py); nothing
// here decides whether a file is valid.

const $ = (id) => document.getElementById(id);
const el = (tag, attrs = {}, ...kids) => {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") n.className = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined && v !== false) n.setAttribute(k, v);
  }
  for (const k of kids.flat()) if (k !== null && k !== undefined) n.append(k);
  return n;
};
const fmt = (v) => typeof v === "number" ? v.toLocaleString("en-US", {maximumFractionDigits: 6}) : String(v);

let S = null;              // the settings being edited
let snapshot = "";         // JSON of S as last loaded or saved
let loadedFrom = "";       // file S was loaded from
let nameTouched = false;   // the user typed a file name
let timer = null;
let seq = 0;

// [path, label, kind, hint]. kind: text | int | int? | num | str? | select:a,b
const SECTIONS = [
  {title: "Run", fields: [
    ["description", "Description", "text", "What this configuration is for. Kept in the file."],
    ["consensus", "Consensus mode", "select:distribution,argmax", "How a lineage's distance from the population consensus is measured. One mode per run: a calibration does not transfer between modes."],
  ]},
  {title: "Population", note: "cal_2 and the experiment phase", fields: [
    ["population.population_size", "Population size N", "int", "Hosts in the initial population."],
    ["population.infected_individuals_at_start", "Infected at start", "int", "Below N."],
    ["population.R", "R", "num", "Reproduction number of standard hosts."],
    ["population.IH_virus_emergence_rate", "IH_virus_emergence_rate", "num", "Intra-host virus emergence rate (kv in output directory names)."],
    ["population.final_time", "Simulated window (days)", "int", "The experiment phase's window. cal_2 calibrates over the same one."],
  ], derived: "population"},
  {title: "Scenarios", custom: "scenarios"},
  {title: "cal_1", note: "long shedders alone: 100% long-shedder population, NSR_long swept", fields: [
    ["cal_1.R", "R", "num", "Reproduction number in the cal_1 population."],
    ["cal_1.final_time", "Base window (days)", "int", "A cohort whose infections need more gets a whole multiple of it."],
    ["cal_1.min_cycles", "Minimum infection cycles", "int", "Complete long-shedder infections each cohort's window must fit."],
    ["cal_1.nsr_long.min", "NSR_long, lowest", "num", ""],
    ["cal_1.nsr_long.max", "NSR_long, highest", "num", "Keep the sweep to one order of magnitude: the fit's model does not follow the curve over two."],
    ["cal_1.nsr_long.steps", "Sweep points", "int", "Geometrically spaced."],
  ], derived: "cal_1"},
  {title: "cal_2", note: "standard hosts alone: 0% long shedders, NSR swept", fields: [
    ["cal_2.nsr_std.min", "NSR, lowest", "num", ""],
    ["cal_2.nsr_std.max", "NSR, highest", "num", "One order of magnitude, as for cal_1."],
    ["cal_2.nsr_std.steps", "Sweep points", "int", "Geometrically spaced."],
  ], derived: "cal_2"},
  {title: "Calibration targets", note: "OSR, substitutions per site per year", fields: [
    ["targets.osr_std", "OSR target, standard", "num", "Global root-to-tip clock cal_2 calibrates NSR to."],
    ["targets.osr_long", "OSR target, long shedders", "num", "Intra-host clock cal_1 calibrates NSR_long to."],
  ]},
  {title: "Seeds", fields: [
    ["seeds.cal", "Seeds per calibration sweep point", "int", "cal_1 and cal_2."],
    ["seeds.exp", "Seeds per scenario", "int", "Experiment phase."],
  ]},
  {title: "Slurm", note: "per task", fields: [
    ["slurm.mem", "Memory", "str", "cal_1 and the experiment phase, e.g. 3G."],
    ["slurm.mem_cal2", "Memory, cal_2", "str?", "Empty: same as Memory."],
    ["slurm.time", "Walltime", "str?", "D-HH:MM:SS. Empty: the pipeline default (SIMPLICITY_SLURM_TIME, 05:00:00)."],
  ]},
];

const SCEN_COLS = [
  ["name", "Name", "str", "Letters and digits joined by single underscores."],
  ["long_shedders_ratio", "Prevalence (fraction)", "num", "long_shedders_ratio: the fraction of hosts who are long shedders. 0 makes the scenario a control."],
  ["inf_duration_long", "Shedding duration (d)", "num?", "inf_duration_long. Null for a control."],
  ["R_long", "R_long", "num?", "Long shedders' reproduction number. Null for a control."],
  ["susceptibility_long", "susceptibility_long", "num", "How much more or less likely a long shedder is to catch an infection."],
];

const get = (o, path) => path.split(".").reduce((a, k) => a == null ? a : a[k], o);
const set = (o, path, v) => {
  const ks = path.split("."), last = ks.pop();
  ks.reduce((a, k) => a[k], o)[last] = v;
};
function parse(kind, raw) {
  const t = raw.trim();
  if (kind === "text") return raw;
  if (t === "") return kind.endsWith("?") ? null : "";
  if (kind.startsWith("str")) return t;
  const n = Number(t);
  return Number.isNaN(n) ? t : n;           // a non-number goes through, for the server to name
}
const show = (v) => v === null || v === undefined ? "" : String(v);

function input(path, kind, hint, obj = S, onchange = changed) {
  const id = "f_" + path.replace(/\W/g, "_");
  let n;
  if (kind.startsWith("select:")) {
    n = el("select", {id}, kind.slice(7).split(",").map(o => el("option", {value: o}, o)));
    n.value = get(obj, path);
    n.addEventListener("change", () => { set(obj, path, n.value); onchange(); });
  } else if (kind === "text") {
    n = el("textarea", {id, rows: 2});
    n.value = show(get(obj, path));
    n.addEventListener("input", () => { set(obj, path, parse(kind, n.value)); onchange(); });
  } else {
    n = el("input", {id, type: "text", spellcheck: "false", autocomplete: "off",
                     inputmode: kind.startsWith("str") ? null : "decimal",
                     placeholder: kind.endsWith("?") ? "null" : null});
    n.value = show(get(obj, path));
    n.addEventListener("input", () => { set(obj, path, parse(kind, n.value)); onchange(); });
  }
  n.dataset.path = path;
  if (hint) n.title = hint;
  return n;
}

function renderForm() {
  const form = $("form");
  form.replaceChildren();
  for (const sec of SECTIONS) {
    const card = el("section", {class: "card"},
      el("h2", {}, sec.title, sec.note ? el("span", {class: "note"}, sec.note) : null));
    if (sec.custom === "scenarios") card.append(scenarioTable());
    else {
      const grid = el("div", {class: "grid"});
      for (const [path, label, kind, hint] of sec.fields) {
        const n = input(path, kind, hint);
        grid.append(el("label", {for: n.id}, label), el("div", {class: "field"}, n,
                       hint ? el("small", {}, hint) : null));
      }
      card.append(grid);
    }
    if (sec.derived) card.append(el("div", {class: "derived", id: "d_" + sec.derived}));
    form.append(card);
  }
}

function scenarioTable() {
  const wrap = el("div", {id: "scen"});
  const head = el("tr", {}, SCEN_COLS.map(([, label, , hint]) => el("th", {title: hint}, label)),
                  el("th", {title: "tau_3_long: shedding duration minus the other infection phases"}, "tau_3_long (d)"),
                  el("th", {title: "round(N × prevalence)"}, "Long shedders at start"), el("th"));
  const body = el("tbody");
  S.scenarios.forEach((x, i) => {
    const ctrl = x.long_shedders_ratio === 0;
    const tr = el("tr", {class: ctrl ? "ctrl" : ""});
    for (const [k, , kind, hint] of SCEN_COLS) {
      const n = input(`scenarios.${i}.${k}`, kind, hint, S, () => {
        const wasCtrl = ctrl, isCtrl = x.long_shedders_ratio === 0;
        if (k === "long_shedders_ratio" && wasCtrl !== isCtrl) {
          if (isCtrl) { x.inf_duration_long = null; x.R_long = null; }
          else { x.inf_duration_long = x.inf_duration_long ?? 63.0;
                 x.R_long = x.R_long ?? S.population.R; }
          rerenderScenarios(`scenarios.${i}.long_shedders_ratio`);
        }
        changed();
      });
      if (ctrl && (k === "inf_duration_long" || k === "R_long")) n.disabled = true;
      tr.append(el("td", {}, n));
    }
    tr.append(el("td", {class: "num", id: `d_tau_${i}`}), el("td", {class: "num", id: `d_long_${i}`}),
              el("td", {}, el("button", {type: "button", class: "x", title: "Remove this scenario",
                 onclick: () => { S.scenarios.splice(i, 1); rerenderScenarios(); changed(); }}, "×")));
    body.append(tr);
  });
  wrap.append(el("div", {class: "tablewrap"}, el("table", {}, el("thead", {}, head), body)),
    el("button", {type: "button", onclick: () => {
      S.scenarios.push({name: `scenario_${S.scenarios.length + 1}`, inf_duration_long: 63.0,
                        long_shedders_ratio: 0.01, R_long: S.population.R, susceptibility_long: 1.0});
      rerenderScenarios(); changed();
    }}, "+ Add scenario"),
    el("small", {class: "block"}, "Scenarios sharing a shedding duration and R_long share one cal_1 cohort."));
  return wrap;
}

function rerenderScenarios(focusPath) {
  $("scen").replaceWith(scenarioTable());
  if (focusPath) document.querySelector(`[data-path="${focusPath}"]`)?.focus();
}

function changed() {
  $("dirty").hidden = JSON.stringify(S) === snapshot;
  $("json").textContent = JSON.stringify(S, null, 2);
  suggestName();
  clearTimeout(timer);
  timer = setTimeout(check, 200);
}

function suggestName() {
  if (nameTouched) return;
  const n = S.population.population_size, t = S.population.final_time;
  const auto = `impact_long_shedders_N${n}_T${t}.json`;
  $("fname").value = loadedFrom && JSON.stringify(S) === snapshot ? loadedFrom : auto;
}

async function api(path, body) {
  const r = await fetch(path, body === undefined ? {} :
    {method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify(body)});
  let data = {};
  try { data = await r.json(); } catch { data = {error: `HTTP ${r.status}`}; }
  return {ok: r.ok, status: r.status, data};
}

const ERR_PATH = /^(?:scenarios\[(\d+)\](?: '[^']*')?(?:\.(\w+))?|([a-z_0-9]+(?:\.[a-z_A-Z0-9]+)*))[:\s]/;
function pathOf(msg) {
  const m = msg.match(ERR_PATH);
  if (!m) return null;
  if (m[1] !== undefined) return m[2] ? `scenarios.${m[1]}.${m[2]}` : `scenarios.${m[1]}`;
  return m[3];
}

async function check() {
  const mine = ++seq;
  const {data} = await api("/api/check", {settings: S});
  if (mine !== seq) return;                       // a newer edit is already being checked
  document.querySelectorAll(".bad").forEach(n => { n.classList.remove("bad"); });
  const box = $("check");
  box.replaceChildren();
  const errs = data.errors || [data.error || "no answer from the server"];
  $("check-card").className = "card " + (errs.length ? "has-err" : data.derived?.warnings.length ? "has-warn" : "ok");
  if (errs.length) {
    box.append(el("p", {class: "lead err"}, `${errs.length} error${errs.length > 1 ? "s" : ""}: the pipeline would refuse this file`),
               el("ul", {class: "msgs err"}, errs.map(e => el("li", {}, e))));
    for (const e of errs) {
      const p = pathOf(e);
      if (!p) continue;
      document.querySelectorAll(`[data-path="${p}"], [data-path^="${p}."]`).forEach(n => {
        n.classList.add("bad"); n.title = e;
      });
    }
    renderDerived(null);
    return;
  }
  const w = data.derived.warnings;
  box.append(el("p", {class: "lead okt"}, "Valid: the pipeline would accept this file."));
  if (w.length) box.append(el("p", {class: "lead warn"}, `${w.length} to look at`),
                           el("ul", {class: "msgs warn"}, w.map(x => el("li", {}, x))));
  renderDerived(data.derived);
}

function kv(rows) {
  return el("dl", {class: "kv"}, rows.flatMap(([k, v]) => [el("dt", {}, k), el("dd", {}, v)]));
}

function renderDerived(d) {
  const clear = (id) => { const n = $(id); if (n) n.replaceChildren(); return n; };
  const pop = clear("d_population"), c1 = clear("d_cal_1"), c2 = clear("d_cal_2"), sims = $("sims");
  S.scenarios.forEach((_, i) => { for (const k of ["tau", "long"]) { const n = $(`d_${k}_${i}`); if (n) n.textContent = ""; } });
  if (!d) { sims.replaceChildren(el("p", {class: "muted"}, "Fix the errors to see the counts.")); return; }
  pop.append(kv([
    ["Reservoir", `${fmt(d.reservoir)} individuals (100 × N)`],
    ["Infected at start", `${d.infected_at_start_pct.toFixed(2)}% of N`],
  ]));
  d.scenarios.forEach((r, i) => {
    const t = $(`d_tau_${i}`), l = $(`d_long_${i}`);
    if (t) { t.textContent = S.scenarios[i].long_shedders_ratio > 0 ? fmt(r.tau_3_long) : "—";
             t.classList.toggle("badtxt", S.scenarios[i].long_shedders_ratio > 0 && r.tau_3_long <= 0); }
    if (l) { l.textContent = fmt(r.initial_long);
             l.classList.toggle("badtxt", S.scenarios[i].long_shedders_ratio > 0 && r.initial_long === 0); }
  });
  c1.append(
    el("h3", {}, `${d.cohorts.length} cohort${d.cohorts.length === 1 ? "" : "s"}`),
    d.cohorts.length ? el("div", {class: "tablewrap"}, el("table", {class: "mini"},
      el("thead", {}, el("tr", {}, ["group", "tau_3_long (d)", "R_long", "window (d)"].map(h => el("th", {}, h)))),
      el("tbody", {}, d.cohorts.map(c => el("tr", {},
        el("td", {}, el("code", {}, c.group)), el("td", {class: "num"}, fmt(c.tau_3_long)),
        el("td", {class: "num"}, fmt(c.R_long)), el("td", {class: "num"}, fmt(c.final_time)))))))
      : el("p", {class: "muted"}, "No long-shedder scenario, so cal_1 has nothing to calibrate."),
    kv([["NSR_long sweep", d.nsr_long_grid.map(fmt).join(" · ")]]));
  c2.append(kv([["Window", `${fmt(d.cal2_final_time)} d, the experiment phase's`],
                ["NSR sweep", d.nsr_std_grid.map(fmt).join(" · ")]]));
  const s = d.simulations;
  sims.replaceChildren(kv([
    ["cal_1", `${fmt(s.cal_1)}  (${d.cohorts.length} cohorts × ${S.cal_1.nsr_long.steps} points × ${S.seeds.cal} seeds)`],
    ["cal_2", `${fmt(s.cal_2)}  (${S.cal_2.nsr_std.steps} points × ${S.seeds.cal} seeds)`],
    ["experiment", `${fmt(s.exp)}  (${S.scenarios.length} scenarios × ${S.seeds.exp} seeds)`],
    ["total", fmt(s.total)],
  ]));
}

async function loadFiles(select) {
  const {data} = await api("/api/files");
  $("dir").textContent = data.dir + "/";
  $("dir2").textContent = data.dir + "/";
  const sel = $("files");
  sel.replaceChildren(...data.files.map(f => el("option", {value: f.name, title: f.description || f.error || ""},
    f.name + (f.name === data.default ? "  (default)" : "") + (f.valid ? "" : "  — invalid"))));
  sel.value = select || data.default;
  return data;
}

async function loadFile(name) {
  const {ok, data} = await api(`/api/file?name=${encodeURIComponent(name)}`);
  if (!ok) { alert(data.error); return; }
  S = data.settings;
  snapshot = JSON.stringify(S);
  loadedFrom = name;
  nameTouched = false;
  $("saved").replaceChildren();
  renderForm();
  changed();
}

async function save() {
  const name = $("fname").value.trim();
  const {ok, status, data} = await api("/api/save", {name, settings: S, overwrite: $("overwrite").checked});
  const out = $("saved");
  if (!ok) {
    out.replaceChildren(el("p", {class: "err"}, data.error),
      ...(data.errors ? [el("ul", {class: "msgs err"}, data.errors.map(e => el("li", {}, e)))] : []));
    if (status === 409) $("overwrite").focus();
    return;
  }
  snapshot = JSON.stringify(S);
  loadedFrom = name;
  nameTouched = false;
  $("overwrite").checked = false;
  changed();
  await loadFiles(name);
  const cmd = el("code", {class: "cmd"}, data.command);
  out.replaceChildren(el("p", {class: "okt"}, "Saved ", el("code", {}, data.saved)),
    el("p", {class: "muted"}, "Launch it on the HPC login node, inside tmux. The experiment is named by the minute it starts:"), cmd,
    el("button", {type: "button", onclick: (e) => {
      navigator.clipboard?.writeText(data.command).then(() => { e.target.textContent = "Copied"; });
    }}, "Copy command"));
}

$("load").addEventListener("click", () => {
  if ($("dirty").hidden || confirm("Discard unsaved changes?")) loadFile($("files").value);
});
$("fname").addEventListener("input", () => { nameTouched = true; });
$("save").addEventListener("click", save);
window.addEventListener("beforeunload", (e) => { if (!$("dirty").hidden) e.preventDefault(); });

loadFiles().then(d => loadFile(d.default));
