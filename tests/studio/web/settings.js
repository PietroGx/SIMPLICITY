/* The one copy of every setting in the studio.

   No tab keeps its own copy and no tab reads localStorage. A tab reads
   `settings`, and calls `setSetting` to change one; localStorage is a mirror so
   a browser reload starts where you left off, never the source of truth.

   Which store a setting reaches decides what changing it costs:

     build parameters   k, mt, win, prom, lo, hi, window   a build must exist
     read parameters    burnin                             applied when read
     view parameters    everything else                    nothing is refetched

   m_T moved from read to build in v2.5.5. Under rule B a lineage founds a
   clade only if its subtree carries at least m_T of the run's transmissions,
   so moving m_T moves the PARTITION: it cannot be applied over a build the way
   the burn-in toggle can, and each value needs its own clade build.

   `settingKind` answers that, so a tab can tell the difference without
   hardcoding a list of names.
*/

const SETTING_DEFAULTS = {
  // build parameters -- these name a clade build on disk
  k: 7, mt: 0.01, win: 15, prom: 0.10, lo: 0.10, hi: 0.90, window: "post",
  // read parameters -- applied over a build when a panel is computed
  mf: 0,
  // read parameters -- applied over a build when a panel is computed
  burnin: 1,
  // view parameters -- never leave the browser
  scenario: "control", seed: 0, originRule: "frac", colourBy: "clade",
  // the clade the Muller plot has selected, or null: see wireCladeSelection
  clade: null,
  treeAxis: "time", episodeView: "grid", episodeCols: 3, showRaw: 1,
  overview: 0, yScale: "shared",
  // which plot each panel tab draws
  focus: "",
  plot1: "clade_points",
  plot3: "rise_latency", axis4: "k",
  // Figure 5 breakouts: peak that makes a sweep big, hidden days, launch ceiling
  brkF: 0.05, brkH: 60, brkLo: 0.001,
  // the per-run points behind the interval, one box per panel -- per panel and
  // not one global switch, because each figure subpanel is captioned alone.
  dist1: 0, dist2: 0, dist3: 0, dist4: 0,
};

/* Settings whose value changes what the rail OFFERS, not just what it shows:
   the plot form decides whether the distribution box can be ticked at all. A rail is built once per tab so
   that dragging a slider cannot lose its focus, so these ask for a rebuild. */
const RAIL_SETTINGS = ["plot1", "plot3", "axis4", "focus"];

const BUILD_SETTINGS = ["k", "mt", "win", "prom", "lo", "hi", "window"];
const READ_SETTINGS = ["burnin", "mf"];

function settingKind(name) {
  if (BUILD_SETTINGS.includes(name)) return "build";
  if (READ_SETTINGS.includes(name)) return "read";
  return "view";
}

const SETTINGS_STORAGE_KEY = "studio";
/* Bump when a stored value changes MEANING, not when one is added: a saved
   blob from an older schema is dropped rather than merged. v2 is m_T moving
   from a percentage to a fraction -- a saved mt=1 meant 1% and would now ask
   for 100%, which has no build. */
const SETTINGS_VERSION = 2;

function loadSettings() {
  let saved = {};
  try {
    const raw = JSON.parse(localStorage.getItem(SETTINGS_STORAGE_KEY) || "{}");
    if (raw && raw.v === SETTINGS_VERSION) saved = raw.s || {};
  }
  catch (e) { saved = {}; }          // private window, blocked storage: defaults
  return Object.assign({}, SETTING_DEFAULTS, saved);
}

const settings = loadSettings();
const settingListeners = [];

function mirrorSettings() {
  try { localStorage.setItem(SETTINGS_STORAGE_KEY,
                             JSON.stringify({v: SETTINGS_VERSION, s: settings})); }
  catch (e) {}                       // a failed mirror must never break a render
}

/* Change one setting and tell whoever is listening. Setting a value it already
   holds notifies nobody, so a slider that fires `input` on every pixel of the
   same value costs nothing. */
function setSetting(name, value) {
  if (settings[name] === value) return false;
  settings[name] = value;
  mirrorSettings();
  const kind = settingKind(name);
  settingListeners.forEach(fn => fn(name, kind));
  return true;
}

function setSettings(patch) {
  const changed = Object.keys(patch).filter(n => setSetting(n, patch[n]));
  return changed.length > 0;
}

function onSettingChange(fn) { settingListeners.push(fn); }

function resetSettings(names) {
  const which = names || Object.keys(SETTING_DEFAULTS);
  setSettings(Object.fromEntries(which.map(n => [n, SETTING_DEFAULTS[n]])));
}

/* The query string for an endpoint that reads a build. One place builds it, so
   a tab cannot accidentally send a different parameter set than the state line
   shows. */
function buildQuery(extra) {
  const q = Object.assign({}, extra || {});
  for (const n of BUILD_SETTINGS.concat(READ_SETTINGS))
    if (!(n in q)) q[n] = settings[n];
  return Object.keys(q).sort()
    .map(n => `${n}=${encodeURIComponent(q[n])}`).join("&");
}

function settingsSummary() {
  return `k=${settings.k} · m_T=${mtPct(settings.mt)}% · win=${settings.win}d`
    + ` · prom=${settings.prom}`
    + ` · ${settings.window === "post" ? "post burn-in" : "whole record"}`;
}
