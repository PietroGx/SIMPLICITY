/* The studio shell: one document, one tab bar, one render loop.

   The studio used to be nine separate HTML documents navigated with <a href>,
   so every tab switch was a browser reload that re-parsed the scripts and
   refetched the API. Here a tab switch shows a different view element and
   fetches nothing it already has.

   A tab registers itself with `registerTab`, so the tab list has one source of
   truth -- this registry -- rather than a PAGES array that each page repeated.

   The render loop has one rule worth knowing: a tab is re-rendered when a
   setting it cares about changes, and only the ACTIVE tab is rendered. The
   others are marked stale and redrawn when you switch to them, which is what
   makes a tab switch instant after the first visit.
*/

const tabs = [];
let activeTabId = null;
const staleTabs = new Set();
let renderTimer = null;

/* Nav groups, in nav order. A tab names its group in registerTab; inside a
   group tabs keep registration order, so a new figure joins by registering. */
const TAB_GROUPS = [
  {id: "explore", label: "Data exploration"},
  {id: "figures", label: "Figures", heading: n => `Figure ${n}: `, label1: true},
  {id: "supp", label: "Supplementary", heading: n => `${n} · `},
  {id: "tools", label: "Tools"},
];

function registerTab(tab) {
  if (!TAB_GROUPS.some(g => g.id === tab.group))
    console.warn(`tab ${tab.id}: unknown group ${tab.group}, shown under Tools`);
  tabs.push(Object.assign({railBuilt: false, everRendered: false}, tab));
}

function tabById(id) { return tabs.find(t => t.id === id); }

function groupOf(t) {
  return TAB_GROUPS.find(g => g.id === t.group) || TAB_GROUPS[TAB_GROUPS.length - 1];
}

/* Tabs in nav order: group order first, registration order inside a group. */
function navOrder() {
  return TAB_GROUPS.flatMap(g => tabs.filter(t => groupOf(t) === g));
}

function tabHeading(t) {
  const g = groupOf(t), pre = g.heading && t.num != null ? g.heading(t.num) : "";
  return pre && !t.title.startsWith(pre.replace(/( ·|:) $/, "")) ? pre + t.title : t.title;
}

/* ---------------------------------------------------------------- chrome */

function buildNav() {
  const nav = document.createElement("nav");
  nav.className = "tabs";
  nav.setAttribute("role", "tablist");
  nav.innerHTML = TAB_GROUPS.map(g => {
    const ts = tabs.filter(t => groupOf(t) === g);
    return !ts.length ? "" :
      `<div class="tgrp" data-group="${g.id}"><span class="glab">${g.label}</span>` +
      ts.map(t => `<button type="button" role="tab" data-tab="${t.id}" aria-selected="false">` +
        /* a figure's tab reads "Figure N: title"; other numbered tabs keep the badge */
        (g.label1 ? tabHeading(t) : (t.num != null && t.num !== "" ? `<span class="num">${t.num}</span>` : "")
          + t.title) + `</button>`).join("") + `</div>`;
  }).join("")
    + `<span class="sp" id="busy" hidden></span>`
    + `<span class="sp" id="stateLine"></span>`;
  nav.addEventListener("click", ev => {
    const b = ev.target.closest("button[data-tab]");
    if (b) showTab(b.dataset.tab);
  });
  document.body.prepend(nav);
}

function paintStateLine() {
  const el = document.getElementById("stateLine");
  if (el) el.textContent = settingsSummary();
}

/* The pending indicator. api.js reports what is in flight and for how long; a
   request that answers inside 150 ms never reaches here, so a cache hit does
   not flash. */
function paintBusy(jobs) {
  const el = document.getElementById("busy");
  if (!el) return;
  if (!jobs.length) { el.hidden = true; el.textContent = ""; return; }
  const secs = Math.round((Date.now() - jobs[0].since) / 1000);
  el.hidden = false;
  el.textContent = jobs.length > 1
    ? `computing ${jobs.length} things…`
    : `${jobs[0].label}…${secs >= 2 ? " " + secs + "s" : ""}`;
}

/* Keeps the elapsed seconds moving while something slow is running, so a long
   wait visibly ticks instead of looking frozen. */
function startBusyTicker() {
  setInterval(() => {
    const el = document.getElementById("busy");
    if (el && !el.hidden) announceBusy();
  }, 1000);
}

/* ---------------------------------------------------------------- views */

function viewElement(tab) {
  let el = document.getElementById(tab.viewId);
  if (el) return el;
  el = document.createElement("div");
  el.id = tab.viewId;
  el.hidden = true;
  el.innerHTML = `<div class="wrap"><div class="cols">
      <div class="maincol"><h1>${tabHeading(tab)}</h1><p class="lede">${tab.lede || ""}</p>
        <div class="main" id="${tab.id}_main"><div class="loading">loading…</div></div></div>
      <div class="rail" id="${tab.id}_rail"></div></div></div>`;
  document.body.append(el);
  return el;
}

async function renderTab(tab) {
  const main = document.getElementById(`${tab.id}_main`);
  const rail = document.getElementById(`${tab.id}_rail`);
  if (!tab.railBuilt && tab.rail) {
    rail.innerHTML = "";          /* a rebuild must not stack a second copy */
    tab.rail(rail);
    tab.railBuilt = true;
  }
  try {
    await tab.render(main, rail);
    tab.everRendered = true;
    staleTabs.delete(tab.id);
  } catch (e) {
    if (isSuperseded(e)) return;          // a newer render is already on its way
    main.innerHTML = "";
    main.append(errorBox(e.message));
    staleTabs.delete(tab.id);
  }
}

function showTab(id) {
  const tab = tabById(id);
  if (!tab) return;
  activeTabId = id;
  tabs.forEach(t => {
    const el = viewElement(t);
    el.hidden = t.id !== id;
  });
  document.querySelectorAll("nav.tabs button[data-tab]").forEach(b =>
    b.setAttribute("aria-selected", String(b.dataset.tab === id)));
  if (location.hash.slice(1) !== id) history.replaceState(null, "", "#" + id);
  if (!tab.everRendered || staleTabs.has(id)) renderTab(tab);
}

/* A setting changed. Mark every tab that cares as stale and redraw the active
   one; debounced, so dragging a slider renders once per frame rather than once
   per pixel. */
function onSettingsChanged(name, kind) {
  paintStateLine();
  /* A control that decides what the rail offers has to rebuild it: the
     distribution box is greyed by the selected plot form, and a rail built once
     keeps whatever state it had when the tab was first opened. */
  if (RAIL_SETTINGS.includes(name)) tabs.forEach(t => { t.railBuilt = false; });
  tabs.forEach(t => {
    if (!t.rerenderOn || t.rerenderOn.includes(kind)) staleTabs.add(t.id);
  });
  clearTimeout(renderTimer);
  renderTimer = setTimeout(() => {
    const tab = tabById(activeTabId);
    if (tab && staleTabs.has(tab.id)) renderTab(tab);
  }, 90);
}

function startStudio() {
  buildNav();
  tabs.forEach(viewElement);
  onSettingChange(onSettingsChanged);
  onBusyChange(paintBusy);
  startBusyTicker();
  paintStateLine();
  const wanted = location.hash.slice(1);
  showTab(tabById(wanted) ? wanted : navOrder()[0].id);
  window.addEventListener("hashchange", () => {
    const id = location.hash.slice(1);
    if (tabById(id) && id !== activeTabId) showTab(id);
  });
  /* Canvases are sized from their container, so a resize needs a redraw. */
  let resizeTimer = null;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => {
      const tab = tabById(activeTabId);
      if (tab && tab.everRendered) renderTab(tab);
    }, 150);
  });
  /* ...and so does a column that narrows without one: a scrollbar arriving
     when a late grid makes the page taller. Width only, so a taller page
     after the redraw does not trigger another. */
  if (typeof ResizeObserver !== "undefined") {
    const seen = new Map();
    const ro = new ResizeObserver(entries => entries.forEach(e => {
      const wNow = Math.round(e.contentRect.width), was = seen.get(e.target);
      seen.set(e.target, wNow);
      if (!was || !wNow || Math.abs(wNow - was) < 2) return;   // 0: a hidden tab being shown
      const tab = tabById(activeTabId);
      if (!tab || e.target.id !== `${tab.id}_main` || !tab.everRendered) return;
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => renderTab(tab), 150);
    }));
    tabs.forEach(t => { const el = document.getElementById(`${t.id}_main`); if (el) ro.observe(el); });
  }
}
