/* The only place the studio talks to the server.

   Everything a tab needs from the server comes through `apiGet`, which gives
   three things no caller has to repeat:

     cached      the same URL is never fetched twice; a tab switch costs nothing
     superseded  a second request in the same SLOT aborts the first, so dragging
                 a slider cannot leave the answer to an abandoned value on screen
     busy        one counter drives the pending indicator, with the label of
                 whatever is actually computing

   Why slots: two requests for different things must run together, but two
   requests for the SAME thing at different parameters must not race. The slot
   names the thing ("panel", "episodes"), not the parameters.

   Every URL is origin-relative, so putting nginx in front needs no change.
*/

const apiCache = new Map();          // url -> parsed response
const apiInflight = new Map();       // url -> {promise, controller, slot}
const apiSlotUrl = new Map();        // slot -> the url it is currently waiting on

/* Shown only after this long, so a 2 ms cache hit never flashes a spinner. */
const BUSY_SHOW_AFTER_MS = 150;

const busyListeners = [];
let busyJobs = [];                   // [{label, since}], newest last

function onBusyChange(fn) { busyListeners.push(fn); }

function announceBusy() {
  const visible = busyJobs.filter(j => Date.now() - j.since >= BUSY_SHOW_AFTER_MS);
  busyListeners.forEach(fn => fn(visible));
}

function busyStart(label) {
  const job = {label, since: Date.now()};
  busyJobs.push(job);
  setTimeout(announceBusy, BUSY_SHOW_AFTER_MS + 10);
  return job;
}

function busyEnd(job) {
  busyJobs = busyJobs.filter(j => j !== job);
  announceBusy();
}

class ApiSuperseded extends Error {
  constructor(url) { super(`superseded: ${url}`); this.superseded = true; }
}

/* True when a render should quietly stop rather than show an error: the answer
   it was waiting for is no longer the one wanted. */
function isSuperseded(e) { return !!(e && (e.superseded || e.name === "AbortError")); }

async function apiGet(url, opts) {
  const {slot, label, fresh} = opts || {};

  if (!fresh && apiCache.has(url)) return apiCache.get(url);

  // the same url already on its way: share it rather than open a second socket
  const running = apiInflight.get(url);
  if (running && !fresh) return running.promise;

  // a different url in the same slot is now stale -- stop waiting for it
  if (slot) {
    const stale = apiSlotUrl.get(slot);
    if (stale && stale !== url) {
      const s = apiInflight.get(stale);
      if (s) { s.controller.abort(); apiInflight.delete(stale); }
    }
    apiSlotUrl.set(slot, url);
  }

  const controller = new AbortController();
  const job = busyStart(label || "loading");

  const promise = (async () => {
    try {
      const r = await fetch(url, {signal: controller.signal,
                                  headers: {"Accept": "application/json"}});
      const body = await r.json().catch(() => ({error: `bad response ${r.status}`}));
      if (!r.ok || body.error) throw new Error(body.error || `HTTP ${r.status}`);
      apiCache.set(url, body);
      return body;
    } catch (e) {
      if (controller.signal.aborted) throw new ApiSuperseded(url);
      throw e;
    } finally {
      busyEnd(job);
      if (apiInflight.get(url) && apiInflight.get(url).controller === controller)
        apiInflight.delete(url);
    }
  })();

  apiInflight.set(url, {promise, controller, slot});
  return promise;
}

/* Drop cached answers whose url starts with `prefix`. Called after a rebuild:
   the server's answers changed, so the browser's copies are wrong. */
function apiInvalidate(prefix) {
  for (const url of [...apiCache.keys()])
    if (!prefix || url.startsWith(prefix)) apiCache.delete(url);
}

function apiCacheSize() { return apiCache.size; }

async function apiPost(url, body) {
  const r = await fetch(url, {method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify(body || {})});
  const j = await r.json().catch(() => ({error: `bad response ${r.status}`}));
  if (!r.ok && !j.error) throw new Error(`HTTP ${r.status}`);
  return j;
}

/* The convenience the tabs use: the slot and the label come from the path, so a
   tab names a URL and gets caching, abort-on-supersede and a pending label
   without repeating the bookkeeping. */
const API_LABELS = {
  "/api/panel": "computing the panel", "/api/across": "reading the clade build",
  "/api/episodes": "measuring episodes", "/api/index": "reading the cache state",
  "/api/jobs": "checking the job", "/api/run": "loading the simulation",
  "/api/runs": "reading the simulation index",
  "/api/run_daily": "reading the daily record",
  "/api/rebuild": "starting the pass",
  "/api/dist": "reading the individual simulations",
  "/api/sens": "walking the build grid",
};
function apiGetSlot(url, opts) {
  const path = String(url).split("?")[0];
  const n = new URLSearchParams(String(url).split("?")[1] || "").get("n");
  return apiGet(url, Object.assign({
    slot: path + (n ? "#" + n : ""),
    label: API_LABELS[path] || "loading",
  }, opts || {}));
}
