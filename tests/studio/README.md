# Figure 2 studio

One local page with eleven tabs over the impact-long-shedders runs. Not an
artifact: it is served from this box, bound to loopback, and reached over an
SSH tunnel.

Vocabulary: `GLOSSARY.md` at the repo root.
Design and the decisions behind each panel: `tests/artifacts/STUDIO_PLAN.md`.

## Run it

```bash
# on the box
python -I tests/studio/server.py            # 127.0.0.1:8600

# from a laptop
ssh -N -L 8600:127.0.0.1:8600 <this box>
# then open http://127.0.0.1:8600/
```

The server refuses to bind anything but loopback. It serves compute endpoints,
and this box also runs the walkingtoursmanager vhosts, so an exposed trigger
would be a production incident rather than a slow tool.

## Tabs

`web/index.html` is the only document. It loads its scripts in dependency
order and each tab registers itself (`id`, `group`, `num`, `title`), so adding
one means adding a file.

| group | tab | file | reads |
|---|---|---|---|
| figures | 1 · Model | (being written) | `fig1/` (`tests/studio/fig1.py`) |
| figures | 2 · Clustering | `web/tab_clades4.js` | `/api/runs`, `/api/run` |
| figures | 3 · Transmission | `web/tab_panel2.js` | `/api/panel?n=2` (`tx/`) |
| figures | 4 · Burden | `web/tab_panel1.js` | `/api/panel?n=1`, `n=1c` |
| figures | 5 · Rise vs latency | `web/tab_panel3.js` | `/api/panel?n=3` |
| supplementary | S1 · Clustering sensitivity | `web/tab_panel4.js` | `/api/panel?n=4`, `/api/sens` |
| data exploration | Clustering | `web/tab_clustering.js` | `/api/runs`, `/api/run` |
| data exploration | m_T effect on clades | `web/tab_across.js` | `/api/across` |
| data exploration | Sweep time | `web/tab_sweep.js` | `/api/run_daily` |
| tools | Method | `web/tab_method.js` | `/api/index` |
| tools | Export | `web/tab_export.js` | `/api/jobs`, `POST /api/rebuild/*` |

One settings object (`web/settings.js`) mirrored to `localStorage`, so a reload
starts where you left off. A setting is one of three kinds: **build** (`k`,
`win`, `prom`, `lo`, `hi`, `window`) needs a clade build, **read** (`mt`,
`burnin`) is applied over a build, and everything else is a view parameter that
never leaves the browser.

## Data

```
Data/figures_studio/<store>/           one store per experiment (sweep.py --print-paths)
  runs/     55 MB   phylogeny + sparse lineage frequency, 250 runs (120 time points)
  tx/               per-day 2x2 transmission counts by host type, and per-day
                    jump tables per transmitter type (export_transmissions.py)
  fig1/             fig1.json + index.json, Figure 1's data (fig1.py build)
  cache/            clade builds + the panel cache, keyed by parameters and
                    a hash of the code that produced them
```

**`runs/` is not deletable.** It is the only source for the phylogeny, the
Muller plot and the episode grid. The artifact-era `clades/` store and
`runs/formation.json` were deleted in v2.5.4; `formation.json` is regenerable
with `tests/clade_formation_grid.py --json`.

A clade build's key is
`k7_w15_p0.1_lo0.1_hi0.9_post__<hash of compute.py AND clustering.py>`. Both,
because a clade build is the output of both: an uncommitted one-line change to
the clade walk once altered what a clade is while `compute.py` was untouched,
and 24 builds carried the new maths under a key claiming the old.

The panel cache sits under `cache/_panels/`, keyed by the request and by a hash
of **both** `compute.py` and `panels.py`. The set of clade builds used to be in
that key, because panel 4 reads every build; panel 4 now caches one ROW per
build instead, so adding a build no longer discards every other panel's answer.
That took panel 4 from 30.4 s to 15 ms warm.

**`runs/` is not good enough for episodes.** It is downsampled to 120 time
points about 9 days apart; episode rise times land at 30-39 days, so a whole
rise would be 3-4 samples. The clade pass therefore reads the raw run output at
daily resolution and clusters with `simplicity.clustering`'s own function
rather than reimplementing it.

The Sweep time tab used to measure on `runs/` anyway, which showed it one
number and panel 3 another -- 275 emergence episodes against 13 on one run.
Since v2.5.4 it reads `/api/run_daily`, which serves the daily record, and the
episodes are detected once in python (`panels.run_daily` calling
`compute.episodes_of` and `compute.rise_of`); the browser draws what it is
given. `runs/` is still what the phylogeny and the Muller plot are drawn from --
those are shapes, not numbers.

## Rebuilding

Both are jobs: one at a time, loopback only, refused rather than queued.

```bash
# the transmission table. Reads each run's phylogeny for the jump tables,
# ~1.5 s a run; 250 runs took 1 min 40 s at --procs 4
python -I tests/export_transmissions.py --out <store>/tx --procs 2

# one clade build -- clades, composition, burden, episodes at one k. ~12 min
python -I tests/studio/compute.py --out Data/figures_studio/cache --procs 2 \
    --k 7 --win 15 --prom 0.10 --window post
```

`fig1.py build --exp-num N --store <store>` writes `<store>/fig1/`: one
process, ~5 min, ~1 GB.

Use `--procs 2`. The box has 4 cores and about 4 GB free with swap already in
use; the exporter costs roughly `0.15 + procs x 0.9` GB.

## Precompute on the HPC

**Staged for the next rebuild:** `pending_compute_ntx_post.patch` adds each
clade's infections on or after the burn-in day (`ntx_post`) to the clade pass.
Not applied yet because it changes `compute.py`'s hash and so every build's
name. Before launching: `git apply -p1 tests/studio/pending_compute_ntx_post.patch`,
then delete the patch. `run_precompute.sh` refuses to submit while it is unapplied.
Figure 4's infections bar uses `ntx_post` once the builds carry it.

`run_precompute.sh` runs the login-node checks, asks, then submits:

| job | script | depends on | resources |
|---|---|---|---|
| grid array, one build a task | `slurm_precompute.sh` | -- | 1 h, 16G, 4 cpu |
| transmission table into `tx/` | `slurm_tx.sh` | -- | 30 min, 8G, 4 cpu |
| Figure 1 data into `fig1/` | `slurm_fig1.sh` | -- | 45 min, 4G, 1 cpu |
| warm array, `--shard I/N` | `slurm_warm.sh` | afterany grid, tx | 1 h, 8G, 4 cpu |
| pack | `pack_cache.sh` | afterany warm, tx, fig1 | 3 h, 16G, 8 cpu |

```bash
EXP_NUM=2 bash tests/studio/run_precompute.sh --dry-run   # print the chain, submit nothing
EXP_NUM=2 bash tests/studio/run_precompute.sh             # checks, ask, submit
#   --no-fig1  --warm-shards 0  --resubmit  --limit N  --calibrate  --check-only
python -I tests/studio/warm_panels.py --exp-num 2 --shard 0/32 --list   # jobs per shard
```

The pack job tars `cache/` (with `_panels/`), `tx/` and `fig1/` into
`<store>/fig2_cache_<hash>_<date>.tar.gz`; a missing `tx/` or `fig1/` is a
warning in the log and the manifest, a partial grid is a refusal (exit 3).
The warm pass skips a memo only if its file exists, and the file name carries
the hash of the code that computes it, so a changed panel re-warms itself.

## Checks

```bash
# no server needed
python -I tests/check_clustering_module.py      # 31  the clade rules against what they replaced
python -I tests/check_manipulation.py           # 28  equal R, equal susceptibility
python -I tests/panel_denominators.py           #     what each panel's denominator is
python -I tests/rise_window_effect.py           #     is the rise the trajectory or the window?
node    tests/artifacts/check_sentinel.js       # 11  "not measurable" survives packing

# with the server up
python -I tests/check_transmissions.py          # 18  the transmission table, against a separate script
python -I tests/check_studio_api.py             #     the API a page actually receives
node    tests/studio/check_app.js               #     every tab rendered headlessly
node    tests/studio/check_episodes.js          # 25  the tab's episodes against the build's
```

`check_app.js` executes every tab against a running server, so a broken render
fails here instead of showing a blank page. It also asserts that no URL is
fetched twice in one render, that revisiting a tab fetches nothing, and that
every endpoint has a pending label. It replaced `check_pages.js`, which tested
the nine per-page documents that no longer exist.

`check_episodes.js` is the one that matters most: it pins the browser's
episode numbers to the build's, clade by clade.

## What the numbers are checked against

Nothing here is trusted because it looks plausible:

- the transmission table is assembled twice per run, from each individual's
  `parent` and from each individual's `new_infections`, and the two must agree
  on every cell and every day
- host counts match `tests/_prevalence_window.csv`, written by a different
  script on a different code path, on all 250 runs
- episodes per clade at the operating point come to 1.51-1.96, against the
  1.78-2.17 `tests/artifacts/rise_lab.html` measured independently
- the browser's emergence rise times match `compute.py`'s clade by clade,
  within 0.010 d (`check_episodes.js`)
- the clade rules in `simplicity/clustering.py` produce one partition from two
  different walks, on real runs (`check_clustering_module.py`)
- clades per run come to 10.1-20.6, against the 10-21 recorded in BACKLOG
- control carries zero burden in every composition class that requires a
  long-born defining mutation, and zero transmission events involving a long
  shedder
- composition classes sum to 100% of each run, because the founding clade is
  kept as the `wt` class rather than dropped

## Publishing later

Add a vhost, a cert, an htpasswd and a reverse proxy to 8600. No code change:
every URL in the pages is origin-relative.

Keep the split. Expose the read endpoints; leave `POST /api/rebuild/*` on
loopback. They are already refused for any non-loopback client, and a compute
trigger on a production host is the one risk here that mitigation does not
remove.
