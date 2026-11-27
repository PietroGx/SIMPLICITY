# figures_gui -- the studio

The only source of the paper's figures. Exploration tabs set the data rules
(and save/load named sets of them), figure tabs set only how each figure looks,
and the Export tab writes PNG, JPEG or PDF at a column width and dpi. What the
studio kept from the retired `scripts/nature_plots` is in `prep/`: the scenario
list, the clade build's episode rules, the burn-in day and Figure 1's modules.

One local page over the impact-long-shedders runs. Not an
artifact: it is served from this box, bound to loopback, and reached over an
SSH tunnel.

Vocabulary: `GLOSSARY.md` at the repo root.
Design and the decisions behind each panel: `tests/artifacts/STUDIO_PLAN.md`.

## Run it

```bash
# on the box
python -I figures_gui/server.py            # 127.0.0.1:8600, experiment 1
STAMP=2 python -I figures_gui/server.py  # another experiment (or --stamp 2)

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
| figures | 1 · Model | `web/tab_fig1.js` | `/api/fig1` (`fig1/`, built here by `fig1.py`) |
| figures | 2 · Clustering | `web/tab_clades4.js` | `/api/runs`, `/api/run` |
| figures | 3 · Transmission | `web/tab_panel2.js` | `/api/panel?n=2` (`tx/`) |
| figures | 4 · Frequency AUC | `web/tab_panel1.js` | `/api/panel?n=1`, `n=1h` |
| figures | 5 · Breakouts | `web/tab_panel3.js` | `/api/explore_panel` |
| supplementary | S1 · All trees | `web/tab_trees.js` | `/api/runs`, `/api/run` |
| supplementary | S2 · Sensitivity to clustering | `web/tab_sa.js` | `/api/across`, `/api/panel?n=1h`, `/api/explore_panel` at each (k, m_T) |
| supplementary | S3 · Sensitivity to outbreak detection | `web/tab_sa.js` | `/api/explore_panel` at each level of `sa_levels.json` |
| data exploration | Clustering | `web/tab_clustering.js` | `/api/runs`, `/api/run` |
| data exploration | m_T effect on clades | `web/tab_across.js` | `/api/across` |
| data exploration | Breakout detection (`#sweep`) | `web/tab_sweep.js` | `/api/explore` |
| data exploration | Cryptic breakouts (`#cryptic`) | `web/tab_sweep.js` | `/api/explore` |
| data exploration | Cryptic by scenario (`#crystats`) | `web/tab_sweep.js` | `/api/explore_panel`, `/api/count_tests` |
| tools | Method | `web/tab_method.js` | `/api/index` |
| tools | Compute | `web/tab_compute.js` | `/api/jobs`, `POST /api/rebuild/*` |
| tools | Export | `web/tab_export.js` | the figure tabs themselves (html2canvas, jsPDF from cdnjs) |

One settings object (`web/settings.js`) mirrored to `localStorage`, so a reload
starts where you left off. A setting is one of three kinds: **build** (`k`,
`win`, `prom`, `lo`, `hi`, `window`) needs a clade build, **read** (`mt`,
`burnin`) is applied over a build, and everything else is a view parameter that
never leaves the browser.

## Data

```
Data/figures/<name>_<stamp>/    one experiment (layout.py, sweep.py --print-paths)
  exploration/                  the store: what panels.py and explore.py read
    cache/<build>/              clade builds, keyed by parameters and a hash of
                                the code that produced them
    cache/_explore/ _panels/    the exploration tabs' caches
    runs/     55 MB             phylogeny + sparse lineage frequency, 250 runs (120 time points),
                                export_clade_explorer.py
    tx/                         per-day 2x2 transmission counts by host type, and per-day
                                jump tables per transmitter type (export_transmissions.py)
  fig1/                         fig1.json + index.json, Figure 1's data (fig1.py build), local only
  fig3/ fig4/ fig5/ cache/      each figure's panel cache (layout.MEMO_HOME)
```

The builds stay at `exploration/cache/<build>`: `panels._build_dir` is in every
panel cache key, so moving them would discard every warmed panel.
`python -I figures_gui/layout.py migrate [--apply] [--prune]` moves a
pre-v2.5.21 store into this layout; `--prune` deletes builds at another code
hash and caches older than the oldest current build.

**`runs/` is not deletable.** It is the only source for the phylogeny, the
Muller plot and the episode grid. The artifact-era `clades/` store and
`runs/formation.json` were deleted in v2.5.4; `formation.json` is regenerable
with `tests/clade_formation_grid.py --json`.

A clade build's key is
`k7_w15_p0.1_lo0.1_hi0.9_post__<hash of compute.py AND clustering.py>`. Both,
because a clade build is the output of both: an uncommitted one-line change to
the clade walk once altered what a clade is while `compute.py` was untouched,
and 24 builds carried the new maths under a key claiming the old.

The panel cache sits in each figure's `cache/` (exploration tabs: `exploration/cache/_panels/`), keyed by the request and by a hash
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

It now reads `/api/explore` (`explore.py`) instead, which measures on the same
daily record with no clade build behind it, so every detection rule moves:
smoothing kernel, window and log space; relative and absolute prominence,
minimum separation and width, censored peaks; the episode start (last trough,
below a fraction of the peak, prominence base); hidden time from emergence in
host or from first transmission; first or largest breakout. Its settings are
its own (`sd*`), never the build parameters. At Figure 5's rules it returns
exactly what `/api/run_daily` does, which `check_explore.py` asserts clade by
clade. `explore.py` is in no cache key, so editing it orphans nothing.

Figure 5 (tab 5, `web/tab_panel3.js`) always measures with the detection tabs'
rules (`/api/explore_panel`, every simulation). The clade build's own episodes
(`/api/panel?n=3b`) are equal to them at Figure 5's rules for every clade of
every run, which `check_explore.py` asserts.

The only slow step is clustering, which depends on k and m_T alone. Per run it
is cached on disk -- the daily record in `cache/_panels/`, the hidden-time
anchors and long-born definer counts in `cache/_explore/` -- by
`warm_explore.py` (about 6 min for 250 runs on 3 local cores), and on the
cluster by `slurm_explore.sh`, which `run_precompute.sh` submits with the grid
and the launcher submits on its own for a stamped experiment. After that a new
rule set costs one detection pass over every run, cached on disk per rule set.

## Rebuilding

Both are jobs: one at a time, loopback only, refused rather than queued.

```bash
# the transmission table. Reads each run's phylogeny for the jump tables,
# ~1.5 s a run; 250 runs took 1 min 40 s at --procs 4
python -I figures_gui/export_transmissions.py --out <store>/tx --procs 2   # <store> = .../exploration

# one clade build -- clades, composition, burden, episodes at one k. ~12 min
python -I figures_gui/compute.py --out Data/figures/impact_long_shedders_1/exploration/cache --procs 2 \
    --k 7 --win 15 --prom 0.10 --window post
```

`fig1.py build --stamp N --store <experiment>` writes `<experiment>/fig1/`: one
process, ~6 min, ~1 GB. The Figure 1 tab starts it on its own (below).

Use `--procs 2`. The box has 4 cores and about 4 GB free with swap already in
use; the exporter costs roughly `0.15 + procs x 0.9` GB.

## Preparing the figure data

One command runs every step after the simulations, on the cluster or here:

```bash
python -I figures_gui/prepare.py --stamp 2 --slurm [--pack] [--dry-run]    # Slurm jobs
python -I figures_gui/prepare.py --stamp 2 --local [--procs 2] [--pack]    # here, in order
#   --local: --only tx,runs  --skip check  --all-builds (the whole grid, ~a day on 4 cores)
#   --slurm: anything after `--` goes to run_precompute.sh (--limit N, --resubmit, ...)
```

A pipeline run does the same once its experiments are sealed:
`run_impact_long_shedders_pipeline.py --prepare-figures [--pack-results]`
(Slurm jobs under `--runner slurm`, `--local` otherwise). `--pack-results`
packs the run's data into `Data_Export/*.tar.xz`, and with `--prepare-figures`
also the figure data into `Data_Export/figures/`. Both are off by default.
`scripts/launchers/run_pipeline_and_grid.sh` runs `prepare.py --slurm` as its stage 2 and
packs both the run data and the figure data unless `--no-compress`.

| step | `--local` runs | `--slurm` job | depends on |
|---|---|---|---|
| check | `run_precompute.sh --check-only` | the same checks, on the login node | -- |
| tx | `export_transmissions.py` | `slurm_tx.sh` (30 min, 8G, 4 cpu) | -- |
| runs | `export_clade_explorer.py` | `slurm_runs.sh` (1 h, 16G, 4 cpu) | -- |
| explore | `warm_explore.py` at the operating point, its Figure 5 panel and S3's levels | `slurm_explore.sh` at every (k, m_T) of the grid (2 h, 16G, 4 cpu), then one task for the Figure 5 panels at each pair and at S4's levels (6 h) | -- |
| grid | `compute.py` at the operating point; `--all-builds`: `sweep.py` | `slurm_precompute.sh`, one build a task (1 h, 16G, 4 cpu) | -- |
| warm | `warm_panels.py` | `slurm_warm.sh`, 32 shards (1 h, 8G, 4 cpu) | grid, tx, explore |
| pack | `pack_cache.sh --here` (`--partial` unless `--all-builds`) | `pack_cache.sh` (3 h, 16G, 8 cpu) | everything, afterany |

`--slurm` is `run_precompute.sh` with `--no-pack` unless `--pack`; calling it
directly still works (`--dry-run`, `--warm-shards 0`, `--resubmit`, `--limit N`,
`--calibrate`, `--check-only`). `--local` skips a step whose prerequisite
failed and exits nonzero if any step failed or was skipped.

**Figure 1 is not a step.** It reads patient data in `Data/RealWorldData/` that
is not on the cluster and cannot be on git. The studio builds it here: opening
Figure 1 with `fig1/` missing, or stale (its code or its input files changed,
`fig1.load`), starts `fig1.py build` as a job (`POST /api/rebuild/fig1`, ~6 min,
~1 GB) and redraws when it ends. `fig1/` is never packed (`layout.SECTIONS`),
because it holds values derived from the patient data.

The pack job (`layout.py pack`) writes `Data_Export/figures/<name>_<stamp>/`:
one `.tar.gz` per section (exploration, fig2..fig5), only
builds at the current hash, split into `_partNN` files of at most 1 GB that each
unpack on their own, plus a manifest with checksums and the extract line. A
missing `tx/` or `runs/` is a warning, a partial grid is a refusal (exit 3).
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
node    figures_gui/check_app.js               #     every tab rendered headlessly
node    figures_gui/check_episodes.js          # 25  the tab's episodes against the build's
python -I figures_gui/check_explore.py         #     /api/explore IS run_daily, explore_panel IS n=3b, at Figure 5's rules
node    figures_gui/test_sa_tabs.js            #     S1-S3: the sensitivity metrics ARE Figure 5's counting, every facet draws
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
