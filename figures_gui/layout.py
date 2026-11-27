"""Where the studio keeps one experiment's data, and how it travels.

    Data/figures/<exp_name>_<stamp>/
        exploration/            the store: panels.py, explore.py, dist.py read it
            cache/<build>/      the clade grid, one build per parameter set
            cache/_explore/     the exploration tabs' cache (explore.py)
            cache/_panels/      the exploration tabs' panel cache (across, daily, episodes)
            tx/                 the transmission table (export_transmissions.py)
            runs/               the per-seed export behind the Muller plots
        fig1/                   Figure 1 (fig1.py), built locally, never packed
        fig3/ fig4/ fig5/       cache/: that figure's panel cache

Packed into Data_Export/figures/<exp_name>_<stamp>/, one archive per section,
none over 1 GB; a section that does not fit is split into parts that each
unpack on their own.

    python -I figures_gui/layout.py migrate [--apply] [--prune]
    python -I figures_gui/layout.py pack --stamp 2 [--exp-name NAME]
"""
import argparse, datetime, gzip, hashlib, json, os, shutil, subprocess, sys, tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BASE = os.path.join(ROOT, "Data", "figures")
EXPORT = os.path.join(ROOT, "Data_Export", "figures")
EXPLORATION = "exploration"
# what pack archives. Not fig1/: it holds values derived from patient data that
# may not leave this machine, and the studio builds it locally (fig1.py)
SECTIONS = (EXPLORATION, "fig2", "fig3", "fig4", "fig5")
# panel-cache tag -> the figure it serves; tags not listed serve the exploration tabs
MEMO_HOME = {"p1": "fig4", "p1c": "fig4", "p1h": "fig4", "d1": "fig4",
             "p2": "fig3", "d2": "fig3",
             "p3": "fig5", "p3b": "fig5", "d3": "fig5"}
LIMIT = 1_000_000_000


def exp_tag(stamp, exp_name):
    return f"{exp_name}_{stamp}"


def root_dir(base, stamp, exp_name):
    return os.path.join(os.path.expanduser(base), exp_tag(stamp, exp_name))


def store_dir(base, stamp, exp_name):
    return os.path.join(root_dir(base, stamp, exp_name), EXPLORATION)


def root_of(store):
    s = os.path.normpath(store)
    return os.path.dirname(s) if os.path.basename(s) == EXPLORATION else s


def memo_dir(store, tag):
    home = MEMO_HOME.get(tag)
    s = os.path.normpath(store)
    if home and os.path.basename(s) == EXPLORATION:
        return os.path.join(os.path.dirname(s), home, "cache")
    return os.path.join(s, "cache", "_panels")


def memo_dirs(store):
    return sorted({memo_dir(store, t) for t in list(MEMO_HOME) + [None]})


def _tag(fname):
    return fname.rsplit("_", 1)[0] if fname.endswith(".json") else None


def _code_hash():
    sys.path.insert(0, HERE)
    import compute
    return compute.code_hash()


def _du(path):
    if os.path.isfile(path):
        return os.path.getsize(path)
    return sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(path) for f in fs)


def _move(src, dst, apply, log):
    """Rename src to dst; into an existing dst, file by file, never overwriting."""
    if not os.path.exists(src):
        return
    if os.path.isdir(src) and not any(fs for _, _, fs in os.walk(src)):
        log.append(f"remove {src}  (empty)")
        if apply:
            shutil.rmtree(src)
        return
    if not os.path.exists(dst):
        log.append(f"move   {src} -> {dst}")
        if apply:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            os.rename(src, dst)
        return
    if os.path.isdir(src) and os.path.isdir(dst):
        for n in sorted(os.listdir(src)):
            _move(os.path.join(src, n), os.path.join(dst, n), apply, log)
        if apply and not os.listdir(src):
            os.rmdir(src)
        return
    log.append(f"KEEP   {src}  ({dst} exists)")


# ---------------------------------------------------------------- migrate

def legacy_stores(base):
    """(old store, exp_name, stamp) for every pre-v2.5.21 store under base."""
    out = []
    if os.path.isdir(os.path.join(base, "cache")):
        sys.path.insert(0, HERE)
        import compute
        out.append((base, compute.EXP_NAME, "1"))
    for n in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        if n.startswith("exp_") and os.path.isdir(os.path.join(base, n)):
            name, stamp = n[4:].rsplit("_", 1)
            out.append((os.path.join(base, n), name, stamp))
    return out


def route_memos(store, apply, log):
    """Panel caches out of exploration/cache/_panels into the figure they serve."""
    src = os.path.join(store, "cache", "_panels")
    if not os.path.isdir(src):
        return
    for f in sorted(os.listdir(src)):
        d = memo_dir(store, _tag(f))
        if d != src:
            if apply:
                os.makedirs(d, exist_ok=True)
                os.rename(os.path.join(src, f), os.path.join(d, f))
            log.append(f"route  {f} -> {os.path.relpath(d, root_of(store))}")


def prune(store, apply, log):
    """Builds at another code hash, and caches older than the oldest current build.

    A cache is keyed on the code that made it, so a stale one is never served;
    this is about disk, not correctness. A cache written before the first
    current build existed cannot have read one.
    """
    h = _code_hash()
    cache = os.path.join(store, "cache")
    names = sorted(os.listdir(cache)) if os.path.isdir(cache) else []
    cur = [os.path.join(cache, n, "index.json") for n in names if n.endswith("__" + h)]
    cur = [f for f in cur if os.path.isfile(f)]
    freed = 0
    for n in names:
        p = os.path.join(cache, n)
        why = ("build at another code hash" if os.path.isdir(p) and not n.startswith("_")
               and not n.endswith("__" + h) else "nothing reads it" if os.path.isfile(p) else None)
        if why:
            freed += _du(p)
            log.append(f"delete {p}  ({why})")
            if apply:
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    t0 = min(os.path.getmtime(f) for f in cur) if cur else float("inf")
    for d in memo_dirs(store) + [os.path.join(cache, "_explore")]:
        old = [f for f in os.listdir(d) if os.path.getmtime(os.path.join(d, f)) < t0] \
            if os.path.isdir(d) else []
        if old:
            freed += sum(os.path.getsize(os.path.join(d, f)) for f in old)
            log.append(f"delete {len(old)} caches in {d}  "
                       + ("(older than the current builds)" if cur else f"(no build at {h})"))
            if apply:
                for f in old:
                    os.remove(os.path.join(d, f))
    return freed


def _roots(base):
    return sorted(os.path.join(base, n) for n in os.listdir(base)
                  if os.path.isdir(os.path.join(base, n, EXPLORATION))) if os.path.isdir(base) else []


def migrate(base, apply, do_prune):
    log, freed = [], 0
    legacy = legacy_stores(base)
    if do_prune:
        for old, _, _ in legacy:
            freed += prune(old, apply, log)
        for root in _roots(base):
            freed += prune(os.path.join(root, EXPLORATION), apply, log)
    for old, name, stamp in legacy:
        root = root_dir(base, stamp, name)
        store = os.path.join(root, EXPLORATION)
        for sub in ("cache", "tx", "runs"):
            _move(os.path.join(old, sub), os.path.join(store, sub), apply, log)
        _move(os.path.join(old, "fig1"), os.path.join(root, "fig1"), apply, log)
        if apply and old != base and os.path.isdir(old) and not os.listdir(old):
            os.rmdir(old)
    for root in _roots(base):
        route_memos(os.path.join(root, EXPLORATION), apply, log)
    return log, freed


# ---------------------------------------------------------------- pack

def section_files(root, section, h):
    """(path on disk, name in the archive) for one section, current builds only."""
    top = os.path.join(root, section)
    base = os.path.dirname(root)
    out = []
    for d, dirs, fs in os.walk(top):
        if section == EXPLORATION and os.path.dirname(d) == top and os.path.basename(d) == "cache":
            dirs[:] = sorted(n for n in dirs if n.startswith("_") or n.endswith("__" + h))
            fs = []
        else:
            dirs.sort()
        for f in sorted(fs):
            if not f.endswith(".tmp"):
                p = os.path.join(d, f)
                out.append((p, os.path.relpath(p, base)))
    return out


class _Part:
    def __init__(self, path):
        self.path, self.raw = path, open(path, "wb")
        self.gz = gzip.GzipFile(fileobj=self.raw, mode="wb", compresslevel=6, mtime=0)
        self.tar, self.n = tarfile.open(fileobj=self.gz, mode="w"), 0

    def close(self):
        self.tar.close(); self.gz.close(); self.raw.close()


def pack(base, stamp, exp_name, out, limit=LIMIT):
    h, root = _code_hash(), root_dir(base, stamp, exp_name)
    tag, day = exp_tag(stamp, exp_name), datetime.date.today().strftime("%Y%m%d")
    dest = os.path.join(out, tag)
    os.makedirs(dest, exist_ok=True)
    made = []
    for section in SECTIONS:
        files = section_files(root, section, h)
        if not files:
            continue
        stem, parts, cur = f"{tag}_{section}_{h}_{day}", [], None
        for p, arc in files:
            size = os.path.getsize(p)
            if cur and cur.n and cur.raw.tell() + size + (1 << 20) > limit:
                cur.close(); cur = None
            if cur is None:
                cur = _Part(os.path.join(dest, f"{stem}_part{len(parts) + 1:02d}.tar.gz"))
                parts.append(cur)
            cur.tar.add(p, arcname=arc, recursive=False)
            cur.n += 1
        cur.close()
        if len(parts) == 1:
            os.replace(parts[0].path, os.path.join(dest, f"{stem}.tar.gz"))
            parts[0].path = os.path.join(dest, f"{stem}.tar.gz")
        for q in parts:
            made.append((section, q.path, q.n))
        print(f"[pack] {section}: {len(files)} files in {len(parts)} archive(s)", flush=True)
    return h, dest, made


def manifest(h, dest, made, stamp, exp_name, extra):
    tag, day = exp_tag(stamp, exp_name), datetime.date.today().strftime("%Y%m%d")
    rows = []
    for section, p, n in made:
        s = hashlib.sha256()
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                s.update(b)
        rows.append(f"  {os.path.basename(p)}  {os.path.getsize(p)/1e6:8.1f} MB  {n:6d} files  sha256 {s.hexdigest()}")
    try:
        commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        commit = "unknown"
    missing = [s for s in SECTIONS if s not in {m[0] for m in made}]
    txt = "\n".join([
        f"studio data, {exp_name} #{stamp}",
        f"packed       {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"build hash   {h}",
        f"commit       {commit}",
        *[f"{k:<12} {v}" for k, v in extra],
        f"sections     {', '.join(sorted({m[0] for m in made}, key=SECTIONS.index))}",
        f"absent       {', '.join(missing) or 'none'}",
        "",
        *rows,
        "",
        "extract on the studio box, from the repo root (every part unpacks on its own):",
        f"    for f in {tag}_*_{h}_{day}*.tar.gz; do tar -xzf \"$f\" -C Data/figures; done",
        "",
    ])
    f = os.path.join(dest, f"{tag}_{h}_{day}.manifest.txt")
    open(f, "w").write(txt)
    return f, txt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("migrate", help="move pre-v2.5.21 stores into this layout")
    m.add_argument("--base", default=BASE)
    m.add_argument("--apply", action="store_true", help="without it, only print the plan")
    m.add_argument("--prune", action="store_true", help="also delete stale builds and caches")
    p = sub.add_parser("pack", help="archive one experiment into Data_Export/figures")
    p.add_argument("--base", default=BASE)
    p.add_argument("--out", default=EXPORT)
    p.add_argument("--stamp", required=True)
    p.add_argument("--exp-name", default=None)
    p.add_argument("--limit-mb", type=int, default=LIMIT // 1_000_000)
    p.add_argument("--note", action="append", default=[], help="KEY=VALUE line for the manifest")
    a = ap.parse_args()
    if a.cmd == "migrate":
        log, freed = migrate(a.base, a.apply, a.prune)
        moves = [l for l in log if not l.startswith("route")]
        routed = len(log) - len(moves)
        print("\n".join(moves) or "nothing to move")
        if routed:
            print(f"route  {routed} figure panel caches out of exploration/cache/_panels")
        if a.prune:
            print(f"frees  {freed/1e9:.2f} GB")
        if not a.apply:
            print("dry run: nothing changed. --apply to do it.")
        return 0
    if a.exp_name is None:
        sys.path.insert(0, HERE)
        import compute
        a.exp_name = compute.EXP_NAME
    h, dest, made = pack(a.base, a.stamp, a.exp_name, a.out, a.limit_mb * 1_000_000)
    if not made:
        print(f"[pack] nothing under {root_dir(a.base, a.stamp, a.exp_name)}", file=sys.stderr)
        return 1
    f, txt = manifest(h, dest, made, a.stamp, a.exp_name, [n.split("=", 1) for n in a.note])
    print(txt)
    print(f"[pack] manifest {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
