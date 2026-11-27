"""layout.py keeps sets apart: prune drops stale builds, caches and stray
sections; pack leaves one set in its folder; unpack checks the manifest and
replaces each section whole, refusing bad copies, a foreign hash or paths
outside <tag>/<section>/. Runs on a throwaway store, never on Data/.

    python -I figures_gui/test_layout_pack.py [--dir DIR]
"""
import argparse, io, os, shutil, sys, tarfile, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import layout

H, TAG = "testhash0001", "test_layout_pack_2"
ok = fail = 0


def check(cond, what):
    global ok, fail
    ok, fail = ok + bool(cond), fail + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {what}")


def put(path, text="{}", age=0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(text)
    if age:
        t = time.time() - age
        os.utime(path, (t, t))


def tree(root):
    return sorted(os.path.relpath(os.path.join(d, f), root)
                  for d, _, fs in os.walk(root) for f in fs) if os.path.isdir(root) else []


def store(base):
    r = os.path.join(base, TAG)
    put(f"{r}/exploration/cache/b1__{H}/index.json", '{"b": 1}', age=100)
    put(f"{r}/exploration/cache/b0__oldhash0000/index.json", age=900)
    put(f"{r}/exploration/cache/_explore/new.json", age=10)
    put(f"{r}/exploration/cache/_explore/old.json", age=500)
    put(f"{r}/exploration/runs/seed.json")
    put(f"{r}/fig3/cache/p2_a.json", age=10)
    put(f"{r}/supplementary/S1/cache/x.json")
    put(f"{r}/fig1/fig1.json", '{"local": true}')
    return r


def refused(fn):
    try:
        fn()
    except SystemExit as e:
        return str(e)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=None)
    a = ap.parse_args()
    top = tempfile.mkdtemp(prefix="test_layout_pack_", dir=a.dir)
    layout._code_hash = lambda: H
    try:
        base, out, base2 = (os.path.join(top, n) for n in ("test_store", "test_export", "test_studio"))
        root = store(base)

        print("prune")
        before = tree(root)
        layout.prune(os.path.join(root, "exploration"), False, [])
        check(tree(root) == before, "dry run changes nothing")
        layout.prune(os.path.join(root, "exploration"), True, [])
        t = tree(root)
        check(not any("oldhash" in p for p in t), "build at another hash deleted")
        check("exploration/cache/_explore/old.json" not in t, "cache older than the builds deleted")
        check("exploration/cache/_explore/new.json" in t and "fig3/cache/p2_a.json" in t, "current caches kept")
        check(not os.path.exists(f"{root}/supplementary"), "stray section deleted")
        check("fig1/fig1.json" in t, "local fig1 kept")
        layout._code_hash = lambda: "otherhash999"
        before = tree(root)
        layout.prune(os.path.join(root, "exploration"), True, log := [])
        check(tree(root) == before and any("nothing deleted" in l for l in log),
              "no build at the current hash: nothing deleted")
        layout._code_hash = lambda: H

        print("pack")
        put(f"{out}/{TAG}/{TAG}_exploration_oldhash0000_20200101.tar.gz", "old")
        put(f"{out}/{TAG}/{TAG}_oldhash0000_20200101.manifest.txt", "old")
        h, dest, made = layout.pack(base, "2", "test_layout_pack", out)
        mf, _ = layout.manifest(h, dest, made, "2", "test_layout_pack", [])
        gone = layout.drop_other_sets(dest, {p for _, p, _ in made} | {mf})
        left = sorted(os.listdir(dest))
        check(len(gone) == 2 and all("oldhash" not in n for n in left), "older set deleted")
        check(not any("fig1" in n for n in left), "fig1 never packed")
        check({s for s, _, _ in made} == {"exploration", "fig3"}, "sections packed: exploration, fig3")

        print("unpack")
        r2 = os.path.join(base2, TAG)
        put(f"{r2}/exploration/cache/b0__oldhash0000/index.json")
        put(f"{r2}/exploration/cache/_explore/stale.json")
        put(f"{r2}/fig2/cache/junk.json")
        put(f"{r2}/supplementary/x.json")
        put(f"{r2}/fig1/fig1.json", '{"mine": true}')
        root2, _, log = layout.unpack(mf, base2)
        exp = sorted(p for p in tree(root) if not p.startswith("fig1/"))
        got = tree(root2)
        check(sorted(p for p in got if not p.startswith("fig1/")) == exp, "store equals what was packed")
        check(open(f"{r2}/fig1/fig1.json").read() == '{"mine": true}', "local fig1 untouched")
        check(not os.path.exists(f"{r2}/fig2") and not os.path.exists(f"{r2}/supplementary"),
              "absent section and stray deleted")
        check(not any(n.startswith(".unpack_") for n in os.listdir(base2)), "staging removed")

        before = tree(root2)
        arch = os.path.join(dest, next(n for n in os.listdir(dest) if "_fig3_" in n))
        good = open(arch, "rb").read()
        open(arch, "ab").write(b"x")
        msg = refused(lambda: layout.unpack(mf, base2))
        check(msg and "checksum" in msg and tree(root2) == before, "altered archive refused, nothing changed")
        open(arch, "wb").write(good)
        layout._code_hash = lambda: "otherhash999"
        msg = refused(lambda: layout.unpack(mf, base2))
        check(msg and "build hash" in msg and tree(root2) == before, "foreign build hash refused")
        layout._code_hash = lambda: H

        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as t:
            for name in (f"{TAG}/fig3/cache/ok.json", f"{TAG}/fig3/../../escape.json"):
                ti = tarfile.TarInfo(name); ti.size = 2
                t.addfile(ti, io.BytesIO(b"{}"))
        open(arch, "wb").write(buf.getvalue())
        txt = open(mf).read()
        lines = [l if "_fig3_" not in l else
                 f"  {os.path.basename(arch)}  0.0 MB  {2:6d} files  sha256 {layout._sha256(arch)}"
                 for l in txt.splitlines()]
        open(mf, "w").write("\n".join(lines) + "\n")
        msg = refused(lambda: layout.unpack(mf, base2))
        check(msg and "outside" in msg and tree(root2) == before
              and not os.path.exists(os.path.join(base2, "escape.json")), "path outside the tag refused")
    finally:
        shutil.rmtree(top, ignore_errors=True)
    print(f"\n{ok} passed, {fail} failed")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
