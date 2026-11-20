"""The AST digest that replaces the whole-file hash, written once and imported
by both compute.py and panels.py.

A build and a panel memo are named after the code that produced them, so that a
change to the maths cannot be served out of a stale answer -- that is on the
record: an uncommitted one-line change to the clade walk altered what a clade
is while compute.py was untouched, and 24 builds carried the new maths under a
key claiming the old.

Hashing the raw bytes also names everything that is NOT the maths. A comment, a
docstring, a blank line or a reformat moved the key and orphaned the lot: 168
clade builds, an HPC array resubmit, and 2 GB of panel memos. The stale
`mt percent` docstring in compute.py sat unfixed for exactly that reason.

So hash what the code DOES: parse it, drop docstrings, dump the tree without
line numbers. Comments never reach the AST, docstrings are removed here, and
`include_attributes=False` drops lineno/col_offset -- which matters most,
because inserting one comment shifts every line number below it.

What this does NOT make free: renaming a local variable, since names are in the
AST. That is the right trade -- a rename is a code edit, and the cheap thing to
get wrong is prose.
"""
import ast
import hashlib
import inspect
import os
import textwrap


def _strip_docstrings(tree):
    """Every docstring out of the tree. Collected first, then mutated: ast.walk
    yields lazily, so editing a body mid-walk changes what it visits."""
    holders = [n for n in ast.walk(tree)
               if isinstance(n, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef))]
    for node in holders:
        body = node.body
        if (body and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            node.body = body[1:] or [ast.Pass()]
    return tree


_digests = {}


def ast_digest(path):
    """One file's behaviour as a string, free of comments, docstrings, blank
    lines and formatting. Cached on (mtime, size): parsing costs ~25 ms and
    panel 4 asks once per build row, 8 s a request uncached."""
    st = os.stat(path)
    key = (os.path.realpath(path), st.st_mtime_ns, st.st_size)
    if key not in _digests:
        _digests[key] = _ast_digest(path)
    return _digests[key]


def _ast_digest(path):
    with open(path, "rb") as fh:
        src = fh.read()
    tree = _strip_docstrings(ast.parse(src))
    return ast.dump(tree, annotate_fields=True, include_attributes=False)


def files_hash(paths, extra=""):
    """A short digest over several files' behaviour, in the order given.
    `extra` is for a caller that wants to version the scheme itself."""
    h = hashlib.sha256()
    for p in paths:
        h.update(ast_digest(p).encode())
        h.update(b"\0")
    if extra:
        h.update(extra.encode())
    return h.hexdigest()[:12]


def source_digest(obj):
    """One function's or class's behaviour, free of comments and docstrings."""
    tree = _strip_docstrings(ast.parse(textwrap.dedent(inspect.getsource(obj))))
    return ast.dump(tree, annotate_fields=True, include_attributes=False)


def names_hash(namespace, names, extra=""):
    """A digest over the behaviour of the NAMED functions only, so a memo is
    invalidated by the code that actually computes it and not by its
    neighbours. `names` is sorted here, so the declaration order of a
    dependency list cannot move a key."""
    h = hashlib.sha256()
    for n in sorted(names):
        h.update(n.encode())
        h.update(source_digest(namespace[n]).encode())
        h.update(b"\0")
    if extra:
        h.update(extra.encode())
    return h.hexdigest()[:12]
