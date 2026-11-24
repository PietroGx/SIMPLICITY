#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import pandas as pd
import simplicity.output_manager as om
import os
from typing import Dict, Any, List

def _serialize_definers_for_clade(def_df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Convert a 'defining mutations for clade' dataframe into a list[dict] with
    stable keys.
    Expected columns:
      - 'mutation'
      - 'emergence_lineage' 
      - 'host_type'         
      - 'time'              
    """
    if def_df is None or def_df.empty:
        return []

    cols = {c.lower(): c for c in def_df.columns}
    get = lambda row, key, default="": row.get(cols.get(key, key), default)

    out: List[Dict[str, Any]] = []
    for row in def_df.to_dict("records"):
        rec = {
            "mutation": get(row, "mutation", ""),
            "emergence_lineage": get(row, "emergence_lineage", ""),
            "host_type": get(row, "host_type", ""),
            "time": float(get(row, "time", 0.0) or 0.0),
        }
        out.append(rec)
    return out

# Rule B's operating point. m_T is a share of the run's own transmissions; see
# GLOSSARY.md. Named here so the pipeline, the figure builder and the studio
# cannot drift to different defaults. A FRACTION since v2.5.15, not a
# percentage: 0.01 is the same floor the old MIN_TX_FRAC = 1.0 computed.
MIN_TX_FRAC = 0.01


def _table_matches(path, shared_mut_threshold, min_tx_frac) -> bool:
    """Was this clustering CSV built by this (k, min_tx_frac)?

    A table written before v2.5.7 carries neither column, so it cannot be shown
    to match and is rebuilt. Any read failure counts as no match too:
    recomputing costs about 1.4s for one run, while serving the wrong clade
    definition costs a figure.
    """
    try:
        head = pd.read_csv(path, nrows=1)
    except Exception:
        return False
    if len(head) == 0:
        return False
    for col, want in (("shared_mut_threshold", shared_mut_threshold),
                      ("min_tx_frac", min_tx_frac)):
        if col not in head.columns:
            return False
        try:
            if float(head[col].iloc[0]) != float(want):
                return False
        except Exception:
            return False
    return True


def get_clustering_for_ssod(
    experiment_name: str,
    seeded_simulation_output_dir: str,
    phylogenetic_df: pd.DataFrame | None = None,
    shared_mut_threshold: int = 5,
    min_tx_frac: float = MIN_TX_FRAC,
    force: bool = False,
) -> str:
    """
    If the clustering CSV on disk was built by THIS definition (and not force),
    returns it. Otherwise computes the clustering and writes a single CSV.

    The table carries rule B: a clade must clear the divergence threshold AND
    have transmitted at least `min_tx_frac` of the run's own transmission
    events. Divergence alone makes a clade of almost every tip -- measured on
    control at k=5, 402 clades per run of which 86% hold one lineage and 78%
    infected nobody, against 14 under rule B.

    The definition is written into the table and checked on the way back out,
    because the filename does not carry it: the path is
    `07_Clustering/{SOD}_seed_{SEED}_clustering.csv`, with no k and no
    `min_tx_frac` in it. Returning the file whenever it merely EXISTS was wrong
    twice over -- a table built at k=5 was served to a caller asking for k=7,
    and a divergence-only table written before rule B was served as though it
    carried rule B. Neither is visible to the caller: both are well-formed CSVs.
    """
    path = om.get_clustering_table_filepath(experiment_name, seeded_simulation_output_dir)
    if (not force) and os.path.exists(path) and _table_matches(
            path, shared_mut_threshold, min_tx_frac):
        return path

    if phylogenetic_df is None:
        phylogenetic_df = om.read_phylogenetic_data(seeded_simulation_output_dir)

    support = min_support = None
    if min_tx_frac > 0:
        ind = om.read_individuals_data(seeded_simulation_output_dir)
        support = subtree_transmission_counts(phylogenetic_df, ind)
        min_support = transmission_floor(ind, min_tx_frac)

    clade_to_lineages, lineage_to_clade, per_clade_mut_df, clade_meta_df = \
        cluster_lin_into_clades_with_meta(phylogenetic_df,
                                         shared_mut_threshold=shared_mut_threshold,
                                         support=support,
                                         min_support=min_support or 0.0)
        
    # Labelled BEFORE the set_index below: label_clades_from_definers reads
    # 'clade' as a COLUMN, so the indexed frame raised KeyError('clade') on
    # every call and the bare handler turned that into no labels at all.
    # Definer-less NON-root clades are withheld from the vote, as
    # the retired tests/export_fig2_builder.py and figures_gui/compute.py do -- at low k
    # the function refuses them by design and they carry n_defining = 0. Every
    # other refusal, an unknown host_type above all, propagates.
    par_of = dict(zip(clade_meta_df["clade"], clade_meta_df["parent_clade"]))
    votable = {c: d for c, d in per_clade_mut_df.items()
               if not ((d is None or len(d) == 0) and isinstance(par_of.get(c), str))}
    labels_s, _ = label_clades_from_definers(votable, clade_meta_df)

    # Set 'clade' as the index
    if isinstance(clade_meta_df, pd.DataFrame) and "clade" in clade_meta_df.columns:
        clade_meta_df = clade_meta_df.set_index("clade")

    rows = []
    meta_cols = set(clade_meta_df.columns) if isinstance(clade_meta_df, pd.DataFrame) else set()

    for clade, lineages in clade_to_lineages.items():
        lineages = sorted(lineages)

        if isinstance(clade_meta_df, pd.DataFrame) and clade in clade_meta_df.index:
            meta = clade_meta_df.loc[clade]
            root_lineage = meta["root_lineage"] if "root_lineage" in meta_cols else meta.get("root_lineage", "")
            parent_clade = meta["parent_clade"] if "parent_clade" in meta_cols else meta.get("parent_clade", "")
            start_time = float(meta["start_time"]) if "start_time" in meta_cols else float(meta.get("start_time", 0.0))
        else:
            root_lineage, parent_clade, start_time = "", "", 0.0

        def_df = per_clade_mut_df.get(clade, pd.DataFrame())
        definers_list = _serialize_definers_for_clade(def_df)

        rows.append({
            # the definition that produced this table; _table_matches reads it
            # back before the file is ever served from cache
            "shared_mut_threshold": int(shared_mut_threshold),
            "min_tx_frac": float(min_tx_frac),
            "clade": clade,
            "root_lineage": root_lineage,
            "parent_clade": parent_clade,
            "start_time": start_time,
            "label": labels_s.get(clade, ""),
            "n_lineages": len(lineages),
            "n_defining": len(definers_list),
            "n_def_standard": sum(1 for d in definers_list if str(d.get("host_type","")).lower().startswith("standard")),
            "n_def_long":   sum(1 for d in definers_list if "long" in str(d.get("host_type","")).lower()),
            "lineages": repr(lineages),
            "defining_mutations": repr(definers_list),
        })

    df = pd.DataFrame(rows).sort_values(["start_time", "clade"]).reset_index(drop=True)
    om.write_clustering_table(experiment_name, seeded_simulation_output_dir, df)
    return path

def parse_genome(genome):
    """Convert a genome ({position: base}) into a frozenset of (position, base).

    No bare except: a malformed genome must raise rather than silently yield an
    empty set, which would collapse every clade into one.
    """
    if not genome:
        return frozenset()
    return frozenset(genome.items())
    
def build_lineage_to_mutation_dict(phylogenetic_data_df):
    """Return {Lineage_name -> frozenset of mutation tuples} from the phylo table."""
    return {
        row["Lineage_name"]: parse_genome(row["Genome"])
        for _, row in phylogenetic_data_df.iterrows()
    }

def subtree_transmission_counts(phylogenetic_data_df, individuals_df):
    """{lineage -> onward infections produced by it and its descendants}.

    The support term of rule B. Cumulative down the tree, so a refusal is
    hereditary: a lineage whose whole subtree transmitted too little cannot
    have a descendant that transmitted more. Ported unchanged from
    tests/clade_tx_minimum.py, which is where the rule was measured.
    """
    par = {l: (None if (not isinstance(p, str) or p == l) else p)
           for l, p in zip(phylogenetic_data_df["Lineage_name"],
                           phylogenetic_data_df["Lineage_parent"])}
    order = list(phylogenetic_data_df
                 .sort_values("Time_emergence", kind="mergesort")["Lineage_name"])
    tx = individuals_df["inherited_lineage"].value_counts().to_dict()
    stx = {l: int(tx.get(l, 0)) for l in order}
    for l in reversed(order):                      # children before parents
        p = par.get(l)
        if p is not None and p in stx:
            stx[p] += stx[l]
    return stx


def transmission_total(individuals_df):
    """The run's own transmission events: inherited_lineage, one row per host."""
    return int(individuals_df["inherited_lineage"].value_counts().sum())


def transmission_floor(individuals_df, min_tx_frac):
    """m_T as an absolute count for one run.

    m_T is a SHARE of the run's own transmissions, not a count: the scenarios
    run 11,846 to 27,510 transmissions, so a fixed count is a different filter
    in each. A FRACTION since v2.5.15 -- 0.01, not 1.0, for one percent.
    """
    return float(min_tx_frac) * transmission_total(individuals_df)


def normalise_host_type(h):
    """Host type as one of "long" or "standard".

    Was written four times -- _host_of in the retired tests/export_fig2_builder.py,
    _is_long in figures_gui/compute.py, _norm_host inside the clade builder,
    and inline in the browser. Every copy tested the same "long" prefix, which
    is the only thing that makes them agree; this is the one that is allowed to
    change.
    """
    if h is None:
        return "standard"
    v = str(h).strip().lower()
    return "long" if v.startswith("long") else "standard"


def resolve_clade_label(clade_labels, clade):
    """A clade's origin label, normalised to long / standard / mixed / founder.

    Restated here from _resolve_label in
    scripts/long_paper_figures/long_shedders_preprocess.py, which left the
    studio importing a clade rule from a figure script. That script still
    carries its own copy: collapsing the two means editing scripts/, which is
    out of scope for figure 2 (see BACKLOG).
    """
    if not clade:
        return "standard"
    raw = str(clade_labels.get(clade, "standard")).lower()
    if "long" in raw:
        return "long"
    if "mix" in raw:
        return "mixed"
    if "founder" in raw:
        return "founder"
    return "standard"


def founds_clade(lineage, clade_root_mutset, lin2mut, shared_mut_threshold,
                 support=None, min_support=0.0):
    """Rule B, as one predicate: does this lineage start a new clade?

    Both walks in this module ask it -- the metadata-building recursion in
    `cluster_lin_into_clades_with_meta` and the flat `designate_clades` -- so
    the rule has one statement rather than one per walk.

      divergence   >= shared_mut_threshold mutations from the clade root,
                   measured on genome SETS, not mutation counts, which disagree
                   on seeds with back-mutations
      support      and, when min_support > 0, at least that much support, which
                   is normally the onward infections of this lineage and its
                   descendants
    """
    if len(lin2mut.get(lineage, frozenset()) - clade_root_mutset) < shared_mut_threshold:
        return False
    if min_support > 0:
        return (support or {}).get(lineage, 0) >= min_support
    return True


def lineage_tree(phylogenetic_data_df, individuals_df=None,
                 order_by="emergence"):
    """Parent map, walk order, and the two cumulative subtree counts.

    One replacement for the three `_tree` functions in tests/ that each built
    this again (clade_tx_minimum, clade_formation_grid, clade_designation_rules).

    Returns (lin2mut, par, order, sub, stx): subtree lineage counts and subtree
    onward infections. `stx` is all zeros without `individuals_df`.

    `order_by` picks the walk order, and the two are not interchangeable:

      "emergence"  sorted by Time_emergence, which is what the clade build and
                   `cluster_lin_into_clades_with_meta` use
      "topology"   down from the roots, so a parent is always placed before its
                   children whatever the timestamps say

    `designate_clades` reads a lineage's parent's clade, so a child placed
    before its parent becomes its own clade root. Under "emergence" that is the
    build's own behaviour and is deliberate; a diagnostic that wants the
    partition the topology implies asks for "topology".
    """
    lin2mut = build_lineage_to_mutation_dict(phylogenetic_data_df)
    par = {l: (None if (not isinstance(p, str) or p == l) else p)
           for l, p in zip(phylogenetic_data_df["Lineage_name"],
                           phylogenetic_data_df["Lineage_parent"])}
    if order_by == "emergence":
        order = list(phylogenetic_data_df
                     .sort_values("Time_emergence", kind="mergesort")["Lineage_name"])
    elif order_by == "topology":
        kids = {}
        roots = []
        for l, p in par.items():
            if p is None:
                roots.append(l)
            else:
                kids.setdefault(p, []).append(l)
        order, stack = [], list(roots)
        while stack:
            v = stack.pop()
            order.append(v)
            stack.extend(kids.get(v, ()))
    else:
        raise ValueError(
            f'order_by must be "emergence" or "topology", got {order_by!r}')
    tx = ({} if individuals_df is None
          else individuals_df["inherited_lineage"].value_counts().to_dict())
    sub = {l: 1 for l in order}
    stx = {l: int(tx.get(l, 0)) for l in order}
    for l in reversed(order):                      # children before parents
        p = par.get(l)
        if p is not None and p in sub:
            sub[p] += sub[l]
            stx[p] += stx[l]
    return lin2mut, par, order, sub, stx


def designate_clades(lin2mut, par, order, shared_mut_threshold,
                     support=None, min_support=0.0):
    """{lineage -> its clade's root lineage}, by rule B.

    The lightweight walk: no clade names, no defining mutations, no metadata.
    `cluster_lin_into_clades_with_meta` is the one that builds those; this is
    for the diagnostics that only need the partition, and both reach the same
    decision because both call `founds_clade`.
    """
    croot = {}
    for v in order:
        p = par.get(v)
        if p is None or p not in croot:
            croot[v] = v
            continue
        cr = croot[p]
        croot[v] = v if founds_clade(v, lin2mut.get(cr, frozenset()), lin2mut,
                                     shared_mut_threshold, support,
                                     min_support) else cr
    return croot


def clade_max_divergence(member_lineages, root_lineage, lin2mut):
    """Mutations between the clade's furthest member and its own root.

    A bound, so it is the clade's worst case, not its typical one. Derived
    independently in the retired tests/export_fig2_builder.py and figures_gui/compute.py
    before this.
    """
    rmut = lin2mut.get(root_lineage, frozenset())
    return max((len(lin2mut.get(m, frozenset()) - rmut) for m in member_lineages),
               default=0)


def keep_clades_by_transmission(clades, n_tx_run, min_tx_frac,
                                ntx=lambda c: c["ntx"],
                                founding=lambda c: c["founding"],
                                n_defining=lambda c: c["ndef"]):
    """Rule B's support term applied to already-built clades.

    The studio stores the divergence-only superset so m_T can be varied without
    re-clustering, then filters here. That filter IS half the clade definition,
    so it belongs beside the other half rather than restated in the panel code.

    The founding clade is kept and labelled separately -- it has no defining
    mutations and no parent, so neither term applies to it. A group with no
    defining mutations and a parent is not a clade.
    """
    floor = float(min_tx_frac) * n_tx_run
    out = []
    for c in clades:
        if founding(c):
            out.append(c)
            continue
        if n_defining(c) <= 0:
            continue
        if ntx(c) <= 0 or ntx(c) < floor:
            continue
        out.append(c)
    return out


def keep_designated_clades(clades,
                           founding=lambda c: c["founding"],
                           n_defining=lambda c: c["ndef"]):
    """Structural filter for clades ALREADY designated under rule B.

    The counterpart to `keep_clades_by_transmission`, and the two must not be
    confused:

      keep_clades_by_transmission   designation used divergence alone, so the
                                    support term is applied here, by DROPPING
                                    clades under m_T and discarding their
                                    lineages
      keep_designated_clades        designation used rule B, so the support
                                    term is already satisfied and every lineage
                                    is already FOLDED into a surviving clade

    Re-applying the floor after rule B designation would be wrong, not merely
    redundant. Designation tests the founding lineage's SUBTREE support, which
    includes descendants that later found clades of their own; a clade's own
    `ntx` counts only the lineages that stayed. So a legitimate rule B clade
    can carry fewer than m_T transmissions itself, and the floor would drop it.

    What is still structural: the founding clade is kept and labelled -- it has
    no defining mutations and no parent, so neither term applies to it -- and a
    group with no defining mutations and a parent is not a clade.
    """
    out = []
    for c in clades:
        if founding(c):
            out.append(c)
            continue
        if n_defining(c) <= 0:
            continue
        out.append(c)
    return out


def cluster_lin_into_clades_with_meta(
    phylogenetic_data_df: pd.DataFrame,
    shared_mut_threshold: int = 5,
    support: dict | None = None,
    min_support: float = 0.0,
):
    """
    Cluster lineages into nested clades and, at the moment a new clade
    is created, record the clade-defining mutations (the mutations accumulated since the previous
    clade root up to the current clade root) together with where they emerged.

    Required columns:
      - 'Lineage_name'
      - 'Lineage_parent'   
      - 'Genome'
      - 'Host_type'
      - 'Time_emergence'   

    Rule B: with `support` and `min_support > 0`, a lineage founds a clade only
    if it has ALSO got at least `min_support` support -- normally the onward
    infections of it and its descendants, from `subtree_transmission_counts`.
    Leaving them out (the default) applies the divergence test alone, which is
    the superset: every rule-B clade is in it, so a caller that wants to vary
    m_T without re-clustering filters that output instead.

    Returns
    -------
    clade_to_lineages : dict
        { clade_name → set(lineage_names) }
    lineage_to_clade  : dict
        { lineage_name → clade_name }
    per_clade_mut_df  : dict
        { clade_name → DataFrame with columns:
            ['mutation', 'emergence_lineage', 'host_type', 'time'] }
        These rows are *only* the clade-defining mutations for that clade.
    clade_meta_df     : pd.DataFrame
        Columns: ['clade','root_lineage','parent_clade','n_defining','start_time']
    """
    required = {
        "Lineage_name", "Lineage_parent", "Genome", "Host_type", "Time_emergence"
    }
    missing = required - set(phylogenetic_data_df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Maps
    lin2mut    = build_lineage_to_mutation_dict(phylogenetic_data_df)
    lin2parent = dict(zip(phylogenetic_data_df["Lineage_name"], phylogenetic_data_df["Lineage_parent"]))
    lin2host   = dict(zip(phylogenetic_data_df["Lineage_name"], phylogenetic_data_df["Host_type"]))
    lin2time   = dict(zip(phylogenetic_data_df["Lineage_name"], phylogenetic_data_df["Time_emergence"]))

    # Root detection
    def _is_root_candidate(lin, par):
        return pd.isna(par) or par in ("", None) or lin == par

    root_candidates = [lin for lin, par in lin2parent.items() if _is_root_candidate(lin, par)]
    if not root_candidates:
        raise ValueError("No root lineage found.")
    root = min(root_candidates, key=lambda lin: (lin2time.get(lin, float("inf")), str(lin)))

    # Children adjacency
    children = defaultdict(list)
    for child, parent in lin2parent.items():
        if _is_root_candidate(child, parent):
            continue
        children[parent].append(child)
    for p in children:
        children[p].sort(key=lambda c: (lin2time.get(c, float("inf")), str(c)))

    # Precompute edge mutations for each lineage
    edge_muts = {}
    for lin, muts in lin2mut.items():
        parent = lin2parent.get(lin)
        parent_muts = lin2mut.get(parent, frozenset()) if parent in lin2mut else frozenset()
        edge_muts[lin] = muts - parent_muts  # typically size 1 in fully resolved trees
        
    # Containers
    clade_to_lineages: dict[str, set] = {}
    lineage_to_clade: dict[str, str] = {}
    clade_branch_counts = defaultdict(int)
    per_clade_mut_df: dict[str, pd.DataFrame] = {}
    clade_meta_rows = []

    def _norm_host(h):
        if pd.isna(h):
            return None
        v = str(h).strip().lower()
        if v.startswith("long"):
            return "long_shedder"
        if v.startswith("norm"):
            return "standard"
        if v in ("", "none", "nan", "unknown"):
            return None
        return v

    def recurse(current_lin: str, clade_root_mutset: frozenset, clade_path: str, path_since_root: list[str]):
        
        current_mut = lin2mut.get(current_lin, frozenset())
        start_new_clade = (clade_path == "") or founds_clade(
            current_lin, clade_root_mutset, lin2mut, shared_mut_threshold,
            support, min_support)

        if start_new_clade:
            parent_clade = clade_path if clade_path != "" else None

            # Build the segment from the previous clade root to the new clade root:
            # nodes on the path since last root (path_since_root) plus current node.
            # path_since_root[0] is the PARENT clade's root. Its edge mutations are
            # already inside clade_root_mutset, so the divergence test never counted
            # them, and they already define the parent clade. Including them here
            # let one mutation define two clades at once; dropped (v2.5.3).
            segment_nodes = (path_since_root[1:] + [current_lin]) if clade_path != "" else [current_lin]

            # Create the new clade name
            if clade_path == "":
                new_clade = "Clade_0"
            else:
                idx = clade_branch_counts[clade_path]
                clade_branch_counts[clade_path] += 1
                new_clade = f"{clade_path}.{idx}"

            # Collect clade-defining mutation records from segment nodes
            records = []
            defining_set = set()
            for node in segment_nodes:
                node_time = lin2time.get(node)
                node_host = _norm_host(lin2host.get(node))
                for m in edge_muts.get(node, ()):
                    defining_set.add(m)
                    records.append({
                        "mutation": m,
                        "emergence_lineage": node,
                        "host_type": node_host if node_host is not None else "unknown",
                        "time": node_time
                    })

            # clade objects
            clade_to_lineages.setdefault(new_clade, set())
            per_clade_mut_df[new_clade] = pd.DataFrame(records, columns=["mutation", "emergence_lineage", "host_type", "time"])
            clade_meta_rows.append({
                "clade": new_clade,
                "root_lineage": current_lin,
                "parent_clade": parent_clade,
                "n_defining": len(defining_set),
                "start_time": lin2time.get(current_lin)
            })

            # reset clade root state
            clade_root_mutset = current_mut
            clade_path = new_clade
            path_since_root = []  # reset path accumulator after starting a new clade

        # Assign current lineage to the current clade
        clade_to_lineages[clade_path].add(current_lin)
        lineage_to_clade[current_lin] = clade_path

        # Descend: accumulate path since current clade root
        next_path_since_root = path_since_root + [current_lin]
        for child in children.get(current_lin, []):
            recurse(child, clade_root_mutset, clade_path, next_path_since_root)

    recurse(root, lin2mut.get(root, frozenset()), "", [])

    clade_meta_df = pd.DataFrame(clade_meta_rows, columns=["clade","root_lineage","parent_clade","n_defining","start_time"])
    return clade_to_lineages, lineage_to_clade, per_clade_mut_df, clade_meta_df


def label_clades_from_definers(per_clade_mut_df: dict[str, pd.DataFrame],
                               clade_meta_df: pd.DataFrame):
    """
    Assign a label to each clade from its defining mutations.

    -----
    - If the per-clade defining-mutation table is empty:
        * If parent_clade is None => label = 'founder' (root clade).

    - Otherwise, every defining mutation must have host_type in {'standard','long_shedder'}.
      If any definer has missing/unknown host_type => raise ValueError.
      
    - Majority vote between 'standard' and 'long_shedder' determines the label.
      If a tie somehow occurs, label 'mixed'.

    Returns
    -------
    clade_labels : pd.Series mapping clade -> {'standard','long_shedder','founder','mixed'}
    clade_summary_df : DataFrame with ['clade','total_standard','total_long_shedder','total_unknown','label']
    """
    labels = {}
    rows = []

    for clade, df in per_clade_mut_df.items():
        # Case 1: no defining mutations
        if df is None or df.empty:
            # Check in clade_meta_df whether this is root
            parent = None
            if clade in set(clade_meta_df["clade"]):
                parent = clade_meta_df.loc[clade_meta_df["clade"] == clade, "parent_clade"].iloc[0]
            if pd.isna(parent) or parent is None:
                labels[clade] = "founder"
                rows.append({
                    "clade": clade,
                    "total_standard": 0,
                    "total_long_shedder": 0,
                    "total_unknown": 0,
                    "label": "founder",
                })
                continue
            else:
                raise ValueError(
                    f"[label_clades_from_definers] Clade {clade} has no definers but is not root "
                    f"(parent={parent})."
                )

        # Case 2: definers present
        tmp = df.copy()
        tmp["host_type"] = tmp["host_type"].astype(str).str.strip().str.lower()
        tmp["host_type"] = tmp["host_type"].replace({
            "long": "long_shedder",
            "long-shedder": "long_shedder",
            "longshedder": "long_shedder",
            "norm": "standard",
            "none": "unknown",
            "nan": "unknown",
            "": "unknown",
        })

        counts = tmp["host_type"].value_counts(dropna=False).to_dict()
        total_standard = int(counts.get("standard", 0))
        total_long   = int(counts.get("long_shedder", 0))
        total_unknown = int(counts.get("unknown", 0))

        if total_unknown > 0:
            bad = tmp[tmp["host_type"] == "unknown"].head(5)
            raise ValueError(
                f"[label_clades_from_definers] Found {total_unknown} 'unknown' definers in clade '{clade}'. "
                f"Expected only 'standard'/'long_shedder'. Examples:\n{bad}"
            )

        if total_long > total_standard:
            label = "long_shedder"
        elif total_standard > total_long:
            label = "standard"
        else:
            # should not happen if threshold is odd
            label = "mixed"

        labels[clade] = label
        rows.append({
            "clade": clade,
            "total_standard": total_standard,
            "total_long_shedder": total_long,
            "total_unknown": total_unknown,
            "label": label
        })

    clade_labels = pd.Series(labels)
    clade_summary_df = pd.DataFrame(
        rows, columns=["clade","total_standard","total_long_shedder","total_unknown","label"]
    )
    return clade_labels, clade_summary_df

