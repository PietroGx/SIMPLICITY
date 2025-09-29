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

def get_clustering_for_ssod(
    experiment_name: str,
    seeded_simulation_output_dir: str,
    phylogenetic_df: pd.DataFrame | None = None,
    shared_mut_threshold: int = 5,
    force: bool = False,
) -> str:
    """
    If the clustering CSV exists (and not force), returns it.
    Otherwise computes clustering and writes a single CSV.
    """
    path = om.get_clustering_table_filepath(experiment_name, seeded_simulation_output_dir)
    if (not force) and os.path.exists(path):
        return path

    if phylogenetic_df is None:
        phylogenetic_df = om.read_phylogenetic_data(seeded_simulation_output_dir)

    clade_to_lineages, lineage_to_clade, per_clade_mut_df, clade_meta_df = \
        cluster_lin_into_clades_with_meta(phylogenetic_df, shared_mut_threshold=shared_mut_threshold)
        
    # Set 'clade' as the index 
    if isinstance(clade_meta_df, pd.DataFrame) and "clade" in clade_meta_df.columns:
        clade_meta_df = clade_meta_df.set_index("clade")

    try:
        labels_s, _ = label_clades_from_definers(per_clade_mut_df, clade_meta_df)
    except Exception:
        labels_s = pd.Series({}, dtype=object)

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
            "clade": clade,
            "root_lineage": root_lineage,
            "parent_clade": parent_clade,
            "start_time": start_time,
            "label": labels_s.get(clade, ""),
            "n_lineages": len(lineages),
            "n_defining": len(definers_list),
            "n_def_normal": sum(1 for d in definers_list if str(d.get("host_type","")).lower().startswith("normal")),
            "n_def_long":   sum(1 for d in definers_list if "long" in str(d.get("host_type","")).lower()),
            "lineages": repr(lineages),
            "defining_mutations": repr(definers_list),
        })

    df = pd.DataFrame(rows).sort_values(["start_time", "clade"]).reset_index(drop=True)
    om.write_clustering_table(experiment_name, seeded_simulation_output_dir, df)
    return path

def parse_genome(genome):
    """Convert a genome (list of mutations) into a frozenset of mutation tuples."""
    if genome is None:
        return frozenset()
    try:
        return frozenset(tuple(x) for x in genome)
    except Exception:
        return frozenset()
    
def build_lineage_to_mutation_dict(phylogenetic_data_df):
    """Return {Lineage_name -> frozenset of mutation tuples} from the phylo table."""
    return {
        row["Lineage_name"]: parse_genome(row["Genome"])
        for _, row in phylogenetic_data_df.iterrows()
    }

def cluster_lin_into_clades_with_meta(
    phylogenetic_data_df: pd.DataFrame,
    shared_mut_threshold: int = 5,
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
            return "normal"
        if v in ("", "none", "nan", "unknown"):
            return None
        return v

    def recurse(current_lin: str, clade_root_mutset: frozenset, clade_path: str, path_since_root: list[str]):
        
        current_mut = lin2mut.get(current_lin, frozenset())
        mut_distance = len(current_mut - clade_root_mutset)
        start_new_clade = (mut_distance >= shared_mut_threshold) or (clade_path == "")

        if start_new_clade:
            parent_clade = clade_path if clade_path != "" else None

            # Build the segment from the previous clade root to the new clade root:
            # nodes on the path since last root (path_since_root) plus current node.
            segment_nodes = (path_since_root + [current_lin]) if clade_path != "" else [current_lin]

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

    - Otherwise, every defining mutation must have host_type in {'normal','long_shedder'}.
      If any definer has missing/unknown host_type => raise ValueError.
      
    - Majority vote between 'normal' and 'long_shedder' determines the label.
      If a tie somehow occurs, label 'mixed'.

    Returns
    -------
    clade_labels : pd.Series mapping clade -> {'normal','long_shedder','founder','mixed'}
    clade_summary_df : DataFrame with ['clade','total_normal','total_long_shedder','total_unknown','label']
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
                    "total_normal": 0,
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
            "norm": "normal",
            "none": "unknown",
            "nan": "unknown",
            "": "unknown",
        })

        counts = tmp["host_type"].value_counts(dropna=False).to_dict()
        total_normal = int(counts.get("normal", 0))
        total_long   = int(counts.get("long_shedder", 0))
        total_unknown = int(counts.get("unknown", 0))

        if total_unknown > 0:
            bad = tmp[tmp["host_type"] == "unknown"].head(5)
            raise ValueError(
                f"[label_clades_from_definers] Found {total_unknown} 'unknown' definers in clade '{clade}'. "
                f"Expected only 'normal'/'long_shedder'. Examples:\n{bad}"
            )

        if total_long > total_normal:
            label = "long_shedder"
        elif total_normal > total_long:
            label = "normal"
        else:
            # should not happen if threshold is odd
            label = "mixed"

        labels[clade] = label
        rows.append({
            "clade": clade,
            "total_normal": total_normal,
            "total_long_shedder": total_long,
            "total_unknown": total_unknown,
            "label": label
        })

    clade_labels = pd.Series(labels)
    clade_summary_df = pd.DataFrame(
        rows, columns=["clade","total_normal","total_long_shedder","total_unknown","label"]
    )
    return clade_labels, clade_summary_df

