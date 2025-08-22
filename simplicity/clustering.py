#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from collections import defaultdict
import pandas as pd

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
    Cluster lineages into nested clades (time-ordered, deterministic) and, at the moment a new clade
    is created, record the *clade-defining mutations* (the mutations accumulated since the previous
    clade root up to the current clade root) together with where they emerged.

    Required columns:
      - 'Lineage_name'     : unique lineage ID
      - 'Lineage_parent'   : parent lineage ID
      - 'Genome'           : list-like of mutations (already parsed upstream)
      - 'Host_type'        : host phenotype where the lineage/edge arose
      - 'Time_emergence'   : numeric time (used for deterministic ordering)

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
    
    
    
    # # --- DEBUG: input host values & genome parsing ---
    # raw_hosts = phylogenetic_data_df["Host_type"].astype(str).str.strip().str.lower()
    # print("[DBG] Host_type unique (raw):", sorted(set(raw_hosts)))
    
    # nonempty_edges = sum(1 for v in lin2mut.values())
    # print(f"[DBG] Total lineages with genomes: {nonempty_edges}/{len(lin2mut)}")
    
    # # Check edge_muts fill after you build it (see next block)

    
    

    # Root detection (earliest Time_emergence among NaN/empty/self-parent)
    def _is_root_candidate(lin, par):
        return pd.isna(par) or par in ("", None) or lin == par

    root_candidates = [lin for lin, par in lin2parent.items() if _is_root_candidate(lin, par)]
    if not root_candidates:
        raise ValueError("No root lineage found.")
    root = min(root_candidates, key=lambda lin: (lin2time.get(lin, float("inf")), str(lin)))

    # Children adjacency, deterministic order: by time then name
    children = defaultdict(list)
    for child, parent in lin2parent.items():
        if _is_root_candidate(child, parent):
            continue
        children[parent].append(child)
    for p in children:
        children[p].sort(key=lambda c: (lin2time.get(c, float("inf")), str(c)))

    # Precompute edge mutations for each lineage: muts(lineage) - muts(parent)
    edge_muts = {}
    for lin, muts in lin2mut.items():
        parent = lin2parent.get(lin)
        parent_muts = lin2mut.get(parent, frozenset()) if parent in lin2mut else frozenset()
        edge_muts[lin] = muts - parent_muts  # typically size 1 in fully resolved trees
        
    
    
    # # --- DEBUG: edge mutation sanity ---
    # non_empty = sum(1 for v in edge_muts.values() if len(v) > 0)
    # print(f"[DBG] edge_muts non-empty: {non_empty}/{len(edge_muts)} ({non_empty/len(edge_muts):.1%})")
    # if non_empty == 0:
    #     print("[DBG] WARNING: No edge mutations found; definers will be empty -> labels become 'unknown'/'mixed'.")

    
    

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
        # NEW mutations since current clade root
        mut_distance = len(current_mut - clade_root_mutset)

        start_new_clade = (mut_distance >= shared_mut_threshold) or (clade_path == "")

        if start_new_clade:
            parent_clade = clade_path if clade_path != "" else None

            # Build the segment from the previous clade root to the *new* clade root:
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
            
            
            # # --- DEBUG: clade-defining segment summary ---
            # host_counts = pd.Series([r["host_type"] for r in records]).value_counts(dropna=False).to_dict()
            # print(
            #     f"[DBG] New {new_clade}: "
            #     f"parent={parent_clade}, segment_nodes={len(segment_nodes)}, "
            #     f"n_defining={len(defining_set)}, start_time={lin2time.get(current_lin)}"
            # )
            # print(f"[DBG]   definers host_counts: {host_counts}")
            # if host_counts.get("unknown", 0) > 0:
            #     print("[DBG]   first few defining rows:", records[:min(3, len(records))])

            
            
            
            
            # Register clade objects
            clade_to_lineages.setdefault(new_clade, set())
            per_clade_mut_df[new_clade] = pd.DataFrame(records, columns=["mutation", "emergence_lineage", "host_type", "time"])
            clade_meta_rows.append({
                "clade": new_clade,
                "root_lineage": current_lin,
                "parent_clade": parent_clade,
                "n_defining": len(defining_set),
                "start_time": lin2time.get(current_lin)
            })

            # Reset clade root state
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


def label_clades_from_definers(per_clade_mut_df: dict[str, pd.DataFrame]):
    """
    Given per-clade mutation tables (only clade-defining mutations), assign each clade a label
    by majority of host types among those mutations and also return a summary table.

    Returns
    -------
    clade_labels : pd.Series
        { clade_name -> 'normal' / 'long_shedder' / 'unknown' }

    clade_summary_df : pd.DataFrame
        columns = ['clade','total_normal','total_long_shedder','total_unknown','label']
    """
    labels = {}
    rows = []

    for clade, df in per_clade_mut_df.items():
        if df is None or df.empty:
            labels[clade] = "unknown"
            rows.append({"clade": clade, "total_normal": 0, "total_long_shedder": 0, "total_unknown": 0, "label": "unknown"})
            continue

        counts = df["host_type"].value_counts(dropna=False).to_dict()
        total_normal = counts.get("normal", 0)
        total_long   = counts.get("long_shedder", 0)
        total_unknown = counts.get("unknown", 0)
        
        # # --- DEBUG: per-clade label inputs ---
        # print(f"[DBG] Labeling {clade}: normal={total_normal}, long={total_long}, unknown={total_unknown}, n_rows={len(df)}")
        # if total_normal == 0 and total_long == 0 and total_unknown > 0:
        #     print("[DBG]   NOTE: all definers are 'unknown' -> will not favor normal/long.")


        if total_long > total_normal:
            label = "long_shedder"
        elif total_normal > total_long:
            label = "normal"
        else:
            # tie or both zero, assign mixed label
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
    clade_summary_df = pd.DataFrame(rows, columns=["clade","total_normal","total_long_shedder","total_unknown","label"])
    return clade_labels, clade_summary_df




