#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare interactive Nextstrain (Auspice v2) datasets from SIMPLICITY outputs.

For each SSOD:
  - Get clustering table (one CSV per SSOD).
  - Get phylogenetic tree JSON.
  - Build Auspice v2 dataset with:
      * node_attrs.num_date (time axis from Time_emergence)
      * node_attrs.cluster  (tip clade/cluster)
      * meta.colorings      (embedded per-clade colors to match Figure 3 palette)
  - Write into 06_Trees/nextstrain/ as:
      SOD_seed_SEEDNR_nextstrain_dataset.json
      SOD_seed_SEEDNR_metadata.tsv
"""
from __future__ import annotations
import argparse
import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Any, List, Tuple
from shutil import which as _which
import subprocess
import pandas as pd

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tree.tree_builder as tb
import simplicity.clustering as clustering
import simplicity.plots_manager as pm

# ---------------- helpers ----------------

def to_decimal_year(dt: datetime) -> float:
    y0 = datetime(dt.year, 1, 1); y1 = datetime(dt.year + 1, 1, 1)
    frac = (dt - y0).total_seconds() / (y1 - y0).total_seconds()
    return dt.year + frac

def decimal_year_to_date_iso(decimal_year: float) -> str:
    year = int(decimal_year)
    y0 = datetime(year, 1, 1); y1 = datetime(year + 1, 1, 1)
    dt = y0 + timedelta(seconds=(decimal_year - year) * (y1 - y0).total_seconds())
    return dt.date().isoformat()

def get_json_tree(ssod: str) -> Path:
    """
    Return path to the phylogenetic JSON for this SSOD; build it if missing.
    """
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)

    tree_path = Path(om.get_tree_file_filepath(
        experiment_name,
        ssod,
        tree_type="phylogenetic",
        tree_subtype="binary",
        file_type="json"
    ))

    if tree_path.exists():
        return tree_path

    tb.get_tree(
        experiment_name=experiment_name,
        seeded_simulation_output_dir=ssod,
        tree_type="phylogenetic",
        tree_subtype="binary",
        coloring="lineage",
        save_plot=False,
        export_filetype="json",
        dashplot=False
    )

    return tree_path

def get_clustering_table(ssod: str, shared_mut_threshold: int) -> Path:
    """
    Return path to the single clustering CSV for this SSOD; compute if missing.
    """        
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)

    csv_path = Path(om.get_clustering_table_filepath(experiment_name, ssod))
    if not csv_path.exists():
        _ = clustering.get_clustering_for_ssod(
            experiment_name=experiment_name,
            seeded_simulation_output_dir=ssod,
            shared_mut_threshold=shared_mut_threshold,
            force=False
        )
        
    return csv_path

def load_lin2clade_and_roots(clustering_csv: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    From single-CSV clustering: build (lineage -> clade, clade -> root_lineage).
    """
    import ast
    df = pd.read_csv(clustering_csv, converters={"lineages": ast.literal_eval})
    lin2clade: Dict[str, str] = {}
    clade2root: Dict[str, str] = {}
    
    WT_LINEAGE_NAME = "wt"

    for _, row in df.iterrows():
        clade = str(row["clade"])
        root_lineage = row.get("root_lineage")

        if pd.isna(root_lineage):
            if clade == "Clade_0":
                clade2root[clade] = WT_LINEAGE_NAME
            else:
                raise ValueError(
                    f"\n\nFATAL: Subclade '{clade}' is missing its 'root_lineage' in the clustering file.\n"
                )
        else:
            clade2root[clade] = str(root_lineage)

        for lin in row["lineages"]:
            lin2clade[str(lin)] = clade
            
    return lin2clade, clade2root

def build_colorings_block(clade2root: Dict[str, str], ssod: str) -> Dict[str, Any]:
    """
    Create Auspice meta.colorings entry for 'cluster'
    """
    cmap_df = pm.make_lineages_colormap(ssod)
    scale = []
    # Sort clades by hierarchy for stable order (root first, then descendants)
    for clade in sorted(clade2root.keys(), key=lambda c: (c.count("."), c)):
        root_lin = clade2root.get(clade, "")
        hexcol = pm.get_lineage_color(root_lin, cmap_df)  
        scale.append([clade, hexcol])
    return {
        "key": "cluster",
        "type": "categorical",
        "title": "cluster",
        "scale": scale
    }

def json_to_auspice(root_obj: Dict[str, Any],
                    start_date: datetime,
                    time_units: str,
                    lineage_to_clade: Dict[str, str]) -> Dict[str, Any]:
    """
    Convert SIMPLICITY tree JSON → Auspice v2 'tree' object.
    """
    def rec(n: Dict[str, Any]) -> Dict[str, Any]:
        name_in = n.get("name") or n.get("label")
        t = float(n.get("time_emergence", 0.0) or 0.0)
        dt = start_date + (timedelta(days=t) if time_units.lower().startswith("day")
                           else timedelta(days=t * 365.2425))
        node = {
            "name": name_in,
            "node_attrs": {"num_date": {"value": to_decimal_year(dt)}}
        }
        kids = n.get("children") or []

        if kids:
            node["children"] = [rec(k) for k in kids]

            # infer the parent's state from its children.
            child_clusters = set()
            for child_node in node["children"]:
                cluster_attr = child_node.get("node_attrs", {}).get("cluster")
                if cluster_attr:
                    child_clusters.add(cluster_attr["value"])

            # If all children have the same, single cluster, assign it to the parent.
            if len(child_clusters) == 1:
                parent_cluster = child_clusters.pop()
                node["node_attrs"]["cluster"] = {"value": parent_cluster}

        else:
            # This is a tip: add cluster from the lookup table.
            cl = lineage_to_clade.get(str(name_in))
            if cl:
                node["node_attrs"]["cluster"] = {"value": str(cl)}

        return node

    return rec(root_obj)

def process_ssod(ssod: str,
                 start_date: datetime,
                 time_units: str,
                 shared_mut_threshold: int) -> Tuple[Path, Path]:
    """
    Build Auspice dataset + metadata for one SSOD.
    """
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)

    # Get inputs
    clustering_csv = get_clustering_table(ssod, shared_mut_threshold)
    lin2clade, clade2root = load_lin2clade_and_roots(clustering_csv)
    json_tree_path = get_json_tree(ssod)

    # Load tree JSON 
    root = json.loads(Path(json_tree_path).read_text())
    if "children" not in root: 
        root = root.get("tree") or root

    # Convert to Auspice v2 tree
    aus_tree = json_to_auspice(root, start_date, time_units, lin2clade)

    # Dataset + metadata paths 
    dataset_path, meta_path, base = om.get_nextstrain_dataset_paths(experiment_name, ssod)

    # Embed color palette for this dataset
    colorings = build_colorings_block(clade2root, ssod)

    # Write dataset JSON
    dataset = {"version": 2, "meta": {"title": base, "colorings": [colorings]}, "tree": aus_tree}
    Path(dataset_path).write_text(json.dumps(dataset, separators=(",", ":"), ensure_ascii=False))

    # Write minimal metadata (strain, date, cluster)
    rows: List[Tuple[str, str, str]] = []
    def collect_tips(n: Dict[str, Any]):
        kids = n.get("children")
        if not kids:
            nm = n["name"]
            nd = float(n["node_attrs"]["num_date"]["value"])
            date_iso = decimal_year_to_date_iso(nd)
            cluster = n["node_attrs"].get("cluster", {}).get("value", "")
            rows.append((nm, date_iso, cluster)); return
        for c in kids:
            collect_tips(c)
    collect_tips(aus_tree)
    with open(meta_path, "w") as f:
        f.write("strain\tdate\tcluster\n")
        for nm, date_iso, cl in rows:
            f.write(f"{nm}\t{date_iso}\t{cl}\n")

    print(f"{Path(dataset_path).name}  (Auspice dataset w/ embedded colors)")
    print(f"{Path(meta_path).name}     (metadata)")
    return Path(dataset_path), Path(meta_path)

def main():

    ap = argparse.ArgumentParser(
        description="Prepare Nextstrain (Auspice) datasets from SIMPLICITY outputs."
    )
    mx = ap.add_mutually_exclusive_group(required=True)
    mx.add_argument("--experiment", type=str, help="Process ALL SSODs in this experiment")
    mx.add_argument("--ssod", type=str, help="Process a single seeded_simulation_output_dir")
    ap.add_argument("--start-date", default="2020-01-01", type=str,
                    help="Simulation start date (YYYY-MM-DD)")
    ap.add_argument("--time-units", default="days", choices=["days", "years"],
                    help="Units of node['time_emergence']")
    ap.add_argument("--shared-mut-threshold", default=5, type=int,
                    help="Clustering shared-mutation threshold")

    args = ap.parse_args()
    start_date = datetime.fromisoformat(args.start_date)

    if args.experiment:
        # process all SSODs under the experiment (no GUI)
        ssods = []
        for sod in dm.get_simulation_output_dirs(args.experiment):
            ssod_list = dm.get_seeded_simulation_output_dirs(sod)
            ssods.extend(ssod_list)
        if not ssods:
            # print("No seeded simulations found for this experiment.")
            return
        # print(ssods)
        for ssod in ssods:
            process_ssod(ssod, start_date, args.time_units, args.shared_mut_threshold)
        return

    # else: single SSOD process and launch Auspice
    ssod = args.ssod
    process_ssod(ssod, start_date, args.time_units, args.shared_mut_threshold)

    # determine dataset directory and launch Auspice
    exp_name = dm.get_experiment_foldername_from_SSOD(ssod)
    dataset_dir = dm.get_nextstrain_dir(exp_name)

    if _which("auspice"):
        cmd = ["auspice", "view", "--datasetDir", str(dataset_dir)]
    elif _which("npx"):
        cmd = ["npx", "auspice", "view", "--datasetDir", str(dataset_dir)]
    else:
        print("Auspice not found. Install with: npm i -g auspice")
        return

    print(f"Launching Auspice... (datasetDir={dataset_dir})  (Ctrl+C to exit)")
    subprocess.run(cmd, check=False)

if __name__ == "__main__":
    main()
