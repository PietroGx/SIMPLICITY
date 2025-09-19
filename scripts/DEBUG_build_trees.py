# Build the simulated *phylogenetic* tree for selected seeded simulations
# and export it in the requested formats (json / newick / nexus).

import argparse
import re
from tqdm import tqdm

import simplicity.tree.tree_builder as tb
import simplicity.output_manager as om
import simplicity.dir_manager as dm

# -----------------------------
# helpers
# -----------------------------

def select_seeded_simulations(experiment_name: str,
                              time_threshold: float,
                              lineage_number_threshold: int):
    """
    Filters seeded simulation output dirs by:
      - final time > time_threshold
      - number of lineages > lineage_number_threshold
    and keeps only seeds common to all simulation runs (like in 03_plot_trees.py).
    """
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    all_ssods = [dm.get_seeded_simulation_output_dirs(d) for d in simulation_output_dirs]

    filtered_dirs = []
    seed_pattern = re.compile(r"seed_(\d+)")
    seeds_per_dir = {}

    for ssod_list in tqdm(all_ssods, desc="Filtering seeded simulations"):
        for ssod in ssod_list:
            final_time = om.read_final_time(ssod)
            if final_time > time_threshold:
                phylogenetic_data = om.read_phylogenetic_data(ssod)
                n_lineages = len(phylogenetic_data["Lineage_name"])
                if n_lineages > lineage_number_threshold:
                    filtered_dirs.append(ssod)
                    m = seed_pattern.search(ssod)
                    if m:
                        seed = m.group(1)
                        seeds_per_dir.setdefault(seed, []).append(ssod)

    common_seeds = [
        seed for seed, dirs in seeds_per_dir.items()
        if len(dirs) == len(simulation_output_dirs)
    ]
    final_dirs = [
        d for seed, dirs in seeds_per_dir.items()
        if seed in common_seeds for d in dirs
    ]
    kept_seeds = sorted(common_seeds)
    return final_dirs, kept_seeds


def ensure_phylo_tree(experiment_name: str,
                      ssod: str,
                      filetype: str,
                      overwrite: bool = False):
    """
    Build/export the phylogenetic tree in the requested filetype if missing
    (or if overwrite=True). Uses the existing tree builder.
    """
    path = om.get_tree_file_filepath(
        experiment_name,
        ssod,
        tree_type="phylogenetic",
        tree_subtype="binary",
        file_type=filetype
    )

    if (not overwrite) and om.os.path.exists(path):
        print(f"✔︎ Exists: {path}")
        return path

    print(f"⧗ Building/exporting phylogenetic tree → {filetype} …")
    tb.get_tree(
        experiment_name,
        ssod,
        tree_type="phylogenetic",
        tree_subtype="binary",
        coloring="lineage",     # coloring has no effect on file export
        save_plot=False,        # ensure no plots are generated
        export_filetype=filetype,
        dashplot=False
    )
    print(f"✓ Wrote: {path}")
    return path


# -----------------------------
# main
# -----------------------------

def main():
    p = argparse.ArgumentParser(
        description="Build phylogenetic trees and export to JSON/Newick/Nexus (no plots)."
    )
    p.add_argument("experiment_name", type=str, help="Experiment name")
    p.add_argument("--time-threshold", type=float, default=365.0,
                   help="Min final time to keep a seeded simulation (default: 365)")
    p.add_argument("--lineage-threshold", type=int, default=30,
                   help="Min lineage count to keep a seeded simulation (default: 30)")
    p.add_argument("--formats", type=str, default="json,newick,nexus",
                   help="Comma-separated list of formats to export (json,newick,nexus)")
    p.add_argument("--overwrite", action="store_true",
                   help="Rebuild even if the target file already exists")

    args = p.parse_args()

    fmts = [f.strip().lower() for f in args.formats.split(",") if f.strip()]
    allowed = {"json", "newick", "nexus"}
    for f in fmts:
        if f not in allowed:
            raise ValueError(f"Unsupported format '{f}'. Use any of: {sorted(allowed)}")

    selected_ssods, kept_seeds = select_seeded_simulations(
        args.experiment_name,
        args.time_threshold,
        args.lineage_threshold
    )

    if not selected_ssods:
        print("No seeded simulations passed the filters. Nothing to do.")
        return

    print(f"Seeds kept: {kept_seeds}")
    print(f"Formats: {fmts}")

    for ssod in tqdm(selected_ssods, desc="Exporting phylogenetic trees"):
        for f in fmts:
            ensure_phylo_tree(
                experiment_name=args.experiment_name,
                ssod=ssod,
                filetype=f,
                overwrite=args.overwrite
            )

    print("Done. Phylogenetic trees exported in requested formats.")

if __name__ == "__main__":
    main()

