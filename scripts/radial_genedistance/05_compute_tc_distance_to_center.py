"""
05_compute_tc_distance_to_center.py
=====================================
Script 1 of 3 for nuclear distance analysis – Distance to nuclear centre for
translocated conditions (T1 and C1).

What this script does
---------------------
For each translocated condition (T1 and C1), this script asks: after a
chromosomal translocation, do the translocated genes move closer to or
further from the centre of the nucleus?

In Chrom3D polymer models, the nucleus is represented as a sphere centred
on the origin (0, 0, 0). The Euclidean distance from a gene's bead position
to the origin is therefore a measure of its radial position in the nucleus —
genes near the centre tend to be in transcriptionally active compartments,
while genes near the periphery tend to be in inactive compartments.

For each .cmm structural model this script:
1. Reads the gene expression TSV produced by script 01 (gene_expression
   extraction), keeping only translocated and control genes.
2. Maps each gene to its nearest genomic bead in the model.
3. Computes the Euclidean distance from that bead to the origin (0, 0, 0).
4. Saves per-model distances and an aggregated summary (mean, std, n_models).

Why only inside_transloc and control?
--------------------------------------
neighbor genes are excluded because the radial repositioning hypothesis
is specifically about whether the translocated genes themselves move toward
or away from the nuclear centre after the translocation — not about what
happens to the genes they land near.

Output of this script is used by script 07.

Usage
-----
    1. Run script 01 (01_extract_gene_expression.py) first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 05_compute_tc_distance_to_center.py

Dependencies
------------
    pandas, numpy
    Install with: pip install pandas numpy

"""

import os
import re
import glob

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    "conditions": {
        "T1": {
            # Gene expression TSV from script 01
            "gene_expression_file": "/path/to/T1_genes_expression.tsv",
            # Directory of Chrom3D .cmm structural model files for T1
            "cmm_dir":              "/path/to/cmm/T1",
            # Where to save the per-model distances (one row per gene per model)
            "output":               "/path/to/T1_distance_to_center.tsv",
            # Where to save the aggregated summary (mean/std across all models)
            "output_agg":           "/path/to/T1_distance_to_center_agg.tsv",
        },
        "C1": {
            "gene_expression_file": "/path/to/C1_genes_expression.tsv",
            "cmm_dir":              "/path/to/cmm/C1",
            "output":               "/path/to/C1_distance_to_center.tsv",
            "output_agg":           "/path/to/C1_distance_to_center_agg.tsv",
        },
    },
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def parse_cmm(file_path: str) -> pd.DataFrame:
    """
    Parse a UCSF Chimera / Chrom3D .cmm marker file and return a DataFrame
    with columns: chr, start, end, x, y, z.

    Each <marker> line with a beadID of the form 'chrN:start-end' represents
    one genomic bead positioned at (x, y, z) in 3D space. <link> elements
    and all other attributes are ignored.
    """
    markers = []
    with open(file_path) as fh:
        for line in fh:
            if "<marker" not in line:
                continue
            attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
            if "beadID" not in attrs:
                continue
            try:
                chrom_part, coords = attrs["beadID"].split(":")
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part
                start, end = map(int, coords.split("-"))
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, start, end, x, y, z])
            except Exception:
                continue  # skip malformed lines silently

    return pd.DataFrame(markers, columns=["chr", "start", "end", "x", "y", "z"])


def build_bead_index(beads_by_chr: dict) -> dict:
    """
    Pre-convert each chromosome's bead DataFrame into numpy arrays for fast
    nearest-bead lookup.

    For each gene we need to find the closest bead by genomic midpoint.
    Pre-extracting midpoints and coordinates as numpy arrays (once per model)
    avoids slow and memory-heavy DataFrame filtering inside the gene loop.

    Returns a dict: chromosome → {'mids': array, 'coords': array}
    """
    index = {}
    for chrom, bdf in beads_by_chr.items():
        index[chrom] = {
            "mids":   ((bdf["start"] + bdf["end"]) / 2).to_numpy(),
            "coords": bdf[["x", "y", "z"]].to_numpy(),
        }
    return index


def map_genes_to_beads(genes_df: pd.DataFrame, bead_index: dict) -> np.ndarray:
    """
    For each gene, find the nearest bead on the same chromosome (by midpoint)
    and return its 3D coordinates.

    Returns
    -------
    coords : (N, 3) float array – x, y, z of the matched bead (NaN if none found)
    """
    n      = len(genes_df)
    coords = np.full((n, 3), np.nan)

    for i, (_, row) in enumerate(genes_df.iterrows()):
        entry = bead_index.get(row["chrom"])
        if entry is None:
            continue
        gene_mid  = (row["start"] + row["end"]) / 2
        nearest   = np.abs(entry["mids"] - gene_mid).argmin()
        coords[i] = entry["coords"][nearest]

    return coords


# =============================================================================
# Main loop over conditions
# =============================================================================

for cond, paths in CONFIG["conditions"].items():
    print(f"\n===== Processing {cond} =====")

    genes = pd.read_csv(paths["gene_expression_file"], sep="\t")

    # Ensure consistent chromosome notation
    genes["chrom"] = genes["chrom"].astype(str).apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )

    # Keep only translocated and control genes.
    # neighbor genes are excluded because the question here is about radial
    # repositioning of the translocated genes themselves.
    genes = genes[genes["region_type"].isin(["inside_transloc", "control"])].copy()
    print(f"Loaded {len(genes)} genes (inside_transloc + control)")

    cmm_files = sorted(glob.glob(os.path.join(paths["cmm_dir"], "*.cmm")))
    print(f"Found {len(cmm_files)} CMM files to process")

    all_results = []

    for f in cmm_files:
        print(f"  Processing model: {os.path.basename(f)}")
        beads = parse_cmm(f)
        if beads.empty:
            print("    No markers found, skipping")
            continue

        # Build numpy bead index and free the raw DataFrame
        beads_by_chr = {c: b for c, b in beads.groupby("chr")}
        bead_index   = build_bead_index(beads_by_chr)
        del beads, beads_by_chr

        # Map each gene to its nearest bead in this model
        coords = map_genes_to_beads(genes, bead_index)
        del bead_index

        # Euclidean distance from each bead to the nuclear centre (0, 0, 0).
        # In Chrom3D models the nucleus is a unit sphere centred on the origin,
        # so this distance is a direct measure of radial nuclear position.
        distances = np.linalg.norm(coords, axis=1)

        tmp = genes[["gene_id", "gene_name", "chrom",
                     "region_type", "transloc_id"]].copy()
        tmp["distance_to_center"] = distances
        tmp["model"]              = os.path.basename(f)
        all_results.append(tmp)

    if not all_results:
        print("  No results — skipping.")
        continue

    df = pd.concat(all_results, ignore_index=True)

    os.makedirs(os.path.dirname(paths["output"]), exist_ok=True)
    df.to_csv(paths["output"], sep="\t", index=False)
    print(f"  Saved per-model distances: {paths['output']}")

    # Aggregate across all models: one row per gene
    df_agg = (
        df.groupby(
            ["gene_name", "chrom", "region_type", "transloc_id"],
            as_index=False,
        )["distance_to_center"]
        .agg(mean_distance="mean", std_distance="std", n_models="count")
    )

    df_agg.to_csv(paths["output_agg"], sep="\t", index=False)
    print(f"  Saved aggregated distances: {paths['output_agg']}")
