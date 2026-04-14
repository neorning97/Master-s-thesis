"""
06_compute_wt_distance_to_center.py
=====================================
Script 2 of 3 for nuclear distance analysis – Distance to nuclear centre for
the wildtype (WT) models.

What this script does
---------------------
This script computes the same distance-to-centre measure as script 05, but
using the wildtype (WT) structural models instead of the translocated ones.

Why do we need WT distances?
----------------------------
To know whether a translocation moves genes toward or away from the nuclear
centre, we need a baseline — where were those same genes located in a normal
(wildtype) cell, before the translocation occurred?

This script maps every gene from the T1 and C1 gene expression tables to its
position in the WT models, computes its distance to the nuclear centre
(origin = 0, 0, 0), and aggregates across all WT models.

The output is then compared with the T1/C1 distances (from script 05) in
script 07 to test whether translocated genes are significantly repositioned.

Design choice
-------------
Gene coordinates in WT models come from the reference genome (GTF), not from
the translocated genome, because WT cells have no rearrangement.

Output of this script is used by script 07.

Usage
-----
    1. Run script 01 (01_extract_gene_expression.py) first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 06_compute_wt_distance_to_center.py

Dependencies
------------
    pandas, numpy
    Install with: pip install pandas numpy

"""

import os
import re
import glob
from collections import defaultdict

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    # Gene expression TSVs from script 01 — used to get the gene list
    "t1_gene_expression_file": "/path/to/T1_genes_expression.tsv",
    "c1_gene_expression_file": "/path/to/C1_genes_expression.tsv",

    # GTF annotation for reference genome coordinates
    # (WT cells are not rearranged, so we use standard gene coordinates)
    "gtf_file": "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",

    # Directory of wildtype .cmm structural model files
    "wt_cmm_dir": "/path/to/cmm/WT",

    # Where to save the aggregated WT distances
    "output_file": "/path/to/WT_distance_to_center_agg.tsv",
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

    Returns a dict: chromosome → {'mids': array, 'coords': array}
    """
    index = {}
    for chrom, bdf in beads_by_chr.items():
        index[chrom] = {
            "mids":   ((bdf["start"] + bdf["end"]) / 2).to_numpy(),
            "coords": bdf[["x", "y", "z"]].to_numpy(),
        }
    return index


# =============================================================================
# Main
# =============================================================================

print("Computing WT distance to nuclear centre...")

# ── Load gene lists from T1 and C1 ───────────────────────────────────────────
# We use the same genes as script 05 so the comparison is fair
df_t1 = pd.read_csv(CONFIG["t1_gene_expression_file"], sep="\t")
df_c1 = pd.read_csv(CONFIG["c1_gene_expression_file"], sep="\t")

# Collect unique gene names and their region types from both conditions.
# Only inside_transloc and control genes are needed (matching script 05).
gene_region: dict[str, str] = {}
for df in (df_t1, df_c1):
    for _, row in df[df["region_type"].isin(["inside_transloc", "control"])].iterrows():
        gene_region[row["gene_name"]] = row["region_type"]

genes_of_interest = list(gene_region.keys())
print(f"Total genes for WT mapping: {len(genes_of_interest)}")

# ── Load GTF to get reference genome coordinates ──────────────────────────────
# Filter the GTF to only the genes we need to keep memory usage down.
# WT cells have no rearrangement, so standard reference coordinates are correct.
gtf_cols = ["chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
gtf_df = pd.read_csv(CONFIG["gtf_file"], sep="\t", comment="#",
                     header=None, names=gtf_cols)
gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')
gtf_df = gtf_df[gtf_df["gene_name"].isin(genes_of_interest)].copy()
gtf_df["chrom"] = gtf_df["chrom"].astype(str).apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"]   = pd.to_numeric(gtf_df["end"],   errors="coerce")
print(f"GTF entries found for these genes: {len(gtf_df)}")

# ── Process each WT structural model ──────────────────────────────────────────
cmm_files = sorted(glob.glob(os.path.join(CONFIG["wt_cmm_dir"], "*.cmm")))
print(f"Found {len(cmm_files)} WT CMM files to process")

# Store one distance value per gene per model
distances: defaultdict[str, list] = defaultdict(list)

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

    # Map every gene to its nearest bead in this WT model
    for _, row in gtf_df.iterrows():
        gene  = row["gene_name"]
        chrom = row["chrom"]
        mid   = (row["start"] + row["end"]) / 2
        entry = bead_index.get(chrom)
        if entry is None:
            continue
        nearest = np.abs(entry["mids"] - mid).argmin()
        xyz     = entry["coords"][nearest]

        # Distance from bead to nuclear centre (0, 0, 0)
        distances[gene].append(np.linalg.norm(xyz))

    del bead_index

# ── Aggregate across all models and save ─────────────────────────────────────
records = []
for gene, dist_list in distances.items():
    records.append({
        "gene_name":   gene,
        "region_type": gene_region.get(gene, "unknown"),
        "mean_WT":     np.mean(dist_list),
        "std_WT":      np.std(dist_list),
        "n_models":    len(dist_list),
    })

wt_df = pd.DataFrame(records)
os.makedirs(os.path.dirname(CONFIG["output_file"]), exist_ok=True)
wt_df.to_csv(CONFIG["output_file"], sep="\t", index=False)

print(f"\nAggregated to {len(wt_df)} unique genes")
print(f"Saved WT distances to: {CONFIG['output_file']}")
