"""
03_compute_wt_distances.py
===========================
Script 3 of 4 for gene pair distance calculations – 3D distance calculation for the wildtype (WT) models.

What this script does
---------------------
This script computes the same gene-pair distances as script 02, but using
the wildtype structural models instead of the translocated ones.

Why do we need WT distances?
----------------------------
The translocated conditions (T1, C1) rearrange parts of the genome, which
may change how close genes are to each other in 3D space. To know whether
any distance change is caused by the translocation, we need a baseline —
what would those distances look like in a normal (wildtype) cell?

By computing distances for the exact same gene pairs in WT models, we can
directly compare:
    WT distance  vs  T1/C1 distance  --> did the translocation bring these
                                        genes closer together or push them apart?

Design choice
-------------
We use the gene pairs already identified in script 02 (from T1 and C1) so
that the comparison is fair — we measure the same pairs in all three
conditions (WT, T1, C1).

Gene coordinates in WT models come from the reference genome (GTF), not from
the translocated genome, because WT cells have no rearrangement.

Output of this script is used by script 04.

Usage
-----
    1. Run scripts 01 and 02 first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 03_compute_wt_distances.py

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
    # Filtered distance TSVs produced by script 02
    "t1_distances_file": "/path/to/T1_distances.tsv",
    "c1_distances_file": "/path/to/C1_distances.tsv",

    # GTF annotation for looking up WT gene coordinates
    "gtf_file": "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",

    # Directory of wildtype .cmm structural model files
    "wt_cmm_dir": "/path/to/cmm/WT",

    # Where to save the aggregated WT distance TSV
    "output_file": "/path/to/WT_distances_agg.tsv",
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def parse_cmm(file_path: str) -> pd.DataFrame:
    """
    Parse a UCSF Chimera .cmm marker file and return a DataFrame of beads
    with columns: chr, start, end, x, y, z.

    Each <marker> line with a beadID of the form 'chrN:start-end' represents
    one genomic bead positioned at (x, y, z) in 3D space.
    """
    markers = []
    with open(file_path, "r") as fh:
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
                start_i, end_i = map(int, coords.split("-"))
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, start_i, end_i, x, y, z])
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

    Returns a dict: chromosome → {'mids': array, 'coords': array, 'chr': str}
    """
    index = {}
    for chrom, bdf in beads_by_chr.items():
        index[chrom] = {
            "mids":   ((bdf["start"] + bdf["end"]) / 2).to_numpy(),
            "coords": bdf[["x", "y", "z"]].to_numpy(),
            "chr":    chrom,
        }
    return index


# =============================================================================
# Main
# =============================================================================

print("Computing WT distances...")

# ── Load gene pairs from T1 and C1 (produced by script 02) ───────────────────
df_t1 = pd.read_csv(CONFIG["t1_distances_file"], sep="\t")
df_c1 = pd.read_csv(CONFIG["c1_distances_file"], sep="\t")

# Collect all unique genes and their region types from both conditions.
# We need genes from both T1 and C1 so the WT baseline covers all conditions.
gene_region_type = {}
for df in (df_t1, df_c1):
    for gene, region in zip(df["gene1"], df["gene1_region_type"]):
        gene_region_type[gene] = region
    for gene, region in zip(df["gene2"], df["gene2_region_type"]):
        gene_region_type[gene] = region

all_genes = list(gene_region_type.keys())
print(f"Total unique genes across T1 + C1: {len(all_genes)}")

# ── Load GTF to get reference genome coordinates for each gene ────────────────
# We filter the GTF to only the genes we need, keeping memory usage down.
# WT cells have no rearrangement, so we use standard reference coordinates.
gtf_cols = ["chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
gtf_df = pd.read_csv(CONFIG["gtf_file"], sep="\t", comment="#", header=None, names=gtf_cols)
gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')
gtf_df = gtf_df[gtf_df["gene_name"].isin(all_genes)].copy()
gtf_df["chrom"] = gtf_df["chrom"].astype(str).apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"]   = pd.to_numeric(gtf_df["end"],   errors="coerce")

# ── Collect the exact gene pairs identified in T1 and C1 ─────────────────────
# Using the same pairs ensures a fair comparison across all three conditions
gene_pairs: set[tuple[str, str]] = set()
for df in (df_t1, df_c1):
    gene_pairs.update(zip(df["gene1"], df["gene2"]))
print(f"Total gene pairs to compute: {len(gene_pairs)}")

# ── Process each WT structural model ─────────────────────────────────────────
cmm_files = sorted(glob.glob(os.path.join(CONFIG["wt_cmm_dir"], "*.cmm")))
print(f"Found {len(cmm_files)} CMM files to process")

# Store all per-model distances in memory (one float per pair per model).
# This is manageable because we only store scalar distances, not matrices.
distances_dict: defaultdict[tuple, list] = defaultdict(list)
chr_dict_cache: dict[str, str] = {}

for f in cmm_files:
    print(f"\nProcessing model: {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No markers found, skipping")
        continue

    # Build numpy bead index and immediately free the raw DataFrame
    beads_by_chr = {c: b for c, b in beads.groupby("chr")}
    bead_index   = build_bead_index(beads_by_chr)
    del beads, beads_by_chr

    # Map every required gene to its nearest bead in this model
    coords_dict: dict[str, np.ndarray] = {}
    chr_dict:    dict[str, str]         = {}

    for _, row in gtf_df.iterrows():
        gene  = row["gene_name"]
        chrom = row["chrom"]
        mid   = (row["start"] + row["end"]) / 2
        entry = bead_index.get(chrom)
        if entry is None:
            continue
        nearest           = np.abs(entry["mids"] - mid).argmin()
        coords_dict[gene] = entry["coords"][nearest]
        chr_dict[gene]    = entry["chr"]

    # Cache chromosome assignments from the first model (consistent across models)
    if not chr_dict_cache:
        chr_dict_cache = chr_dict

    # Compute the Euclidean distance for each required gene pair
    for gene1, gene2 in gene_pairs:
        if gene1 in coords_dict and gene2 in coords_dict:
            dist = np.linalg.norm(coords_dict[gene1] - coords_dict[gene2])
            distances_dict[(gene1, gene2)].append(dist)

    del bead_index

# ── Aggregate across all models and save ─────────────────────────────────────
records = []
for (gene1, gene2), dist_list in distances_dict.items():
    records.append({
        "gene1":             gene1,
        "gene1_region_type": gene_region_type.get(gene1, "unknown"),
        "chr1":              chr_dict_cache.get(gene1, ""),
        "gene2":             gene2,
        "gene2_region_type": gene_region_type.get(gene2, "unknown"),
        "chr2":              chr_dict_cache.get(gene2, ""),
        "mean_distance":     np.mean(dist_list),
        "std_distance":      np.std(dist_list),
        "n_models":          len(dist_list),
    })

summary_df = pd.DataFrame(records)
os.makedirs(os.path.dirname(CONFIG["output_file"]), exist_ok=True)
summary_df.to_csv(CONFIG["output_file"], sep="\t", index=False)

print(f"\nAggregated to {len(summary_df)} unique gene pairs")
print(f"Saved aggregated distances to: {CONFIG['output_file']}")
