"""
03_compute_wt_distances.py
===========================
Script 3 of 4 – 3D distance calculation for wildtype (WT) models.

Computes the same gene-pair distances as script 02, but using wildtype
models instead of translocated ones. This gives us a baseline to compare
against: did the translocation bring genes closer or push them apart?

I reuse the exact gene pairs from script 02 so the comparison is fair.

Usage:  Edit the paths below, then run: python 03_compute_wt_distances.py
Dependencies:  pip install pandas numpy
"""

import os
import re
import glob
from collections import defaultdict

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG – edit these paths before running
# =============================================================================

T1_DISTANCES_FILE = "/path/to/T1_distances.tsv"
C1_DISTANCES_FILE = "/path/to/C1_distances.tsv"
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
WT_CMM_DIR = "/path/to/cmm/WT"
OUTPUT_FILE = "/path/to/WT_distances_agg.tsv"


# =============================================================================
# Load gene pairs from T1 and C1 (produced by script 02)
# =============================================================================

print("Computing WT distances...")

# Load distance data for translocated conditions
df_t1 = pd.read_csv(T1_DISTANCES_FILE, sep="\t")
df_c1 = pd.read_csv(C1_DISTANCES_FILE, sep="\t")

# Dictionary: gene --> region type (translocated, neighbor, control)
gene_region_type = {}
# Loop through both datasets
for df in (df_t1, df_c1):
    # Extract gene1 info
    for gene, region in zip(df["gene1"], df["gene1_region_type"]):
        gene_region_type[gene] = region
    # Extract gene2 info
    for gene, region in zip(df["gene2"], df["gene2_region_type"]):
        gene_region_type[gene] = region

# List of all unique genes involved in any pair
all_genes = list(gene_region_type.keys())
print(f"Total unique genes across T1 + C1: {len(all_genes)}")

# Collect all unique gene pairs (gene1, gene2)
gene_pairs = set()
for df in (df_t1, df_c1):
    gene_pairs.update(zip(df["gene1"], df["gene2"]))
print(f"Total gene pairs to compute: {len(gene_pairs)}")


# =============================================================================
# Load GTF to get reference genome coordinates
# =============================================================================

# WT cells have no rearrangement, so we use standard reference coordinates
gtf_cols = ["chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"]
# Load GTF file (ignore comment lines starting with '#')
gtf_df = pd.read_csv(GTF_FILE, sep="\t", comment="#", header=None, names=gtf_cols)
# Keep only gene entries (ignore transcripts, exons, etc.)
gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()
# Extract gene_name from attribute column
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')

# Only keep genes that appear in the dataset
gtf_df = gtf_df[gtf_df["gene_name"].isin(all_genes)].copy()
# Ensure chromosome names have "chr" prefix
gtf_df["chrom"] = gtf_df["chrom"].astype(str).apply(
    lambda x: x if x.startswith("chr") else "chr" + x
)
# Convert coordinates to numeric
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")


# =============================================================================
# Process each WT structural model
# =============================================================================

# Find all .cmm files (each is one 3D genome model)
cmm_files = sorted(glob.glob(os.path.join(WT_CMM_DIR, "*.cmm")))
print(f"Found {len(cmm_files)} CMM files to process")

# Dictionary: (gene1, geen2) --> list of distances across models
distances_dict = defaultdict(list)
# Cache chromosome assignments (same across models)
chr_dict_cache = {}

# Loop over each WT model
for f in cmm_files:
    print(f"\nProcessing model: {os.path.basename(f)}")

    # --------------------------------------------------
    # Parse .cmm file (extract bead positions)
    # --------------------------------------------------
    markers = []
    with open(f, "r") as fh:
        for line in fh:
            # Skip lines that are not bead markers
            if "<marker" not in line:
                continue
            # Extract attributes from XML-like structures
            attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
            # Skip if no beadID
            if "beadID" not in attrs:
                continue
            try:
                # beadID format: chr:start-end
                chrom_part, coords = attrs["beadID"].split(":")
                # Ensure chromosome format
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part
                # Extract genomic coordinates
                start_i, end_i = map(int, coords.split("-"))
                # Extract 3D coordinates (x, y, z)
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                # Store bead info
                markers.append([chrom_part, start_i, end_i, x, y, z])
            except Exception:
                continue    # Skip malformed entries

    # Convert markers to DataFrame
    beads = pd.DataFrame(markers, columns=["chr", "start", "end", "x", "y", "z"])
    # Skip if no beads found
    if beads.empty:
        print("  No markers found, skipping")
        continue
    
    # --------------------------------------------------
    # Build bead index (fast lookup)
    # --------------------------------------------------
    bead_index = {}
    # Group beads by chromosome
    for chrom, bdf in beads.groupby("chr"):
        bead_index[chrom] = {
            "mids": ((bdf["start"] + bdf["end"]) / 2).to_numpy(), # Genomic midpoints
            "coords": bdf[["x", "y", "z"]].to_numpy(),            # 3D coordinates
        }
    # Free memory
    del beads

    # --------------------------------------------------
    # Map genes to nearest bead
    # --------------------------------------------------
    coords_dict = {}    # Gene --> 3D coordinate
    chr_dict = {}       # Gene --> chromosome

    # Loop over each gene
    for _, row in gtf_df.iterrows():
        gene = row["gene_name"]
        chrom = row["chrom"]
        # Compute gene midpoint
        mid = (row["start"] + row["end"]) / 2
        # Get beads for this chromosome
        entry = bead_index.get(chrom)
        if entry is None:
            continue
        # Find closest bead
        nearest = np.abs(entry["mids"] - mid).argmin()
        # Assign coordinates
        coords_dict[gene] = entry["coords"][nearest]
        chr_dict[gene] = chrom

    # Save chromosome mapping (only need once)
    if not chr_dict_cache:
        chr_dict_cache = chr_dict

    # --------------------------------------------------
    # Compute distances for all gene pairs
    # --------------------------------------------------
    for gene1, gene2 in gene_pairs:
        # Only compute if both genes were mapped
        if gene1 in coords_dict and gene2 in coords_dict:
            # Euclidean distance in 3D
            dist = np.linalg.norm(coords_dict[gene1] - coords_dict[gene2])
            # Store distance for this pair
            distances_dict[(gene1, gene2)].append(dist)

    # Free memory
    del bead_index


# =============================================================================
# Aggregate across all models and save
# =============================================================================

records = []
# Loop over each gene pair
for (gene1, gene2), dist_list in distances_dict.items():
    records.append({
        "gene1": gene1,
        "gene1_region_type": gene_region_type.get(gene1, "unknown"),
        "chr1": chr_dict_cache.get(gene1, ""),
        "gene2": gene2,
        "gene2_region_type": gene_region_type.get(gene2, "unknown"),
        "chr2": chr_dict_cache.get(gene2, ""),
        "mean_distance": np.mean(dist_list),    # Average distance across models
        "std_distance": np.std(dist_list),      # Variability
        "n_models": len(dist_list),             # Number of models used
    })

# Convert to DataFrame
summary_df = pd.DataFrame(records)
# Ensure output directory exists
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
# Save results
summary_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"\nAggregated to {len(summary_df)} unique gene pairs")
print(f"Saved to: {OUTPUT_FILE}")
print("Done!")
