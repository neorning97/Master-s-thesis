"""
01_extract_gene_expression.py
==============================
Script 1 of 4 for gene pair distance calculations – Gene expression extraction.

For each condition (T1 and C1), this script:
1. Reads BED files defining translocated and neighbor regions.
2. Finds all genes in those regions using a GTF annotation.
3. Samples matching control genes from uninvolved chromosomes.
4. Attaches RNA-seq TPM expression values to every gene.
5. Saves the result as a TSV file.

Usage:  Edit the paths below, then run: python 01_extract_gene_expression.py
Dependencies:  pip install pandas bioframe
"""

import os
import pandas as pd
import bioframe


# =============================================================================
# CONFIG – edit these paths before running
# =============================================================================

#GTF_FILE  = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
#EXPR_FILE = "/path/to/GSE246689_gene_tpm.tsv"
#RESULTS_DIR = "/path/to/results"
GTF_FILE  = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
EXPR_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/GSE246689_gene_tpm.tsv"
RESULTS_DIR = "/Users/nadiaorning/Desktop/UiO/Vår2026/results"



# Dictionary defining experimental conditions and their input files:
# - transloc_bed: genomic regions with translocations
# - neighbor_bed: nearby regions on receiving chromosome
CONDITIONS = {
    "T1": {
        "transloc_bed": "/path/to/T1_translocations.bed",
        "neighbor_bed": "/path/to/T1_neighbor.bed",
    },
    "C1": {
        "transloc_bed": "/path/to/C1_translocations.bed",
        "neighbor_bed": "/path/to/C1_neighbor.bed",
    },
}

# =============================================================================
# Load GTF genes
# =============================================================================

print("Loading GTF annotation...")
# Reads GTF file into DataFrame, and keeps only rows describing genes
gtf = bioframe.read_table(GTF_FILE, schema="gtf")
gtf = gtf[gtf["feature"] == "gene"].copy()

# GTF uses 1-based coordinates; convert start to 0-based to match BED files
gtf["start"] -= 1
gtf["chrom"] = gtf["chrom"].astype(str).str.strip()

# Pull gene metadata out of the attributes column using regex
gtf["gene_id"]   = gtf["attributes"].str.extract(r'gene_id "([^"]+)"')
gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')
gtf["gene_type"] = gtf["attributes"].str.extract(r'gene_type "([^"]+)"')

# Strip version suffixes (e.g. ENSG00000001234.5 -> ENSG00000001234)
gtf["gene_id"] = gtf["gene_id"].str.split(".").str[0]


# =============================================================================
# Load expression data
# =============================================================================

print("Loading expression data...")
expr = pd.read_csv(EXPR_FILE, sep="\t")

# Checks if there is a gene column, and if there is -> renames to gene_id
if "gene" in expr.columns:
    expr = expr.rename(columns={"gene": "gene_id"})

# Convert gene_id to string, if not already, and strip version suffixes
expr["gene_id"] = expr["gene_id"].astype(str).str.split(".").str[0]

# Drop gene_name if it exists -> we'll get it from the GTF instead
if "gene_name" in expr.columns:
    expr = expr.drop(columns=["gene_name"])


# =============================================================================
# Process each condition
# =============================================================================

os.makedirs(RESULTS_DIR, exist_ok=True)

# Loop through each condition (T1, C1)
for cond, paths in CONDITIONS.items():
    print(f"\nProcessing condition: {cond}")

    #---------------------------------
    # Load translocation regions
    #---------------------------------
    transloc_df = pd.read_csv(paths["transloc_bed"], sep="\t")

    # Clean and standardize columns
    transloc_df["chrom"] = transloc_df["chrom"].astype(str).str.strip()
    transloc_df["start"] = pd.to_numeric(transloc_df["start"], errors="coerce").astype(int)
    transloc_df["end"]   = pd.to_numeric(transloc_df["end"],   errors="coerce").astype(int)

    #---------------------------------
    # Load neighboring regions
    #---------------------------------
    neighbor_df = pd.read_csv(paths["neighbor_bed"], sep="\t")

    # Clean and standardize columns
    neighbor_df["chrom"] = neighbor_df["chrom"].astype(str).str.strip() if "chrom" in neighbor_df.columns else None
    neighbor_df["start"] = pd.to_numeric(neighbor_df["start"], errors="coerce").astype(int)
    neighbor_df["end"]   = pd.to_numeric(neighbor_df["end"],   errors="coerce").astype(int)

    # Container for all gene hits
    all_genes = []

    #-------------------------------------------
    # Find genes inside translocated regions 
    #-------------------------------------------
    for _, row in transloc_df.iterrows():
        # Select genes overlapping this region
        hits = gtf[
            (gtf["chrom"] == row["chrom"]) &
            (gtf["start"] < row["end"]) &
            (gtf["end"]   > row["start"])
        ].copy()

        # If any genes overlap, annotate and store them
        if not hits.empty:
            hits["region_type"] = "inside_transloc"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"] = row["transloc_id"]
            all_genes.append(hits)

    #-------------------------------------------
    # Find genes in neighboring regions 
    #-------------------------------------------
    # Neighbor BED uses 'neighbor_chr' because it's on the receiving chromosome
    for _, row in neighbor_df.iterrows():
        hits = gtf[
            (gtf["chrom"] == row["neighbor_chr"]) &
            (gtf["start"] < row["end"]) &
            (gtf["end"]   > row["start"])
        ].copy()

        if not hits.empty:
            hits["region_type"] = "neighbor"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"] = row["transloc_id"]
            all_genes.append(hits)

    # If no genes found at all --> skip condition
    if not all_genes:
        print(f"  No genes found – skipping {cond}")
        continue

    # Combine all found genes into one DataFrame
    all_genes_df = pd.concat(all_genes, ignore_index=True)

    #-------------------------------------------------------------------------
    # Sample control genes from chromosomes not involved in the translocation 
    #-------------------------------------------------------------------------
    # Identify chromosomes involved in translocations
    trans_chroms = set(transloc_df["chrom"])

    # Control genes = genes NOT on those chromosomes
    control_pool = gtf[~gtf["chrom"].isin(trans_chroms)]

    control_rows = []

    # Loop per translocation ID
    for tid in all_genes_df["transloc_id"].unique():
        # Count how many translocated genes this translocation has
        mask = (all_genes_df["transloc_id"] == tid) & (all_genes_df["region_type"] == "inside_transloc")
        n_trans = mask.sum()

        if n_trans == 0:
            continue

        # Sample same number of control genes (reproducible via seed)
        seed = abs(hash(tid)) % (2**32)
        sampled = control_pool.sample(n=n_trans, random_state=seed).copy()
        sampled["region_type"] = "control"
        sampled["transloc_id"] = tid
        sampled["translocation_label"] = "control"

        control_rows.append(sampled)

    # Add control genes to dataset
    if control_rows:
        all_genes_df = pd.concat([all_genes_df] + control_rows, ignore_index=True)

    #-------------------------------------------
    # Merge gene annotations with expression data 
    #-------------------------------------------
    merged = all_genes_df.merge(expr, on="gene_id", how="left")

    # Put gene metadata columns first, then all MCF10A expression columns
    base_cols = [
        "chrom", "start", "end",
        "gene_id", "gene_name", "gene_type",
        "region_type", "translocation_label", "transloc_id",
        "source", "feature", "score", "strand", "frame",
    ]

    # Select expression columns
    expr_cols = [c for c in merged.columns if c.startswith("MCF10A")]
    merged = merged[base_cols + expr_cols]

    #-------------------------------------------
    # Save results 
    #-------------------------------------------
    output_file = os.path.join(RESULTS_DIR, f"{cond}_genes_expression.tsv")
    merged.to_csv(output_file, sep="\t", index=False)
    print(f"  Saved: {output_file}")

print("\nDone!")
