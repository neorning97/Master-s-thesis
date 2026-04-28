"""
01_extract_gene_expression.py
==============================
This script builds a table of genes for each condition (T1 and C1),
labelling each gene as either "inside_transloc", "neighbor", or "control",
and adding their RNA-seq expression values.

What this gives us:
- A table per condition with one row per gene.
- Each gene is tagged with what KIND of gene it is:
  * "inside_transloc": the gene sits inside a translocated region.
  * "neighbor":        the gene is in a region near where the translocation
                       landed on the receiving chromosome.
  * "control":         a randomly picked gene from a chromosome that wasn't
                       involved in the translocation. Used as a baseline
                       to compare against.
- And the gene's expression values are attached.

Control genes:
- Control genes are sampled at random from chromosomes that were NOT
  involved in the translocation.
- We pick the SAME NUMBER of controls as we have translocated genes
  (so the comparison is fair).

Edit the file paths below and run the script.

Required libraries: pandas, bioframe
Install them with: pip install pandas bioframe
"""

import os            
import pandas as pd
import bioframe     


# =============================================================================
# CONFIG SECTION - Edit these paths before running
# =============================================================================

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------
# GTF = standard genome annotation file (tells us where every gene is)
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# TPM file = gene expression values (TPM = Transcripts Per Million)
EXPR_FILE = "/path/to/GSE246689_gene_tpm.tsv"

# Where to save the output TSV files
RESULTS_DIR = "/path/to/results"

# -----------------------------------------------------------------------------
# BED files for translocations and their neighbors, per condition
# -----------------------------------------------------------------------------
# A BED file lists genomic regions: chromosome + start + end + extra columns.
# - transloc_bed: regions that got translocated
# - neighbor_bed: regions on the receiving chromosome NEXT TO the translocation
T1_TRANSLOC_BED = "/path/to/T1_translocations.bed"
T1_NEIGHBOR_BED = "/path/to/T1_neighbor.bed"
C1_TRANSLOC_BED = "/path/to/C1_translocations.bed"
C1_NEIGHBOR_BED = "/path/to/C1_neighbor.bed"

# Group these in a dictionary so we can loop through both conditions
conditions = {
    "T1": {
        "transloc_bed": T1_TRANSLOC_BED,
        "neighbor_bed": T1_NEIGHBOR_BED,
    },
    "C1": {
        "transloc_bed": C1_TRANSLOC_BED,
        "neighbor_bed": C1_NEIGHBOR_BED,
    },
}


# =============================================================================
# Step 1: Load the GTF file (do this once, then reuse for both conditions)
# =============================================================================

print("Loading GTF annotation...")

# bioframe.read_table can read GTFs and parse the columns automatically.
# schema="gtf" tells it the file is in GTF format.
gtf = bioframe.read_table(GTF_FILE, schema="gtf")

# Keep only rows where feature is "gene"
# (a GTF can also have rows for transcripts, exons, etc.)
gtf = gtf[gtf["feature"] == "gene"].copy()

# -----------------------------------------------------------------------------
# Convert GTF coordinates to match BED file convention
# -----------------------------------------------------------------------------
# GTF files use 1-based coordinates (the first base is position 1).
# BED files use 0-based coordinates (the first base is position 0).
# So we subtract 1 from each start position to convert.
# We don't need to change "end" because BED is half-open ([start, end))
# and GTF is closed ([start, end]). The math works out so end stays the same.
gtf["start"] = gtf["start"] - 1

# Make sure chromosome names are clean strings (no leading/trailing spaces)
gtf["chrom"] = gtf["chrom"].astype(str).str.strip()

# -----------------------------------------------------------------------------
# Pull gene_id, gene_name, gene_type out of the "attributes" column
# -----------------------------------------------------------------------------
# The attributes column looks like:
#   gene_id "ENSG00..."; gene_name "ABC1"; gene_type "protein_coding"
# str.extract with a regex grabs the part inside the parentheses ( ).
gtf["gene_id"]   = gtf["attributes"].str.extract(r'gene_id "([^"]+)"')
gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')
gtf["gene_type"] = gtf["attributes"].str.extract(r'gene_type "([^"]+)"')

# -----------------------------------------------------------------------------
# Strip version numbers from gene IDs
# -----------------------------------------------------------------------------
# Ensembl IDs sometimes have versions like "ENSG00000001234.5".
# We want just "ENSG00000001234" so we can match with the expression file.
# .str.split(".") splits on the dot, .str[0] takes the first part.
gtf["gene_id"] = gtf["gene_id"].str.split(".").str[0]


# =============================================================================
# Step 2: Load the expression data
# =============================================================================

print("Loading expression data...")

# Read the TPM file (tab-separated)
expr = pd.read_csv(EXPR_FILE, sep="\t")

# -----------------------------------------------------------------------------
# Some expression files use "gene" as the column name, others use "gene_id".
# Standardize to "gene_id" so the merge later works.
# -----------------------------------------------------------------------------
if "gene" in expr.columns:
    expr = expr.rename(columns={"gene": "gene_id"})

# -----------------------------------------------------------------------------
# Same trick as above: strip version numbers from gene IDs
# -----------------------------------------------------------------------------
expr["gene_id"] = expr["gene_id"].astype(str).str.split(".").str[0]

# -----------------------------------------------------------------------------
# Drop the gene_name column if it's there
# -----------------------------------------------------------------------------
# We'll get gene names from the GTF file instead, they're more reliable.
if "gene_name" in expr.columns:
    expr = expr.drop(columns=["gene_name"])


# =============================================================================
# Step 3: Process each condition
# =============================================================================

# Make sure the output folder exists
os.makedirs(RESULTS_DIR, exist_ok=True)

# Loop through each condition (T1, C1)
# .items() gives both the key (cond) and value (paths dict)
for cond, paths in conditions.items():

    print("")
    print("Processing condition: " + cond)

    # -------------------------------------------------------------------------
    # Step 3a: Load this condition's BED files
    # -------------------------------------------------------------------------

    # Load the translocation BED file
    transloc_df = pd.read_csv(paths["transloc_bed"], sep="\t")

    # Clean up the chromosome and coordinate columns.
    # pd.to_numeric converts text to numbers.
    # errors="coerce" means: if something can't be converted, make it NaN
    # (instead of crashing). Then .astype(int) converts to integers.
    transloc_df["chrom"] = transloc_df["chrom"].astype(str).str.strip()
    transloc_df["start"] = pd.to_numeric(transloc_df["start"], errors="coerce").astype(int)
    transloc_df["end"] = pd.to_numeric(transloc_df["end"], errors="coerce").astype(int)

    # Load the neighbor BED file
    neighbor_df = pd.read_csv(paths["neighbor_bed"], sep="\t")

    # Clean up neighbor coordinate columns
    # Note: the chromosome column for neighbors is called "neighbor_chr"
    # (not "chrom") because it's a different chromosome than the original.
    if "chrom" in neighbor_df.columns:
        neighbor_df["chrom"] = neighbor_df["chrom"].astype(str).str.strip()

    neighbor_df["start"] = pd.to_numeric(neighbor_df["start"], errors="coerce").astype(int)
    neighbor_df["end"] = pd.to_numeric(neighbor_df["end"], errors="coerce").astype(int)

    # We'll collect all the gene-hit DataFrames in this list,
    # then combine them at the end with pd.concat.
    all_genes = []

    # -------------------------------------------------------------------------
    # Step 3b: Find genes inside the translocated regions
    # -------------------------------------------------------------------------

    # Loop through each translocated region
    for _, row in transloc_df.iterrows():

        # Find all genes overlapping this region.
        # The overlap rule: same chromosome AND
        # gene's start < region's end AND
        # gene's end > region's start
        same_chrom = (gtf["chrom"] == row["chrom"])
        starts_before_end = (gtf["start"] < row["end"])
        ends_after_start = (gtf["end"] > row["start"])

        # Combine with &
        mask = same_chrom & starts_before_end & ends_after_start

        # Get the matching genes and make a copy
        hits = gtf[mask].copy()

        # If we found any genes, add labels and store them
        if not hits.empty:
            hits["region_type"] = "inside_transloc"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"] = row["transloc_id"]
            all_genes.append(hits)

    # -------------------------------------------------------------------------
    # Step 3c: Find genes in the neighbor regions
    # -------------------------------------------------------------------------
    # Same idea, but the chromosome column is called "neighbor_chr" here
    # (because neighbors are on the receiving chromosome).

    for _, row in neighbor_df.iterrows():

        # Same overlap rule as above
        same_chrom = (gtf["chrom"] == row["neighbor_chr"])
        starts_before_end = (gtf["start"] < row["end"])
        ends_after_start = (gtf["end"] > row["start"])
        mask = same_chrom & starts_before_end & ends_after_start

        hits = gtf[mask].copy()

        if not hits.empty:
            hits["region_type"] = "neighbor"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"] = row["transloc_id"]
            all_genes.append(hits)

    # -------------------------------------------------------------------------
    # If we didn't find any genes at all, skip to the next condition
    # -------------------------------------------------------------------------
    if len(all_genes) == 0:
        print("  No genes found, skipping " + cond)
        continue

    # Combine all the gene hits into one DataFrame
    # ignore_index=True gives fresh row numbers (instead of duplicates)
    all_genes_df = pd.concat(all_genes, ignore_index=True)

    # -------------------------------------------------------------------------
    # Step 3d: Sample control genes from uninvolved chromosomes
    # -------------------------------------------------------------------------
    # Control genes are random genes from chromosomes that were NOT
    # involved in any translocation. They serve as a "baseline" to
    # compare against.

    # Get the set of chromosomes involved in translocations
    trans_chroms = set(transloc_df["chrom"])

    # The "control pool" = all genes on chromosomes NOT in trans_chroms.
    # ~ inverts a boolean (NOT in trans_chroms).
    # .isin() returns True/False per row.
    control_pool = gtf[~gtf["chrom"].isin(trans_chroms)]

    # We'll collect each batch of control genes here
    control_rows = []

    # Loop through each unique transloc_id (one batch of controls per translocation)
    unique_tids = all_genes_df["transloc_id"].unique()

    for tid in unique_tids:

        # ---------------------------------------------------------------------
        # Count how many translocated genes this transloc_id has
        # ---------------------------------------------------------------------
        # We want to sample the SAME NUMBER of control genes for a fair comparison.
        is_this_tid = (all_genes_df["transloc_id"] == tid)
        is_inside = (all_genes_df["region_type"] == "inside_transloc")
        mask = is_this_tid & is_inside

        # .sum() on a boolean Series counts the Trues
        n_trans = mask.sum()

        # If this translocation has no inside genes, skip
        if n_trans == 0:
            continue

        # ---------------------------------------------------------------------
        # Sample n_trans random genes from the control pool
        # ---------------------------------------------------------------------
        # We use a "seed" to make the random sampling REPRODUCIBLE.
        # If you re-run the script, you'll get the same control genes.
        # The seed is derived from the translocation ID so each translocation
        # gets a different (but stable) sample.
        # abs(hash(tid)) % (2**32) just turns the tid string into a number
        # that fits in 32 bits (a typical seed format).
        seed = abs(hash(tid)) % (2**32)

        # .sample(n=...) picks n random rows
        # random_state=seed makes it reproducible
        sampled = control_pool.sample(n=n_trans, random_state=seed).copy()

        # Add the labels
        sampled["region_type"] = "control"
        sampled["transloc_id"] = tid
        sampled["translocation_label"] = "control"

        control_rows.append(sampled)

    # -------------------------------------------------------------------------
    # Add the control genes to our gene table
    # -------------------------------------------------------------------------
    if len(control_rows) > 0:
        # The list [all_genes_df] + control_rows builds a list of DataFrames
        # to concatenate. all_genes_df becomes the first item.
        all_genes_df = pd.concat(
            [all_genes_df] + control_rows,
            ignore_index=True
        )

    # -------------------------------------------------------------------------
    # Step 3e: Add expression data (TPM values) to each gene
    # -------------------------------------------------------------------------
    # merge() joins two tables on a common column.
    # how="left" keeps all rows from all_genes_df (the left table),
    # adding NaN for genes that don't have expression data.
    merged = all_genes_df.merge(expr, on="gene_id", how="left")

    # -------------------------------------------------------------------------
    # Step 3f: Reorder the columns
    # -------------------------------------------------------------------------
    # Put the gene metadata columns first, then the expression columns.

    # The "base columns" we always want first
    base_cols = [
        "chrom", "start", "end",
        "gene_id", "gene_name", "gene_type",
        "region_type", "translocation_label", "transloc_id",
        "source", "feature", "score", "strand", "frame",
    ]

    # The expression columns all start with "MCF10A".
    # We pick them out with a list comprehension.
    # (A list comprehension is a one-line for loop that builds a list.)
    expr_cols = []
    for c in merged.columns:
        if c.startswith("MCF10A"):
            expr_cols.append(c)

    # Reorder the columns: base columns first, then expression columns
    merged = merged[base_cols + expr_cols]

    # -------------------------------------------------------------------------
    # Step 3g: Save the result
    # -------------------------------------------------------------------------
    output_filename = cond + "_genes_expression.tsv"
    output_file = os.path.join(RESULTS_DIR, output_filename)

    # Save as TSV (tab-separated)
    # index=False means: don't write the row numbers as a column
    merged.to_csv(output_file, sep="\t", index=False)

    print("  Saved: " + output_file)


# Done!
print("")
print("Done!")
