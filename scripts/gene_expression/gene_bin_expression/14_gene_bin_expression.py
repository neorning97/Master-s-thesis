"""
14_gene_bin_expression.py
==========================
This script connects the bin-classification analysis (script 10) with the
differential expression analysis (script 13).

The big question:
- Script 10 told us how each 100 kb bin behaved after translocation:
  retained / adopted / other.
- Script 13 told us which genes went UP, DOWN, or didn't change (ns).
- This script asks: do bins that ADOPTED a new subcompartment contain
  more differentially expressed genes than bins that just RETAINED
  their original subcompartment?

If the answer is yes, that supports the idea that compartment switching
DRIVES the gene expression changes.

Step by step:
1. Load the bin classifications (from script 10).
2. Load gene coordinates from a GTF file.
3. Map each gene to the bin it sits in (so each gene inherits the
   bin's classification labels).
4. Load TPM expression values, calculate log2 fold changes vs WT.
5. Load DESeq2 results (from script 13) and label each gene as
   "up" / "down" / "ns" based on log2 fold change and padj thresholds.
6. For each condition (T1, C1) and each classification scheme, run:
   - Kruskal-Wallis test: do the log2 fold changes differ between
     bin categories (retained vs adopted vs other)?
   - Chi-square test: do the up/down/ns proportions differ between
     bin categories?
7. Make stacked bar plots.

Why these specific stat tests?
- Kruskal-Wallis is like a t-test but doesn't assume normal distribution.
  Used to compare CONTINUOUS values (log2 fold changes) across groups.
- Chi-square is for COUNTS in categories. We use it to test if the
  up/down/ns counts differ between bin categories.

Run scripts 10 and 13 first.
Then edit the file paths below and run this script.

Required libraries: pandas, numpy, matplotlib, scipy
Install them with: pip install pandas numpy matplotlib scipy
"""

import os                                          
import numpy as np                                  
import pandas as pd                                 
import matplotlib.pyplot as plt                     
from scipy.stats import kruskal, chi2_contingency   


# =============================================================================
# CONFIG SECTION - Edit these paths and settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# Bin classification table (from script 10)
# -----------------------------------------------------------------------------
BINS_FILE = "/path/to/bins_annotation.csv"

# -----------------------------------------------------------------------------
# GTF gene annotation file
# -----------------------------------------------------------------------------
# Standard genome annotation file. We use it to get gene coordinates.
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# -----------------------------------------------------------------------------
# RNA-seq TPM expression file
# -----------------------------------------------------------------------------
# TPM = Transcripts Per Million, a unit of gene expression.
# Higher TPM means a gene is more strongly expressed.
TPM_FILE = "/path/to/GSE246689_gene_tpm.tsv"

# -----------------------------------------------------------------------------
# DESeq2 result files (from script 13)
# -----------------------------------------------------------------------------
T1_DE_FILE = "/path/to/DE_WT_vs_T1.tsv"
C1_DE_FILE = "/path/to/DE_WT_vs_C1.tsv"

# Group them in a dict for looping
de_files = {
    "T1": T1_DE_FILE,
    "C1": C1_DE_FILE,
}

# -----------------------------------------------------------------------------
# Replicate column names in the TPM file
# -----------------------------------------------------------------------------
# Each condition has 3 biological replicates. We'll average them.
WT_REPLICATES = ["MCF10AWT_REP1", "MCF10AWT_REP2", "MCF10AWT_REP3"]
T1_REPLICATES = ["MCF10AT1_REP1", "MCF10AT1_REP2", "MCF10AT1_REP3"]
C1_REPLICATES = ["MCF10AC1_REP1", "MCF10AC1_REP2", "MCF10AC1_REP3"]

# Group them for looping
tpm_replicates = {
    "WT": WT_REPLICATES,
    "T1": T1_REPLICATES,
    "C1": C1_REPLICATES,
}

# -----------------------------------------------------------------------------
# Thresholds for calling a gene differentially expressed
# -----------------------------------------------------------------------------
# A gene is called "up" or "down" if BOTH:
#   - the adjusted p-value is below 0.05 (statistically significant)
#   - the absolute log2 fold change is above 1.0 (i.e. at least 2x change)
# Otherwise it's "ns" (not significant).
LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

# -----------------------------------------------------------------------------
# Output folder
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/results/plots/gene_bin_expression"

# Conditions we'll process (we don't process WT, it's the reference)
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Load gene coordinates from a GTF file
# =============================================================================
# Same idea as in script 06, read the GTF, keep only "gene" rows, pull
# out the gene_id and gene_name from the "attributes" column.

def load_gtf_genes(gtf_file):
    """
    Load a GTF file and return a DataFrame with one row per gene.
    Columns: chrom, start, end, gene_id, gene_name
    """
    # GTF files have 9 columns, we provide names since GTF files don't
    # have a header row.
    column_names = [
        "chrom", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes"
    ]

    # Read the GTF file. comment="#" skips comment lines starting with #.
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=column_names
    )

    # Keep only rows where feature is "gene"
    # (a GTF can have rows for genes, transcripts, exons, etc.)
    gtf = gtf[gtf["feature"] == "gene"].copy()

    # Extract gene_id and gene_name from the attributes column.
    # The attributes column looks like:
    #   gene_id "ENSG00..."; gene_name "ABC1"; gene_biotype "protein_coding"
    # str.extract grabs the part inside the parentheses ( ) of the regex.
    gtf["gene_id"] = gtf["attributes"].str.extract(r'gene_id "([^"]+)"')
    gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')

    # Keep only the columns we need
    keep_cols = ["chrom", "start", "end", "gene_id", "gene_name"]
    genes = gtf[keep_cols].copy()

    # Make sure start and end are integers
    genes["start"] = genes["start"].astype(int)
    genes["end"] = genes["end"].astype(int)

    # Reset row numbers so they're 0, 1, 2, ... again
    genes = genes.reset_index(drop=True)

    return genes


# =============================================================================
# Helper function: Find which bin a gene sits in
# =============================================================================
# For each gene, we want to find the first bin that overlaps with it,
# and inherit its classification labels.

def assign_bin(gene_row, bins_df):
    """
    Find the first bin that overlaps the given gene.
    Returns a tuple of 5 labels from that bin:
    (transloc_id, condition, behavior, behavior_collapsed, change_direction)
    Returns NaNs if no bin overlaps.
    """
    # Same overlap idea as in earlier scripts:
    # two intervals [a, b) and [c, d) overlap if a < d AND b > c
    same_chrom = (bins_df["chrom"] == gene_row["chrom"])
    starts_before_end = (bins_df["start"] < gene_row["end"])
    ends_after_start = (bins_df["end"] > gene_row["start"])

    # Combine all three with & (AND)
    mask = same_chrom & starts_before_end & ends_after_start

    # Get the rows that overlap
    overlap = bins_df[mask]

    # If no overlap, return all NaNs
    if overlap.empty:
        # This is a tuple of 5 NaN values
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

    # Take the first overlapping bin
    # iloc[0] gets the first row by position
    first_bin = overlap.iloc[0]

    # Build the tuple of labels and return it
    labels = (
        first_bin["transloc_id"],
        first_bin["condition"],
        first_bin["behavior"],
        first_bin["behavior_collapsed"],
        first_bin["change_direction"],
    )
    return labels


# =============================================================================
# Helper function: Classify a gene as up / down / ns
# =============================================================================

def de_call(log2fc, padj):
    """
    Classify a gene based on its log2 fold change and adjusted p-value.
    Returns "up", "down", or "ns" (not significant).
    """
    # If either value is missing (NaN), call it "ns"
    # pd.notna returns True if the value is NOT NaN
    if pd.isna(log2fc) or pd.isna(padj):
        return "ns"

    # If the p-value isn't small enough, call it "ns"
    if padj >= PADJ_THRESHOLD:
        return "ns"

    # The padj is low enough, now check the fold change direction
    if log2fc > LFC_THRESHOLD:
        return "up"
    elif log2fc < -LFC_THRESHOLD:
        return "down"
    else:
        # padj is significant but the fold change is too small
        return "ns"


# =============================================================================
# Helper function: Run stats and make a plot for one condition + scheme
# =============================================================================
# We use this 3 times per condition (once per classification scheme).

def run_analysis(df_bins, df_genome, cond, category_col,
                 category_order, label):
    """
    Run Kruskal-Wallis and chi-square tests for one condition and one
    classification scheme, then save a stacked bar plot.

    - df_bins:       genes that overlap a translocated bin (with bin labels)
    - df_genome:     all genes (used for the genome-wide control bar)
    - cond:          "T1" or "C1"
    - category_col:  which column holds the classification (e.g. "behavior")
    - category_order: list of expected category values in display order
    - label:         readable label for the plot title
    """
    # -------------------------------------------------------------------------
    # Filter to just this condition's data
    # -------------------------------------------------------------------------
    df = df_bins[df_bins["condition"] == cond].copy()

    # Build the column names that contain the data we need
    # f-strings would be cleaner but we use string concatenation here for clarity
    de_col = "DE_" + cond
    fc_col = "log2FC_" + cond + "_vs_WT"

    # -------------------------------------------------------------------------
    # Step 1: Kruskal-Wallis test
    # -------------------------------------------------------------------------
    # Question: do the log2 fold changes differ between bin categories?
    #
    # Kruskal-Wallis takes several lists of numbers (one per group) and
    # tests whether they have the same median.

    # Build a list of fold change values, one list per category.
    # We skip categories that aren't actually present in this data.
    groups = []
    for cat in category_order:

        # Check if this category exists in the data
        if cat not in df[category_col].values:
            continue

        # Get the fold changes for this category
        cat_rows = df[df[category_col] == cat]
        fold_changes = cat_rows[fc_col]

        # Drop any NaN values (Kruskal-Wallis can't handle them)
        fold_changes = fold_changes.dropna()

        # Add this category's fold changes to our list
        groups.append(fold_changes)

    # Run the test if we have at least 2 groups
    # The * unpacks the list into separate arguments (kruskal expects
    # each group as a separate argument, not one list)
    if len(groups) > 1:
        kw_result = kruskal(*groups)
        p_kw = kw_result.pvalue
    else:
        # Not enough groups to run the test
        p_kw = np.nan

    # -------------------------------------------------------------------------
    # Step 2: Chi-square test
    # -------------------------------------------------------------------------
    # Question: do the up/down/ns proportions differ between bin categories?
    #
    # We need a "contingency table" with one row per bin category and one
    # column per DE status. pd.crosstab does this for us.

    ct_counts = pd.crosstab(df[category_col], df[de_col])

    # Run the test if we have at least 2 rows (categories)
    if ct_counts.shape[0] > 1:
        # chi2_contingency returns a tuple (chi2, p, dof, expected)
        # We just want the p-value (index 1)
        chi_result = chi2_contingency(ct_counts)
        p_chi = chi_result[1]
    else:
        p_chi = np.nan

    # -------------------------------------------------------------------------
    # Step 3: Print results
    # -------------------------------------------------------------------------
    print("")
    print(cond + " — " + label)
    # {:.3e} = scientific notation with 3 decimals
    print("  Kruskal-Wallis p = {:.3e}".format(p_kw))
    print("  Chi-square p     = {:.3e}".format(p_chi))
    print(ct_counts)

    # -------------------------------------------------------------------------
    # Step 4: Build the genome-wide control bar
    # -------------------------------------------------------------------------
    # The genome-wide control shows the background rate of differential
    # expression across all genes (not just genes in translocated bins).
    # We make a fake "category" column where every gene is in the
    # "genome-wide" group, then build a contingency table from that.

    # Make a Series with "genome-wide" repeated for every row
    fake_category = pd.Series(
        ["genome-wide"] * len(df_genome),
        name=category_col   # name has to match for crosstab
    )

    # Make the contingency table
    ct_genome = pd.crosstab(fake_category, df_genome[de_col])

    print("")
    print("  Genome-wide:")
    print(ct_genome)

    # -------------------------------------------------------------------------
    # Step 5: Combine and normalize for plotting
    # -------------------------------------------------------------------------
    # Stack the two contingency tables on top of each other
    ct_plot = pd.concat([ct_counts, ct_genome])

    # Convert counts to fractions (each row should add up to 1)
    # ct_plot.sum(axis=1) gives the row totals
    # .div(... , axis=0) divides each row by its total
    row_sums = ct_plot.sum(axis=1)
    ct_frac = ct_plot.div(row_sums, axis=0)

    # Make sure the columns are in the order ns, up, down (so the plot
    # always has the same colour stacking).
    # fill_value=0 means: if a column is missing, fill it with zeros.
    ct_frac = ct_frac.reindex(columns=["ns", "up", "down"], fill_value=0)

    # -------------------------------------------------------------------------
    # Step 6: Make and save the stacked bar plot
    # -------------------------------------------------------------------------
    # Create the figure
    fig, ax = plt.subplots(figsize=(6, 4))

    # Define colours for each DE status
    color_dict = {
        "ns": "#d9d9d9",   # light grey
        "up": "green",
        "down": "red",
    }

    # Make the stacked bar plot using pandas built-in plotting
    ct_frac.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=color_dict
    )

    # Add labels and title
    ax.set_ylabel("Fraction of genes")
    ax.set_title(cond + ": " + label)

    plt.tight_layout()

    # Build the output filename and save
    filename = cond + "_" + category_col + "_DE_fraction.png"
    plot_path = os.path.join(OUTPUT_FOLDER, filename)
    plt.savefig(plot_path, dpi=300)

    # Close the figure to free memory
    plt.close(fig)


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)


# -----------------------------------------------------------------------------
# Step 1: Load the bin classifications (from script 10)
# -----------------------------------------------------------------------------
print("Loading bin classifications...")

bins = pd.read_csv(BINS_FILE)
bins["start"] = bins["start"].astype(int)
bins["end"] = bins["end"].astype(int)


# -----------------------------------------------------------------------------
# Step 2: Load gene coordinates from the GTF file
# -----------------------------------------------------------------------------
print("Loading GTF gene annotations...")
genes = load_gtf_genes(GTF_FILE)


# -----------------------------------------------------------------------------
# Step 3: Load the TPM expression file and average replicates
# -----------------------------------------------------------------------------
print("Loading TPM expression...")
expr = pd.read_csv(TPM_FILE, sep="\t")

# For each condition, average its 3 replicates into a single column
# Loop through the dict: cond is the key, reps is the list of column names
for cond, reps in tpm_replicates.items():
    # mean(axis=1) means "take the mean across columns, for each row"
    expr[cond] = expr[reps].mean(axis=1)


# -----------------------------------------------------------------------------
# Step 4: Build the genome-wide gene expression table
# -----------------------------------------------------------------------------
# This table has ALL genes (not just genes in translocated bins).
# It's used for the genome-wide control bar in the plots.
print("Building genome-wide gene expression table...")

# Merge gene coordinates with expression values.
# how="left" keeps all genes even if they're missing from the expression file.
gene_expr = genes.merge(
    expr[["gene_id", "WT", "T1", "C1"]],
    on="gene_id",
    how="left"
)

# Calculate log2(TPM + 1) for each condition.
# We add 1 first so log2(0) doesn't blow up to negative infinity.
for cond in ["WT", "T1", "C1"]:
    log_col_name = "log2_" + cond
    gene_expr[log_col_name] = np.log2(gene_expr[cond] + 1)

# Calculate log2 fold change vs WT
for cond in ["T1", "C1"]:
    fc_col_name = "log2FC_" + cond + "_vs_WT"
    log_cond_col = "log2_" + cond
    log_wt_col = "log2_WT"
    gene_expr[fc_col_name] = gene_expr[log_cond_col] - gene_expr[log_wt_col]


# -----------------------------------------------------------------------------
# Step 5: Load DESeq2 results and call up/down/ns for each gene
# -----------------------------------------------------------------------------
print("Loading DESeq2 results...")

# Loop through the conditions (T1 and C1)
for cond, de_path in de_files.items():

    # Read the DESeq2 results file
    de_full = pd.read_csv(de_path, sep="\t")

    # Keep only the columns we need
    de = de_full[["gene_id", "log2FoldChange", "padj"]].copy()

    # Rename columns so they have the condition name in them
    # (otherwise we'd have name collisions when we merge T1 and C1)
    de = de.rename(columns={
        "log2FoldChange": "log2FoldChange_" + cond,
        "padj": "padj_" + cond,
    })

    # Merge into the genome-wide table
    gene_expr = gene_expr.merge(de, on="gene_id", how="left")

    # Call up/down/ns using our helper function.
    # We use .apply() to run de_call on each row.
    # Note: we use a plain for loop instead of apply for clarity.
    de_calls = []
    log2fc_col = "log2FoldChange_" + cond
    padj_col = "padj_" + cond

    for _, row in gene_expr.iterrows():
        call = de_call(row[log2fc_col], row[padj_col])
        de_calls.append(call)

    # Add the calls as a new column
    gene_expr["DE_" + cond] = de_calls


# Save the genome-wide table to a file
genomewide_path = os.path.join(OUTPUT_FOLDER, "gene_expression_genomewide.csv")
gene_expr.to_csv(genomewide_path, index=False)
print("Saved genome-wide gene expression table.")


# -----------------------------------------------------------------------------
# Step 6: Map genes to bins, separately per condition
# -----------------------------------------------------------------------------
# IMPORTANT: We must do this per condition!
#
# Why? If we mix T1 and C1 bins together, a gene on a chromosome shared
# between both conditions (e.g. chr3) will only get assigned to ONE of them
# (whichever comes first in the bins table). So if T1 rows happen to come
# first, the C1 entries on those chromosomes are essentially ignored.
#
# To avoid this, we map genes to T1 bins, then map genes to C1 bins,
# and combine at the end.

print("")
print("Mapping " + str(len(genes)) + " genes to bins per condition...")

# These are the columns that come back from assign_bin()
bin_label_cols = [
    "transloc_id",
    "condition",
    "behavior",
    "behavior_collapsed",
    "change_direction"
]

# We'll collect the per-condition tables here
genes_bins_per_cond = []

# Loop through each condition
for cond in CONDITIONS:

    # Get only this condition's bins
    bins_cond = bins[bins["condition"] == cond].copy()

    # Make a fresh copy of the genes table to add bin labels to
    genes_cond = genes.copy()

    # ---------------------------------------------------------------------
    # For each gene, find which bin it sits in
    # ---------------------------------------------------------------------
    # We loop through each gene and call assign_bin().
    # The result is a tuple of 5 labels.
    # We collect them all and add them as new columns at the end.

    # Empty lists for each label column
    transloc_ids = []
    conditions_col = []
    behaviors = []
    behaviors_collapsed = []
    change_directions = []

    for _, gene_row in genes_cond.iterrows():
        labels = assign_bin(gene_row, bins_cond)
        # Unpack the tuple into 5 separate values
        tid, c, beh, beh_col, cd = labels
        transloc_ids.append(tid)
        conditions_col.append(c)
        behaviors.append(beh)
        behaviors_collapsed.append(beh_col)
        change_directions.append(cd)

    # Add the lists as new columns
    genes_cond["transloc_id"] = transloc_ids
    genes_cond["condition"] = conditions_col
    genes_cond["behavior"] = behaviors
    genes_cond["behavior_collapsed"] = behaviors_collapsed
    genes_cond["change_direction"] = change_directions

    # Drop genes that didn't overlap any bin (their condition will be NaN)
    genes_cond = genes_cond.dropna(subset=["condition"]).copy()

    print("  " + cond + ": " + str(len(genes_cond))
          + " genes mapped to translocated bins.")

    # ---------------------------------------------------------------------
    # Add expression data to these bin-annotated genes
    # ---------------------------------------------------------------------
    # Merge with the TPM table
    genes_cond = genes_cond.merge(
        expr[["gene_id", "WT", "T1", "C1"]],
        on="gene_id",
        how="left"
    )

    # Calculate log2(TPM + 1) for each condition
    for c in ["WT", "T1", "C1"]:
        log_col = "log2_" + c
        genes_cond[log_col] = np.log2(genes_cond[c] + 1)

    # Calculate log2 fold changes vs WT
    for c in ["T1", "C1"]:
        fc_col = "log2FC_" + c + "_vs_WT"
        log_cond_col = "log2_" + c
        log_wt_col = "log2_WT"
        genes_cond[fc_col] = genes_cond[log_cond_col] - genes_cond[log_wt_col]

    # ---------------------------------------------------------------------
    # Add DESeq2 results and DE calls
    # ---------------------------------------------------------------------
    for c, de_path in de_files.items():

        # Read DESeq2 file (same as before)
        de_full = pd.read_csv(de_path, sep="\t")
        de = de_full[["gene_id", "log2FoldChange", "padj"]].copy()
        de = de.rename(columns={
            "log2FoldChange": "log2FoldChange_" + c,
            "padj": "padj_" + c,
        })

        # Merge into genes_cond
        genes_cond = genes_cond.merge(de, on="gene_id", how="left")

        # Call up/down/ns
        log2fc_col = "log2FoldChange_" + c
        padj_col = "padj_" + c

        de_calls = []
        for _, row in genes_cond.iterrows():
            call = de_call(row[log2fc_col], row[padj_col])
            de_calls.append(call)

        genes_cond["DE_" + c] = de_calls

    # Add this condition's bin-annotated table to our collection
    genes_bins_per_cond.append(genes_cond)


# -----------------------------------------------------------------------------
# Step 7: Combine the per-condition bin-annotated tables and save
# -----------------------------------------------------------------------------
gene_expr_bins = pd.concat(genes_bins_per_cond, ignore_index=True)

# Save the combined table
annotated_path = os.path.join(OUTPUT_FOLDER, "gene_bin_expression_annotation.csv")
gene_expr_bins.to_csv(annotated_path, index=False)
print("Saved annotated gene tables.")

# Print how many genes we have per condition
print("")
print("Gene counts per condition in bin-annotated table:")
print(gene_expr_bins["condition"].value_counts())


# -----------------------------------------------------------------------------
# Step 8: Run the analyses for each condition x classification scheme
# -----------------------------------------------------------------------------
print("")
print("Running analyses...")

# We have 3 classification schemes:
# 1. behavior_collapsed (retained / adopted / other) at A/B level
# 2. change_direction (A_to_B / B_to_A / retained)
# 3. behavior (retained / adopted / other) at full subcompartment level
# We run each one for both T1 and C1.

for cond in CONDITIONS:

    # Get only this condition's bin-annotated genes
    gene_expr_bins_cond = gene_expr_bins[
        gene_expr_bins["condition"] == cond
    ].copy()

    # Get the genome-wide control set: all genes that have a valid DE call
    # for this condition (i.e. not missing).
    de_col = "DE_" + cond
    df_genome_cond = gene_expr[gene_expr[de_col].notna()].copy()

    # ---------------------------------------------------------------------
    # Run the three analyses for this condition
    # ---------------------------------------------------------------------

    # Analysis 1: collapsed A/B
    run_analysis(
        df_bins=gene_expr_bins_cond,
        df_genome=df_genome_cond,
        cond=cond,
        category_col="behavior_collapsed",
        category_order=["retained", "adopted", "other"],
        label="Collapsed A/B compartments"
    )

    # Analysis 2: A/B switching direction
    run_analysis(
        df_bins=gene_expr_bins_cond,
        df_genome=df_genome_cond,
        cond=cond,
        category_col="change_direction",
        category_order=["A_to_B", "B_to_A", "retained"],
        label="A/B direction"
    )

    # Analysis 3: full subcompartment behaviour
    run_analysis(
        df_bins=gene_expr_bins_cond,
        df_genome=df_genome_cond,
        cond=cond,
        category_col="behavior",
        category_order=["retained", "adopted", "other"],
        label="Subcompartment behaviour"
    )


# Done!
print("")
print("Done!")
print("Plots saved to: " + OUTPUT_FOLDER)
