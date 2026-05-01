"""
16_volcano_plots.py
====================
This script makes "volcano plots" of the differential expression results,
with translocated genes highlighted in colour.

What's a volcano plot?
- It's a scatter plot where each dot is one gene.
- X-axis: log2 fold change (how much a gene went up or down between
  conditions). Genes that didn't change much are near the middle.
- Y-axis: -log10(adjusted p-value). Genes with very small p-values
  (statistically significant) appear HIGH up. Genes with large p-values
  (not significant) are near the bottom.
- The shape often looks like an exploding volcano: most dots cluster at
  the bottom-middle, with a few significant genes shooting up and out
  to the sides.

What this script does:
- Highlights the translocated genes in different colours so we can see
  if they cluster among the most differentially expressed genes.
- Exports ranked tables of top up- and down-regulated translocated genes
  for downstream analysis (like GO enrichment in script 18).

For each condition (T1, C1), the script:
1. Loads the DESeq2 results (from script 13).
2. Loads the distance table (from scripts 01-04) to get the list of
   translocated genes.
3. Classifies each gene as Up / Down / No change.
4. Plots all genes as small grey dots.
5. Overlays the translocated genes in colour (green=Up, red=Down, blue=No change).
6. Saves the plot and ranked TSV tables.

Run scripts 01-04 and script 13 first.
Then edit the file paths below and run this script.

Required libraries: pandas, numpy, matplotlib
Install them with: pip install pandas numpy matplotlib
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Increase font sizes for all plots
plt.rcParams.update({
    "font.size": 22,         # default text size
    "axes.titlesize": 30,    # plot title
    "axes.labelsize": 30,    # x and y axis labels
    "xtick.labelsize": 24,   # x tick labels
    "ytick.labelsize": 20,   # y tick labels
    "legend.fontsize": 24,   # legend text
    "legend.title_fontsize": 21  # legend title
})

# =============================================================================
# CONFIG SECTION - Edit these paths and settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# DESeq2 result files (from script 13)
# -----------------------------------------------------------------------------
T1_DE_FILE = "/path/to/DE_WT_vs_T1.tsv"
C1_DE_FILE = "/path/to/DE_WT_vs_C1.tsv"

de_files = {
    "T1": T1_DE_FILE,
    "C1": C1_DE_FILE,
}

# -----------------------------------------------------------------------------
# Distance tables (from scripts 01-04)
# -----------------------------------------------------------------------------
# Used only to figure out which genes are inside translocated regions.
T1_DIST_FILE = "/path/to/T1_distances_agg.tsv"
C1_DIST_FILE = "/path/to/C1_distances_agg.tsv"

dist_files = {
    "T1": T1_DIST_FILE,
    "C1": C1_DIST_FILE,
}

# -----------------------------------------------------------------------------
# Thresholds for calling a gene differentially expressed
# -----------------------------------------------------------------------------
# Same as scripts 14 and 15.
LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

# -----------------------------------------------------------------------------
# How many top genes to export in the ranked TSV tables
# -----------------------------------------------------------------------------
TOP_N_TABLE = 20

# -----------------------------------------------------------------------------
# Colours for translocated genes by DE status
# -----------------------------------------------------------------------------
de_colors = {
    "Up":        "green",
    "Down":      "red",
    "No change": "blue",
}

# -----------------------------------------------------------------------------
# Output folder for plots and tables
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/plots/volcano_plots"

# Conditions we'll process
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Classify each gene as Up / Down / No change
# =============================================================================

def classify_de(df):
    """
    Add a DE_status column (Up/Down/No change) and a -log10(padj) column
    to a DESeq2 results DataFrame.
    """
    # Make a copy so we don't change the original DataFrame
    df = df.copy()

    # -------------------------------------------------------------------------
    # Handle missing or zero p-values before taking log
    # -------------------------------------------------------------------------
    # If padj is missing (NaN), set it to 1 (= "definitely not significant").
    # If padj is exactly 0, replace with a tiny number to avoid log(0)
    # which would give negative infinity.
    df["padj"] = df["padj"].fillna(1)
    df["padj"] = df["padj"].replace(0, 1e-300)

    # Calculate -log10(padj) for the y-axis of the volcano plot.
    # Smaller padj -> bigger -log10(padj) -> dot appears HIGHER in the plot.
    df["-log10(padj)"] = -np.log10(df["padj"])

    # -------------------------------------------------------------------------
    # Classify each gene as Up / Down / No change
    # -------------------------------------------------------------------------
    de_statuses = []

    for _, row in df.iterrows():
        log2fc = row["log2FoldChange"]
        padj = row["padj"]

        # If padj isn't significant, call it "No change"
        if padj >= PADJ_THRESHOLD:
            de_statuses.append("No change")
            continue

        # padj is significant, check the fold change direction
        if log2fc > LFC_THRESHOLD:
            de_statuses.append("Up")
        elif log2fc < -LFC_THRESHOLD:
            de_statuses.append("Down")
        else:
            # padj is significant but the fold change is too small
            de_statuses.append("No change")

    df["DE_status"] = de_statuses

    return df


# =============================================================================
# Helper function: Load translocated gene names from a distance table
# =============================================================================

def load_transloc_genes(dist_file):
    """
    Read a distance table and return a Python set of gene names that are
    inside translocated regions.
    """
    # Read the file
    dist_df = pd.read_csv(dist_file, sep="\t")

    # Keep only rows where gene1 is inside a translocation
    inside_rows = dist_df[dist_df["gene1_region_type"] == "inside_transloc"]

    # Get the unique gene names and put them in a set for fast lookup
    unique_names = inside_rows["gene1"].unique()
    transloc_set = set(unique_names)

    return transloc_set


# =============================================================================
# Helper function: Make the volcano plot for one condition
# =============================================================================

def plot_volcano(de_df, transloc_genes, cond):
    """
    Make the volcano plot for one condition (T1 or C1) and export
    the ranked gene tables.
    """
    # -------------------------------------------------------------------------
    # Step 1: Mark which genes are translocated
    # -------------------------------------------------------------------------
    de_df = de_df.copy()

    # .isin(some_set) returns True/False for each row
    de_df["is_translocated"] = de_df["gene_name"].isin(transloc_genes)

    # Split into translocated and non-translocated
    trans_df = de_df[de_df["is_translocated"]]

    # The ~ (tilde) inverts a boolean, so this gets rows where
    # is_translocated is False
    bg_df = de_df[~de_df["is_translocated"]]

    # -------------------------------------------------------------------------
    # Step 2: Create the figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8))

    # -------------------------------------------------------------------------
    # Step 3: Plot all the non-translocated genes as gray dots
    # -------------------------------------------------------------------------
    ax.scatter(
        bg_df["log2FoldChange"],
        bg_df["-log10(padj)"],
        c="gray",
        alpha=0.4,
        s=8,
        label="All other genes"
    )

    # -------------------------------------------------------------------------
    # Step 4: Overlay translocated genes coloured by DE status
    # -------------------------------------------------------------------------
    for status, color in de_colors.items():

        # Get just the rows with this DE status
        subset = trans_df[trans_df["DE_status"] == status]

        # Plot them. We use bigger dots (s=30) and add a black outline
        # (edgecolor) so they stand out from the background gray dots.
        ax.scatter(
            subset["log2FoldChange"],
            subset["-log10(padj)"],
            c=color,
            alpha=0.9,
            s=30,
            edgecolor="black",
            linewidth=0.3,
            label="Translocated — " + status
        )

    # -------------------------------------------------------------------------
    # Step 5: Add dashed lines showing the DE thresholds
    # -------------------------------------------------------------------------
    # axvline = vertical line at given x value
    # axhline = horizontal line at given y value
    ax.axvline(x=LFC_THRESHOLD, color="black", linestyle="--", lw=1)
    ax.axvline(x=-LFC_THRESHOLD, color="black", linestyle="--", lw=1)
    ax.axhline(y=-np.log10(PADJ_THRESHOLD), color="black", linestyle="--", lw=1)

    # -------------------------------------------------------------------------
    # Step 6: Add labels and title
    # -------------------------------------------------------------------------
    # The $...$ syntax in matplotlib labels uses LaTeX math formatting.
    # This is how we get nice subscripts like log_2 and log_10.
    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$(adjusted p-value)")
    ax.set_title("WT vs " + cond + " — translocated genes highlighted")

    # frameon=False removes the box around the legend
    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    # -------------------------------------------------------------------------
    # Step 7: Save the plot
    # -------------------------------------------------------------------------
    fig.tight_layout()

    plot_filename = "volcano_WT_vs_" + cond + "_highlighted.png"
    plot_path = os.path.join(OUTPUT_FOLDER, plot_filename)
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")

    plt.close(fig)

    print("  Saved plot: " + plot_path)

    # -------------------------------------------------------------------------
    # Step 8: Find the significant translocated genes (for export)
    # -------------------------------------------------------------------------
    is_up = (trans_df["DE_status"] == "Up")
    is_down = (trans_df["DE_status"] == "Down")
    trans_sig = trans_df[is_up | is_down]   # | is OR

    # -------------------------------------------------------------------------
    # Step 9: Export ranked top-N TSV tables
    # -------------------------------------------------------------------------
    # We rank by FOLD CHANGE (not p-value) to get the biggest movers up
    # and down.
    # nlargest returns the rows with the biggest log2FC (most upregulated).
    # nsmallest returns the rows with the smallest log2FC (most downregulated).

    only_up = trans_sig[trans_sig["DE_status"] == "Up"]
    top_up_table = only_up.nlargest(TOP_N_TABLE, "log2FoldChange")

    only_down = trans_sig[trans_sig["DE_status"] == "Down"]
    top_down_table = only_down.nsmallest(TOP_N_TABLE, "log2FoldChange")

    # Build paths and save
    up_filename = "top_up_transloc_genes_WT_vs_" + cond + ".tsv"
    down_filename = "top_down_transloc_genes_WT_vs_" + cond + ".tsv"

    up_path = os.path.join(OUTPUT_FOLDER, up_filename)
    down_path = os.path.join(OUTPUT_FOLDER, down_filename)

    top_up_table.to_csv(up_path, sep="\t", index=False)
    top_down_table.to_csv(down_path, sep="\t", index=False)

    print("")
    print("  Exported top " + str(TOP_N_TABLE) + " upregulated genes:   " + up_path)
    print("  Exported top " + str(TOP_N_TABLE) + " downregulated genes: " + down_path)

    # -------------------------------------------------------------------------
    # Step 10: Export ALL significant translocated genes (no cap)
    # -------------------------------------------------------------------------
    # These full lists are used by script 18 for GO enrichment.
    # GO enrichment needs many genes to be statistically powerful, so we
    # don't limit to top-N here.

    # Sort all upregulated genes by log2FC (largest first)
    all_up_table = trans_sig[trans_sig["DE_status"] == "Up"]
    all_up_table = all_up_table.sort_values("log2FoldChange", ascending=False)

    # Sort all downregulated genes by log2FC (most negative first)
    all_down_table = trans_sig[trans_sig["DE_status"] == "Down"]
    all_down_table = all_down_table.sort_values("log2FoldChange", ascending=True)

    # Build paths and save
    all_up_filename = "all_up_transloc_genes_WT_vs_" + cond + ".tsv"
    all_down_filename = "all_down_transloc_genes_WT_vs_" + cond + ".tsv"

    all_up_path = os.path.join(OUTPUT_FOLDER, all_up_filename)
    all_down_path = os.path.join(OUTPUT_FOLDER, all_down_filename)

    all_up_table.to_csv(all_up_path, sep="\t", index=False)
    all_down_table.to_csv(all_down_path, sep="\t", index=False)

    print("  Exported all " + str(len(all_up_table))
          + " upregulated genes:   " + all_up_path)
    print("  Exported all " + str(len(all_down_table))
          + " downregulated genes: " + all_down_path)


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Loop through each condition (T1 and C1)
for cond in CONDITIONS:

    print("")
    print("===== Processing " + cond + " =====")

    # -------------------------------------------------------------------------
    # Step 1: Load and classify the DESeq2 results
    # -------------------------------------------------------------------------
    de_df = pd.read_csv(de_files[cond], sep="\t")
    de_df = classify_de(de_df)

    # -------------------------------------------------------------------------
    # Step 2: Load the translocated gene names
    # -------------------------------------------------------------------------
    transloc_genes = load_transloc_genes(dist_files[cond])

    print("  Translocated genes identified: " + str(len(transloc_genes)))

    # Count how many translocated genes appear in the DE table.
    # .isin returns True/False per row, .sum counts the Trues.
    n_in_de_table = de_df["gene_name"].isin(transloc_genes).sum()
    print("  Translocated genes in DE table: " + str(n_in_de_table))

    # -------------------------------------------------------------------------
    # Step 3: Make the volcano plot and export tables
    # -------------------------------------------------------------------------
    plot_volcano(de_df, transloc_genes, cond)


# Done!
print("")
print("Done!")
