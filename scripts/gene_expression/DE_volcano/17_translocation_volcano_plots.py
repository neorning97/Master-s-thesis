"""
17_translocation_volcano_plots.py
===================================
This script makes volcano plots highlighting genes from SPECIFIC named
translocations (like Der(17)t(3;17)), and runs statistical tests comparing
those genes against the rest of the genome.

How is this different from script 16?
- Script 16 highlights translocated genes that came from the distance pipeline.
- This script lets you define translocation regions DIRECTLY by their
  genomic coordinates (chromosome + start + end), and finds all genes inside.
- This is useful for analysing each named translocation separately
  (e.g. only Der(17)t(3;17), only Der(3)t(3;17)).

What's a volcano plot?
- See script 16's docstring for a full explanation.
- Quick recap: each dot = one gene. X-axis = log2 fold change.
  Y-axis = -log10(padj). Significant changes shoot up and out to the sides.

What this script does for each translocation x condition:
1. Finds all genes whose coordinates overlap the translocation regions.
2. Loads the DESeq2 results and labels each gene as Up / Down / No change.
3. Makes a volcano plot with translocation genes highlighted in colour.
4. Runs two Fisher's exact tests:
   a. "DE enrichment": are translocation genes more often DE than other genes?
   b. "Up/Down bias": among DE genes, do translocation genes lean up vs down?
5. Exports ranked TSV tables of the top up- and down-regulated genes.

Run script 13 first to make the DESeq2 result files.
Then edit the file paths and translocation regions below.

Required libraries: pandas, numpy, matplotlib, scipy
Install them with: pip install pandas numpy matplotlib scipy
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

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
# GTF file with gene coordinates
# -----------------------------------------------------------------------------
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# -----------------------------------------------------------------------------
# Translocation regions to analyse
# -----------------------------------------------------------------------------
# This is a dictionary where:
# - The KEY is a readable name for the translocation (used in plot titles
#   and filenames).
# - The VALUE is a list of (chromosome, start, end) tuples, one tuple
#   per genomic segment that's part of this translocation.
#
# Each translocation involves TWO chromosome segments (a fragment + a body),
# so each list has 2 entries.
#
# The underscore in numbers like 58_600_000 is just for readability
# (Python ignores it). It's like writing 58,600,000.

translocations = {
    "Der(17)t(3;17)": [
        ("chr3",  0,          58_600_000),    # the chr3 fragment
        ("chr17", 22_700_001, 83_257_441),    # the chr17 body
    ],
    "Der(3)t(3;17)": [
        ("chr3",  58_600_001, 198_295_559),   # the chr3 body
        ("chr17", 0,          22_700_000),    # the chr17 fragment
    ],
}

# -----------------------------------------------------------------------------
# Thresholds for calling a gene differentially expressed
# -----------------------------------------------------------------------------
# Same as in scripts 14, 15, 16.
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
# Output folder
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/plots/translocation_volcano_plots"

# Conditions we'll process
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Load gene coordinates from a GTF file
# =============================================================================

def load_gtf_genes(gtf_file):
    """
    Load a GTF file and return a DataFrame with one row per gene.
    Columns: chr, start, end, gene_name
    """
    column_names = [
        "chr", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes"
    ]

    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=column_names
    )

    # Keep only rows where feature is "gene"
    gtf = gtf[gtf["feature"] == "gene"].copy()

    # Extract gene_name from the attributes column.
    gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')

    # Keep only the columns we need
    keep_cols = ["chr", "start", "end", "gene_name"]
    genes = gtf[keep_cols].copy()

    return genes


# =============================================================================
# Helper function: Find all gene names that overlap a list of regions
# =============================================================================

def get_transloc_genes(gtf, regions):
    """
    Find all gene names whose coordinates overlap any of the given
    (chromosome, start, end) regions. Returns a Python set of names.
    """
    genes = set()

    for chrom, region_start, region_end in regions:

        # A gene overlaps if it's on the same chromosome AND
        # the gene's end is at/after the region's start AND
        # the gene's start is at/before the region's end.
        same_chrom = (gtf["chr"] == chrom)
        ends_after_region_starts = (gtf["end"] >= region_start)
        starts_before_region_ends = (gtf["start"] <= region_end)

        mask = same_chrom & ends_after_region_starts & starts_before_region_ends

        # Get the matching gene names and add to set
        matching_names = gtf[mask]["gene_name"]
        genes.update(matching_names)

    return genes


# =============================================================================
# Helper function: Classify each gene as Up / Down / No change
# =============================================================================

def classify_de(df):
    """
    Add DE_status (Up/Down/No change) and -log10(padj) columns to a
    DESeq2 results DataFrame.
    """
    df = df.copy()

    # Handle missing or zero p-values before taking log
    df["padj"] = df["padj"].fillna(1)
    df["padj"] = df["padj"].replace(0, 1e-300)

    # Calculate -log10(padj) for the y-axis of the volcano plot
    df["-log10(padj)"] = -np.log10(df["padj"])

    de_statuses = []

    for _, row in df.iterrows():
        log2fc = row["log2FoldChange"]
        padj = row["padj"]

        # Not significant if padj is too high
        if padj >= PADJ_THRESHOLD:
            de_statuses.append("No change")
            continue

        # padj is significant, check fold change direction
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
# Helper function: Run two Fisher's exact tests
# =============================================================================

def run_statistics(df, transloc_label, cond):
    """
    Run two Fisher's exact tests comparing translocation genes vs the
    rest of the genome:
    1. DE vs No change (are translocation genes more often DE?)
    2. Up vs Down (among DE genes, are translocation genes biased toward Up?)
    Print the results.
    """
    trans_df = df[df["is_translocated"]]
    genome_df = df[~df["is_translocated"]]

    n_trans = len(trans_df)

    if n_trans == 0:
        print("  No translocation genes found for "
              + transloc_label + " in " + cond)
        return

    # -------------------------------------------------------------------------
    # Print DE fractions in the translocation group
    # -------------------------------------------------------------------------
    print("")
    print("  " + transloc_label + " | WT vs " + cond)
    print("  Translocation genes: " + str(n_trans))

    for status in ["Up", "Down", "No change"]:
        n = (trans_df["DE_status"] == status).sum()
        fraction = n / n_trans
        print("    " + status + ": " + str(n)
              + " ({:.1%})".format(fraction))

    # -------------------------------------------------------------------------
    # Test 1: DE enrichment (DE vs No change)
    # -------------------------------------------------------------------------
    trans_de = (trans_df["DE_status"] != "No change").sum()
    trans_no = (trans_df["DE_status"] == "No change").sum()
    genome_de = (genome_df["DE_status"] != "No change").sum()
    genome_no = (genome_df["DE_status"] == "No change").sum()

    de_table = [
        [trans_de, trans_no],
        [genome_de, genome_no]
    ]

    or_de, p_de = fisher_exact(de_table)

    print("")
    print("  DE enrichment Fisher test:")
    print("    Translocation: DE=" + str(trans_de)
          + ", No change=" + str(trans_no))
    print("    Genome-wide:   DE=" + str(genome_de)
          + ", No change=" + str(genome_no))
    print("    OR={:.3f},  p={:.4g}".format(or_de, p_de))

    # -------------------------------------------------------------------------
    # Test 2: Up/Down bias (among DE genes only)
    # -------------------------------------------------------------------------
    trans_up = (trans_df["DE_status"] == "Up").sum()
    trans_down = (trans_df["DE_status"] == "Down").sum()
    genome_up = (genome_df["DE_status"] == "Up").sum()
    genome_down = (genome_df["DE_status"] == "Down").sum()

    ud_table = [
        [trans_up, trans_down],
        [genome_up, genome_down]
    ]

    or_ud, p_ud = fisher_exact(ud_table)

    print("")
    print("  Up/Down bias Fisher test (among DE genes only):")
    print("    Translocation: Up=" + str(trans_up)
          + ", Down=" + str(trans_down))
    print("    Genome-wide:   Up=" + str(genome_up)
          + ", Down=" + str(genome_down))
    print("    OR={:.3f},  p={:.4g}".format(or_ud, p_ud))


# =============================================================================
# Helper function: Make the volcano plot for one translocation x condition
# =============================================================================

def plot_volcano(df, transloc_genes, transloc_label, cond):
    """
    Make a volcano plot for one translocation in one condition.
    Also runs the statistical tests and exports TSV tables.
    """
    # -------------------------------------------------------------------------
    # Step 1: Mark which genes are in this translocation
    # -------------------------------------------------------------------------
    df = df.copy()
    df["is_translocated"] = df["gene_name"].isin(transloc_genes)

    trans_df = df[df["is_translocated"]]
    bg_df = df[~df["is_translocated"]]

    # -------------------------------------------------------------------------
    # Step 2: Create the figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8))

    # -------------------------------------------------------------------------
    # Step 3: Draw the background (non-translocation) genes as gray dots
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
    # Step 4: Overlay translocation genes coloured by DE status
    # -------------------------------------------------------------------------
    for status, color in de_colors.items():

        subset = trans_df[trans_df["DE_status"] == status]

        # Only plot if there are any genes in this category
        if len(subset) > 0:
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
    # Step 5: Add dashed threshold lines
    # -------------------------------------------------------------------------
    ax.axvline(x=LFC_THRESHOLD, color="black", linestyle="--", lw=1)
    ax.axvline(x=-LFC_THRESHOLD, color="black", linestyle="--", lw=1)
    ax.axhline(y=-np.log10(PADJ_THRESHOLD), color="black", linestyle="--", lw=1)

    # -------------------------------------------------------------------------
    # Step 6: Labels and title
    # -------------------------------------------------------------------------
    # The $...$ syntax uses LaTeX math formatting for nice subscripts
    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$(adjusted p-value)")
    ax.set_title("WT vs " + cond + " — " + transloc_label
                 + " genes highlighted")

    # Move legend outside the plot to the right
    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    # -------------------------------------------------------------------------
    # Step 7: Save the plot
    # -------------------------------------------------------------------------
    fig.tight_layout()

    # Make a filename-safe version of the translocation label.
    safe_label = transloc_label.replace("(", "")
    safe_label = safe_label.replace(")", "")
    safe_label = safe_label.replace(";", "_")

    plot_filename = "volcano_WT_vs_" + cond + "_" + safe_label + ".png"
    plot_path = os.path.join(OUTPUT_FOLDER, plot_filename)
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("  Saved plot: " + plot_path)

    # -------------------------------------------------------------------------
    # Step 8: Find the significant translocation genes (for export)
    # -------------------------------------------------------------------------
    is_up = (trans_df["DE_status"] == "Up")
    is_down = (trans_df["DE_status"] == "Down")
    trans_sig = trans_df[is_up | is_down]   # | is OR

    # -------------------------------------------------------------------------
    # Step 9: Export ranked TSV tables
    # -------------------------------------------------------------------------
    # We rank by FOLD CHANGE (not p-value).
    # nlargest gives us the most upregulated genes.
    # nsmallest gives us the most downregulated (smallest = most negative).

    only_up = trans_sig[trans_sig["DE_status"] == "Up"]
    top_up_table = only_up.nlargest(TOP_N_TABLE, "log2FoldChange")

    only_down = trans_sig[trans_sig["DE_status"] == "Down"]
    top_down_table = only_down.nsmallest(TOP_N_TABLE, "log2FoldChange")

    # Build the file paths
    up_filename = "top_up_" + safe_label + "_WT_vs_" + cond + ".tsv"
    down_filename = "top_down_" + safe_label + "_WT_vs_" + cond + ".tsv"

    up_path = os.path.join(OUTPUT_FOLDER, up_filename)
    down_path = os.path.join(OUTPUT_FOLDER, down_filename)

    # Save the tables
    top_up_table.to_csv(up_path, sep="\t", index=False)
    top_down_table.to_csv(down_path, sep="\t", index=False)

    print("")
    print("  Exported top " + str(TOP_N_TABLE)
          + " upregulated genes:   " + up_path)
    print("  Exported top " + str(TOP_N_TABLE)
          + " downregulated genes: " + down_path)

    # -------------------------------------------------------------------------
    # Step 10: Run the statistical tests
    # -------------------------------------------------------------------------
    run_statistics(df, transloc_label, cond)


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)


# -----------------------------------------------------------------------------
# Step 1: Load the GTF file (only needs to be done once)
# -----------------------------------------------------------------------------
print("Loading GTF gene annotations...")
gtf = load_gtf_genes(GTF_FILE)
print("  " + str(len(gtf)) + " protein-coding genes loaded.")


# -----------------------------------------------------------------------------
# Step 2: Load and classify the DESeq2 results for each condition
# -----------------------------------------------------------------------------
print("Loading DESeq2 results...")

de_results = {}

for cond in CONDITIONS:

    de_path = de_files[cond]
    de_df = pd.read_csv(de_path, sep="\t")

    # Add the DE_status and -log10(padj) columns
    de_df = classify_de(de_df)

    de_results[cond] = de_df

    print("  " + cond + ": " + str(len(de_df)) + " genes loaded.")


# -----------------------------------------------------------------------------
# Step 3: Process each translocation x condition combination
# -----------------------------------------------------------------------------
for transloc_label, regions in translocations.items():

    print("")
    print("=" * 55)
    print("Translocation: " + transloc_label)
    print("=" * 55)

    # Find all genes whose coordinates overlap this translocation's regions
    transloc_genes = get_transloc_genes(gtf, regions)
    print("  Genes in translocation region: " + str(len(transloc_genes)))

    # For each condition, make a plot
    for cond in CONDITIONS:
        plot_volcano(
            df=de_results[cond],
            transloc_genes=transloc_genes,
            transloc_label=transloc_label,
            cond=cond
        )

# Done!
print("")
print("Done!")
