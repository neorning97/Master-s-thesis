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
# How many top genes to label and export
# -----------------------------------------------------------------------------
TOP_N_LABELS = 5    # how many to label on each plot
TOP_N_TABLE = 20    # how many to export per TSV file

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
OUTPUT_FOLDER = "/path/to/results/plots/translocation_volcano_plots"

# Conditions we'll process
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Load gene coordinates from a GTF file
# =============================================================================
# Same as in script 14, read the GTF, keep only "gene" rows, pull the
# gene_name out of the attributes column.

def load_gtf_genes(gtf_file):
    """
    Load a GTF file and return a DataFrame with one row per gene.
    Columns: chr, start, end, gene_name
    """
    # GTF files have 9 columns, we provide names ourselves.
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
    # str.extract with a regex grabs the part inside the parentheses.
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
    # Empty set to start, we'll add gene names to it as we find them
    genes = set()

    # Loop through each region tuple
    for chrom, region_start, region_end in regions:

        # Same overlap idea as in earlier scripts:
        # A gene overlaps if it's on the same chromosome AND
        # the gene's end is at/after the region's start AND
        # the gene's start is at/before the region's end.
        same_chrom = (gtf["chr"] == chrom)
        ends_after_region_starts = (gtf["end"] >= region_start)
        starts_before_region_ends = (gtf["start"] <= region_end)

        # Combine with & (AND)
        mask = same_chrom & ends_after_region_starts & starts_before_region_ends

        # Get the matching gene names
        matching_names = gtf[mask]["gene_name"]

        # Add them to our set.
        # set.update() adds all values from a collection.
        # Since it's a set, duplicates are automatically ignored.
        genes.update(matching_names)

    return genes


# =============================================================================
# Helper function: Classify each gene as Up / Down / No change
# =============================================================================
# Same as in script 16.

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

    # Classify each gene with an explicit for loop (instead of np.select)
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
    # Split into translocation and genome-wide groups.
    # ~ inverts a boolean (so ~True = False).
    trans_df = df[df["is_translocated"]]
    genome_df = df[~df["is_translocated"]]

    n_trans = len(trans_df)
    n_genome = len(genome_df)

    # If there are no translocation genes, there's nothing to test
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
        # Count rows with this DE status
        n = (trans_df["DE_status"] == status).sum()

        # {:.1%} formats a number as a percentage with 1 decimal
        # (e.g. 0.156 becomes "15.6%")
        fraction = n / n_trans
        print("    " + status + ": " + str(n)
              + " ({:.1%})".format(fraction))

    # -------------------------------------------------------------------------
    # Test 1: DE enrichment (DE vs No change)
    # -------------------------------------------------------------------------
    # Build the 2x2 table:
    #                       DE      Not DE
    # Translocation         a         b
    # Genome-wide           c         d

    # Count Trues with .sum() (True counts as 1, False as 0)
    trans_de = (trans_df["DE_status"] != "No change").sum()
    trans_no = (trans_df["DE_status"] == "No change").sum()
    genome_de = (genome_df["DE_status"] != "No change").sum()
    genome_no = (genome_df["DE_status"] == "No change").sum()

    # Build the table as a list of lists
    de_table = [
        [trans_de, trans_no],
        [genome_de, genome_no]
    ]

    # Run Fisher's exact test
    or_de, p_de = fisher_exact(de_table)

    # Print results
    print("")
    print("  DE enrichment Fisher test:")
    print("    Translocation: DE=" + str(trans_de)
          + ", No change=" + str(trans_no))
    print("    Genome-wide:   DE=" + str(genome_de)
          + ", No change=" + str(genome_no))
    # {:.3f} = 3 decimal places, {:.4g} = up to 4 significant digits
    print("    OR={:.3f},  p={:.4g}".format(or_de, p_de))

    # -------------------------------------------------------------------------
    # Test 2: Up/Down bias (among DE genes only)
    # -------------------------------------------------------------------------
    # Build the 2x2 table:
    #                       Up      Down
    # Translocation         a         b
    # Genome-wide           c         d

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

    # .isin returns True/False per row
    df["is_translocated"] = df["gene_name"].isin(transloc_genes)

    # Split into translocation and background groups
    trans_df = df[df["is_translocated"]]
    bg_df = df[~df["is_translocated"]]

    # -------------------------------------------------------------------------
    # Step 2: Create the figure
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    # -------------------------------------------------------------------------
    # Step 3: Draw the background (non-translocation) genes as light grey dots
    # -------------------------------------------------------------------------
    ax.scatter(
        bg_df["log2FoldChange"],
        bg_df["-log10(padj)"],
        c="lightgray",
        alpha=0.4,
        s=8,
        label="All other genes"
    )

    # -------------------------------------------------------------------------
    # Step 4: Overlay translocation genes coloured by DE status
    # -------------------------------------------------------------------------
    # We loop through the colour dict so each status gets its own legend entry
    for status, color in de_colors.items():

        # Get just the rows with this DE status
        subset = trans_df[trans_df["DE_status"] == status]

        # Only plot if there are any genes in this category
        # (otherwise we'd add an empty legend entry)
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
    ax.legend(frameon=False, fontsize=8)

    # -------------------------------------------------------------------------
    # Step 7: Find the top significant translocation genes to label
    # -------------------------------------------------------------------------
    # Get only the genes that are differentially expressed
    is_up = (trans_df["DE_status"] == "Up")
    is_down = (trans_df["DE_status"] == "Down")
    trans_sig = trans_df[is_up | is_down]   # | is OR

    # nlargest(N, "col") returns N rows with the largest values in col
    top_up = trans_sig[trans_sig["DE_status"] == "Up"]
    top_up = top_up.nlargest(TOP_N_LABELS, "-log10(padj)")

    top_down = trans_sig[trans_sig["DE_status"] == "Down"]
    top_down = top_down.nlargest(TOP_N_LABELS, "-log10(padj)")

    # Stick them together for labelling
    top_genes = pd.concat([top_up, top_down])

    # -------------------------------------------------------------------------
    # Step 8: Add the gene name labels on the plot
    # -------------------------------------------------------------------------
    # We loop through each top gene and add a text label.
    # i is the row index (0, 1, 2, ...), we use it to stagger labels vertically.
    for i, (_, row) in enumerate(top_genes.iterrows()):

        x_pos = row["log2FoldChange"]
        y_pos = row["-log10(padj)"]
        gene_name = row["gene_name"]

        # np.sign returns +1 if positive, -1 if negative.
        # x_offset is +0.05 for upregulated (label to the right)
        # and -0.05 for downregulated (label to the left).
        x_offset = 0.05 * np.sign(x_pos)

        # y_offset increases with i so labels stagger upward, not on top of each other
        y_offset = 0.1 + 0.05 * i

        # Choose horizontal alignment based on direction
        if x_pos > 0:
            horizontal_alignment = "left"
        else:
            horizontal_alignment = "right"

        # Add the label
        ax.text(
            x_pos + x_offset,
            y_pos + y_offset,
            gene_name,
            fontsize=6,
            rotation=30,
            ha=horizontal_alignment,
            va="bottom"
        )

    # -------------------------------------------------------------------------
    # Step 9: Save the plot
    # -------------------------------------------------------------------------
    fig.tight_layout()

    # Make a filename-safe version of the translocation label.
    # We can't have parentheses or semicolons in filenames on some systems.
    safe_label = transloc_label.replace("(", "")
    safe_label = safe_label.replace(")", "")
    safe_label = safe_label.replace(";", "_")

    plot_filename = "volcano_WT_vs_" + cond + "_" + safe_label + ".png"
    plot_path = os.path.join(OUTPUT_FOLDER, plot_filename)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    print("  Saved plot: " + plot_path)

    # -------------------------------------------------------------------------
    # Step 10: Print top genes to the terminal
    # -------------------------------------------------------------------------
    print("")
    print("  Top " + str(TOP_N_LABELS) + " upregulated genes ("
          + transloc_label + ", " + cond + "):")
    for gene in top_up["gene_name"].values:
        print("    " + gene)

    print("")
    print("  Top " + str(TOP_N_LABELS) + " downregulated genes ("
          + transloc_label + ", " + cond + "):")
    for gene in top_down["gene_name"].values:
        print("    " + gene)

    # -------------------------------------------------------------------------
    # Step 11: Export ranked TSV tables
    # -------------------------------------------------------------------------
    # Note: here we rank by FOLD CHANGE (not p-value).
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
    # Step 12: Run the statistical tests
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

# Empty dictionary, we'll fill it with one DataFrame per condition
de_results = {}

for cond in CONDITIONS:

    # Read the file
    de_path = de_files[cond]
    de_df = pd.read_csv(de_path, sep="\t")

    # Add the DE_status and -log10(padj) columns
    de_df = classify_de(de_df)

    # Store in our dict
    de_results[cond] = de_df

    print("  " + cond + ": " + str(len(de_df)) + " genes loaded.")


# -----------------------------------------------------------------------------
# Step 3: Process each translocation x condition combination
# -----------------------------------------------------------------------------
# Loop through each translocation
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
