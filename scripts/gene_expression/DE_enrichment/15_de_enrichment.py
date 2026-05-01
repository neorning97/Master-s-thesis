"""
15_de_enrichment.py
====================
This script tests whether translocated genes are MORE LIKELY to be
differentially expressed than other genes in the genome.

The big questions:
1. Are translocated genes more often differentially expressed (Up or Down)
   than genes elsewhere in the genome?
2. Among the genes that ARE differentially expressed, do translocated genes
   tend to go UP more often than DOWN compared to the rest of the genome?

How we answer this:
- Split all genes into two groups:
  - "translocation genes": genes inside translocated regions
  - "genome background": all other genes
- Compare the two groups using Fisher's exact test.

What's Fisher's exact test?
- It's a statistical test for 2x2 tables.
- A 2x2 table here looks like:
                       DE       Not DE
  Translocation        a         b
  Genome background    c         d
- The test asks: "Are the proportions in the two rows significantly
  different?"
- We use Fisher's exact (instead of chi-square) because it works well
  even when some of the cells have small counts. Translocation genes
  are usually a small group (~100-500 genes), so this safer test is
  the right choice.

Why two separate tests?
- Test 1 ("DE vs No change"): Are translocated genes more affected overall?
- Test 2 ("Up vs Down"): Among DE genes, are translocated genes biased
  toward going UP rather than DOWN?
- These are different biological questions, so we test them separately.

Run scripts 01-04 and script 13 first.
Then edit the file paths below and run this script.

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
    "legend.fontsize": 22,   # legend text
    "legend.title_fontsize": 21  # legend title
})

# =============================================================================
# CONFIG SECTION - Edit these paths and settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# Distance table from scripts 01-04
# -----------------------------------------------------------------------------
# We only use this file to get the LIST of translocated gene names.
# We don't actually use the distance values here.
DIST_FILE = "/path/to/WT_distances_agg.tsv"

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
# Thresholds for calling a gene differentially expressed
# -----------------------------------------------------------------------------
# Same idea as in script 14: we need BOTH a small adjusted p-value AND
# a big enough fold change to call a gene "Up" or "Down".
LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

# -----------------------------------------------------------------------------
# Output folder for plots
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/plots/DE_enrichment"

# Conditions we'll process
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Classify each gene as Up / Down / No change
# =============================================================================
# This adds a "DE_status" column to a DESeq2 results DataFrame.
# A gene is "Up" or "Down" only if BOTH thresholds are met.

def classify_de(df):
    """
    Add a DE_status column to a DESeq2 results DataFrame.
    Each gene becomes "Up", "Down", or "No change".
    """
    # We'll build the DE_status as a list, one entry per row,
    # then add it as a new column at the end.
    de_statuses = []

    # Loop through each row
    for _, row in df.iterrows():

        log2fc = row["log2FoldChange"]
        padj = row["padj"]

        # If either value is missing (NaN), call it "No change"
        if pd.isna(log2fc) or pd.isna(padj):
            de_statuses.append("No change")
            continue

        # Check if padj is small enough
        if padj >= PADJ_THRESHOLD:
            de_statuses.append("No change")
            continue

        # padj is significant, now check the fold change direction
        if log2fc > LFC_THRESHOLD:
            de_statuses.append("Up")
        elif log2fc < -LFC_THRESHOLD:
            de_statuses.append("Down")
        else:
            # padj is significant but the fold change is too small
            de_statuses.append("No change")

    # Add the list as a new column
    df["DE_status"] = de_statuses

    return df


# =============================================================================
# Helper function: Count Up / Down / No change in a set of genes
# =============================================================================

def get_de_counts(gene_set, de_df):
    """
    For the genes in gene_set, count how many are Up, Down, and No change.
    Returns a pandas Series with three numbers (one per category).
    """
    # Get only rows where the gene name is in our set.
    # .isin() returns True/False for each row.
    subset = de_df[de_df["gene_name"].isin(gene_set)]

    # value_counts() counts how many times each value appears
    counts = subset["DE_status"].value_counts()

    # reindex makes sure we always have all 3 categories in the same order.
    # If a category is missing, fill_value=0 puts a 0 there.
    counts = counts.reindex(["Up", "Down", "No change"], fill_value=0)

    return counts


# =============================================================================
# Helper function: Run a Fisher's exact test
# =============================================================================
# We use this for both:
# - DE vs No change (test 1)
# - Up vs Down      (test 2)

def fisher_test(counts_a, counts_b, test_type):
    """
    Run a two-sided Fisher's exact test comparing two groups.

    counts_a, counts_b: pandas Series with Up / Down / No change counts.
    test_type: either "DE_vs_NoChange" or "Up_vs_Down".

    Returns a tuple: (the 2x2 table as a list, odds ratio, p-value).
    """
    # Build the 2x2 table based on which test we're running
    if test_type == "DE_vs_NoChange":

        # We combine Up + Down into "DE", and put "No change" in the other column
        # int() converts to a regular Python integer (instead of numpy int)
        a = int(counts_a["Up"] + counts_a["Down"])    # group A, DE
        b = int(counts_a["No change"])                # group A, not DE
        c = int(counts_b["Up"] + counts_b["Down"])    # group B, DE
        d = int(counts_b["No change"])                # group B, not DE

    elif test_type == "Up_vs_Down":

        # Compare Up to Down (ignoring No change)
        a = int(counts_a["Up"])     # group A, Up
        b = int(counts_a["Down"])   # group A, Down
        c = int(counts_b["Up"])     # group B, Up
        d = int(counts_b["Down"])   # group B, Down

    else:
        # Defensive check, raise an error if test_type is something unexpected.
        # raise stops the script and shows the error message.
        raise ValueError("Unknown test_type '" + test_type
                         + "'. Use 'DE_vs_NoChange' or 'Up_vs_Down'.")

    # Build the 2x2 contingency table as a list of lists
    table = [[a, b], [c, d]]

    # Run Fisher's exact test.
    # It returns the odds ratio and the p-value.
    oddsratio, p = fisher_exact(table)

    return table, oddsratio, p


# =============================================================================
# Helper function: Make a grouped bar plot
# =============================================================================
# Side-by-side bars showing the proportions of Up/Down/No change in each
# of the two groups (translocation vs genome background).

def plot_grouped_bar(counts_transloc, counts_control, cond):
    """
    Make a grouped bar plot comparing the translocation group to the
    genome background group.
    """
    # The three categories shown on the x-axis
    labels = ["Up", "Down", "No change"]

    # Get the values in the right order
    # .values turns a pandas Series into a numpy array
    trans_vals = counts_transloc[labels].values
    control_vals = counts_control[labels].values

    # -------------------------------------------------------------------------
    # Convert raw counts to proportions
    # -------------------------------------------------------------------------
    # The two groups have very different sizes (e.g. 200 vs 20,000 genes).
    # Comparing raw counts wouldn't be meaningful, we'd just see that one
    # bar is way taller than the other. Proportions let us compare them fairly.
    trans_prop = trans_vals / trans_vals.sum()
    control_prop = control_vals / control_vals.sum()

    # -------------------------------------------------------------------------
    # Set up bar positions
    # -------------------------------------------------------------------------
    # x is the position of each category on the x-axis: [0, 1, 2]
    # np.arange(3) gives us array([0, 1, 2])
    x = np.arange(len(labels))

    # width is how wide each bar should be
    width = 0.35

    # -------------------------------------------------------------------------
    # Create the figure and draw the bars
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8))

    # Draw translocation bars slightly to the LEFT of each x position
    # (x - width/2 puts them centered just to the left)
    ax.bar(
        x - width / 2,
        trans_prop,
        width,
        label="Translocation genes",
        color="steelblue"
    )

    # Draw control bars slightly to the RIGHT of each x position
    ax.bar(
        x + width / 2,
        control_prop,
        width,
        label="Genome background",
        color="orange"
    )

    # -------------------------------------------------------------------------
    # Add labels, title, and legend
    # -------------------------------------------------------------------------
    # Set x-axis tick positions and labels
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    ax.set_ylabel("Fraction of genes")
    ax.set_title(cond + ": DE categories")

    # Show the legend (using the labels we provided when calling ax.bar)
    ax.legend()

    # Adjust spacing
    fig.tight_layout()

    # -------------------------------------------------------------------------
    # Save the plot to a file
    # -------------------------------------------------------------------------
    filename = cond + "_grouped_DE_proportions.png"
    out_path = os.path.join(OUTPUT_FOLDER, filename)
    fig.savefig(out_path, dpi=300)

    # Close the figure to free memory
    plt.close(fig)

    print("  Saved: " + out_path)


# =============================================================================
# Helper function: Run the full analysis for one condition
# =============================================================================

def run_analysis(transloc_genes, de_df, cond):
    """
    For one condition, split genes into the two groups, run both Fisher's
    tests, and make the bar plot.
    """
    # -------------------------------------------------------------------------
    # Step 1: Build the gene sets
    # -------------------------------------------------------------------------
    # The "universe" is all genes that have a DE result for this condition.
    # We use Python sets here because we'll do set operations like
    # intersection (&) and difference (-) below.
    all_genes = set(de_df["gene_name"])

    # Translocated genes that are also in the DE table
    # & is "set intersection", genes in BOTH sets
    trans_genes = all_genes & transloc_genes

    # Everything else, the genome background,
    # is "set difference", genes in all_genes but NOT in trans_genes
    control_genes = all_genes - trans_genes

    print("")
    print(cond)
    print("  Total genes in DE table:  " + str(len(all_genes)))
    print("  Translocated genes:       " + str(len(trans_genes)))
    print("  Genome background genes:  " + str(len(control_genes)))

    # -------------------------------------------------------------------------
    # Step 2: Count Up/Down/No change in each group
    # -------------------------------------------------------------------------
    counts_transloc = get_de_counts(trans_genes, de_df)
    counts_control = get_de_counts(control_genes, de_df)

    print("")
    print("  Counts — translocation: " + str(counts_transloc.to_dict()))
    print("  Counts — genome background: " + str(counts_control.to_dict()))

    # -------------------------------------------------------------------------
    # Step 3: Test 1 - DE vs No change
    # -------------------------------------------------------------------------
    # Are translocated genes more often differentially expressed?
    table_de, OR_de, p_de = fisher_test(
        counts_transloc,
        counts_control,
        "DE_vs_NoChange"
    )

    print("")
    print("  DE enrichment test (translocation vs genome background):")
    print("    Table: " + str(table_de))
    # {:.4f} = 4 decimal places
    # {:.4g} = up to 4 significant digits (used for p-values that may be tiny)
    print("    Odds ratio: {:.4f},  p = {:.4g}".format(OR_de, p_de))

    # -------------------------------------------------------------------------
    # Step 4: Test 2 - Up vs Down
    # -------------------------------------------------------------------------
    # Among DE genes, do translocated genes lean more toward Up than Down?
    table_ud, OR_ud, p_ud = fisher_test(
        counts_transloc,
        counts_control,
        "Up_vs_Down"
    )

    print("")
    print("  Up/Down bias test (among DE genes only):")
    print("    Table: " + str(table_ud))
    print("    Odds ratio: {:.4f},  p = {:.4g}".format(OR_ud, p_ud))

    # -------------------------------------------------------------------------
    # Step 5: Make the plot
    # -------------------------------------------------------------------------
    plot_grouped_bar(counts_transloc, counts_control, cond)


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)


# -----------------------------------------------------------------------------
# Step 1: Load the distance table and pull out translocated gene names
# -----------------------------------------------------------------------------
# We only use this file to figure out which genes are inside translocations.
# The actual DE analysis uses the DESeq2 results.

print("Loading distance table...")
dist = pd.read_csv(DIST_FILE, sep="\t")

# Keep only rows where gene1 is inside a translocation
inside_rows = dist[dist["gene1_region_type"] == "inside_transloc"]

# Get the unique gene names from those rows.
# .unique() gives a numpy array of unique values.
# We turn it into a Python set for fast lookups later.
transloc_gene_array = inside_rows["gene1"].unique()
transloc_genes = set(transloc_gene_array)

print("Translocated genes identified: " + str(len(transloc_genes)))


# -----------------------------------------------------------------------------
# Step 2: For each condition, load the DE file and run the analysis
# -----------------------------------------------------------------------------
for cond in CONDITIONS:

    # Load the DESeq2 results for this condition
    de_path = de_files[cond]
    de_df = pd.read_csv(de_path, sep="\t")

    # ---------------------------------------------------------------------
    # Handle duplicate gene names
    # ---------------------------------------------------------------------
    # Sometimes the same gene_name appears multiple times in DESeq2 output
    # (e.g. due to duplicate Ensembl IDs mapping to the same name).
    # We keep just the row with the most significant p-value (smallest padj).
    #
    # Sort by padj ascending, smallest padj at the top.
    # Then drop duplicates, keeping the first (= smallest padj).
    de_df = de_df.sort_values("padj")
    de_df = de_df.drop_duplicates("gene_name", keep="first")

    # ---------------------------------------------------------------------
    # Add the DE_status column (Up / Down / No change)
    # ---------------------------------------------------------------------
    de_df = classify_de(de_df)

    # ---------------------------------------------------------------------
    # Run the full analysis (both tests + plot)
    # ---------------------------------------------------------------------
    run_analysis(transloc_genes, de_df, cond)


# Done!
print("")
print("Done!")
