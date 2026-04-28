"""
04_plot_results.py
===================
This script makes scatter plots and runs statistical tests comparing
the gene-pair distances from the translocated cells (T1, C1) against
the wildtype (WT) baseline.

The big question:
- After a translocation, did gene pairs get CLOSER together or FURTHER
  apart in 3D space, compared to normal (WT) cells?

How we answer it:
- We have distances for the same gene pairs in WT, T1, and C1
  (from scripts 02 and 03).
- For each comparison (WT vs T1, WT vs C1), we make:
  1. A scatter plot: one dot per gene pair.
     - X-axis: the distance in WT.
     - Y-axis: the distance in the translocated condition.
     - The red dashed line is y = x (no change).
     - Dots BELOW the line = pair is CLOSER in the condition.
     - Dots ABOVE the line = pair is FURTHER apart in the condition.
  2. A Wilcoxon signed-rank test:
     - Tells us if the change is statistically significant.
     - We use this (not a t-test) because the distances aren't necessarily
       normally distributed.

We split the comparisons into translocated-neighbor pairs vs control-neighbor
pairs, so we can see if the translocation has a SPECIFIC effect (translocated
genes get closer to neighbors) or a general effect (everything moves).

Run scripts 02 and 03 first.
Then edit the file paths below and run this script.

Required libraries: pandas, matplotlib, seaborn, scipy
Install them with: pip install pandas matplotlib seaborn scipy
"""

import os                             
import pandas as pd                    
import matplotlib.pyplot as plt        
import seaborn as sns                  
from scipy.stats import wilcoxon       


# =============================================================================
# CONFIG SECTION - Edit these paths before running
# =============================================================================

# Folder containing the aggregated distance files (T1, C1, WT)
RESULTS_DIR = "/path/to/results"

# Folder where we'll save the plots
PLOTS_DIR = "/path/to/results/plots/genepair"


# =============================================================================
# Helper function: Add "chr" prefix to a chromosome name if it's missing
# =============================================================================
# Same helper used in earlier scripts.

def fix_chrom_name(c):
    """If c starts with 'chr', return as-is. Otherwise prepend 'chr'."""
    c = str(c).strip()
    if c.startswith("chr"):
        return c
    else:
        return "chr" + c


# =============================================================================
# Main code starts here
# =============================================================================

print("Loading distance tables...")


# -----------------------------------------------------------------------------
# Step 1: Load the three distance tables
# -----------------------------------------------------------------------------
# Build the file paths
t1_path = os.path.join(RESULTS_DIR, "T1_distances_agg.tsv")
c1_path = os.path.join(RESULTS_DIR, "C1_distances_agg.tsv")
wt_path = os.path.join(RESULTS_DIR, "WT_distances_agg.tsv")

# Read each file
t1_df = pd.read_csv(t1_path, sep="\t")
c1_df = pd.read_csv(c1_path, sep="\t")
wt_df = pd.read_csv(wt_path, sep="\t")


# -----------------------------------------------------------------------------
# Step 2: Make sure all chromosome names start with "chr"
# -----------------------------------------------------------------------------
# We do this for chr1 and chr2 columns in all three tables.

for df in [wt_df, t1_df, c1_df]:
    for col in ["chr1", "chr2"]:
        # List comprehension (one-line for loop) to fix every value
        df[col] = [fix_chrom_name(c) for c in df[col]]


# -----------------------------------------------------------------------------
# Step 3: Print a summary of how many genes are in each category
# -----------------------------------------------------------------------------
# This is just a sanity check, we make sure we have the genes we expect.

# We'll loop through the three tables, with their labels
table_label_pairs = [("T1", t1_df), ("C1", c1_df), ("WT", wt_df)]

for label, df in table_label_pairs:

    # ---------------------------------------------------------------------
    # Get the unique gene names in each region category
    # ---------------------------------------------------------------------
    # A gene can appear as gene1 OR gene2 in our pair table.
    # So we collect from both columns and combine with set union (|).

    # Translocated (inside_transloc) genes
    g1_trans = df.loc[df["gene1_region_type"] == "inside_transloc", "gene1"]
    g2_trans = df.loc[df["gene2_region_type"] == "inside_transloc", "gene2"]
    trans = set(g1_trans) | set(g2_trans)   # | is set union

    # Neighbor genes
    g1_neigh = df.loc[df["gene1_region_type"] == "neighbor", "gene1"]
    g2_neigh = df.loc[df["gene2_region_type"] == "neighbor", "gene2"]
    neighbor = set(g1_neigh) | set(g2_neigh)

    # Control genes
    g1_ctrl = df.loc[df["gene1_region_type"] == "control", "gene1"]
    g2_ctrl = df.loc[df["gene2_region_type"] == "control", "gene2"]
    control = set(g1_ctrl) | set(g2_ctrl)

    # Total unique genes across all three categories
    all_unique = trans | neighbor | control

    print("")
    print(label + ": "
          + str(len(trans)) + " transloc, "
          + str(len(neighbor)) + " neighbor, "
          + str(len(control)) + " control, "
          + str(len(all_unique)) + " total unique")


# -----------------------------------------------------------------------------
# Step 4: Rename distance columns so we can tell datasets apart after merging
# -----------------------------------------------------------------------------
# Each table currently has "mean_distance" and "std_distance" columns.
# After we merge them, they'd collide. So we rename them to be specific.

# Columns we want to KEEP from each table
keep_cols = [
    "gene1", "gene1_region_type", "chr1",
    "gene2", "gene2_region_type", "chr2",
    "mean_distance", "std_distance"
]

# Subset and rename for T1
t1_df = t1_df[keep_cols].rename(columns={
    "mean_distance": "mean_T1",
    "std_distance": "std_T1"
})

# Subset and rename for C1
c1_df = c1_df[keep_cols].rename(columns={
    "mean_distance": "mean_C1",
    "std_distance": "std_C1"
})

# Subset and rename for WT
wt_df = wt_df[keep_cols].rename(columns={
    "mean_distance": "mean_WT",
    "std_distance": "std_WT"
})


# -----------------------------------------------------------------------------
# Step 5: Split WT into translocated-neighbor pairs vs control-neighbor pairs
# -----------------------------------------------------------------------------
# We want to compare these two groups separately.
# A "translocated-neighbor" pair has one gene that's inside_transloc and
# the other is neighbor. (gene1 and gene2 can be in either order.)

is_trans_g1 = (wt_df["gene1_region_type"] == "inside_transloc") & \
              (wt_df["gene2_region_type"] == "neighbor")
is_trans_g2 = (wt_df["gene1_region_type"] == "neighbor") & \
              (wt_df["gene2_region_type"] == "inside_transloc")

# A pair is "trans-neighbor" if either combination is true
is_trans_neigh = is_trans_g1 | is_trans_g2

# Split WT into two groups
wt_transloc = wt_df[is_trans_neigh]      # translocated-neighbor pairs
wt_control = wt_df[~is_trans_neigh]      # everything else (control-neighbor)


# -----------------------------------------------------------------------------
# Step 6: Make sure the output folder exists
# -----------------------------------------------------------------------------
os.makedirs(PLOTS_DIR, exist_ok=True)


# -----------------------------------------------------------------------------
# Step 7: Set up the four comparisons we want to make
# -----------------------------------------------------------------------------
# We want to compare:
# 1. WT trans-neighbor pairs vs T1 (these are the translocated genes' distances in T1)
# 2. WT control-neighbor pairs vs T1
# 3. WT trans-neighbor pairs vs C1
# 4. WT control-neighbor pairs vs C1
#
# Each "comparison" is a tuple: (name, WT subset, condition df, condition column)

comparisons = [
    ("T1_transloc", wt_transloc, t1_df, "mean_T1"),
    ("T1_control",  wt_control,  t1_df, "mean_T1"),
    ("C1_transloc", wt_transloc, c1_df, "mean_C1"),
    ("C1_control",  wt_control,  c1_df, "mean_C1"),
]

# Columns that uniquely identify a gene pair
# (we use these for merging, they should match exactly between WT and condition)
pair_cols = [
    "gene1", "gene2",
    "gene1_region_type", "gene2_region_type",
    "chr1", "chr2"
]


# -----------------------------------------------------------------------------
# Step 8: Process each comparison
# -----------------------------------------------------------------------------
for plot_name, wt_sub, cond_df, cond_col in comparisons:

    print("")
    print("── " + plot_name + " ──")

    # ---------------------------------------------------------------------
    # Step 8a: Filter the condition table to matching region types
    # ---------------------------------------------------------------------
    # We only want condition rows whose gene1_region_type and gene2_region_type
    # are present in the WT subset.
    wt_g1_types = wt_sub["gene1_region_type"].unique()
    wt_g2_types = wt_sub["gene2_region_type"].unique()

    g1_match = cond_df["gene1_region_type"].isin(wt_g1_types)
    g2_match = cond_df["gene2_region_type"].isin(wt_g2_types)

    cond_sub = cond_df[g1_match & g2_match]

    # ---------------------------------------------------------------------
    # Step 8b: Find gene pairs that exist in BOTH WT and the condition
    # ---------------------------------------------------------------------
    # We want to compare the SAME gene pairs in both conditions.
    # An "inner" merge keeps only rows that match between two tables.

    # Get unique pairs from WT
    wt_pairs = wt_sub[pair_cols].drop_duplicates()

    # Get unique pairs from the condition
    cond_pairs = cond_sub[pair_cols].drop_duplicates()

    # Find pairs that appear in both
    common = wt_pairs.merge(cond_pairs, on=pair_cols, how="inner")

    # Filter both tables down to just these shared pairs
    wt_eq = wt_sub.merge(common, on=pair_cols, how="inner")
    cond_eq = cond_sub.merge(common, on=pair_cols, how="inner")

    # Now merge the two filtered tables together side-by-side
    merged = wt_eq.merge(cond_eq, on=pair_cols, how="inner")

    # If no shared pairs, skip this comparison
    if merged.empty:
        print("  No shared pairs to plot, skipping")
        continue

    # ---------------------------------------------------------------------
    # Step 8c: Set up colours and plot limits
    # ---------------------------------------------------------------------
    # Colour palette by gene type
    palette = {
        "inside_transloc": "blue",
        "control": "green",
        "neighbor": "orange"
    }

    # Find the smallest and biggest distances to set the plot range
    # (so the y=x line spans the whole plot)
    wt_min = merged["mean_WT"].min()
    cond_min_val = merged[cond_col].min()
    wt_max = merged["mean_WT"].max()
    cond_max_val = merged[cond_col].max()

    dist_min = min(wt_min, cond_min_val)
    dist_max = max(wt_max, cond_max_val)

    # ---------------------------------------------------------------------
    # Step 8d: Make the scatter plot
    # ---------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    # Use seaborn to draw the scatter plot
    # hue colors the dots by gene1_region_type
    # alpha=0.6 makes them slightly transparent so we can see overlaps
    # s=20 sets the dot size
    sns.scatterplot(
        data=merged,
        x="mean_WT",
        y=cond_col,
        hue="gene1_region_type",
        palette=palette,
        ax=ax,
        alpha=0.6,
        s=20
    )

    # Draw the y = x reference line in red dashed
    # Just two points (smallest and biggest) connected by a straight line
    ax.plot(
        [dist_min, dist_max],
        [dist_min, dist_max],
        color="red",
        linestyle="--",
        label="y = x"
    )

    # Labels and title
    ax.set_xlabel("WT mean distance")
    ax.set_ylabel(plot_name + " mean distance")
    ax.set_title("WT vs " + plot_name)
    ax.legend(title="Gene type")

    # Adjust spacing
    fig.tight_layout()

    # Save the plot
    plot_filename = plot_name + "_scatter.png"
    plot_path = os.path.join(PLOTS_DIR, plot_filename)
    fig.savefig(plot_path)

    # Close the figure to free memory
    plt.close(fig)

    # ---------------------------------------------------------------------
    # Step 8e: Run the Wilcoxon signed-rank test
    # ---------------------------------------------------------------------
    # Tests if the median distance changed significantly between WT and
    # the condition. A small p-value (e.g. < 0.05) means the change is
    # statistically significant.
    stat, p = wilcoxon(merged["mean_WT"], merged[cond_col])

    # Calculate the difference (WT - condition)
    # Positive diff = pair was further in WT (closer in the condition)
    diff = merged["mean_WT"] - merged[cond_col]

    # Count how many pairs became closer in the condition
    # (cond distance is SMALLER than WT distance)
    is_closer = (merged[cond_col] < merged["mean_WT"])
    n_below = is_closer.sum()

    # Total number of pairs we compared
    n_total = len(merged)

    # ---------------------------------------------------------------------
    # Step 8f: Print the results
    # ---------------------------------------------------------------------
    # Format codes:
    #   {:.2f} = 2 decimal places
    #   {:.4g} = up to 4 significant digits
    #   {:.1%} = percentage with 1 decimal
    print("  Wilcoxon: W={:.2f}, p={:.4g}".format(stat, p))

    print("  Distance difference (WT - " + plot_name + "):")
    # .describe() gives count, mean, std, min, quartiles, max
    print(diff.describe().to_string())

    fraction_closer = n_below / n_total
    print("  Pairs closer in " + plot_name + " than WT: "
          + str(n_below) + "/" + str(n_total)
          + " ({:.1%})".format(fraction_closer))


# Done!
print("")
print("All plots saved to: " + PLOTS_DIR)
print("Done!")
