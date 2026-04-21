"""
04_plot_results.py
===================
Script 4 of 4 – Visualisation and statistical testing.

Loads the aggregated distance tables from scripts 02 and 03, then for each
condition (T1, C1) makes scatter plots of WT vs condition distances and
runs Wilcoxon signed-rank tests to check if distances changed significantly.

Usage:  Edit the paths below, then run: python 04_plot_results.py
Dependencies:  pip install pandas matplotlib seaborn scipy
"""

import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon


# =============================================================================
# CONFIG – edit these paths before running
# =============================================================================

RESULTS_DIR = "/path/to/results"
PLOTS_DIR = "/path/to/plots"


# =============================================================================
# Load the aggregated distance files from scripts 02 and 03
# =============================================================================

print("Loading distance tables...")

# Load aggregated distances for each condition
t1_df = pd.read_csv(os.path.join(RESULTS_DIR, "T1_distances_agg.tsv"), sep="\t")
c1_df = pd.read_csv(os.path.join(RESULTS_DIR, "C1_distances_agg.tsv"), sep="\t")
wt_df = pd.read_csv(os.path.join(RESULTS_DIR, "WT_distances_agg.tsv"), sep="\t")

# Make sure chromosome columns all have the 'chr' prefix (e.g., "chr1")
for df in (wt_df, t1_df, c1_df):
    for col in ("chr1", "chr2"):
        # Convert to string and ensure "chr" prefix exists
        df[col] = df[col].astype(str).apply(
            lambda x: x if x.startswith("chr") else "chr" + x
        )

# Print summary of gene categories
for label, df in [("T1", t1_df), ("C1", c1_df), ("WT", wt_df)]:
    # Collect genes that are inside translocations
    trans = set(df.loc[df["gene1_region_type"] == "inside_transloc", "gene1"]) | \
            set(df.loc[df["gene2_region_type"] == "inside_transloc", "gene2"])
    # Collect neighbor genes
    neighbor = set(df.loc[df["gene1_region_type"] == "neighbor", "gene1"]) | \
               set(df.loc[df["gene2_region_type"] == "neighbor", "gene2"])
    # Collect control genes
    control = set(df.loc[df["gene1_region_type"] == "control", "gene1"]) | \
              set(df.loc[df["gene2_region_type"] == "control", "gene2"])
    # Print counts
    print(f"\n{label}: {len(trans)} transloc, {len(neighbor)} neighbor, "
          f"{len(control)} control, {len(trans | neighbor | control)} total unique")

# Rename distance columns to distinguish datasets after merging
cols = ["gene1", "gene1_region_type", "chr1",
        "gene2", "gene2_region_type", "chr2",
        "mean_distance", "std_distance"]
# Rename columns for each dataset
t1_df = t1_df[cols].rename(columns={"mean_distance": "mean_T1", "std_distance": "std_T1"})
c1_df = c1_df[cols].rename(columns={"mean_distance": "mean_C1", "std_distance": "std_C1"})
wt_df = wt_df[cols].rename(columns={"mean_distance": "mean_WT", "std_distance": "std_WT"})


# =============================================================================
# Split WT into transloc-neighbor pairs vs control-neighbor pairs
# =============================================================================

# Identify WT pair where:
# one gene is inside_transloc AND the other is neighbor
is_trans_neigh = (
    ((wt_df["gene1_region_type"] == "inside_transloc") & (wt_df["gene2_region_type"] == "neighbor")) |
    ((wt_df["gene1_region_type"] == "neighbor") & (wt_df["gene2_region_type"] == "inside_transloc"))
)
# SUbset WT data
wt_transloc = wt_df[is_trans_neigh]     # Translocated gene pairs
wt_control = wt_df[~is_trans_neigh]     # Control gene pairs


# =============================================================================
# Make scatter plots and run stats
# =============================================================================

os.makedirs(PLOTS_DIR, exist_ok=True)

# Need to do 4 comparisons:
# T1 transloc vs WT,  T1 control vs WT,  C1 transloc vs WT,  C1 control vs WT
comparisons = [
    ("T1_transloc", wt_transloc, t1_df, "mean_T1"),
    ("T1_control",  wt_control,  t1_df, "mean_T1"),
    ("C1_transloc", wt_transloc, c1_df, "mean_C1"),
    ("C1_control",  wt_control,  c1_df, "mean_C1"),
]

# Columns that uniquely define a gene pair
pair_cols = ["gene1", "gene2", "gene1_region_type", "gene2_region_type", "chr1", "chr2"]

# Loop over each comparison
for plot_name, wt_sub, cond_df, cond_col in comparisons:
    print(f"\n── {plot_name} ──")

    # Only keep matching region types between WT and condition
    cond_sub = cond_df[
        (cond_df["gene1_region_type"].isin(wt_sub["gene1_region_type"].unique())) &
        (cond_df["gene2_region_type"].isin(wt_sub["gene2_region_type"].unique()))
    ]

    # Restrict to gene pairs that appear in BOTH WT and condition
    # so we're comparing the exact same pairs
    common = wt_sub[pair_cols].drop_duplicates().merge(
        cond_sub[pair_cols].drop_duplicates(), on=pair_cols, how="inner"
    )
    # Filter WT and condition to only shared pairs
    wt_eq = wt_sub.merge(common, on=pair_cols, how="inner")
    cond_eq = cond_sub.merge(common, on=pair_cols, how="inner")

    # Merged WT and condition distances into one table
    merged = wt_eq.merge(
        cond_eq,
        on=pair_cols,
        how="inner",
    )

    if merged.empty:
        print(f"  No shared pairs to plot, skipping")
        continue

    # --------------------------
    # Create scatter plots
    # --------------------------
    # Points below the red line = closer in condition than in WT
    # Color palette by gene type
    palette = {"inside_transloc": "blue", "control": "green", "neighbor": "orange"}
    # Determine axis limits
    dist_min = min(merged["mean_WT"].min(), merged[cond_col].min())
    dist_max = max(merged["mean_WT"].max(), merged[cond_col].max())

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=merged,
        x="mean_WT", y=cond_col, 
        hue="gene1_region_type",
        palette=palette, ax=ax,
        alpha=0.6, s=20,
    )
    # Add diagonal reference line (y = x)
    ax.plot([dist_min, dist_max], [dist_min, dist_max],
            color="red", linestyle="--", label="y = x")
    ax.set_xlabel("WT mean distance")
    ax.set_ylabel(f"{plot_name} mean distance")
    ax.set_title(f"WT vs {plot_name}")
    ax.legend(title="Gene type")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS_DIR, f"{plot_name}_scatter.png"))
    # Close figure to free memory
    plt.close(fig)

    # --------------------------------------------
    # Statistical test: Wilcoxon signed-rank test
    # --------------------------------------------
    # Non-parametric test: did the median distance change between WT and condition?
    stat, p = wilcoxon(merged["mean_WT"], merged[cond_col])
    # Compute difference (WT - condition)
    diff = merged["mean_WT"] - merged[cond_col]
    # Count how many pairs became closer in condition
    n_below = (merged[cond_col] < merged["mean_WT"]).sum()

    print(f"  Wilcoxon: W={stat:.2f}, p={p:.4g}")
    print(f"  Distance difference (WT - {plot_name}):")
    print(diff.describe().to_string())
    print(f"  Pairs closer in {plot_name} than WT: "
          f"{n_below}/{len(merged)} ({n_below / len(merged):.1%})")

print(f"\nAll plots saved to: {PLOTS_DIR}")
print("Done!")
