"""
04_plot_results.py
===================
Script 4 of 4 for gene pair distance calculations – Visualisation and statistical testing.

What this script does
---------------------
Loads the aggregated distance tables produced by scripts 02 and 03, then
for each condition (T1, C1) and each gene-pair type (translocated–neighbour
and control–neighbour) produces:

  1. A scatter plot of WT mean distance vs condition mean distance.
     Each point is one gene pair. The red dashed line is y = x (equal distance).
     Points below the line mean the pair is CLOSER in the condition than in WT.
     Points above the line mean the pair is FURTHER APART.

  2. A Wilcoxon signed-rank test comparing the WT and condition distances.
     This non-parametric test is used because distance distributions are
     non-normal. It tests whether the median distance is significantly
     shifted between WT and the translocated condition.

Why enforce equal pairs?
------------------------
Before comparing WT and condition distances we restrict both datasets to the
gene pairs that appear in both. This ensures we are always comparing the same
pairs like-for-like, rather than comparing different subsets.

Usage
-----
    1. Run scripts 01, 02, and 03 first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 04_plot_results.py

Dependencies
------------
    pandas, matplotlib, seaborn, scipy
    Install with: pip install pandas matplotlib seaborn scipy

"""

import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    # Results directory from scripts 02 and 03
    "results_dir": "/path/to/results",

    # Where to save the output plots (created automatically)
    "plots_dir": "/path/to/plots",
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def count_gene_types(df: pd.DataFrame, label: str) -> None:
    """Print the number of unique genes in each region-type category."""
    trans    = set(df.loc[df["gene1_region_type"] == "inside_transloc", "gene1"]) | \
               set(df.loc[df["gene2_region_type"] == "inside_transloc", "gene2"])
    neighbor = set(df.loc[df["gene1_region_type"] == "neighbor", "gene1"]) | \
               set(df.loc[df["gene2_region_type"] == "neighbor", "gene2"])
    control  = set(df.loc[df["gene1_region_type"] == "control", "gene1"]) | \
               set(df.loc[df["gene2_region_type"] == "control", "gene2"])
    print(f"\nGene counts in {label}:")
    print(f"  inside_transloc : {len(trans)}")
    print(f"  neighbor        : {len(neighbor)}")
    print(f"  control         : {len(control)}")
    print(f"  total unique    : {len(trans | neighbor | control)}")


def enforce_equal_pairs(
    wt_df: pd.DataFrame,
    cond_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Restrict both DataFrames to the gene pairs that appear in both.

    This ensures the WT and condition distances are compared on exactly the
    same set of pairs. Pairs that could not be mapped in one condition (e.g.
    because a gene had no bead on that chromosome) are excluded from both.
    """
    pair_cols = ["gene1", "gene2", "gene1_region_type", "gene2_region_type", "chr1", "chr2"]
    common    = wt_df[pair_cols].drop_duplicates().merge(
                    cond_df[pair_cols].drop_duplicates(), on=pair_cols, how="inner"
                )
    return (
        wt_df.merge(common, on=pair_cols, how="inner"),
        cond_df.merge(common, on=pair_cols, how="inner"),
    )


# =============================================================================
# Main plotting function
# =============================================================================

def plot_wt_vs_condition(
    wt_sub: pd.DataFrame,
    cond_df: pd.DataFrame,
    cond_col: str,
    plot_name: str,
    output_dir: str,
) -> pd.DataFrame:
    """
    Generate a scatter plot comparing WT and one condition, then run a
    Wilcoxon signed-rank test.

    Parameters
    ----------
    wt_sub     : WT aggregated distances (filtered to relevant region types).
    cond_df    : Condition (T1 or C1) aggregated distances.
    cond_col   : Column name for condition mean distances (e.g. 'mean_T1').
    plot_name  : Label used in plot titles and output filenames.
    output_dir : Directory where the plot is saved.

    Returns
    -------
    merged : DataFrame of matched WT/condition pairs (used for the plot).
    """
    # Keep only region-type combinations that appear in the WT subset
    cond_sub = cond_df[
        (cond_df["gene1_region_type"].isin(wt_sub["gene1_region_type"].unique())) &
        (cond_df["gene2_region_type"].isin(wt_sub["gene2_region_type"].unique()))
    ]

    # Restrict to pairs present in both WT and condition
    wt_eq, cond_eq = enforce_equal_pairs(wt_sub, cond_sub)

    merged = wt_eq.merge(
        cond_eq,
        on=["gene1", "gene2", "gene1_region_type", "gene2_region_type", "chr1", "chr2"],
        how="inner",
    )

    if merged.empty:
        print(f"  No shared pairs to plot for {plot_name}")
        return merged

    # ── Scatter plot ──────────────────────────────────────────────────────────
    palette  = {"inside_transloc": "blue", "control": "green", "neighbor": "orange"}
    dist_min = min(merged["mean_WT"].min(), merged[cond_col].min())
    dist_max = max(merged["mean_WT"].max(), merged[cond_col].max())

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=merged,
        x="mean_WT", y=cond_col,
        hue="gene1_region_type",
        palette=palette, ax=ax,
        alpha=0.6,
        s=20
    )
    # y = x reference line: points below = closer in condition than WT
    ax.plot([dist_min, dist_max], [dist_min, dist_max],
            color="red", linestyle="--", label="y = x")
    ax.set_xlabel("WT mean distance")
    ax.set_ylabel(f"{plot_name} mean distance")
    ax.set_title(f"WT vs {plot_name}")
    ax.legend(title="Gene type")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{plot_name}_scatter.png"))
    plt.close(fig)

    # ── Wilcoxon signed-rank test ─────────────────────────────────────────────
    # Tests whether the median WT distance differs significantly from the
    # condition distance across all gene pairs
    stat, p = wilcoxon(merged["mean_WT"], merged[cond_col])
    diff    = merged["mean_WT"] - merged[cond_col]
    n_below = (merged[cond_col] < merged["mean_WT"]).sum()

    print(f"\n  {plot_name} Wilcoxon signed-rank test: W={stat:.2f}, p={p:.4g}")
    print(f"  Distance difference (WT − {plot_name}) summary:")
    print(diff.describe().to_string())
    print(f"  Gene pairs closer in {plot_name} than WT: "
          f"{n_below}/{len(merged)} ({n_below / len(merged):.1%})")

    return merged


def plot_results(
    t1_agg_file: str,
    c1_agg_file: str,
    wt_agg_file: str,
    output_dir: str,
) -> None:
    """
    Load aggregated distance tables and produce scatter plots + Wilcoxon tests
    for all condition–WT comparisons.

    Four comparisons are made:
      T1_transloc : translocated–neighbour pairs in T1 vs WT
      T1_control  : control–neighbour pairs in T1 vs WT
      C1_transloc : translocated–neighbour pairs in C1 vs WT
      C1_control  : control–neighbour pairs in C1 vs WT

    Parameters
    ----------
    t1_agg_file : Aggregated T1 distances (from script 02).
    c1_agg_file : Aggregated C1 distances (from script 02).
    wt_agg_file : Aggregated WT distances (from script 03).
    output_dir  : Directory where plots are saved.
    """
    os.makedirs(output_dir, exist_ok=True)

    t1_df = pd.read_csv(t1_agg_file, sep="\t")
    c1_df = pd.read_csv(c1_agg_file, sep="\t")
    wt_df = pd.read_csv(wt_agg_file, sep="\t")

    # Ensure consistent 'chr' prefix on all chromosome columns
    for df in (wt_df, t1_df, c1_df):
        for col in ("chr1", "chr2"):
            df[col] = df[col].astype(str).apply(
                lambda x: x if x.startswith("chr") else "chr" + x
            )

    count_gene_types(t1_df, "T1")
    count_gene_types(c1_df, "C1")
    count_gene_types(wt_df, "WT")

    # Rename distance columns so all three conditions can be merged later
    cols  = ["gene1", "gene1_region_type", "chr1",
             "gene2", "gene2_region_type", "chr2",
             "mean_distance", "std_distance"]
    t1_df = t1_df[cols].rename(columns={"mean_distance": "mean_T1", "std_distance": "std_T1"})
    c1_df = c1_df[cols].rename(columns={"mean_distance": "mean_C1", "std_distance": "std_C1"})
    wt_df = wt_df[cols].rename(columns={"mean_distance": "mean_WT", "std_distance": "std_WT"})

    # Split WT into transloc–neighbour and control–neighbour subsets so each
    # comparison is made against the correct WT baseline
    def is_transloc_neighbor(df):
        return (
            ((df["gene1_region_type"] == "inside_transloc") & (df["gene2_region_type"] == "neighbor")) |
            ((df["gene1_region_type"] == "neighbor")        & (df["gene2_region_type"] == "inside_transloc"))
        )

    wt_transloc = wt_df[is_transloc_neighbor(wt_df)]
    wt_control  = wt_df[~is_transloc_neighbor(wt_df)]

    # ── T1 comparisons ────────────────────────────────────────────────────────
    print("\n── T1 vs WT ──")
    plot_wt_vs_condition(wt_transloc, t1_df, "mean_T1", "T1_transloc", output_dir)
    plot_wt_vs_condition(wt_control,  t1_df, "mean_T1", "T1_control",  output_dir)

    # ── C1 comparisons ────────────────────────────────────────────────────────
    print("\n── C1 vs WT ──")
    plot_wt_vs_condition(wt_transloc, c1_df, "mean_C1", "C1_transloc", output_dir)
    plot_wt_vs_condition(wt_control,  c1_df, "mean_C1", "C1_control",  output_dir)

    print(f"\n  All plots saved to: {output_dir}")


# =============================================================================
# Run
# =============================================================================

if __name__ == "__main__":
    results_dir = CONFIG["results_dir"]
    plots_dir   = CONFIG["plots_dir"]

    print("Generating plots and running statistical tests...")
    plot_results(
        t1_agg_file=os.path.join(results_dir, "T1_distances_agg.tsv"),
        c1_agg_file=os.path.join(results_dir, "C1_distances_agg.tsv"),
        wt_agg_file=os.path.join(results_dir, "WT_distances_agg.tsv"),
        output_dir=plots_dir,
    )

    print("\nDone. Output plots:")
    for name in ["T1_transloc", "T1_control", "C1_transloc", "C1_control"]:
        print(f"  {plots_dir}/{name}_scatter.png")
