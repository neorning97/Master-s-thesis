"""
15_de_enrichment.py
====================
Script for testing whether genes in translocated regions are enriched for
differential expression compared to all other genes in the genome.

What this script does
---------------------
Scripts 01–04 computed 3D distances between translocated genes and their
neighbors and identified which genes are in translocated regions
(inside_transloc). Script 13 identified which genes are differentially
expressed between WT and T1/C1 using DESeq2.

This script connects those two results by asking two questions:

  1. DE enrichment: Are translocated genes more likely to be differentially
     expressed (Up or Down) than all other genes in the genome?

  2. Up/Down bias: Among differentially expressed genes, do translocated
     genes show a bias toward upregulation compared to the genome background?

For each condition (T1, C1), the script:
1. Loads the aggregated WT distance table (from scripts 01–04) to identify
   which genes are in translocated regions.
2. Loads DESeq2 results (from script 13) and classifies each gene as
   Up, Down, or No change.
3. Splits all genes into two groups:
     - translocation genes : genes found in translocated regions
     - genome background   : all other genes (the whole-genome control)
4. Runs two Fisher's exact tests:
     - DE vs No change: are translocated genes more often DE?
     - Up vs Down:      among DE genes, are translocated genes more often up?
5. Produces a grouped bar plot showing proportions side by side.

Why Fisher's exact test?
------------------------
We have a 2x2 contingency table:
  rows    = group (translocation vs genome background)
  columns = DE status (e.g. DE vs not DE)

Fisher's exact test is appropriate here because it makes no assumption about
the minimum cell count, unlike chi-square, which requires sufficiently large
counts in all cells. The translocation group is typically small (~100–500
genes), making Fisher's the safer choice.

Why two separate tests?
-----------------------
The first test (DE vs No change) asks whether translocated genes are
generally more affected by the translocation. The second test (Up vs Down)
asks a more specific question: conditional on being differentially expressed,
are translocated genes more likely to be activated than repressed? These
are biologically distinct questions and are answered separately.

Why use the whole genome as a control?
---------------------------------------
The control is all genes not in any translocated region. This is the most
natural comparison and makes the result interpretable: a significant result
means translocated genes have a different DE rate than the rest of the genome.

Usage
-----
    1. Run scripts 01–04 and script 13 first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 15_de_enrichment.py

Dependencies
------------
    pandas, numpy, matplotlib, scipy
    Install with: pip install pandas numpy matplotlib scipy

"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact


# =============================================================================
# CONFIG – edit all paths and thresholds here before running
# =============================================================================
CONFIG = {
    # Aggregated WT distance table from scripts 01–04.
    # Used only to extract the list of translocated gene names
    # (gene1 where gene1_region_type == 'inside_transloc').
    "dist_file": "/path/to/WT_distances_agg.tsv",

    # DESeq2 results from script 13
    "de_files": {
        "T1": "/path/to/DE_WT_vs_T1.tsv",
        "C1": "/path/to/DE_WT_vs_C1.tsv",
    },

    # Thresholds for calling a gene differentially expressed.
    # A gene is called Up/Down if padj < padj_threshold AND
    # |log2FoldChange| > lfc_threshold.
    "lfc_threshold":  1.0,
    "padj_threshold": 0.05,

    # Where to save output plots (created automatically)
    "output_dir": "/path/to/plots",
}

CONDITIONS = ["T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def classify_de(df: pd.DataFrame, lfc_thresh: float, padj_thresh: float) -> pd.DataFrame:
    """
    Add a 'DE_status' column classifying each gene as 'Up', 'Down', or
    'No change' based on its log2 fold change and adjusted p-value.

    A gene is only called Up or Down if BOTH thresholds are met:
    padj < padj_threshold AND |log2FoldChange| > lfc_threshold.
    """
    conditions = [
        (df["log2FoldChange"] >  lfc_thresh) & (df["padj"] < padj_thresh),
        (df["log2FoldChange"] < -lfc_thresh) & (df["padj"] < padj_thresh),
    ]
    df["DE_status"] = np.select(conditions, ["Up", "Down"], default="No change")
    return df


def get_de_counts(gene_set: set, de_df: pd.DataFrame) -> pd.Series:
    """
    Count how many genes in a given set fall into each DE category
    (Up, Down, No change).

    Only genes present in both the gene set and the DE table are counted.
    """
    subset = de_df[de_df["gene_name"].isin(gene_set)]
    return subset["DE_status"].value_counts().reindex(
        ["Up", "Down", "No change"], fill_value=0
    )


def fisher_test(
    counts_a: pd.Series,
    counts_b: pd.Series,
    test_type: str,
) -> tuple[list, float, float]:
    """
    Run a two-sided Fisher's exact test comparing two gene groups.

    Two test types are supported:

    'DE_vs_NoChange'
        Tests whether group A has a higher rate of differential expression
        (Up or Down combined) than group B.
        Table:
                        DE      Not DE
          Group A   [a]       [b]
          Group B   [c]       [d]

    'Up_vs_Down'
        Tests whether group A has a different Up/Down ratio than group B.
        Only genes that are DE (Up or Down) are included.
        Table:
                        Up      Down
          Group A   [a]       [b]
          Group B   [c]       [d]

    Returns
    -------
    table     : 2x2 contingency table as a list of lists
    oddsratio : Odds ratio from the Fisher test
    p_value   : Two-sided p-value
    """
    if test_type == "DE_vs_NoChange":
        a = int(counts_a["Up"] + counts_a["Down"])
        b = int(counts_a["No change"])
        c = int(counts_b["Up"] + counts_b["Down"])
        d = int(counts_b["No change"])
    elif test_type == "Up_vs_Down":
        a = int(counts_a["Up"])
        b = int(counts_a["Down"])
        c = int(counts_b["Up"])
        d = int(counts_b["Down"])
    else:
        raise ValueError(f"Unknown test_type '{test_type}'. "
                         "Use 'DE_vs_NoChange' or 'Up_vs_Down'.")

    table = [[a, b], [c, d]]
    oddsratio, p = fisher_exact(table)
    return table, oddsratio, p


def plot_grouped_bar(
    counts_transloc: pd.Series,
    counts_control: pd.Series,
    cond: str,
    output_dir: str,
) -> None:
    """
    Produce a grouped bar plot showing the fraction of Up, Down, and
    No change genes in the translocation and genome background groups.

    Proportions rather than raw counts are shown so the two groups —
    which differ greatly in size — can be compared visually.
    """
    labels       = ["Up", "Down", "No change"]
    trans_vals   = counts_transloc[labels].values
    control_vals = counts_control[labels].values

    trans_prop   = trans_vals   / trans_vals.sum()
    control_prop = control_vals / control_vals.sum()

    x     = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - width / 2, trans_prop,   width,
           label="Translocation genes", color="steelblue")
    ax.bar(x + width / 2, control_prop, width,
           label="Genome background",   color="orange")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of genes")
    ax.set_title(f"{cond}: DE categories — translocation vs genome background")
    ax.legend()

    fig.tight_layout()
    out_path = os.path.join(output_dir, f"{cond}_grouped_DE_proportions.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def run_analysis(
    transloc_genes: set,
    de_df: pd.DataFrame,
    cond: str,
    output_dir: str,
) -> None:
    """
    Split all genes into translocation and genome background groups, run
    both Fisher tests, and produce a grouped bar plot.

    Parameters
    ----------
    transloc_genes : Set of gene names identified as inside translocated regions.
    de_df          : DESeq2 results with DE_status already assigned.
    cond           : Condition label ('T1' or 'C1').
    output_dir     : Where to save the plot.
    """
    # All genes with a DE result are the universe for this analysis
    all_genes     = set(de_df["gene_name"])
    trans_genes   = all_genes & transloc_genes          # intersection
    control_genes = all_genes - trans_genes             # everything else

    print(f"\n{cond}")
    print(f"  Total genes in DE table:  {len(all_genes)}")
    print(f"  Translocated genes:       {len(trans_genes)}")
    print(f"  Genome background genes:  {len(control_genes)}")

    counts_transloc = get_de_counts(trans_genes,   de_df)
    counts_control  = get_de_counts(control_genes, de_df)

    print(f"\n  Counts — translocation: {counts_transloc.to_dict()}")
    print(f"  Counts — genome background: {counts_control.to_dict()}")

    # ── Test 1: DE enrichment (DE vs No change) ───────────────────────────────
    table, OR, p = fisher_test(counts_transloc, counts_control, "DE_vs_NoChange")
    print(f"\n  DE enrichment test (translocation vs genome background):")
    print(f"    Table: {table}")
    print(f"    Odds ratio: {OR:.4f},  p = {p:.4g}")

    # ── Test 2: Directional bias (Up vs Down) ────────────────────────────────
    table_ud, OR_ud, p_ud = fisher_test(counts_transloc, counts_control, "Up_vs_Down")
    print(f"\n  Up/Down bias test (among DE genes only):")
    print(f"    Table: {table_ud}")
    print(f"    Odds ratio: {OR_ud:.4f},  p = {p_ud:.4g}")

    plot_grouped_bar(counts_transloc, counts_control, cond, output_dir)


# =============================================================================
# Main
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Load distance table and extract translocated gene names ──────────────────
# The distance table is used only to identify which genes were in translocated
# regions. The DE analysis itself uses the DESeq2 results directly.
print("Loading distance table...")
dist = pd.read_csv(CONFIG["dist_file"], sep="\t")

transloc_genes = set(
    dist.loc[dist["gene1_region_type"] == "inside_transloc", "gene1"].unique()
)
print(f"Translocated genes identified: {len(transloc_genes)}")

# ── Run analysis for each condition ───────────────────────────────────────────
for cond in CONDITIONS:
    de_df = pd.read_csv(CONFIG["de_files"][cond], sep="\t")

    # If the same gene appears multiple times (e.g. duplicate Ensembl entries),
    # keep only the row with the most significant p-value
    de_df = de_df.sort_values("padj").drop_duplicates("gene_name", keep="first")

    de_df = classify_de(de_df, CONFIG["lfc_threshold"], CONFIG["padj_threshold"])

    run_analysis(transloc_genes, de_df, cond, CONFIG["output_dir"])

print("\nDone.")
