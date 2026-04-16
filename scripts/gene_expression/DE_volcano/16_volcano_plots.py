"""
16_volcano_plots.py
====================
Script for visualising differential expression results as volcano plots,
with translocated genes highlighted and top significant genes labelled.

What this script does
---------------------
A volcano plot displays every gene in a DESeq2 result as a single point,
with log2 fold change on the x-axis and -log10(adjusted p-value) on the
y-axis. Genes that are both strongly changed and statistically significant
appear in the upper left and upper right corners of the plot.

This script extends a standard volcano plot by overlaying the genes that
sit inside translocated genomic regions in a different colour, making it
easy to see at a glance whether translocated genes cluster among the most
differentially expressed genes or are scattered across the plot.

For each condition (T1, C1), the script:
1. Loads DESeq2 results (from script 13).
2. Loads the aggregated distance table (from scripts 01–04) to identify
   which genes are in translocated regions.
3. Classifies every gene as Up, Down, or No change using the configured
   thresholds.
4. Plots all genes as small grey points, then overlays translocated genes
   coloured by their DE status (green = Up, red = Down, blue = No change).
5. Annotates the top N most significant translocated up- and
   downregulated genes with their gene names.
6. Saves the plot as a PNG file.
7. Exports ranked tables of the top translocated up- and downregulated genes.

Why highlight translocated genes on a volcano plot?
----------------------------------------------------
Volcano plots are a standard way to communicate DESeq2 results.
Overlaying the translocated genes in colour immediately
shows whether they are concentrated among the most significantly changed
genes, or whether their expression is largely unchanged despite the
genomic rearrangement. This is a more intuitive visualisation than tables
of numbers alone.

Usage
-----
    1. Run scripts 01–04 (gene-pair distance pipeline) and script 13
       (differential expression) first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 16_volcano_plots.py

Dependencies
------------
    pandas, numpy, matplotlib
    Install with: pip install pandas numpy matplotlib

"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# CONFIG – edit all paths and settings here before running
# =============================================================================
CONFIG = {
    # DESeq2 results from script 13
    "de_files": {
        "T1": "/path/to/DE_WT_vs_T1.tsv",
        "C1": "/path/to/DE_WT_vs_C1.tsv",
    },

    # Aggregated distance tables from scripts 01–04.
    # Used to identify which genes are in translocated regions.
    "dist_files": {
        "T1": "/path/to/T1_distances_agg.tsv",
        "C1": "/path/to/C1_distances_agg.tsv",
    },

    # Thresholds for calling a gene differentially expressed.
    # A gene is called Up/Down if padj < padj_threshold AND
    # |log2FoldChange| > lfc_threshold.
    "lfc_threshold":  1.0,
    "padj_threshold": 0.05,

    # How many top translocated genes to label on the plot
    "top_n_labels": 5,

    # How many top translocated genes to export to the ranked TSV tables
    "top_n_table": 20,

    # Colours for translocated genes by DE status
    "colors": {
        "Up":        "green",
        "Down":      "red",
        "No change": "blue",
    },

    # Where to save output plots and tables (created automatically)
    "output_dir": "/path/to/results/volcano_plots",
}

CONDITIONS = ["T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def classify_de(
    df: pd.DataFrame,
    lfc_thresh: float,
    padj_thresh: float,
) -> pd.DataFrame:
    """
    Add a 'DE_status' column classifying each gene as 'Up', 'Down', or
    'No change' based on its log2 fold change and adjusted p-value.

    Also adds a '-log10(padj)' column for the y-axis of the volcano plot.
    Missing padj values are set to 1 (not significant) and exact zeros are
    replaced with a very small number to avoid log(0) errors.
    """
    df = df.copy()
    df["padj"] = df["padj"].fillna(1).replace(0, 1e-300)
    df["-log10(padj)"] = -np.log10(df["padj"])

    conditions = [
        (df["log2FoldChange"] >  lfc_thresh) & (df["padj"] < padj_thresh),
        (df["log2FoldChange"] < -lfc_thresh) & (df["padj"] < padj_thresh),
    ]
    df["DE_status"] = np.select(conditions, ["Up", "Down"], default="No change")
    return df


def load_transloc_genes(dist_file: str) -> set:
    """
    Load an aggregated distance table and return the set of gene names
    that are inside translocated regions (gene1_region_type == 'inside_transloc').
    """
    dist_df = pd.read_csv(dist_file, sep="\t")
    return set(
        dist_df.loc[dist_df["gene1_region_type"] == "inside_transloc", "gene1"].unique()
    )


def annotate_top_genes(
    ax,
    top_genes: pd.DataFrame,
    n_labels: int,
) -> None:
    """
    Add gene name labels to the top N most significant translocated genes
    on the volcano plot.

    Labels are offset slightly from the point and staggered vertically
    so they do not overlap each other. The horizontal alignment is chosen
    based on whether the gene is up- or downregulated.
    """
    for i, (_, row) in enumerate(top_genes.head(n_labels * 2).iterrows()):
        x_offset = 0.05 * np.sign(row["log2FoldChange"])
        y_offset = 0.1 + 0.05 * i
        ax.text(
            row["log2FoldChange"] + x_offset,
            row["-log10(padj)"]   + y_offset,
            row["gene_name"],
            fontsize=6,
            rotation=30,
            ha="left" if row["log2FoldChange"] > 0 else "right",
            va="bottom",
        )


def plot_volcano(
    de_df: pd.DataFrame,
    transloc_genes: set,
    cond: str,
    lfc_thresh: float,
    padj_thresh: float,
    colors: dict,
    top_n_labels: int,
    top_n_table: int,
    output_dir: str,
) -> None:
    """
    Draw a volcano plot for one condition with translocated genes highlighted,
    annotate the top significant genes, and export ranked gene tables.

    Parameters
    ----------
    de_df          : DESeq2 results with DE_status and -log10(padj) columns.
    transloc_genes : Set of gene names in translocated regions.
    cond           : Condition label ('T1' or 'C1') for titles and filenames.
    lfc_thresh     : log2FC threshold for DE calling.
    padj_thresh    : Adjusted p-value threshold for DE calling.
    colors         : Dict mapping DE status to plot colour.
    top_n_labels   : Number of top genes to label on the plot.
    top_n_table    : Number of top genes to export to TSV tables.
    output_dir     : Where to save the plot and tables.
    """
    # Mark which genes are in translocated regions
    de_df = de_df.copy()
    de_df["is_translocated"] = de_df["gene_name"].isin(transloc_genes)

    trans_df = de_df[de_df["is_translocated"]]
    bg_df    = de_df[~de_df["is_translocated"]]

    # ── Draw plot ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # All background genes in light grey
    ax.scatter(
        bg_df["log2FoldChange"], bg_df["-log10(padj)"],
        c="lightgray", alpha=0.4, s=8,
        label="All other genes",
    )

    # Translocated genes coloured by DE status
    for status, color in colors.items():
        subset = trans_df[trans_df["DE_status"] == status]
        ax.scatter(
            subset["log2FoldChange"], subset["-log10(padj)"],
            c=color, alpha=0.9, s=30,
            edgecolor="black", linewidth=0.3,
            label=f"Translocated — {status}",
        )

    # Dashed threshold lines marking the DE cutoffs
    ax.axvline(x= lfc_thresh,            color="black", linestyle="--", lw=1)
    ax.axvline(x=-lfc_thresh,            color="black", linestyle="--", lw=1)
    ax.axhline(y=-np.log10(padj_thresh), color="black", linestyle="--", lw=1)

    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$(adjusted p-value)")
    ax.set_title(f"WT vs {cond} — translocated genes highlighted")
    ax.legend(frameon=False, fontsize=8)

    # ── Annotate top significant translocated genes ───────────────────────────
    trans_sig = trans_df[
        (trans_df["DE_status"].isin(["Up", "Down"]))
    ]
    top_up   = trans_sig[trans_sig["DE_status"] == "Up"].nlargest(top_n_labels, "-log10(padj)")
    top_down = trans_sig[trans_sig["DE_status"] == "Down"].nlargest(top_n_labels, "-log10(padj)")
    top_genes = pd.concat([top_up, top_down])

    annotate_top_genes(ax, top_genes, top_n_labels)

    fig.tight_layout()
    plot_path = os.path.join(output_dir, f"volcano_WT_vs_{cond}_highlighted.png")
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"  Saved plot: {plot_path}")

    # ── Print top genes to terminal ───────────────────────────────────────────
    print(f"\n  Top {top_n_labels} upregulated translocated genes in {cond}:")
    for gene in top_up["gene_name"].values:
        print(f"    {gene}")
    print(f"\n  Top {top_n_labels} downregulated translocated genes in {cond}:")
    for gene in top_down["gene_name"].values:
        print(f"    {gene}")

    # ── Export ranked gene tables ─────────────────────────────────────────────
    # Top upregulated: ranked by largest log2FC
    top_up_table = (
        trans_sig[trans_sig["DE_status"] == "Up"]
        .nlargest(top_n_table, "log2FoldChange")
    )
    # Top downregulated: ranked by most negative log2FC
    top_down_table = (
        trans_sig[trans_sig["DE_status"] == "Down"]
        .nsmallest(top_n_table, "log2FoldChange")
    )

    up_path   = os.path.join(output_dir, f"top_up_transloc_genes_WT_vs_{cond}.tsv")
    down_path = os.path.join(output_dir, f"top_down_transloc_genes_WT_vs_{cond}.tsv")

    top_up_table.to_csv(up_path,   sep="\t", index=False)
    top_down_table.to_csv(down_path, sep="\t", index=False)

    print(f"\n  Exported top {top_n_table} upregulated genes:   {up_path}")
    print(f"  Exported top {top_n_table} downregulated genes: {down_path}")


# =============================================================================
# Main
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

for cond in CONDITIONS:
    print(f"\n===== Processing {cond} =====")

    # Load and classify DESeq2 results
    de_df = pd.read_csv(CONFIG["de_files"][cond], sep="\t")
    de_df = classify_de(de_df, CONFIG["lfc_threshold"], CONFIG["padj_threshold"])

    # Load translocated gene names from the distance table
    transloc_genes = load_transloc_genes(CONFIG["dist_files"][cond])
    print(f"  Translocated genes identified: {len(transloc_genes)}")
    print(f"  Translocated genes in DE table: "
          f"{de_df['gene_name'].isin(transloc_genes).sum()}")

    plot_volcano(
        de_df          = de_df,
        transloc_genes = transloc_genes,
        cond           = cond,
        lfc_thresh     = CONFIG["lfc_threshold"],
        padj_thresh    = CONFIG["padj_threshold"],
        colors         = CONFIG["colors"],
        top_n_labels   = CONFIG["top_n_labels"],
        top_n_table    = CONFIG["top_n_table"],
        output_dir     = CONFIG["output_dir"],
    )

print("\nDone.")
