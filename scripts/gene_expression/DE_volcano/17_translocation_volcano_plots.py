"""
17_translocation_volcano_plots.py
===================================
Script for producing volcano plots with genes from specific translocation
regions highlighted, and running statistical tests comparing those genes to
the genome-wide background.

What this script does
---------------------
Script 16 highlighted translocated genes identified from the gene-pair
distance pipeline. This script takes a complementary approach: it defines
translocated regions directly by their genomic coordinates (from the
translocation BED files) and uses the GTF annotation to find all genes
overlapping those regions. This allows analysis of specific named
translocations (e.g. Der(17)t(3;17), Der(3)t(3;17)) independently.

For each translocation and each condition (T1, C1), the script:
1. Finds all genes whose coordinates overlap the translocation regions
   using the GTF annotation.
2. Loads DESeq2 results and classifies genes as Up, Down, or No change.
3. Produces a volcano plot with translocation genes highlighted in colour.
4. Runs two statistical tests:
     a. DE enrichment Fisher test: are genes in this translocation more
        often DE than the genome background?
     b. Up/Down bias Fisher test: among DE genes, are translocation genes
        more often upregulated?
5. Exports ranked tables of the top up- and downregulated genes.

Difference from script 16
--------------------------
Script 16 identifies translocated genes from the distance table produced
by scripts 01–04. This script identifies them directly from user-defined
genomic coordinates, making it useful for analysing specific translocation
events by name rather than relying on the distance pipeline output.

Translocations analysed
------------------------
The translocations are defined in CONFIG as a dictionary mapping each
translocation name to a list of (chromosome, start, end) regions. Each
translocation can span multiple chromosomal segments (e.g. Der(17)t(3;17)
involves both the chr3 fragment and the chr17 body).

Usage
-----
    1. Run script 13 (differential expression) first.
    2. Edit the CONFIG section below to point to your files and define
       the genomic coordinates of the translocations to analyse.
    3. Run: python 17_translocation_volcano_plots.py

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
# CONFIG – edit all paths and settings here before running
# =============================================================================
CONFIG = {
    # DESeq2 results from script 13
    "de_files": {
        "T1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results1/DE/DE_WT_vs_T1.tsv",
        "C1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results1/DE/DE_WT_vs_C1.tsv",
    },

    # GTF gene annotation file (Ensembl GRCh38 protein-coding genes)
    "gtf_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",

    # Translocation regions to analyse.
    # Each key is a human-readable translocation name used in plot titles
    # and filenames. Each value is a list of (chromosome, start, end) tuples
    # defining the genomic regions involved. A gene is included if any part
    # of it overlaps any of the listed regions.
    "translocations": {
        "Der(17)t(3;17)": [
            ("chr3",  0,        58_600_000),   # chr3 fragment
            ("chr17", 22_700_001, 83_257_441), # chr17 body
        ],
        "Der(3)t(3;17)": [
            ("chr3",  58_600_001, 198_295_559), # chr3 body
            ("chr17", 0,          22_700_000),  # chr17 fragment
        ],
    },

    # Thresholds for calling a gene differentially expressed
    "lfc_threshold":  1.0,
    "padj_threshold": 0.05,

    # Number of top genes to label on the plot and export to TSV tables
    "top_n_labels": 5,
    "top_n_table":  20,

    # Colours for translocated genes by DE status
    "colors": {
        "Up":        "green",
        "Down":      "red",
        "No change": "blue",
    },

    # Where to save all output files (created automatically)
    "output_dir": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results1/translocation_volcano_plots",
}

CONDITIONS = ["T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def load_gtf_genes(gtf_file: str) -> pd.DataFrame:
    """
    Load a GTF annotation and return one row per gene with columns:
    chr, start, end, gene_name.

    Only 'gene' feature rows are kept.
    """
    gtf_cols = ["chr", "source", "feature", "start", "end",
                "score", "strand", "frame", "attributes"]
    gtf = pd.read_csv(
        gtf_file, sep="\t", comment="#", header=None, names=gtf_cols
    )
    gtf = gtf[gtf["feature"] == "gene"].copy()
    gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')
    return gtf[["chr", "start", "end", "gene_name"]].copy()


def get_transloc_genes(gtf: pd.DataFrame, regions: list[tuple]) -> set:
    """
    Return the set of gene names whose coordinates overlap any of the
    given genomic regions.

    A gene overlaps a region if any part of it falls within [start, end].
    This is the standard BED-style overlap check.

    Parameters
    ----------
    gtf     : GTF DataFrame with columns chr, start, end, gene_name.
    regions : List of (chromosome, start, end) tuples.
    """
    genes = set()
    for chrom, start, end in regions:
        overlapping = gtf[
            (gtf["chr"] == chrom) &
            (gtf["end"]   >= start) &
            (gtf["start"] <= end)
        ]["gene_name"]
        genes.update(overlapping)
    return genes


def classify_de(
    df: pd.DataFrame,
    lfc_thresh: float,
    padj_thresh: float,
) -> pd.DataFrame:
    """
    Add 'DE_status' and '-log10(padj)' columns to a DESeq2 results DataFrame.

    Missing padj values are set to 1 (not significant) and exact zeros are
    replaced with 1e-300 to avoid log(0) errors.
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


def run_statistics(
    df: pd.DataFrame,
    transloc_label: str,
    cond: str,
) -> None:
    """
    Run two statistical tests comparing genes in the translocation region
    against the genome-wide background, and print the results.

    Tests:
      1. DE enrichment (Fisher): are translocation genes more often DE?
      2. Up/Down bias (Fisher): among DE genes, are translocation genes
         more often upregulated?

    Parameters
    ----------
    df             : Full DESeq2 DataFrame with 'is_translocated' and
                     'DE_status' columns already assigned.
    transloc_label : Name of the translocation (for print output).
    cond           : Condition label ('T1' or 'C1').
    """
    trans_df  = df[df["is_translocated"]]
    genome_df = df[~df["is_translocated"]]

    n_trans  = len(trans_df)
    n_genome = len(genome_df)

    if n_trans == 0:
        print(f"  No translocation genes found for {transloc_label} in {cond}")
        return

    # ── Print DE fractions ────────────────────────────────────────────────────
    print(f"\n  {transloc_label} | WT vs {cond}")
    print(f"  Translocation genes: {n_trans}")
    for status in ["Up", "Down", "No change"]:
        n = (trans_df["DE_status"] == status).sum()
        print(f"    {status}: {n} ({n/n_trans:.1%})")

    # ── Fisher test 1: DE enrichment ─────────────────────────────────────────
    trans_de  = (trans_df["DE_status"]  != "No change").sum()
    trans_no  = (trans_df["DE_status"]  == "No change").sum()
    genome_de = (genome_df["DE_status"] != "No change").sum()
    genome_no = (genome_df["DE_status"] == "No change").sum()

    or_de, p_de = fisher_exact(
        [[trans_de, trans_no], [genome_de, genome_no]]
    )
    print(f"\n  DE enrichment Fisher test:")
    print(f"    Translocation: DE={trans_de}, No change={trans_no}")
    print(f"    Genome-wide:   DE={genome_de}, No change={genome_no}")
    print(f"    OR={or_de:.3f},  p={p_de:.4g}")

    # ── Fisher test 2: Up/Down bias ───────────────────────────────────────────
    # Restrict to DE genes only for this test
    trans_up    = (trans_df["DE_status"]  == "Up").sum()
    trans_down  = (trans_df["DE_status"]  == "Down").sum()
    genome_up   = (genome_df["DE_status"] == "Up").sum()
    genome_down = (genome_df["DE_status"] == "Down").sum()

    or_ud, p_ud = fisher_exact(
        [[trans_up, trans_down], [genome_up, genome_down]]
    )
    print(f"\n  Up/Down bias Fisher test (among DE genes only):")
    print(f"    Translocation: Up={trans_up}, Down={trans_down}")
    print(f"    Genome-wide:   Up={genome_up}, Down={genome_down}")
    print(f"    OR={or_ud:.3f},  p={p_ud:.4g}")


def plot_volcano(
    df: pd.DataFrame,
    transloc_genes: set,
    transloc_label: str,
    cond: str,
    lfc_thresh: float,
    padj_thresh: float,
    colors: dict,
    top_n_labels: int,
    top_n_table: int,
    output_dir: str,
) -> None:
    """
    Draw a volcano plot for one condition and one translocation, run
    statistical tests, and export ranked gene tables.

    Parameters
    ----------
    df             : DESeq2 results with DE_status and -log10(padj) columns.
    transloc_genes : Set of gene names in the translocation region.
    transloc_label : Human-readable translocation name.
    cond           : Condition label ('T1' or 'C1').
    lfc_thresh     : log2FC threshold used for DE calling.
    padj_thresh    : Adjusted p-value threshold used for DE calling.
    colors         : Dict mapping DE status to plot colour.
    top_n_labels   : Number of top genes to label on the plot.
    top_n_table    : Number of top genes to export to TSV.
    output_dir     : Where to save outputs.
    """
    df = df.copy()
    df["is_translocated"] = df["gene_name"].isin(transloc_genes)

    trans_df = df[df["is_translocated"]]
    bg_df    = df[~df["is_translocated"]]

    # ── Draw volcano plot ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # All background genes in light grey
    ax.scatter(
        bg_df["log2FoldChange"], bg_df["-log10(padj)"],
        c="lightgray", alpha=0.4, s=8,
        label="All other genes",
    )

    # Translocation genes coloured by DE status
    for status, color in colors.items():
        subset = trans_df[trans_df["DE_status"] == status]
        if len(subset) > 0:
            ax.scatter(
                subset["log2FoldChange"], subset["-log10(padj)"],
                c=color, alpha=0.9, s=30,
                edgecolor="black", linewidth=0.3,
                label=f"Translocated — {status}",
            )

    # Dashed threshold lines
    ax.axvline(x= lfc_thresh,            color="black", linestyle="--", lw=1)
    ax.axvline(x=-lfc_thresh,            color="black", linestyle="--", lw=1)
    ax.axhline(y=-np.log10(padj_thresh), color="black", linestyle="--", lw=1)

    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$(adjusted p-value)")
    ax.set_title(f"WT vs {cond} — {transloc_label} genes highlighted")
    ax.legend(frameon=False, fontsize=8)

    # ── Annotate top significant translocation genes ──────────────────────────
    trans_sig  = trans_df[trans_df["DE_status"].isin(["Up", "Down"])]
    top_up     = trans_sig[trans_sig["DE_status"] == "Up"].nlargest(top_n_labels, "-log10(padj)")
    top_down   = trans_sig[trans_sig["DE_status"] == "Down"].nlargest(top_n_labels, "-log10(padj)")
    top_genes  = pd.concat([top_up, top_down])

    for i, (_, row) in enumerate(top_genes.iterrows()):
        x_offset = 0.05 * np.sign(row["log2FoldChange"])
        y_offset = 0.1 + 0.05 * i
        ax.text(
            row["log2FoldChange"] + x_offset,
            row["-log10(padj)"]   + y_offset,
            row["gene_name"],
            fontsize=6, rotation=30,
            ha="left" if row["log2FoldChange"] > 0 else "right",
            va="bottom",
        )

    fig.tight_layout()

    # Make a filename-safe version of the translocation label
    safe_label = transloc_label.replace("(", "").replace(")", "").replace(";", "_")
    plot_path  = os.path.join(output_dir,
                              f"volcano_WT_vs_{cond}_{safe_label}.png")
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"  Saved plot: {plot_path}")

    # ── Print top genes ───────────────────────────────────────────────────────
    print(f"\n  Top {top_n_labels} upregulated genes ({transloc_label}, {cond}):")
    for gene in top_up["gene_name"].values:
        print(f"    {gene}")
    print(f"\n  Top {top_n_labels} downregulated genes ({transloc_label}, {cond}):")
    for gene in top_down["gene_name"].values:
        print(f"    {gene}")

    # ── Export ranked gene tables ─────────────────────────────────────────────
    up_path   = os.path.join(output_dir, f"top_up_{safe_label}_WT_vs_{cond}.tsv")
    down_path = os.path.join(output_dir, f"top_down_{safe_label}_WT_vs_{cond}.tsv")

    (trans_sig[trans_sig["DE_status"] == "Up"]
     .nlargest(top_n_table, "log2FoldChange")
     .to_csv(up_path, sep="\t", index=False))

    (trans_sig[trans_sig["DE_status"] == "Down"]
     .nsmallest(top_n_table, "log2FoldChange")
     .to_csv(down_path, sep="\t", index=False))

    print(f"\n  Exported top {top_n_table} upregulated genes:   {up_path}")
    print(f"  Exported top {top_n_table} downregulated genes: {down_path}")

    # ── Run statistical tests ─────────────────────────────────────────────────
    run_statistics(df, transloc_label, cond)


# =============================================================================
# Main
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Load GTF once — used to map genomic coordinates to gene names ─────────────
print("Loading GTF gene annotations...")
gtf = load_gtf_genes(CONFIG["gtf_file"])
print(f"  {len(gtf)} protein-coding genes loaded.")

# ── Load and classify DESeq2 results ─────────────────────────────────────────
print("Loading DESeq2 results...")
de_results = {}
for cond in CONDITIONS:
    de_results[cond] = classify_de(
        pd.read_csv(CONFIG["de_files"][cond], sep="\t"),
        CONFIG["lfc_threshold"],
        CONFIG["padj_threshold"],
    )
    print(f"  {cond}: {len(de_results[cond])} genes loaded.")

# ── Run analysis for each translocation × condition combination ───────────────
for transloc_label, regions in CONFIG["translocations"].items():
    print(f"\n{'='*55}")
    print(f"Translocation: {transloc_label}")
    print(f"{'='*55}")

    # Find all genes overlapping this translocation's genomic regions
    transloc_genes = get_transloc_genes(gtf, regions)
    print(f"  Genes in translocation region: {len(transloc_genes)}")

    for cond in CONDITIONS:
        plot_volcano(
            df             = de_results[cond],
            transloc_genes = transloc_genes,
            transloc_label = transloc_label,
            cond           = cond,
            lfc_thresh     = CONFIG["lfc_threshold"],
            padj_thresh    = CONFIG["padj_threshold"],
            colors         = CONFIG["colors"],
            top_n_labels   = CONFIG["top_n_labels"],
            top_n_table    = CONFIG["top_n_table"],
            output_dir     = CONFIG["output_dir"],
        )

print("\nDone.")
