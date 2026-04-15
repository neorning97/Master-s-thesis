"""
14_gene_bin_expression.py
==========================
Script for linking subcompartment behaviour of translocated genomic bins to
gene expression changes, and testing whether compartment switching is
associated with differential expression.

What this script does
---------------------
Script 10 classified each 100 kb bin in translocated regions as retained,
adopted, or other. Script 13 identified which genes are differentially
expressed between WT and T1/C1. This script connects those two analyses:
it asks whether genes sitting in bins that adopted a new subcompartment are
more likely to be differentially expressed than genes in retained bins.

Step by step:
1. Loads the bin classification table (from script 10).
2. Loads gene coordinates from the GTF annotation and maps each gene to its
   overlapping 100 kb bin, inheriting the bin's subcompartment behaviour
   labels (retained/adopted/other, A_to_B/B_to_A/retained).
3. Loads RNA-seq TPM expression values and computes log2(TPM+1) and
   log2 fold changes vs WT for each condition.
4. Loads DESeq2 results (from script 13) and classifies each gene as
   upregulated, downregulated, or not significant (ns) based on the
   configured thresholds.
5. For each condition (T1, C1) and each classification scheme
   (collapsed A/B, A↔B direction, full subcompartment), runs:
     - A Kruskal-Wallis test: does log2FC differ across bin categories?
     - A chi-square test: does the DE rate (up/down/ns) differ across
       bin categories?
6. Produces stacked bar plots showing the fraction of up/down/ns genes
   per bin category, with a genome-wide control bar for comparison.

Why Kruskal-Wallis?
-------------------
Log2 fold changes are not normally distributed across genes, so a
non-parametric test is used instead of ANOVA. The Kruskal-Wallis test
asks whether the median log2FC differs significantly across bin categories
(retained, adopted, other).

Why chi-square?
---------------
The chi-square test of independence asks whether the proportions of
upregulated, downregulated, and non-significant genes differ across bin
categories. This is complementary to the Kruskal-Wallis test: it captures
discrete DE status rather than the continuous fold change.

Why a genome-wide control?
--------------------------
Adding a genome-wide bar (all genes regardless of bin classification) shows
the background rate of differential expression. If adopted bins show a higher
fraction of upregulated genes than the genome-wide rate, it suggests the
compartment switch has a specific effect on expression.

Output
------
  gene_expression_genomewide.csv    : All genes with expression + DE status
  gene_bin_expression_annotation.csv: Genes in translocated bins with all
                                       classification labels
  {cond}_{col}_DE_fraction.png      : Stacked bar plots per condition and
                                       classification scheme

Usage
-----
    1. Run scripts 10 and 13 first to generate the bin annotation and DE files.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 14_gene_bin_expression.py

Dependencies
------------
    pandas, numpy, matplotlib, scipy
    Install with: pip install pandas numpy matplotlib scipy
    
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal, chi2_contingency


# =============================================================================
# CONFIG – edit all paths and thresholds here before running
# =============================================================================
CONFIG = {
    # Bin classification table produced by script 10
    "bins_file": "/path/to/bins_annotation.csv",

    # GTF annotation file (Ensembl GRCh38 protein-coding genes)
    "gtf_file": "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",

    # RNA-seq TPM expression table (GEO: GSE246689)
    "tpm_file": "/path/to/GSE246689_gene_tpm.tsv",

    # DESeq2 results from script 13
    "de_files": {
        "T1": "/path/to/DE_WT_vs_T1.tsv",
        "C1": "/path/to/DE_WT_vs_C1.tsv",
    },

    # Column names for each replicate in the TPM file
    "tpm_replicates": {
        "WT": ["MCF10AWT_REP1", "MCF10AWT_REP2", "MCF10AWT_REP3"],
        "T1": ["MCF10AT1_REP1", "MCF10AT1_REP2", "MCF10AT1_REP3"],
        "C1": ["MCF10AC1_REP1", "MCF10AC1_REP2", "MCF10AC1_REP3"],
    },

    # Thresholds for calling a gene differentially expressed.
    # A gene is called up/down if padj < padj_threshold AND
    # |log2FC| > lfc_threshold. Otherwise it is labelled 'ns'.
    "lfc_threshold":  1.0,
    "padj_threshold": 0.05,

    # Where to save output plots and tables (created automatically)
    "output_dir": "/path/to/results/gene_bin_expression",
}

CONDITIONS = ["T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def load_gtf_genes(gtf_file: str) -> pd.DataFrame:
    """
    Load a GTF annotation file and return one row per gene with columns:
    chrom, start, end, gene_id, gene_name.

    Only 'gene' feature rows are kept. Coordinates remain 1-based as in GTF.
    """
    gtf = pd.read_csv(
        gtf_file, sep="\t", comment="#", header=None,
        names=["chrom", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"],
    )
    gtf = gtf[gtf["feature"] == "gene"].copy()
    gtf["gene_id"]   = gtf["attributes"].str.extract(r'gene_id "([^"]+)"')
    gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')
    genes = gtf[["chrom", "start", "end", "gene_id", "gene_name"]].copy()
    genes[["start", "end"]] = genes[["start", "end"]].astype(int)
    return genes


def assign_bin(gene_row: pd.Series, bins_df: pd.DataFrame) -> tuple:
    """
    Find the first 100 kb bin that overlaps a gene and return its
    subcompartment classification labels.

    A bin overlaps a gene if any part of it falls within the gene's
    genomic coordinates. When multiple bins overlap, the first one
    (by row order in the bins table) is used.

    Returns a 5-tuple: (transloc_id, condition, behavior,
                         behavior_collapsed, change_direction).
    Returns (NaN, NaN, NaN, NaN, NaN) if no overlapping bin is found.
    """
    overlap = bins_df[
        (bins_df["chrom"] == gene_row["chrom"]) &
        (bins_df["start"]  < gene_row["end"]) &
        (bins_df["end"]    > gene_row["start"])
    ]
    if not overlap.empty:
        row = overlap.iloc[0]
        return (
            row["transloc_id"],
            row["condition"],
            row["behavior"],
            row["behavior_collapsed"],
            row["change_direction"],
        )
    return (np.nan,) * 5


def de_call(log2fc: float, padj: float,
            lfc_thresh: float, padj_thresh: float) -> str:
    """
    Classify a gene as upregulated ('up'), downregulated ('down'), or
    not significant ('ns') based on its fold change and adjusted p-value.

    A gene is only called up or down if BOTH thresholds are met:
    padj < padj_threshold AND |log2FC| > lfc_threshold. This two-threshold
    approach avoids calling statistically significant but biologically
    negligible changes.
    """
    if pd.notna(log2fc) and pd.notna(padj) and padj < padj_thresh:
        if   log2fc >  lfc_thresh: return "up"
        elif log2fc < -lfc_thresh: return "down"
    return "ns"


def run_analysis(
    df_bins: pd.DataFrame,
    df_genome: pd.DataFrame,
    cond: str,
    category_col: str,
    category_order: list[str],
    label: str,
    lfc_thresh: float,
    padj_thresh: float,
    output_dir: str,
) -> None:
    """
    Run Kruskal-Wallis and chi-square tests for one condition and one
    classification scheme, then produce a stacked bar plot.

    Parameters
    ----------
    df_bins    : Genes in translocated bins (with bin labels).
    df_genome  : All genes (for the genome-wide control bar).
    cond       : Condition name ('T1' or 'C1').
    category_col : Column to group by ('behavior', 'behavior_collapsed',
                   or 'change_direction').
    category_order : Expected category values, in display order.
    label      : Human-readable label used in the plot title.
    lfc_thresh, padj_thresh : Thresholds for DE calling.
    output_dir : Where to save the plot.
    """
    df = df_bins[df_bins["condition"] == cond].copy()
    de_col = f"DE_{cond}"
    fc_col = f"log2FC_{cond}_vs_WT"

    # ── Kruskal-Wallis: does log2FC differ across bin categories? ────────────
    groups = [
        df[df[category_col] == cat][fc_col].dropna()
        for cat in category_order if cat in df[category_col].values
    ]
    p_kw = kruskal(*groups).pvalue if len(groups) > 1 else np.nan

    # ── Chi-square: does DE rate differ across bin categories? ───────────────
    ct_counts = pd.crosstab(df[category_col], df[de_col])
    p_chi = chi2_contingency(ct_counts)[1] if ct_counts.shape[0] > 1 else np.nan

    print(f"\n{cond} — {label}")
    print(f"  Kruskal-Wallis p = {p_kw:.3e}")
    print(f"  Chi-square p     = {p_chi:.3e}")
    print(ct_counts)

    # ── Genome-wide control bar ───────────────────────────────────────────────
    ct_genome = pd.crosstab(
        pd.Series(["genome-wide"] * len(df_genome), name=category_col),
        df_genome[de_col],
    )

    print(f"\n  Genome-wide:")
    print(ct_genome)

    # ── Combine, normalise, and plot ──────────────────────────────────────────
    ct_plot = pd.concat([ct_counts, ct_genome])
    ct_frac = ct_plot.div(ct_plot.sum(axis=1), axis=0)
    ct_frac = ct_frac.reindex(columns=["ns", "up", "down"], fill_value=0)

    fig, ax = plt.subplots(figsize=(6, 4))
    ct_frac.plot(
        kind="bar", stacked=True, ax=ax,
        color={"ns": "#d9d9d9", "up": "green", "down": "red"},
    )
    ax.set_ylabel("Fraction of genes")
    ax.set_title(f"{cond}: {label}")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{cond}_{category_col}_DE_fraction.png"),
        dpi=300,
    )
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Load bin classifications ──────────────────────────────────────────────────
print("Loading bin classifications...")
bins = pd.read_csv(CONFIG["bins_file"])
bins[["start", "end"]] = bins[["start", "end"]].astype(int)

# ── Load gene coordinates from GTF ───────────────────────────────────────────
print("Loading GTF gene annotations...")
genes = load_gtf_genes(CONFIG["gtf_file"])

# ── Load TPM expression and compute mean per condition ────────────────────────
print("Loading TPM expression...")
expr = pd.read_csv(CONFIG["tpm_file"], sep="\t")
for cond, reps in CONFIG["tpm_replicates"].items():
    expr[cond] = expr[reps].mean(axis=1)

# Genome-wide gene expression table (all genes, no bin filtering)
gene_expr = genes.merge(expr[["gene_id", "WT", "T1", "C1"]],
                         on="gene_id", how="left")
for cond in ["WT", "T1", "C1"]:
    gene_expr[f"log2_{cond}"] = np.log2(gene_expr[cond] + 1)
for cond in ["T1", "C1"]:
    gene_expr[f"log2FC_{cond}_vs_WT"] = (gene_expr[f"log2_{cond}"]
                                          - gene_expr["log2_WT"])

# ── Load DESeq2 results and classify genes as up/down/ns ─────────────────────
print("Loading DESeq2 results...")
for cond, de_path in CONFIG["de_files"].items():
    de = pd.read_csv(de_path, sep="\t")[["gene_id", "log2FoldChange", "padj"]]
    de = de.rename(columns={
        "log2FoldChange": f"log2FoldChange_{cond}",
        "padj":           f"padj_{cond}",
    })
    gene_expr = gene_expr.merge(de, on="gene_id", how="left")
    gene_expr[f"DE_{cond}"] = gene_expr.apply(
        lambda r, c=cond: de_call(
            r[f"log2FoldChange_{c}"], r[f"padj_{c}"],
            CONFIG["lfc_threshold"], CONFIG["padj_threshold"],
        ),
        axis=1,
    )

gene_expr.to_csv(
    os.path.join(CONFIG["output_dir"], "gene_expression_genomewide.csv"),
    index=False,
)

# ── Map genes to bins separately per condition ────────────────────────────────
# IMPORTANT: the mapping must be done per condition, not across all bins at
# once. If all bins are combined first, genes on chromosomes shared between
# T1 and C1 (chr3, chr6, chr17, chr19) will always be assigned to whichever
# condition's rows appear first in the file, leaving the other condition with
# almost no genes on those chromosomes. Mapping per condition ensures every
# gene is evaluated against that condition's bins only.
print(f"\nMapping {len(genes)} genes to bins per condition...")

bin_labels = ["transloc_id", "condition", "behavior",
              "behavior_collapsed", "change_direction"]

genes_bins_per_cond = []

for cond in CONDITIONS:
    bins_cond = bins[bins["condition"] == cond].copy()
    genes_cond = genes.copy()
    genes_cond[bin_labels] = genes_cond.apply(
        lambda r: pd.Series(assign_bin(r, bins_cond)), axis=1
    )
    genes_cond = genes_cond.dropna(subset=["condition"]).copy()
    print(f"  {cond}: {len(genes_cond)} genes mapped to translocated bins.")

    # Merge expression
    genes_cond = genes_cond.merge(
        expr[["gene_id", "WT", "T1", "C1"]], on="gene_id", how="left"
    )
    for c in ["WT", "T1", "C1"]:
        genes_cond[f"log2_{c}"] = np.log2(genes_cond[c] + 1)
    for c in ["T1", "C1"]:
        genes_cond[f"log2FC_{c}_vs_WT"] = (genes_cond[f"log2_{c}"]
                                             - genes_cond["log2_WT"])

    # Merge DE results
    for c, de_path in CONFIG["de_files"].items():
        de = pd.read_csv(de_path, sep="\t")[["gene_id", "log2FoldChange", "padj"]]
        de = de.rename(columns={
            "log2FoldChange": f"log2FoldChange_{c}",
            "padj":           f"padj_{c}",
        })
        genes_cond = genes_cond.merge(de, on="gene_id", how="left")
        genes_cond[f"DE_{c}"] = genes_cond.apply(
            lambda r, cc=c: de_call(
                r[f"log2FoldChange_{cc}"], r[f"padj_{cc}"],
                CONFIG["lfc_threshold"], CONFIG["padj_threshold"],
            ),
            axis=1,
        )

    genes_bins_per_cond.append(genes_cond)

# Combined bin-annotated table (for saving)
gene_expr_bins = pd.concat(genes_bins_per_cond, ignore_index=True)
gene_expr_bins.to_csv(
    os.path.join(CONFIG["output_dir"], "gene_bin_expression_annotation.csv"),
    index=False,
)
print("Saved annotated gene tables.")

# ── Run analyses and produce plots ────────────────────────────────────────────
print("\nGene counts per condition in bin-annotated table:")
print(gene_expr_bins["condition"].value_counts())

print("\nRunning analyses...")

for cond in CONDITIONS:
    # Use only the bins for this condition
    gene_expr_bins_cond = gene_expr_bins[
        gene_expr_bins["condition"] == cond
    ].copy()

    # Genome-wide control: all genes with valid DE calls for this condition
    df_genome_cond = gene_expr[gene_expr[f"DE_{cond}"].notna()].copy()

    run_analysis(gene_expr_bins_cond, df_genome_cond, cond,
                 "behavior_collapsed", ["retained", "adopted", "other"],
                 "Collapsed A/B compartments",
                 CONFIG["lfc_threshold"], CONFIG["padj_threshold"],
                 CONFIG["output_dir"])

    run_analysis(gene_expr_bins_cond, df_genome_cond, cond,
                 "change_direction", ["A_to_B", "B_to_A", "retained"],
                 "A/B direction",
                 CONFIG["lfc_threshold"], CONFIG["padj_threshold"],
                 CONFIG["output_dir"])

    run_analysis(gene_expr_bins_cond, df_genome_cond, cond,
                 "behavior", ["retained", "adopted", "other"],
                 "Subcompartment behaviour",
                 CONFIG["lfc_threshold"], CONFIG["padj_threshold"],
                 CONFIG["output_dir"])

print("\nAnalysis complete.")
print(f"Plots saved to: {CONFIG['output_dir']}")
