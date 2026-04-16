# Differential Expression and Gene–Bin Expression Analysis

Two scripts for identifying differentially expressed genes and linking
expression changes to chromatin subcompartment behaviour in translocated
genomic regions.

---

## Biological background

Chromosomal translocations may alter gene expression by repositioning genes
into new chromatin environments. Script 10 classified each 100 kb bin in
translocated regions as retained, adopted, or other, reflecting whether the
bin kept its wildtype subcompartment identity or switched to match its new
neighborhood. This pipeline tests whether those compartment switches are
associated with changes in gene expression.

The analysis proceeds in two steps:

1. **Script 13** runs DESeq2 to identify which genes are significantly
   differentially expressed between WT and each translocated cell line
   (T1, C1).

2. **Script 14** maps genes to their overlapping 100 kb bins, inherits
   the bin's subcompartment classification, and tests whether genes in
   adopted bins are more likely to be differentially expressed than genes
   in retained or other bins.

---

## Overview

| Script | Language | What it does | Input | Output |
|--------|----------|-------------|-------|--------|
| `13_differential_expression.R` | R | Runs DESeq2 to compare WT vs T1 and WT vs C1 | Raw gene counts, sample metadata | `DE_WT_vs_T1.tsv`, `DE_WT_vs_C1.tsv` |
| `14_gene_bin_expression.py` | Python | Maps genes to bins per condition, loads DE results, runs Kruskal-Wallis and chi-square tests, produces bar plots | Bin annotation (script 10), GTF, TPM file, DE results (script 13) | Annotated gene tables + stacked bar plots |

Run script 13 before script 14.

---

## Requirements

### Script 13 — R ≥ 4.1

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Script 14 — Python ≥ 3.10

```bash
pip install pandas numpy matplotlib scipy
```

---

## Input files

| File | Description | Source |
|------|-------------|--------|
| `GSE246689_gene_counts.tsv` | Raw gene count matrix (genes x samples) | GEO: GSE246689 |
| `coldata.tsv` | Sample metadata with `condition` column (WT, T1, C1) | This study |
| `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` | Ensembl GRCh38 gene annotation | Ensembl |
| `GSE246689_gene_tpm.tsv` | RNA-seq TPM expression table | GEO: GSE246689 |
| `bins_annotation.csv` | Bin classifications from script 10 | Script 10 |
| `DE_WT_vs_T1.tsv` | DESeq2 results from script 13 | Script 13 |
| `DE_WT_vs_C1.tsv` | DESeq2 results from script 13 | Script 13 |

### coldata.tsv format

Tab-separated, sample IDs as row names, must include a `condition` column:

```
                condition
MCF10AWT_REP1   WT
MCF10AWT_REP2   WT
MCF10AWT_REP3   WT
MCF10AT1_REP1   T1
MCF10AT1_REP2   T1
MCF10AT1_REP3   T1
MCF10AC1_REP1   C1
MCF10AC1_REP2   C1
MCF10AC1_REP3   C1
```

---

## How to run

### Script 13

Edit the CONFIG section at the top of `13_differential_expression.R`:

```r
COUNTS_FILE  <- "/path/to/GSE246689_gene_counts.tsv"
COLDATA_FILE <- "/path/to/coldata.tsv"
OUTPUT_DIR   <- "/path/to/results/DE"
MIN_COUNT    <- 10   # minimum total counts to include a gene
```

Then run:

```bash
Rscript 13_differential_expression.R
```

### Script 14

Edit the CONFIG section at the top of `14_gene_bin_expression.py`:

```python
CONFIG = {
    "bins_file":      "/path/to/results/bins_annotation.csv",
    "gtf_file":       "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",
    "tpm_file":       "/path/to/GSE246689_gene_tpm.tsv",
    "de_files": {
        "T1": "/path/to/results/DE_WT_vs_T1.tsv",
        "C1": "/path/to/results/DE_WT_vs_C1.tsv",
    },
    "lfc_threshold":  1.0,    # |log2FC| must exceed this to be called DE
    "padj_threshold": 0.05,   # adjusted p-value threshold
    "output_dir":     "/path/to/results/gene_bin_expression",
}
```

Then run:

```bash
python 14_gene_bin_expression.py
```

---

## Output files

### Script 13

| File | Description |
|------|-------------|
| `DE_WT_vs_T1.tsv` | DESeq2 results for WT vs T1, one row per gene |
| `DE_WT_vs_C1.tsv` | DESeq2 results for WT vs C1 |

Columns: `gene_id`, `gene_name`, `baseMean`, `log2FoldChange`
(apeglm-shrunken), `lfcSE`, `svalue`, `padj`.

### Script 14

| File | Description |
|------|-------------|
| `gene_expression_genomewide.csv` | All genes with expression values, log2FC, and DE status — used as the genome-wide control |
| `gene_bin_expression_annotation.csv` | Genes in translocated bins with bin classification labels and DE status, one row per gene per condition |
| `T1_behavior_collapsed_DE_fraction.png` | Fraction up/down/ns per collapsed A/B category in T1 |
| `T1_change_direction_DE_fraction.png` | Fraction up/down/ns per A→B / B→A direction in T1 |
| `T1_behavior_DE_fraction.png` | Fraction up/down/ns per full subcompartment category in T1 |
| `C1_*.png` | Same three plots for C1 |

Each bar plot includes a **genome-wide control bar** showing the background
rate of differential expression across all genes, providing a reference for
whether the rates in adopted/retained/other bins are elevated or suppressed.

---

## Design decisions

**Why DESeq2 with apeglm shrinkage?**
Raw log2 fold changes from DESeq2 are noisy for genes with low counts, a
gene with 3 counts in one condition and 1 in another will have a large raw
fold change but very little statistical support. apeglm shrinkage pulls
these unreliable estimates toward zero, producing more conservative fold
changes better suited for downstream analysis and visualisation.

**Why a dual threshold for DE calling (|log2FC| > 1 AND padj < 0.05)?**
Using only a p-value threshold would call statistically significant but
biologically negligible changes as differentially expressed. The fold change
threshold ensures that only genes with at least a twofold expression change
are classified as up or down. Both thresholds can be adjusted in the CONFIG
section of script 14.

**Why is the gene-to-bin mapping done per condition?**
The `bins_annotation.csv` file contains rows for both T1 and C1. If genes
are mapped against all bins at once, genes on chromosomes shared between T1
and C1 (chr3, chr6, chr17, chr19) will always be assigned to whichever
condition's rows appear first in the file, leaving the other condition with
almost no genes on those chromosomes. By mapping separately for each
condition, filtering the bins table to T1 rows only, then C1 rows only,
every gene is evaluated independently for each condition, producing fair and
comparable gene counts.

**Why Kruskal-Wallis rather than ANOVA?**
Log2 fold changes are not normally distributed across genes, violating the
assumptions of ANOVA. The Kruskal-Wallis test is a non-parametric alternative
that tests whether the median fold change differs significantly across bin
categories (retained, adopted, other) without assuming normality.

**Why chi-square in addition to Kruskal-Wallis?**
The two tests capture different aspects of the data. Kruskal-Wallis compares
the continuous distribution of fold changes. Chi-square compares the discrete
proportions of up/down/ns genes. It is possible for fold change distributions
to differ significantly while proportions remain similar, or vice versa, so
both are reported.

**Why a genome-wide control bar?**
The genome-wide bar shows the background rate of differential expression
across all detected genes, regardless of which bin they fall in. If adopted
bins show a higher fraction of upregulated genes than the genome-wide average,
it suggests compartment adoption specifically promotes gene activation rather
than the elevation being a general property of T1 or C1.

**Only verified translocations are analysed**
Script 14 inherits its translocation regions from `bins_annotation.csv`,
which was produced by script 10 using the verified translocation BED files
(`verify_T1_translocations.bed`, `verify_C1_translocations.bed`). Only genes
overlapping bins in those verified regions are included in the bin-stratified
analysis. The genome-wide control covers all genes regardless.

**First-overlap bin assignment**
When a gene overlaps multiple 100 kb bins within a condition, the first
overlapping bin (by row order in the bins table) is used. This is a
simplification that works well in practice because most genes are smaller
than one 100 kb bin.
