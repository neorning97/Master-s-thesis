# Gene Ontology Enrichment Analysis

R script for running Gene Ontology (GO) enrichment analysis on
differentially expressed translocated genes and producing dot plots showing
the most significantly enriched biological processes.

---

## Biological background

Differential expression analysis identifies *which* genes change in
expression, but not *what biological processes* those changes affect.
GO enrichment analysis answers that question: given a list of differentially
expressed genes, which biological processes, molecular functions, or cellular
components are represented more often than expected by chance?

This script focuses on **Biological Process (BP)**, the GO ontology that
describes higher-level biological pathways such as cell cycle regulation,
apoptosis, or immune response. This is the most commonly used ontology in
transcriptomic studies because it is most directly interpretable in a
biological context.

---

## Overview

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `18_go_enrichment.R` | Runs GO enrichment on up- and downregulated translocated genes; produces dot plots | Ranked gene tables (scripts 16) | One dot plot per condition + printed top GO terms |

---

## Dependencies on earlier scripts

```
scripts 16        -->  all_up_transloc_genes_WT_vs_T1.tsv   ---
                  -->  all_down_transloc_genes_WT_vs_T1.tsv ---|
                  -->  all_up_transloc_genes_WT_vs_C1.tsv   ---|-->  script 18
                  -->  all_down_transloc_genes_WT_vs_C1.tsv ---
```

---

## Requirements

**R ≥ 4.1**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
install.packages("ggplot2")
```

---

## Input files

| File | Description | Produced by |
|------|-------------|-------------|
| `all_up_transloc_genes_WT_vs_T1.tsv` | Upregulated translocated genes in T1 | Scripts 16 |
| `all_down_transloc_genes_WT_vs_T1.tsv` | Downregulated translocated genes in T1 | Scripts 16  |
| `all_up_transloc_genes_WT_vs_C1.tsv` | Same for C1 | Scripts 16 |
| `all_down_transloc_genes_WT_vs_C1.tsv` | Same for C1 | Scripts 16 |

Each input file must contain a `gene_name` column with HGNC gene symbols.

---

## How to run

### Step 1 — Edit the CONFIG section

Edit the variables in the `CONFIG SECTION` near the top of
`18_go_enrichment.R`:

```r
# Input gene lists (from scripts 16 or 17)
T1_UP_FILE   <- "/path/to/all_up_transloc_genes_WT_vs_T1.tsv"
T1_DOWN_FILE <- "/path/to/all_down_transloc_genes_WT_vs_T1.tsv"
C1_UP_FILE   <- "/path/to/all_up_transloc_genes_WT_vs_C1.tsv"
C1_DOWN_FILE <- "/path/to/all_down_transloc_genes_WT_vs_C1.tsv"

# Output folder for the dot plots
OUTPUT_DIR <- "/path/to/results/go_enrichment"

# How many top GO terms to show in each plot
TOP_N_CATEGORIES <- 10

# Which GO ontology to use: "BP", "MF", or "CC"
GO_ONTOLOGY <- "BP"

# Adjusted p-value cutoff for calling a GO term enriched
PADJ_CUTOFF <- 0.05
```

### Step 2 — Run

From RStudio: open the script and click run.

---

## Output files

| File | Description |
|------|-------------|
| `T1_GO_BP_dotplot.png` | Dot plot of top 10 enriched GO Biological Process terms in T1 |
| `C1_GO_BP_dotplot.png` | Same for C1 |

### Reading the dot plot

Each row in the dot plot is one GO term. The columns show:

- **x-axis (GeneRatio)**: the fraction of your input genes that belong to
  this GO term. A higher ratio means more of your genes are involved in
  this process.
- **Dot size**: the number of input genes in the GO term.
- **Dot colour**: the adjusted p-value (padj). Darker/more saturated
  colour = more significant.

GO terms are ordered by gene ratio by default, so the most represented
processes appear at the top.

---

## Design decisions

**Why combine up and downregulated genes?**
Combining up- and downregulated genes gives an overview of which biological
processes are affected by the translocation overall.

**Why Biological Process (BP)?**
The three GO ontologies are Biological Process (BP), Molecular Function
(MF), and Cellular Component (CC). BP describes higher-level pathways
(e.g. cell cycle, apoptosis) which are more interpretable in the context
of cancer biology. MF and CC can be selected by changing `GO_ONTOLOGY`
in CONFIG.

**Why Benjamini-Hochberg (BH) correction?**
Thousands of GO terms are tested simultaneously, so many would appear
significant by chance without correction. BH correction controls the false
discovery rate (FDR) at 5%, meaning at most 5% of called enrichments are
expected to be false positives. This is the standard approach for GO
enrichment in transcriptomic studies.

**Why convert to Entrez IDs?**
clusterProfiler's `enrichGO` function requires Entrez IDs rather than
HGNC gene symbols. The `bitr()` function handles the conversion using
`org.Hs.eg.db`, the human gene annotation database. A small fraction of
gene symbols may not map to Entrez IDs, this is normal and those genes
are silently dropped.

**What is the background gene set?**
The background (universe) for the hypergeometric test is the full set of
human genes in `org.Hs.eg.db`. This is the default in clusterProfiler
and is appropriate when your gene list comes from a genome-wide analysis.
