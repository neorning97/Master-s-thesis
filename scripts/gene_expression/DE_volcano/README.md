# Volcano Plots with Translocated Genes Highlighted

Script for visualising DESeq2 differential expression results as volcano
plots, with genes in translocated regions highlighted in colour and the
most significant translocated genes labelled by name.

---

## Biological background

A volcano plot is a standard way to visualise the results of a differential
expression analysis. Every gene is shown as a point, with its log2 fold
change on the x-axis and the statistical significance (-log10 adjusted
p-value) on the y-axis. Genes that are strongly and significantly changed
appear in the upper corners of the plot.

By overlaying the genes from translocated regions in colour on top of all
other genes shown in grey, the plot immediately shows whether translocated
genes are enriched among the most significantly differentially expressed
genes, or whether their expression is largely unchanged despite the genomic
rearrangement.

---

## Overview

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `16_volcano_plots.py` | Plots all genes as a volcano, highlights translocated genes by DE status, labels top significant genes, exports ranked gene tables | DESeq2 results (script 13), distance tables (scripts 01–04) | One PNG per condition + two ranked TSV tables per condition |

---

## Dependencies on earlier scripts

```
scripts 01–04  -->  T1_distances_agg.tsv, C1_distances_agg.tsv  --
                                                                     |-->  script 16
script 13      -->  DE_WT_vs_T1.tsv, DE_WT_vs_C1.tsv            --
```

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib
```

---

## Input files

| File | Description | Produced by |
|------|-------------|-------------|
| `DE_WT_vs_T1.tsv` | DESeq2 results for WT vs T1 | Script 13 |
| `DE_WT_vs_C1.tsv` | DESeq2 results for WT vs C1 | Script 13 |
| `T1_distances_agg.tsv` | Aggregated gene-pair distances for T1, used to identify translocated genes | Scripts 01–04 |
| `C1_distances_agg.tsv` | Same for C1 | Scripts 01–04 |

---

## How to run

### Step 1 — Edit the CONFIG section

```python
CONFIG = {
    "de_files": {
        "T1": "/path/to/DE_WT_vs_T1.tsv",
        "C1": "/path/to/DE_WT_vs_C1.tsv",
    },
    "dist_files": {
        "T1": "/path/to/T1_distances_agg.tsv",
        "C1": "/path/to/C1_distances_agg.tsv",
    },
    "lfc_threshold":  1.0,   # |log2FC| threshold for DE calling
    "padj_threshold": 0.05,  # adjusted p-value threshold
    "top_n_labels":   5,     # number of top genes to label on the plot
    "top_n_table":    20,    # number of top genes to export to TSV
    "colors": {
        "Up":        "green",
        "Down":      "red",
        "No change": "blue",
    },
    "output_dir": "/path/to/results/volcano_plots",
}
```

### Step 2 — Run

```bash
python 16_volcano_plots.py
```

---

## Output files

| File | Description |
|------|-------------|
| `volcano_WT_vs_T1_highlighted.png` | Volcano plot for T1 with translocated genes highlighted |
| `volcano_WT_vs_C1_highlighted.png` | Volcano plot for C1 with translocated genes highlighted |
| `top_up_transloc_genes_WT_vs_T1.tsv` | Top 20 upregulated translocated genes in T1, ranked by log2FC |
| `top_down_transloc_genes_WT_vs_T1.tsv` | Top 20 downregulated translocated genes in T1, ranked by log2FC |
| `top_up_transloc_genes_WT_vs_C1.tsv` | Same for C1 |
| `top_down_transloc_genes_WT_vs_C1.tsv` | Same for C1 |

### Plot description

Each volcano plot shows:
- **Light grey points** — all genes not in translocated regions
- **Coloured points** — translocated genes, coloured by DE status:
  - Green = upregulated (log2FC > 1, padj < 0.05)
  - Red = downregulated (log2FC < −1, padj < 0.05)
  - Blue = not significantly changed
- **Dashed lines** — the log2FC and padj thresholds
- **Gene name labels** — the top 5 most significant upregulated and
  downregulated translocated genes

---

## Design decisions

**Why show all genes in grey rather than only translocated genes?**
Showing the full genome-wide distribution in the background provides
essential context. Without it, you cannot judge whether the translocated
genes are unusually significant or simply typical of what you would expect
from genes with similar expression levels.

**Why label by significance (-log10 padj) rather than by fold change?**
The most statistically robust genes are labelled, not the most extremely
changed ones. A gene with a very large fold change but a high p-value
(e.g. because it has low counts) would be unreliable to highlight. Ranking
by -log10(padj) prioritises the most reproducible signals.

**Why export separate ranked tables?**
The top labelled genes on the plot (top_n_labels = 5) are intentionally few
to keep the plot readable. The exported TSV tables allow you to inspect a
larger set (top_n_table = 20) without cluttering the figure. The tables are
ranked by fold change magnitude rather than significance, which is more
useful for biological follow-up.

**Handling of missing or zero padj values**
DESeq2 sometimes outputs NA for padj (when a gene is filtered by independent
filtering) or exact zeros for extremely significant genes. NA values are
replaced with 1 (not significant) and zeros are replaced with 1x10^{300} to
avoid log(0) errors while still placing these genes at the top of the plot.
