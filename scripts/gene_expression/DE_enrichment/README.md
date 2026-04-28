# Differential Expression Enrichment Analysis

Script for testing whether genes in translocated genomic regions are more
likely to be differentially expressed than the rest of the genome, and
whether they show a directional bias toward upregulation.

---

## Biological background

Chromosomal translocations physically move segments of DNA into new genomic
neighborhoods, which may alter their gene expression by placing them in
new regulatory environments. If the translocation has a specific effect on
expression, we would expect the translocated genes to show higher rates of
differential expression than genes from elsewhere in the genome.

This script tests two related hypotheses:

1. **DE enrichment**: Are genes in translocated regions more often
   differentially expressed (up or down) than all other genes in the genome?

2. **Up/Down bias**: Among differentially expressed genes, are translocated
   genes more often upregulated than the genome background? This would be
   consistent with translocation into a more active chromatin environment.

The control group is simply all genes in the DESeq2 results table that are
not in any translocated region, i.e. the whole-genome background.

---

## Overview

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `15_de_enrichment.py` | Compares DE rates between translocated genes and the genome background; runs two Fisher's exact tests; produces grouped bar plots | Distance table (scripts 01–04), DESeq2 results (script 13) | One bar plot per condition + printed statistics |

---

## Dependencies on earlier scripts

```
scripts 01–04  -->  WT_distances_agg.tsv  ----
                                              |-->  script 15
script 13      -->  DE_WT_vs_T1.tsv           |
               -->  DE_WT_vs_C1.tsv       ----
```

Run scripts 01–04 and script 13 before running this script.

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib scipy
```

---

## Input files

| File | Description | Produced by |
|------|-------------|-------------|
| `WT_distances_agg.tsv` | Aggregated gene-pair distances with region type labels, used only to extract translocated gene names | Scripts 01–04 |
| `DE_WT_vs_T1.tsv` | DESeq2 results for WT vs T1 | Script 13 |
| `DE_WT_vs_C1.tsv` | DESeq2 results for WT vs C1 | Script 13 |

### Distance table format

Tab-separated. Must include:

```
gene1   gene1_region_type    chr1  	gene2  	gene2_region_type  	chr2  	mean_distance  	std_distance  	n_models
```

Only rows where `gene1_region_type == "inside_transloc"` are used, to
build the list of translocated gene names.

---

## How to run

### Step 1 — Edit the CONFIG section

Edit the variables in the `CONFIG SECTION` near the top of
`15_de_enrichment.py`:

```python
# Distance table from scripts 01-04 (used only to get translocated gene names)
DIST_FILE = "/path/to/results/WT_distances_agg.tsv"

# DESeq2 result files (from script 13)
T1_DE_FILE = "/path/to/results/DE_WT_vs_T1.tsv"
C1_DE_FILE = "/path/to/results/DE_WT_vs_C1.tsv"

# Thresholds for calling a gene differentially expressed
LFC_THRESHOLD = 1.0       # |log2FoldChange| must exceed this
PADJ_THRESHOLD = 0.05     # adjusted p-value threshold

# Output folder
OUTPUT_FOLDER = "/path/to/results/de_enrichment"
```

### Step 2 — Run

```bash
python 15_de_enrichment.py
```

---

## Output

### Plots

| File | Description |
|------|-------------|
| `T1_grouped_DE_proportions.png` | Fraction of Up/Down/No change genes, translocated vs genome background in T1 |
| `C1_grouped_DE_proportions.png` | Same for C1 |

Each plot shows two bars per category (Up, Down, No change): blue for
translocated genes and orange for the genome background. Proportions are
shown rather than raw counts so the two groups, which differ greatly in
size, can be compared visually on the same scale.

### Printed statistics

For each condition:

```
T1
  Total genes in DE table:  ....
  Translocated genes:       ....
  Genome background genes:  ....

  Counts — translocation:       {'Up': ..., 'Down': ..., 'No change': ...}
  Counts — genome background:   {'Up': ..., 'Down': ..., 'No change': ...}

  DE enrichment test (translocation vs genome background):
    Table: [[DE_trans, notDE_trans], [DE_bg, notDE_bg]]
    Odds ratio: ...,  p = ...

  Up/Down bias test (among DE genes only):
    Table: [[Up_trans, Down_trans], [Up_bg, Down_bg]]
    Odds ratio: ...,  p = ...
```

---

## Design decisions

**Why Fisher's exact test rather than chi-square, Wilcoxon, or Mann-Whitney U?**
The analysis compares two gene groups using a 2x2 contingency table of
counts (DE vs not DE, or Up vs Down). This is fundamentally count data in
categories, so Wilcoxon and Mann-Whitney U don't apply, those tests work
on continuous distributions, not on counts of items falling into discrete
categories. Chi-square would also be appropriate in principle, but it
requires sufficiently large expected counts in all cells, and the
translocated gene group is typically small. Fisher's exact
test is exact (computes p-values without that approximation) and gives
trustworthy results even when one of the groups is small.

**Why two tests?**
The DE enrichment test (DE vs No change) and the Up/Down bias test answer
different biological questions. A significant DE enrichment result means
translocated genes are generally more affected by the translocation. A
significant Up/Down bias means that among the affected genes, there is a
directional shift, e.g., more activation than repression, which would be
consistent with translocation into a more active chromatin environment.

**Why the whole genome as a control?**
The control is all genes in the DESeq2 results table that are not in any
translocated region. This makes the result directly interpretable: a
significant result means translocated genes behave differently from the
rest of the genome, not just differently from a specific random sample.

**Duplicate gene names**
If the same gene name appears multiple times in the DESeq2 results (e.g. due
to duplicate Ensembl ID mappings), only the row with the most significant
adjusted p-value is kept. This ensures each gene is counted only once.

**DE classification thresholds**
A gene is called Up or Down only if both `|log2FoldChange| > 1.0` AND
`padj < 0.05` are satisfied. The fold change threshold ensures that only
biologically meaningful changes are counted. Both thresholds can be adjusted
in the CONFIG section.
