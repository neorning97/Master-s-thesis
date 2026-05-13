# Volcano Plot Scripts

Two scripts for visualising DESeq2 differential expression results as
volcano plots with translocated genes highlighted, and for testing whether
those genes show enriched or directionally biased differential expression
compared to the genome-wide background.

---

## Biological background

A volcano plot is a standard way to visualise differential expression
results. Every gene is shown as a point with its log2 fold change on the
x-axis and -log10(adjusted p-value) on the y-axis. Genes that are both
strongly changed and statistically significant appear in the upper corners.

By overlaying genes from translocated regions in colour on top of all other
genes shown in gray, the plot immediately shows whether translocated genes
cluster among the most significantly changed genes or are scattered
throughout the background distribution.

---

## Two scripts, two approaches

| Script | How genes are identified | What it analyses | Statistical tests |
|--------|-------------------------|-----------------|-------------------|
| `16_volcano_plots.py` | From the gene-pair distance table (scripts 01–04) | All translocations combined | None (visualisation only) |
| `17_translocation_volcano_plots.py` | Directly from user-defined genomic coordinates | Each named translocation separately (Der(17)t(3;17), Der(3)t(3;17)) | Fisher DE enrichment + Fisher Up/Down bias |

Use script 16 for a genome-wide overview of where translocated genes sit
in the DE landscape, and to generate the full gene lists used in GO
enrichment (script 18). Use script 17 when you want to analyse specific
named translocations independently and run statistical tests.

---

## Dependencies on earlier scripts

```
scripts 01–04  -->  T1_distances_agg.tsv, C1_distances_agg.tsv  -->  script 16
script 13      -->  DE_WT_vs_T1.tsv, DE_WT_vs_C1.tsv            -->  scripts 16 and 17
GTF annotation                                                  -->  script 17
script 16      -->  all_up/down_transloc_genes_*.tsv            -->  script 18
```

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib scipy
```

Note: `scipy` is only needed for script 17. Script 16 does not run any
statistical tests.

---

## Input files

### Script 16

| File | Description | Produced by |
|------|-------------|-------------|
| `DE_WT_vs_T1.tsv` | DESeq2 results for WT vs T1 | Script 13 |
| `DE_WT_vs_C1.tsv` | DESeq2 results for WT vs C1 | Script 13 |
| `T1_distances_agg.tsv` | Aggregated gene-pair distances for T1, used to identify translocated genes | Scripts 01–04 |
| `C1_distances_agg.tsv` | Same for C1 | Scripts 01–04 |

### Script 17

| File | Description | Produced by |
|------|-------------|-------------|
| `DE_WT_vs_T1.tsv` | DESeq2 results for WT vs T1 | Script 13 |
| `DE_WT_vs_C1.tsv` | DESeq2 results for WT vs C1 | Script 13 |
| `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` | Ensembl GRCh38 gene annotation | Ensembl |

---

## How to run

### Script 16

Edit the variables in the `CONFIG SECTION` near the top of
`16_volcano_plots.py`:

```python
# DESeq2 result files (from script 13)
T1_DE_FILE = "/path/to/DE_WT_vs_T1.tsv"
C1_DE_FILE = "/path/to/DE_WT_vs_C1.tsv"

# Distance tables (from scripts 01-04), used to identify translocated genes
T1_DIST_FILE = "/path/to/T1_distances_agg.tsv"
C1_DIST_FILE = "/path/to/C1_distances_agg.tsv"

# Thresholds for calling a gene differentially expressed
LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

# How many top genes to label/export
TOP_N_LABELS = 5     # genes labelled on the plot
TOP_N_TABLE = 20     # genes in the top_* TSV tables

# Colours for translocated genes by DE status
de_colors = {
    "Up":        "green",
    "Down":      "red",
    "No change": "blue",
}

# Output folder
OUTPUT_FOLDER = "/path/to/results/volcano_plots"
```

Then run:

```bash
python 16_volcano_plots.py
```

### Script 17

Edit the variables in the `CONFIG SECTION` near the top of
`17_translocation_volcano_plots.py`. The most important setting is
`translocations`, define each translocation as a list of
(chromosome, start, end) coordinate tuples. A gene is included if it
overlaps any of the listed regions:

```python
# DESeq2 result files (from script 13)
T1_DE_FILE = "/path/to/DE_WT_vs_T1.tsv"
C1_DE_FILE = "/path/to/DE_WT_vs_C1.tsv"

# GTF file with gene coordinates
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# Translocation regions to analyse
# Each translocation is a list of (chromosome, start, end) tuples,
# one tuple per genomic segment
translocations = {
    "Der(17)t(3;17)": [
        ("chr3",  0,          58_600_000),
        ("chr17", 22_700_001, 83_257_441),
    ],
    "Der(3)t(3;17)": [
        ("chr3",  58_600_001, 198_295_559),
        ("chr17", 0,          22_700_000),
    ],
}

# Thresholds for calling a gene differentially expressed
LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

# How many top genes to label/export
TOP_N_LABELS = 5
TOP_N_TABLE = 20

# Colours for translocated genes by DE status
de_colors = {
    "Up":        "green",
    "Down":      "red",
    "No change": "blue",
}

# Output folder
OUTPUT_FOLDER = "/path/to/results/translocation_volcano_plots"
```

To add another translocation (e.g. t(6;19)), add a new entry to
`translocations` with the correct chromosomal coordinates. No other
changes are needed.

Then run:

```bash
python 17_translocation_volcano_plots.py
```

---

## Output files

### Script 16

Script 16 exports **four TSV files per condition**, two capped at
`TOP_N_TABLE` rows for quick inspection, and two containing all significant
DE genes for use in GO enrichment analysis (script 18):

| File | Description |
|------|-------------|
| `volcano_WT_vs_T1_highlighted.png` | Volcano plot for T1 with translocated genes highlighted |
| `volcano_WT_vs_C1_highlighted.png` | Volcano plot for C1 |
| `top_up_transloc_genes_WT_vs_T1.tsv` | Top 20 upregulated translocated genes in T1, ranked by log2FC |
| `top_down_transloc_genes_WT_vs_T1.tsv` | Top 20 downregulated translocated genes in T1 |
| `top_up_transloc_genes_WT_vs_C1.tsv` | Same for C1 |
| `top_down_transloc_genes_WT_vs_C1.tsv` | Same for C1 |
| `all_up_transloc_genes_WT_vs_T1.tsv` | **All** upregulated translocated genes in T1 (no cap), input for script 18 |
| `all_down_transloc_genes_WT_vs_T1.tsv` | **All** downregulated translocated genes in T1 (no cap), input for script 18 |
| `all_up_transloc_genes_WT_vs_C1.tsv` | Same for C1 |
| `all_down_transloc_genes_WT_vs_C1.tsv` | Same for C1 |

### Script 17

One PNG and two TSV files per translocation x condition combination. With
two translocations and two conditions, that is 12 files total:

| File | Description |
|------|-------------|
| `volcano_WT_vs_T1_Der17t3_17.png` | Volcano plot for Der(17)t(3;17) in T1 |
| `volcano_WT_vs_C1_Der17t3_17.png` | Same for C1 |
| `volcano_WT_vs_T1_Der3t3_17.png` | Volcano plot for Der(3)t(3;17) in T1 |
| `volcano_WT_vs_C1_Der3t3_17.png` | Same for C1 |
| `top_up_Der17t3_17_WT_vs_T1.tsv` | Top 20 upregulated Der(17) genes in T1 |
| `top_down_Der17t3_17_WT_vs_T1.tsv` | Top 20 downregulated Der(17) genes in T1 |
| ... | Same pattern for Der(3) and for C1 |

### Plot description (both scripts)

Each volcano plot shows:
- **Light gray points**: all genes not in translocated regions
- **Coloured points**: translocated genes by DE status:
  green = upregulated, red = downregulated, blue = not significantly changed
- **Dashed lines**: the log2FC and padj thresholds
- **Gene name labels**: the top 5 most significant up- and downregulated
  translocated genes

### Statistical output (script 17 only)

For each translocation x condition, the script prints two Fisher's exact
test results:

```
Der(17)t(3;17) | WT vs C1
  Translocation genes: ...
    Up: ...
    Down: ...
    No change: ...

  DE enrichment Fisher test:
    Translocation: DE=..., No change=...
    Genome-wide:   DE=..., No change=...
    OR=...,  p=...

  Up/Down bias Fisher test (among DE genes only):
    Translocation: Up=..., Down=...
    Genome-wide:   Up=...,, Down=...
    OR=...,  p=...
```

---

## Design decisions

**Why two separate scripts?**
Script 16 uses the distance table from scripts 01–04 to identify translocated
genes, and is designed as a quick visualisation companion to that pipeline.
It also exports full gene lists for GO enrichment. Script 17 identifies genes
directly from genomic coordinates and runs statistical tests, making it more
suitable for reporting specific translocation results. Keeping them separate
avoids each script becoming too complex.

**Why export both top-N and all-genes tables in script 16?**
The top-N tables (`TOP_N_TABLE = 20` by default) are convenient for quick
inspection and match the genes labelled on the plot. The all-genes tables
contain every significant DE translocated gene with no cap, which is
required for GO enrichment analysis (script 18), GO enrichment needs a
sufficiently large gene list to detect statistically significant
overrepresentation after multiple testing correction. Using only 20 genes
per list is typically too few.

**Why show all genes in gray rather than only translocated genes?**
The full genome-wide distribution in the background provides essential
context. Without it, you cannot judge whether the translocated genes are
unusually significant or simply typical of what you would expect from genes
with similar expression levels.

**Why label by significance rather than by fold change?**
The most statistically robust genes are labelled, not the most extremely
changed ones. A gene with a very large fold change but a high p-value
(e.g. because it has low counts) would be unreliable to highlight. Ranking
by -log10(padj) prioritises the most reproducible signals.

**Why use GTF coordinates in script 17 rather than the distance table?**
Using genomic coordinates is more transparent, any examiner can verify
which genes are included by checking the coordinates against the BED files.
It also allows the two reciprocal derivatives of t(3;17) to be analysed
independently, which would not be possible with the combined distance table.

**Why two Fisher tests in script 17?**
The DE enrichment test asks whether the overall rate of differential
expression is elevated. The Up/Down bias test asks, among DE genes, whether
there is a directional shift toward activation or repression, a more
specific question directly relevant to the compartment adoption hypothesis.

**Handling of missing or zero padj values**
DESeq2 sometimes outputs NA for padj (when a gene is filtered by independent
filtering) or exact zeros for extremely significant genes. NA values are
replaced with 1 (not significant) and zeros are replaced with 1x10^{-300} to
avoid log(0) errors while still placing these genes at the top of the plot.

