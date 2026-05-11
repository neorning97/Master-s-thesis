# 3D Gene-Pair Distance Analysis

Scripts for analysing 3D genomic distances between translocated genes and
their neighbors in Chrom3D polymer structural models (.cmm files).

This analysis is part of a master's thesis investigating whether chromosomal
translocations in MCF10A breast epithelial cells cause spatial reorganisation
of the genome. Specifically, whether translocated genes are brought into
closer 3D proximity to genes on the receiving/origin (inter/intra) chromosome.

---

## Biological background

Chromosomal translocations physically move segments of DNA from one chromosome
to another. This may change which genomic regions a translocated gene sits
near in 3D nuclear space, potentially altering its gene expression by placing
it in a new regulatory environment (a new "compartment neighborhood").

This pipeline quantifies that spatial change by computing Euclidean distances
between gene pairs in 3D structural models, comparing two translocated cell
lines (T1 and C1) against a wildtype (WT) baseline.

---

## Overview

The analysis runs as four scripts in order. Each script saves its output to
disk so the next script can read it as input.

| Script | What it does | Key input | Key output |
|--------|-------------|-----------|------------|
| `01_extract_gene_expression.py` | Finds genes in translocated and neighboring regions, samples matched control genes, attaches RNA-seq TPM expression | Translocation BED, neighbor BED, GTF, TPM table | `T1_genes_expression.tsv`, `C1_genes_expression.tsv` |
| `02_compute_tc_distances.py` | Computes 3D gene-pair distances in T1 and C1 structural models | Output of 01, neighbor BED, CMM files | `T1_distances.tsv`, `T1_distances_agg.tsv`, `C1_distances.tsv`, `C1_distances_agg.tsv` |
| `03_compute_wt_distances.py` | Computes the same distances in wildtype models (the baseline) | Output of 02, GTF, WT CMM files | `WT_distances_agg.tsv` |
| `04_plot_results.py` | Scatter plots + Wilcoxon signed-rank tests comparing WT vs T1/C1 | Output of 02 and 03 | Scatter plots (.png) |

---

## Requirements

**Python ≥ 3.10**

Install all dependencies at once:

```bash
pip install pandas numpy matplotlib seaborn scipy bioframe
```

---

## Input files

The following files are required. Download the publicly available files from
GEO accession **GSE246689** and place them alongside your BED and CMM files.

| File | Description | Source |
|------|-------------|--------|
| `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` | Ensembl GRCh38 protein-coding gene annotation | Ensembl |
| `GSE246689_gene_tpm.tsv` | RNA-seq TPM expression table for MCF10A cell lines | GEO: GSE246689 |
| `T1_translocations.bed` | Genomic regions translocated in condition T1 | This study |
| `T1_neighbor.bed` | neighboring regions on the receiving chromosome for T1 | This study |
| `T1_neighbor_originchr.bed` | neighboring regions on the original chromosome for T1 | This study |
| `C1_translocations.bed` | Genomic regions translocated in condition C1 | This study |
| `C1_neighbor.bed` | neighboring regions on the receiving chromosome for C1 | This study |
| `C1_neighbor_originchr.bed` | neighboring regions on the original chromosome for C1 | This study |
| `cmm/T1/*.cmm` | Chrom3D 3D structural models for condition T1 | GEO: GSE246947 |
| `cmm/C1/*.cmm` | Chrom3D 3D structural models for condition C1 | GEO: GSE246947 |
| `cmm/WT/*.cmm` | Chrom3D 3D structural models for the wildtype | GEO: GSE246947 |

### BED file formats

**Translocation BED** (tab-separated):

```
chrom   start   end   label   transloc_id
```

`chrom` is the chromosome the segment was translocated *from*.

**neighbor BED** (tab-separated):

```
chrom   start   end   label  type   transloc_id   neighbor_chr
```

`neighbor_chr` is the chromosome the segment was translocated *to*
(the receiving chromosome). `chrom` is the origin chromosome.

### CMM file format

Standard UCSF Chimera / Chrom3D marker file. The file contains `<marker>`
elements (one per genomic bead) connected by `<link>` elements (representing
the polymer chain). Only `<marker>` elements are used by these scripts;
`<link>` elements are ignored.

Each `<marker>` must have a `beadID` attribute in the format `chrN:start-end`
and `x`, `y`, `z` coordinates giving its position in 3D space:

```xml
<marker_set name="chrom3d_model">
  <marker id="0" x="-2.08585" y="-1.03324" z="-2.36597" radius="0.2"
          r="0.9" g="0.1" b="0.1" chrID="chr1" beadID="chr1:0-1950000"/>
  <marker id="1" x="-2.148" y="-0.858654" z="-2.72046" radius="0.2"
          r="0.9" g="0.1" b="0.1" chrID="chr1" beadID="chr1:1950000-3450000"/>
  <link id1="0" id2="1" r="0.50996" g="0.516893" b="0.627882" radius="0.1"/>
  ...
</marker_set>
```

---

## How to run

### Step 1 — Edit the CONFIG section in each script

Each script has a clearly marked `CONFIG SECTION` near the top with flat
top-level variables (e.g. `T1_GENE_FILE`, `C1_CMM_FOLDER`, `OUTPUT_FILE`,
etc.). Replace all `"/path/to/..."` placeholders with the actual paths to
your files. Output directories are created automatically.

### Step 2 — Run the scripts in order

```bash
python 01_extract_gene_expression.py
python 02_compute_tc_distances.py
python 03_compute_wt_distances.py
python 04_plot_results.py
```

Each script prints progress to the terminal and reports where output files
have been saved.

### Scenario setting (inter vs intra)

Script 02 has a `SCENARIO` variable at the top of the CONFIG section:

- `"inter"`: keep only gene pairs where the translocated gene and its
  neighbor are on **different chromosomes** (interchromosomal). Use this
  for translocations between chromosomes.
- `"intra"`: keep only pairs on the **same chromosome** (intrachromosomal).
  Use this when analysing neighbors on the origin chromosome.

Change the neighbor BED file accordingly:
- For `"inter"`: use `T1_neighbor.bed` (neighbors on the receiving chr)
- For `"intra"`: use `T1_neighbor_originchr.bed` (neighbors on the origin chr)

---

## Output files

### Results directory

| File | Description |
|------|-------------|
| `T1_genes_expression.tsv` | Genes in T1 regions with TPM expression and region labels |
| `C1_genes_expression.tsv` | Same for C1 |
| `T1_distances.tsv` | Per-model raw distances for all T1 gene pairs (large file) |
| `T1_distances_agg.tsv` | Distances aggregated across all T1 models (mean, std, n_models) |
| `C1_distances.tsv` | Per-model raw distances for C1 |
| `C1_distances_agg.tsv` | Aggregated distances for C1 |
| `WT_distances_agg.tsv` | Aggregated distances for WT (baseline) |

The large files can be accessed [here](https://doi.org/10.5281/zenodo.20126785).

### Plots directory

| File | Description |
|------|-------------|
| `T1_transloc_scatter.png` | WT vs T1 distances — translocated–neighbor pairs |
| `T1_control_scatter.png` | WT vs T1 distances — control–neighbor pairs |
| `C1_transloc_scatter.png` | WT vs C1 distances — translocated–neighbor pairs |
| `C1_control_scatter.png` | WT vs C1 distances — control–neighbor pairs |

In each scatter plot, each point is one gene pair. The red dashed line is
y = x (equal distance in WT and condition). Points **below** the line are
closer in the translocated condition than in WT; points **above** are
further apart.

---

## Design decisions

**Control gene sampling**: For each translocation event, the same number of
control genes as translocated genes is sampled at random from chromosomes not
involved in the translocation. This creates a size-matched background
distribution. The random seed is derived deterministically from the
translocation ID so sampling is reproducible across runs.

**Bead mapping**: Each gene is mapped to the Chrom3D bead whose genomic
midpoint is closest to the gene's midpoint on the same chromosome.

**Memory management**: Script 02 writes each model's results to disk one
model at a time, instead of keeping all results in memory at once. Each
gene's nearest-bead lookup uses simple per-chromosome filtering inside the
loop. `scipy.spatial.distance.cdist` is used to compute all pairwise
distances at once for each translocation, which is much faster than
nested Python loops.

**Statistical test**: The Wilcoxon signed-rank test is used rather than a
paired t-test because 3D distance distributions are not normally distributed.
It tests whether the median pairwise distance is significantly different
between the wildtype and the translocated condition.

**Pair matching**: Before plotting, both the WT and condition distance tables
are restricted to the gene pairs that appear in both (using inner merges in
script 04). This ensures the comparison is always like-for-like.
