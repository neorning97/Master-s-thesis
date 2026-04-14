# 3D Gene-Pair Distance Analysis

Scripts for analysing 3D genomic distances between translocated genes and
their neighbors in Hi-C / polymer structural models (.cmm files).

Written as part of a master's thesis.  

---

## Overview

The analysis is split into four scripts that must be run **in order**. Each
script produces output files that the next script reads as input.

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `01_extract_gene_expression.py` | Finds genes in translocated/neighbor regions, samples controls, attaches TPM expression | BED files, GTF, expression TSV | `T1_genes_expression.tsv`, `C1_genes_expression.tsv` |
| `02_compute_tc_distances.py` | Computes 3D distances between gene pairs in T1/C1 structural models | Output of 01, neighbor BED, CMM files | `T1_distances.tsv`, `T1_distances_agg.tsv`, (same for C1) |
| `03_compute_wt_distances.py` | Computes the same distances in wildtype models (the baseline) | Output of 02, GTF, WT CMM files | `WT_distances_agg.tsv` |
| `04_plot_results.py` | Scatter plots + Wilcoxon tests comparing WT vs T1/C1 | Output of 02 and 03 | Scatter plots (.png) |

---

## Why separate scripts?

The 3D structural models contain tens of thousands of genomic beads, and with
~4000–5000 translocated genes per condition, the distance matrices are large.
Running all four steps in a single Python process would require more RAM than
is available on a standard laptop (~8 GB). By running each script separately,
Python starts with a fresh memory allocation each time.

---

## Requirements

Python ≥ 3.10

Install all dependencies at once:

```bash
pip install pandas numpy matplotlib seaborn scipy bioframe
```

---

## Input files

| File | Description |
|------|-------------|
| `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` | Ensembl protein-coding gene annotation (GRCh38) |
| `GSE246689_gene_tpm.tsv` | RNA-seq TPM expression table (GEO: GSE246689) |
| `T1_translocations.bed` | Genomic regions translocated in condition T1 |
| `T1_neighbor1.bed` | Neighboring regions on the receiving chromosome for T1 |
| `C1_translocations.bed` | Genomic regions translocated in condition C1 |
| `C1_neighbor1.bed` | Neighboring regions on the receiving chromosome for C1 |
| `cmm/T1/*.cmm` | 3D structural models for condition T1 |
| `cmm/C1/*.cmm` | 3D structural models for condition C1 |
| `cmm/WT/*.cmm` | 3D structural models for WT |

### BED file format

**Translocation BED** – tab-separated, must include at least:

```
chrom   start   end   label   transloc_id
```

**Neighbor BED** – tab-separated, must include at least:

```
neighbor_chr   start   end   label   transloc_id   chrom
```

Note: `neighbor_chr` is the chromosome the neighbor sits on (the receiving
chromosome), while `chrom` is the origin chromosome of the translocation.

### CMM file format

Standard UCSF Chimera marker file. Each `<marker>` element must have:
- A `beadID` attribute in the format `chrN:start-end`
- `x`, `y`, `z` coordinate attributes

Example:
```xml
<marker id="1" beadID="chr1:0-200000" x="1.23" y="4.56" z="7.89" />
```

---

## How to run

### Step 1 – Edit the CONFIG section in each script

Each script has a `CONFIG` dictionary at the top. Replace the
`"/path/to/..."` placeholders with the actual paths to your files.
The `results_dir` and `plots_dir` settings are where output files will be
saved — these are created automatically if they do not exist.

### Step 2 – Run the scripts in order

```bash
python 01_extract_gene_expression.py
python 02_compute_tc_distances.py
python 03_compute_wt_distances.py
python 04_plot_results.py
```

Each script prints progress to the terminal and reports where output files
have been saved.

---

## Output files

### Results directory

| File | Description |
|------|-------------|
| `T1_genes_expression.tsv` | Genes in T1 regions with TPM expression values |
| `C1_genes_expression.tsv` | Genes in C1 regions with TPM expression values |
| `T1_distances.tsv` | Per-model gene-pair distances for T1 (one row per model per pair) |
| `T1_distances_agg.tsv` | Distances aggregated across all T1 models (mean, std, n) |
| `C1_distances.tsv` | Per-model gene-pair distances for C1 |
| `C1_distances_agg.tsv` | Distances aggregated across all C1 models |
| `WT_distances_agg.tsv` | Distances aggregated across all WT models |

### Plots directory

| File | Description |
|------|-------------|
| `T1_transloc_scatter.png` | WT vs T1 distances for translocated–neighbor pairs |
| `T1_control_scatter.png` | WT vs T1 distances for control–neighbor pairs |
| `C1_transloc_scatter.png` | WT vs C1 distances for translocated–neighbor pairs |
| `C1_control_scatter.png` | WT vs C1 distances for control–neighbor pairs |

In each scatter plot, each point is one gene pair. The red dashed line is
y = x (equal distance in WT and condition). Points **below** the line are
closer in the translocated condition than in WT; points **above** are further
apart.

---

## Design decisions

**Control gene sampling** — For each translocation event, the same number of
control genes as translocated genes is sampled at random from chromosomes not
involved in the translocation. This creates a size-matched background
distribution for fair statistical comparison. The random seed is derived
deterministically from the translocation ID, so results are reproducible.

**Bead mapping** — Each gene is mapped to the structural bead whose genomic
midpoint is closest to the gene's midpoint on the same chromosome.

**Scenario filtering** — After computing all distances, the script keeps only
interchromosomal pairs (translocated chr ≠ neighbor chr) because the
translocations in this study are between chromosomes. Change `"scenario"` to
`"intra"` in the CONFIG to analyse intrachromosomal contacts instead.

**Statistical test** — The Wilcoxon signed-rank test is used rather than a
paired t-test because 3D distance distributions are not normally distributed.
It tests whether the median pairwise distance is significantly different
between the wildtype and the translocated condition.

**Memory management** — Script 02 writes results to disk one model at a time
rather than accumulating all models in memory. This keeps RAM usage
proportional to the size of one model, not all models combined.

---

## Reproducing the analysis

1. Clone or download this folder.
2. Install dependencies: `pip install pandas numpy matplotlib seaborn scipy bioframe`
3. Download input files from GEO accession **GSE246689** and place them as
   described in the *Input files* section above.
4. Edit the `CONFIG` block at the top of each script.
5. Run scripts 01 → 02 → 03 → 04 in order.
