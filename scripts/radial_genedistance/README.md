# Nuclear Distance Analysis

Scripts for computing and comparing the **radial nuclear position** of
translocated genes in Chrom3D polymer structural models — i.e. how far each
gene is from the nuclear centre before and after a chromosomal translocation.

This analysis is a companion to the gene-pair distance pipeline
(`01_extract_gene_expression.py` --> `04_plot_results.py`). It uses the same
input files but asks a different question: rather than measuring how close
translocated genes are to their new neighbors, it measures whether the
translocation moves genes toward or away from the centre of the nucleus.

---

## Biological background

In Chrom3D polymer models, the nucleus is modelled as a unit sphere centred
on the origin (0, 0, 0). Genes near the centre tend to be in transcriptionally
active A-compartments, while genes near the nuclear periphery tend to be in
inactive B-compartments.

A chromosomal translocation physically moves a segment of DNA from one
chromosome to another, potentially repositioning those genes in 3D nuclear
space. If translocated genes move inward (toward the centre), they may be
entering a more active nuclear environment, which could explain
transcriptional changes observed after the translocation.

This pipeline tests that hypothesis by comparing the distance-to-centre of
translocated genes (inside_transloc) vs. matched control genes in:
- **T1 and C1** models (the translocated cell lines)
- **WT** models (the wildtype baseline)

---

## Overview

| Script | What it does | Key input | Key output |
|--------|-------------|-----------|------------|
| `05_compute_tc_distance_to_center.py` | Computes distance to nuclear centre for translocated and control genes in T1/C1 models | Gene expression TSV (from script 01), CMM files | `T1_distance_to_center_agg.tsv`, `C1_distance_to_center_agg.tsv` |
| `06_compute_wt_distance_to_center.py` | Computes the same distances in wildtype models (the baseline) | Gene expression TSV (from script 01), GTF, WT CMM files | `WT_distance_to_center_agg.tsv` |
| `07_plot_distance_to_center.py` | Scatter plots + Wilcoxon tests comparing WT vs T1/C1 radial positions | Output of 05 and 06 | Scatter plots (.png) |

These scripts are numbered 05–07 to indicate they follow the gene-pair
distance pipeline (scripts 01–04). Script 01
(`01_extract_gene_expression.py`) from that pipeline must be run first, as
its output (`T1_genes_expression.tsv`, `C1_genes_expression.tsv`) is used
here as the gene list.

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib seaborn scipy
```

---

## Input files

| File | Description | Produced by |
|------|-------------|-------------|
| `T1_genes_expression.tsv` | Genes in T1 regions with region type labels | Script 01 |
| `C1_genes_expression.tsv` | Same for C1 | Script 01 |
| `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` | Ensembl gene annotation (for WT coordinates) | Ensembl |
| `cmm/T1/*.cmm` | Chrom3D structural models for T1 | GEO: GSE246947 |
| `cmm/C1/*.cmm` | Chrom3D structural models for C1 | GEO: GSE246947 |
| `cmm/WT/*.cmm` | Chrom3D structural models for the wildtype | GEO: GSE246947 |

### CMM file format

Standard UCSF Chimera / Chrom3D marker file. Each `<marker>` element
represents one genomic bead. Only `beadID`, `x`, `y`, and `z` are used;
all other attributes and `<link>` elements are ignored.

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

### Step 1 — Run the gene expression extraction script first

This pipeline depends on the gene expression TSV files produced by script 01
of the gene-pair distance pipeline:

```bash
python 01_extract_gene_expression.py
```

### Step 2 — Edit the CONFIG block in each script

Each script has a `CONFIG` dictionary at the top. Replace all
`"/path/to/..."` placeholders with the actual paths to your files.

### Step 3 — Run the scripts in order

```bash
python 05_compute_tc_distance_to_center.py
python 06_compute_wt_distance_to_center.py
python 07_plot_distance_to_center.py
```

---

## Output files

### Results

| File | Description |
|------|-------------|
| `T1_distance_to_center.tsv` | Per-model distance to centre for each gene in T1 (large file) |
| `T1_distance_to_center_agg.tsv` | Distances aggregated across all T1 models (mean, std, n_models) |
| `C1_distance_to_center.tsv` | Per-model distances for C1 |
| `C1_distance_to_center_agg.tsv` | Aggregated distances for C1 |
| `WT_distance_to_center_agg.tsv` | Aggregated distances for WT (baseline) |

### Plots

| File | Description |
|------|-------------|
| `T1_inside_transloc_scatter.png` | WT vs T1 radial position — translocated genes |
| `T1_control_scatter.png` | WT vs T1 radial position — control genes |
| `C1_inside_transloc_scatter.png` | WT vs C1 radial position — translocated genes |
| `C1_control_scatter.png` | WT vs C1 radial position — control genes |

In each scatter plot, each point is one gene. The red dashed line is y = x
(equal distance in WT and condition). Points **below** the line moved closer
to the nuclear centre after the translocation; points **above** moved further
away.

---

## Design decisions

**What does distance to origin represent?**
In Chrom3D, the nucleus is modelled as a sphere of radius 5 µm centred on
the origin (0, 0, 0). The Euclidean distance from a gene's bead to the
origin is therefore a measure of its radial position in the nucleus. Smaller
values indicate a more central location, larger values indicate a position
closer to the nuclear envelope. A distance of 5 µm corresponds to a gene
sitting on the nuclear periphery.

**Control genes**
Control genes are randomly sampled from non-translocated chromosomes (in
script 01) and should show no systematic inward or outward movement between
WT and the translocated condition. If they do, it would suggest a global
structural difference between the models rather than a translocation-specific
effect.

**Statistical test**
The Wilcoxon signed-rank test is used because distance distributions are not
normally distributed. It tests whether the median distance to the nuclear
centre is significantly different between WT and the translocated condition
for each gene group.

**Bead mapping**
Each gene is mapped to the Chrom3D bead whose genomic midpoint is closest
to the gene's midpoint on the same chromosome. The bead index is pre-built
as numpy arrays once per model to avoid slow DataFrame filtering in the inner
loop.
