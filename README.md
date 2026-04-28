# Chromosomal Translocations and Their Impact on 3D Genome Organization and Gene Regulation in Breast Cancer Progression

**Author:** Nadia Ørning
**Institution:** University of Oslo
**Year:** 2026

---

## Overview

This repository contains all code, BED files, and Chrom3D pipeline files associated with the master's thesis:

> *Chromosomal Translocations and Their Impact on 3D Genome Organization and Gene Regulation in Breast Cancer Progression*

The thesis investigates whether chromosomal translocations in breast epithelial cells alter the spatial organization of the genome in the nucleus, and whether any such changes are reflected in gene expression. Using Hi-C-derived 3D structural models (Chrom3D), subcompartment annotations, and RNA-seq data, the study follows an in vitro model of breast cancer progression across three cell lines:

- **WT**: wildtype MCF10A breast epithelial cells
- **T1**: premalignant cells
- **C1**: malignant cells

The analyses cover: spatial distances between translocated genomic regions, radial nuclear positioning of translocated genes, chromatin subcompartment identity, and differential gene expression.

---

## Biological background

Breast cancer is the most prevalent cancer in Norwegian women, with 4215 new cases identified in 2024. A key driver of cancer progression is genomic instability, which produces both point mutations and large-scale structural rearrangements such as chromosomal translocations, the exchange of genetic material between two non-homologous chromosomes.

Beyond their classical effects (oncogenic fusion proteins or enhancer hijacking), translocations may disrupt the three-dimensional organization of the genome in the nucleus. The human genome is arranged hierarchically: chromatin loops, topologically associating domains (TADs), large-scale A/B compartments (active and inactive), and chromosome territories. These structures coordinate gene regulation by controlling which regulatory elements can physically contact which genes. Translocation-induced disruptions to this hierarchy could rewire long-range regulatory relationships and drive aberrant gene expression in cancer.

This thesis uses 3D genome models reconstructed from Hi-C data to quantify how translocations in MCF10A cell lines affect nuclear spatial organization and gene expression across wildtype, premalignant, and malignant states.

---

## Repository structure

```
Master-s-thesis/
|-- scripts/                    All analysis scripts (Python and R)
|   |-- 3D_models/              Script 08: chromosome centroid distance
|   |-- 3D_properties/          Scripts 10–12: subcompartment classification,
|   |   |                       tracks, and alluvial plots
|   |   |-- subcompartment_tracks/   Script 11
|   |   |-- alluvial/                Script 12
|   |-- gene_expression/        Scripts 13–18: differential expression,
|   |   |                       GO enrichment, volcano plots
|   |   |-- DE_enrichment/           Script 15
|   |   |-- DE_volcano/              Scripts 16 and 17
|   |   |-- GO_analysis/             Script 18
|   |   |-- gene_bin_expression/     Scripts 13 and 14
|   |-- genepair/               Scripts 01–04: 3D gene-pair distances
|   |-- radial_chromosome/      Script 09: nuclear radial positioning
|   |-- radial_distance/        Scripts 05–07: distance to nuclear centre
|
|-- BED/                        Translocation and neighbor BED files,
|                               sample metadata (coldata.tsv)
|
|-- chrom3d/                    Nextflow pipeline for running Chrom3D
    |-- bin/
    |-- conf/
    |-- containers/
    |-- samplesheets/
    |-- subworkflows/
    |   |-- local/
    |   |   |-- chrom3d/
    |   |   |-- preprocessing/
    |   |   |-- samplesheet/
    |   |-- nchg/
    |-- custom_resources.config
    |-- main.nf
    |-- nextflow.config
```

---

## Analysis scripts

Each subfolder in `scripts/` contains one or more numbered Python or R
scripts, a README file explaining how to run them, and a `plots/`
subfolder containing the output figures. Scripts are numbered in the order
they should be run.

| Script | Language | What it does |
|--------|----------|-------------|
| 01 | Python | Extract gene expression for translocated regions |
| 02 | Python | Compute 3D distances between translocated and neighbor genes (T1, C1) |
| 03 | Python | Compute the same distances in wildtype models |
| 04 | Python | Plot and test gene-pair distance results |
| 05 | Python | Compute distance to nuclear centre for translocated genes (T1, C1) |
| 06 | Python | Compute distance to nuclear centre in wildtype models |
| 07 | Python | Plot and test nuclear distance results |
| 08 | Python | Compute chromosome centroid distances across WT, T1, C1 |
| 09 | Python | Nuclear radial positioning of translocated segments |
| 10 | Python | Subcompartment classification (retained / adopted / other) |
| 11 | Python | Subcompartment track visualisation |
| 12 | R      | Subcompartment alluvial plots |
| 13 | R      | Differential expression analysis (DESeq2) |
| 14 | Python | Link subcompartment behaviour to gene expression |
| 15 | Python | DE enrichment test (Fisher's exact), translocated vs genome |
| 16 | Python | Genome-wide volcano plots with translocated genes highlighted |
| 17 | Python | Translocation-specific volcano plots (Der(17), Der(3)) |
| 18 | R      | Gene Ontology enrichment analysis |

Detailed documentation for each script is in the README file inside its
subfolder.

---

## Data sources

All raw data files are publicly available and should be downloaded before
running the scripts.

| Dataset | Description | Source |
|---------|-------------|--------|
| Hi-C contact matrices | Merged .mcool files for WT, T1, C1 (`GSE246947_hg38_MCF10A_WT_merged.mcool`, `_T1_merged.mcool`, `_C1_merged.mcool`) | [GSE246947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246947) |
| Genomic domain annotations | 50 kb resolution BED files for Chrom3D modeling (`*.domains.50000.bed.gz`) | [GSE246947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246947) |
| Subcompartment annotations | 100 kb resolution bedGraph for compartment analysis (`GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph.gz`) | [GSE246947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246947) |
| RNA-seq counts | Raw gene count matrix (`GSE246689_gene_counts.tsv`) | [GSE246689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246689) |
| RNA-seq TPM | TPM expression table (`GSE246689_gene_tpm.tsv`) | [GSE246689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246689) |
| Gene annotation | Protein-coding gene GTF derived from GENCODE GRCh38.p13. Download the full GTF from [GENCODE]((https://www.gencodegenes.org/human/release_43.html)) and filter to protein-coding genes (see below) | [GENCODE](https://www.gencodegenes.org/human/release_43.html) |
| 3D structural models | Chrom3D .cmm model files for WT, T1, C1 (10 models each) | [WT](https://doi.org/10.5281/zenodo.18299145), [T1](https://doi.org/10.5281/zenodo.18299102), [C1](https://doi.org/10.5281/zenodo.18299039) |

### Generating the protein-coding gene GTF

The file `Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf` used in
scripts 03, 06, and 17 is not a direct download, it was derived from
the GENCODE GRCh38.p13 annotation by filtering to protein-coding genes
only. To reproduce it, download the full annotation and run:

```bash
# Download the full GENCODE GRCh38.p13 GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip gencode.v43.annotation.gtf.gz

# Filter to protein-coding genes only
grep -v "^#" gencode.v43.annotation.gtf \
  | awk '$3 == "gene"' \
  | grep 'gene_type "protein_coding"' \
  > Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf
```

---

## Dependencies

### Python scripts (scripts 01–11, 14–17)

```bash
pip install pandas numpy matplotlib seaborn scipy bioframe
```

Python ≥ 3.10 is required.

### R scripts (scripts 12, 13, 18)

```r
install.packages(c("tidyverse", "ggalluvial", "ggplot2"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
```

R ≥ 4.1 is required.

### Chrom3D pipeline

The `chrom3d/` folder contains a Nextflow pipeline for running Chrom3D,
based on the [chrom3d-nf](https://github.com/robomics/chrom3d-nf) workflow
by Rossini R. Full setup and configuration instructions are available in
that repository's README.

The pipeline requires either Docker or Apptainer/Singularity to run
containers. The following input files from GEO accession
[GSE246947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246947)
must be downloaded and listed in the samplesheet at
`chrom3d/samplesheets/samplesheet.tsv`:

- Hi-C contact matrices: `GSE246947_hg38_MCF10A_WT_merged.mcool`,
  `_T1_merged.mcool`, `_C1_merged.mcool`
- Genomic domain annotations at 50 kb resolution: `*.domains.50000.bed.gz`
  (one file per condition)

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Run the pipeline (from inside the chrom3d/ folder)
nextflow run . \
  --sample_sheet samplesheets/samplesheet.tsv \
  --output-dir results/ \
  --number_of_models 10 \
  -resume \
  -c custom_resources.config \
  -with-apptainer    # or -with-docker
```

The pipeline was run on the SAGA HPC cluster (NRIS). A custom resource
configuration file (`custom_resources.config`) is included in the
`chrom3d/` folder to increase memory allocation for NCHG filtering steps
and extend runtime for Chrom3D simulations. All other parameters were
left at workflow defaults.

Computed models are available on Zenodo and do not need to be regenerated
unless you want to modify the pipeline parameters:
- WT models: https://doi.org/10.5281/zenodo.18299145
- T1 models: https://doi.org/10.5281/zenodo.18299102
- C1 models: https://doi.org/10.5281/zenodo.18299039

---

## How to reproduce the analysis

1. Clone the repository:
   ```bash
   git clone https://github.com/neorning97/Master-s-thesis.git
   cd Master-s-thesis
   ```

2. Download the raw data from GEO (accessions above) and place the files
   in a local data directory.

3. Edit the `CONFIG` block at the top of each script to point to your
   local file paths. All scripts have a clearly labelled CONFIG section,
   no changes are needed outside of that section.

4. Run the scripts in numerical order (01 --> 18), as earlier scripts produce
   output files that are used as input by later scripts. Each script's
   README documents which earlier scripts must be run first.
