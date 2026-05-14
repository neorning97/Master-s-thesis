# Chrom3D Pipeline

This folder contains the Nextflow pipeline used to generate the 3D
structural models analysed in this thesis. The pipeline is based on
[chrom3d-nf](https://github.com/robomics/chrom3d-nf) by Rossini, R.

## What is in this folder

All files are inherited from the upstream `chrom3d-nf` pipeline, with
one customization: `custom_resources.config`. This file increases memory
allocation for the NCHG filtering steps (32 GB) and extends the maximum
runtime for the `CHROM3D:SIMULATE` step to 7 days. All other parameters
were left at workflow defaults.

## How to run

The pipeline was run on the SAGA HPC cluster (NRIS) with:

\`\`\`bash
nextflow run . \
  --sample_sheet samplesheets/samplesheet.tsv \
  --output-dir results/ \
  --number_of_models 10 \
  -resume \
  -c custom_resources.config \
  -with-apptainer
\`\`\`

Required inputs (download from GEO accession
[GSE246947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246947)
and list in `samplesheets/samplesheet.tsv`):

- Hi-C contact matrices at 50 kb resolution (`*_merged.mcool`)
- Domain annotations (`*.domains.50000.bed.gz`)

In practice, the three conditions were run sequentially rather than in
a single batch, with the samplesheet edited to contain one row at a time
to manage HPC resource allocation. The combined samplesheet committed
here produces equivalent results.

## Output

The pipeline produces 10 `.cmm` model files per condition. The exact
models used in this thesis are also available on Zenodo:

- WT: https://doi.org/10.5281/zenodo.18299145
- T1: https://doi.org/10.5281/zenodo.18299102
- C1: https://doi.org/10.5281/zenodo.18299039

## Reference

See the Chrom3D section of the thesis methods chapter for the full
description of the pipeline, the NCHG model, and the modeling parameters.
