# Chromosome Centroid Distance Analysis

Script for measuring and comparing the **3D distance between two chromosome
centroids** across wildtype (WT), T1, and C1 Chrom3D structural models.

This is a standalone script (`08_chromosome_centroid_distance.py`) that does
not depend on any of the other analysis scripts in this repository. It only
requires the raw `.cmm` structural model files.


---

## Biological background

Chromosomal translocations physically join segments from two different
chromosomes. For a translocation to occur, the two chromosomes involved must
have been spatially close in the nucleus at the time of the rearrangement.
After the translocation, the two chromosomes may remain in closer spatial
proximity than they would be in a normal cell, because they are now physically
connected.

This script tests that hypothesis by computing the distance between the
**centroids** of two chromosomes (e.g. chr3 and chr17) in each Chrom3D
structural model, and comparing that distance between:
- **WT**: wildtype models, where no translocation has occurred
- **T1** and **C1**: models of the two translocated cell lines

If the translocation brings the two chromosomes closer together, we expect
the centroid distance to be significantly smaller in T1 and/or C1 than in WT.

---

## What is a chromosome centroid?

Each chromosome is represented in Chrom3D as a chain of genomic beads, each
positioned at (x, y, z) coordinates in 3D space. The centroid is simply the
mean x, y, z position across all beads on a chromosome, a single point that
summarises where that chromosome sits in the nucleus on average.

The distance between the centroids of two chromosomes is then the Euclidean
distance between those two summary points, giving a single number per model
that represents how spatially separated the two chromosomes are overall.

---

## Overview

This is a single self-contained script:

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `08_chromosome_centroid_distance.py` | Computes centroid-to-centroid distance for two chromosomes across all WT, T1, and C1 models; produces a box plot and Mann-Whitney U tests | CMM files for WT, T1, C1 | Box plot (.png), printed statistics |

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib seaborn scipy
```

---

## Input files

| File | Description | Source |
|------|-------------|--------|
| `cmm/WT/*.cmm` | Chrom3D structural models for the wildtype | GEO: GSE246947 |
| `cmm/T1/*.cmm` | Chrom3D structural models for condition T1 | GEO: GSE246947 |
| `cmm/C1/*.cmm` | Chrom3D structural models for condition C1 | GEO: GSE246947 |

### CMM file format

Standard UCSF Chimera / Chrom3D marker file. Each `<marker>` element
represents one genomic bead. Only `beadID`, `x`, `y`, and `z` are used;
start/end positions, `<link>` elements, and all other attributes are ignored.

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

### Step 1: Edit the CONFIG section

Open `08_chromosome_centroid_distance.py` and edit the variables in the
`CONFIG SECTION` near the top of the script:

```python
# The two chromosomes to compare, set to the chromosomes involved
# in the translocation you are studying
CHR_A = "chr3"
CHR_B = "chr17"

# Folders containing .cmm files for each condition
WT_FOLDER = "/path/to/cmm/WT"
T1_FOLDER = "/path/to/cmm/T1"
C1_FOLDER = "/path/to/cmm/C1"

# Where to save the output plot
OUTPUT_FOLDER = "/path/to/results/chromosome_centroid_distance"
```

To analyse a different chromosome pair, simply change `CHR_A` and `CHR_B`.

### Step 2: Run the script

```bash
python 08_chromosome_centroid_distance.py
```

---

## Output

### Plot

`chromosome_distance_boxplot.png`: a box plot showing the distribution of
centroid distances across all structural models, one box per condition (WT,
T1, C1). Individual model distances are overlaid as points so the full
distribution is visible.

### Printed statistics

The script prints a summary table of distances per condition and the results
of three pairwise Mann-Whitney U tests:

```
Statistical tests (Mann-Whitney U, two-sided):
  WT vs T1: U=..., p=...  (median WT=..., median T1=...: T1 closer than WT)
  WT vs C1: U=..., p=...
  T1 vs C1: U=..., p=...
```

---

## Design decisions

**Why chromosome centroids rather than minimum distances?**
The centroid gives a single stable summary of where each chromosome sits in
the nucleus, averaged across all its beads. Minimum bead-to-bead distances
would be noisier and more sensitive to individual bead positions at the
chromosome ends.

**Why Mann-Whitney U instead of a t-test or Wilcoxon?**
The distribution of centroid distances across structural models is not
assumed to be normally distributed, which rules out a t-test. Wilcoxon
signed-rank would also be inappropriate because it requires *paired* data,
and there is no natural pairing between an individual WT model and an
individual T1 or C1 model, they are independent structural ensembles. The
Mann-Whitney U test is the correct non-parametric choice for comparing two
unpaired groups, testing whether one condition tends to produce
systematically larger or smaller distances than another without assuming a
particular distribution shape.

**Why are some models skipped?**
If either chromosome has no beads in a particular model (which can happen for
very small chromosomal regions or due to model-specific constraints), that
model is skipped for the distance calculation. The script prints a message
whenever this happens so you can track how many models were used.

**Coordinate scale**
The Chrom3D models for this dataset used a nuclear radius of 5 µm, so
distances are in units of micrometres. A centroid distance of 5 µm would
correspond to the two chromosome centroids sitting on opposite sides of the
nuclear envelope.
