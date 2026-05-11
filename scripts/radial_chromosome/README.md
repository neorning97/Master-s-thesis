# Nuclear Radial Positioning Analysis

Script for tracking the **radial nuclear position** of translocated genomic
segments across wildtype (WT), T1, and C1 Chrom3D structural models, and
testing whether translocation moves those segments toward or away from the
nuclear centre.

---

## Biological background

In Chrom3D polymer models, the nucleus is modelled as a sphere of radius
5 µm centred on the origin (0, 0, 0). Genomic regions near the centre tend
to be in transcriptionally active A-compartments, while regions near the
nuclear periphery tend to be in inactive B-compartments.

When a chromosomal translocation moves a segment of DNA onto a new
chromosome, it may also reposition that segment radially in the nucleus,
either by pulling it into the new chromosome's nuclear neighbourhood, or by
dragging it away from its original environment.

This script tests two specific predictions arising from compartment analysis
of the t(3;17) translocation:

- **Der(17)** carries the chr3-derived segment and is expected to move toward
  the nuclear **interior** (smaller radial distance) from WT → T1 → C1,
  consistent with adoption of an active A-compartment environment.
- **Der(3)** carries the chr17-derived segment and is expected to move toward
  the nuclear **periphery** (larger radial distance) from WT → T1 → C1,
  consistent with adoption of an inactive B-compartment environment.

---

## Overview

This is a single self-contained script that produces two plots and a
statistics table.

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `09_nuclear_radial_positioning.py` | Tracks radial position of translocated segments across WT, T1, C1; produces trajectory plots and runs Mann-Whitney U tests | CMM files (WT, T1, C1), translocation BED files | 2 plots + `radial_positioning_stats.tsv` |

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas numpy matplotlib scipy
```

---

## Input files

| File | Description | Source |
|------|-------------|--------|
| `cmm/WT/*.cmm` | Chrom3D structural models for the wildtype | GEO: GSE246947 |
| `cmm/T1/*.cmm` | Chrom3D structural models for condition T1 | GEO: GSE246947 |
| `cmm/C1/*.cmm` | Chrom3D structural models for condition C1 | GEO: GSE246947 |
| `verify_T1_translocations.bed` | Genomic coordinates of translocated segments in T1 | This study |
| `verify_C1_translocations.bed` | Genomic coordinates of translocated segments in C1 | This study |

### Translocation BED format

Tab-separated, must include at least:

```
chrom   start   end   label   transloc_id
```

`transloc_id` values (e.g. T1, T2, T3) must match the keys in the
`T1_NAME_MAP` and `C1_NAME_MAP` settings.

### CMM file format

Standard UCSF Chimera / Chrom3D marker file. Each `<marker>` element
represents one genomic bead. Only `beadID`, `x`, `y`, and `z` are used;
`<link>` elements and all other attributes are ignored.

```xml
<marker_set name="chrom3d_model">
  <marker id="0" x="-2.08585" y="-1.03324" z="-2.36597" radius="0.2"
          r="0.9" g="0.1" b="0.1" chrID="chr1" beadID="chr1:0-1950000"/>
  <link id1="0" id2="1" r="0.50996" g="0.516893" b="0.627882" radius="0.1"/>
  ...
</marker_set>
```

---

## How to run

### Step 1 — Edit the CONFIG section

Open `09_nuclear_radial_positioning.py` and edit the variables in the
`CONFIG SECTION` near the top of the script:

```python
# Folders containing the 3D model files (.cmm) for each condition
WT_CMM_FOLDER = "/path/to/cmm/WT"
T1_CMM_FOLDER = "/path/to/cmm/T1"
C1_CMM_FOLDER = "/path/to/cmm/C1"

# BED files describing the translocated segments
T1_BED_FILE = "/path/to/verify_T1_translocations.bed"
C1_BED_FILE = "/path/to/verify_C1_translocations.bed"

# Maps transloc_id values in the BED files to readable segment names
T1_NAME_MAP = {
    "T1": "Der(17)t(3;17)",
    "T2": "Der(3)t(3;17)",
    "T3": "t(6;19)",
}
C1_NAME_MAP = {
    "T1": "t(2;10)",
    "T2": "Der(17)t(3;17)",
    "T3": "Der(3)t(3;17)",
    "T4": "t(6;19)",
}

# One colour per segment label (hex codes)
segment_colors = {
    "Der(3)t(3;17)":             "#e41a1c",
    "Der(17)t(3;17)":            "#2171b5",
    "t(6;19)":                   "#756bb1",
    "t(2;10)":                   "#31a354",
    "chr3":                      "#e41a1c",
    "chr17":                     "#2171b5",
    "chr3 fragment in Der(17)":  "#fd8d3c",
    "chr17 fragment in Der(3)":  "#31a354",
}

# Where to save all the output files
OUTPUT_FOLDER = "/path/to/results/plots/nuclear_positioning"
```

The keys in `T1_NAME_MAP` and `C1_NAME_MAP` must match the `transloc_id`
column values in your BED files. Check these by opening the BED file and
reading the `transloc_id` column.

### Step 2 — Run the script

```bash
python 09_nuclear_radial_positioning.py
```

The script prints progress and a sanity check showing which chromosomes and
coordinate ranges are present in the first CMM file per condition, so you can
verify that the BED coordinates match the model content before waiting for
all models to be processed.

---

## Output files

### Plots

| File | Description |
|------|-------------|
| `trajectory_radial_distance.png` | Mean radial distance per segment across WT → T1 → C1, with ±SEM error bars. One line per translocation segment. |
| `trajectory_chr3_chr17_and_fragments.png` | Same as the trajectory plot but showing whole chr3 and chr17 alongside the small translocated fragments, to contextualise the segment-level shifts. |

### Data files

| File | Description |
|------|-------------|
| `per_bead_radial_distances.csv` | Raw per-bead radial distances for every bead, model, and condition (large file) |
| `radial_positioning_stats.tsv` | Mann-Whitney U test results for all segment x condition-pair comparisons |

---

## Design decisions

**What does radial distance represent?**
Each bead's radial distance is its Euclidean distance from the origin (0, 0, 0),
which is the nuclear centre in Chrom3D models. The nuclear radius was set to
5 µm for this dataset, so distances are in micrometres. Smaller values indicate
a more central (typically active) nuclear position; larger values indicate a
more peripheral (typically inactive) position.

**Why use per-bead distances and then average per model?**
A translocated segment may span many genomic beads (e.g. a 10 Mb segment
corresponds to roughly 50–100 beads at 200 kb resolution). Averaging across
all beads in a segment gives one representative radial position per model,
which is then used for statistical testing. This avoids pseudoreplication
(treating each bead as an independent observation when they are all from the
same model).

**Why Mann-Whitney U instead of a t-test or Wilcoxon signed-rank?**
The radial distance distributions are not assumed to be normal, which rules
out a t-test. Mann-Whitney U is the appropriate non-parametric test for
comparing two independent samples without assuming normality. We treat the
conditions as independent because each Chrom3D model is generated as an
independent Monte Carlo simulation with a unique random seed, so model 0 in
WT has no natural correspondence to model 0 in T1, even though they share
a file index. Wilcoxon signed-rank would require a meaningful pairing
between the groups, which is not the case here. This is also consistent with
the test choice used in `08_chromosome_centroid_distance.py` for comparing
chromosome centroid distances across the same three conditions.

**Whole chromosomes in plot 2**
Including whole-chromosome radial positions alongside the small translocated
fragments allows you to see whether the fragment's radial shift is larger or
smaller than the background shift of the entire chromosome. If a fragment
moves more than the full chromosome, the effect is likely driven by the
translocation specifically.

**Fragment coordinates come from the BED files**
The chr3 and chr17 fragment segments (`chr3 fragment in Der(17)` and
`chr17 fragment in Der(3)`) are extracted directly from the translocation
BED files rather than being defined manually in CONFIG. This ensures there
is a single source of truth for the fragment coordinates and avoids
duplication.
