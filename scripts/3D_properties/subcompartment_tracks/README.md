# Subcompartment Track Visualisation

Script for visualising chromatin subcompartment states along the derivative
chromosomes formed by all translocations in the T1 and C1 cell lines.

---

## Biological background

A chromosomal translocation fuses portions of two chromosomes into a new
derivative chromosome. To understand how the translocated segment's chromatin
environment changes, this script draws horizontal subcompartment tracks
showing the subcompartment state along the full length of each derivative
chromosome, displaying both the translocated fragment and the chromosomal
body it joined, side by side.

Subcompartments (A0, A1, A2, A3, B0, B1, B2, B3) reflect the transcriptional
and epigenetic state of 100 kb genomic bins. A-compartment regions (red shades)
are transcriptionally active; B-compartment regions (blue shades) are
inactive. By comparing the tracks across conditions (WT, T1, C1), you can
see whether the subcompartment identity of a translocated segment changes
after the translocation and whether it begins to resemble its new neighborhood.

---

## Overview

This is a single self-contained script that produces one figure per
translocation per applicable condition.

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `11_subcompartment_tracks.py` | Draws subcompartment tracks for the full original chromosomes and the assembled derivative for each translocation x condition | Subcompartment bedGraph, translocation BED files | One PNG per translocation per condition |

---

## Translocations covered

| Translocation | Conditions plotted | Fragment | Body |
|---------------|--------------------|---------|------|
| Der(17)t(3;17) | WT, T1, C1 | chr3: 0–58.6 Mb | chr17: 22.7–83.3 Mb |
| Der(3)t(3;17)  | WT, T1, C1 | chr17: 0–22.7 Mb | chr3: 58.6–198.3 Mb |
| t(6;19)        | WT, T1, C1 | chr19: 0–31.9 Mb | chr6: 0–170.8 Mb |
| t(2;10)        | WT, C1 only | chr2: 0–131.7 Mb | chr10: 39.8–133.8 Mb |

t(2;10) is only plotted for WT and C1 because it does not occur in T1.
Breakpoint coordinates are read directly from the translocation BED files,
so they do not need to be hardcoded in the script.

---

## Requirements

**Python ≥ 3.10**

```bash
pip install pandas matplotlib
```

---

## Input files

| File | Description | Source |
|------|-------------|--------|
| `GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph` | 100 kb subcompartment annotations for WT, T1, and C1 | GEO: GSE246947 |
| `verify_T1_translocations.bed` | Genomic coordinates of translocated segments in T1 | This study |
| `verify_C1_translocations.bed` | Genomic coordinates of translocated segments in C1 | This study |

### Translocation BED format

Tab-separated with one row per segment. Paired rows share a `transloc_id`.
The shorter row in each pair is automatically identified as the fragment;
the longer row is the body. Example:

```
chrom   start       end         label                   transloc_id
chr3    0           58600000    T1_transloc_1           T1
chr17   22700001    83257441    T1_transloc_1_partner   T1
chr3    58600001    198295559   T1_transloc_2           T2
chr17   0           22700000    T1_transloc_2_partner   T2
chr6    0           170805979   T1_transloc_3           T3
chr19   0           31900000    T1_transloc_3_partner   T3
```

### Subcompartment bedGraph format

Tab-separated with one row per 100 kb bin. Must include columns named
`MCF10A_WT.state`, `MCF10A_T1.state`, and `MCF10A_C1.state` containing
subcompartment labels such as A0, A1, A2, A3, B0, B1, B2, B3.

---

## How to run

### Step 1: Edit the CONFIG section

Open `11_subcompartment_tracks.py` and edit the variables in the
`CONFIG SECTION` near the top of the script:

```python
# Subcompartment file
SUBCOMPARTMENT_FILE = "/path/to/subcompartments.bedGraph"

# Translocation BED files
T1_TRANSLOC_BED = "/path/to/verify_T1_translocations.bed"
C1_TRANSLOC_BED = "/path/to/verify_C1_translocations.bed"

# Maps transloc_id values in the BED to human-readable names
# used in plot titles and output filenames
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

# Which conditions to plot for each translocation name
transloc_conditions = {
    "t(2;10)":        ["WT", "C1"],
    "Der(17)t(3;17)": ["WT", "T1", "C1"],
    "Der(3)t(3;17)":  ["WT", "T1", "C1"],
    "t(6;19)":        ["WT", "T1", "C1"],
}

# Output folder
OUTPUT_FOLDER = "/path/to/results/subcompartment_tracks"
```

### Step 2: Run the script

```bash
python 11_subcompartment_tracks.py
```

The script prints each translocation and condition as it is processed, so
you can track progress.

---

## Output files

One PNG per translocation per condition, named automatically:

| File | Description |
|------|-------------|
| `subcompartment_Der17t3_17_WT.png` | Der(17)t(3;17) in WT |
| `subcompartment_Der17t3_17_T1.png` | Der(17)t(3;17) in T1 |
| `subcompartment_Der17t3_17_C1.png` | Der(17)t(3;17) in C1 |
| `subcompartment_Der3t3_17_*.png` | Same for the reciprocal Der(3) |
| `subcompartment_t6_19_*.png` | t(6;19) for WT, T1, C1 |
| `subcompartment_t2_10_WT.png` | t(2;10) in WT |
| `subcompartment_t2_10_C1.png` | t(2;10) in C1 |

Each figure contains four panels:

1. **Fragment chromosome (original)**: full chromosome coloured by
   subcompartment state. Shows the pre-translocation environment of the
   segment that was moved.
2. **Body chromosome (original)**: full chromosome coloured by
   subcompartment state. Shows the pre-translocation environment of the
   region the fragment lands next to.
3. **Derivative chromosome**: the fragment segment joined to the body
   segment. A solid green vertical line marks the breakpoint, the junction
   between the two chromosomal regions.
4. **Legend**: colour key for all subcompartment states.

---

## Design decisions

**Why show the full original chromosomes alongside the derivative?**
The derivative track alone is not enough to interpret subcompartment changes.
Showing the full original chromosomes lets you see what the pre-translocation
environment looked like for both the translocated fragment and its new
neighbors, making it easy to judge whether any changes in T1/C1 are driven
by the translocation.

**Why is t(2;10) only plotted for WT and C1?**
The t(2;10) translocation does not occur in T1. Plotting the T1 track for
t(2;10) coordinates would show the unaffected T1 genome at those positions
and would be misleading. The `transloc_conditions` setting in CONFIG controls
which conditions are plotted per translocation.

**Colour scheme**
A-compartment subcompartments are shown in shades of red (A0 = light pink,
A3 = dark red) and B-compartment subcompartments in shades of blue (B0 =
light blue, B3 = dark navy). Any unrecognised state is shown in grey.
