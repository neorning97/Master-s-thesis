# Subcompartment Alluvial Plots

R script for visualising chromatin subcompartment transitions across
WT --> T1 --> C1 (or WT --> C1 for t(2;10)) as alluvial diagrams.

Two types of plot are produced for each translocation: one per individual
genomic segment (fragment or body), and one combining both segments of the
translocation into a single overview plot.

---

## Biological background

A chromosomal translocation may cause genomic bins to switch chromatin
subcompartment, moving from an active A-compartment state to an inactive
B-compartment state, or vice versa. A bar chart can only show the marginal
distribution of subcompartment states per condition, but cannot reveal which
bins changed state or where they went.

An alluvial diagram makes transitions explicit. Each vertical block (stratum)
represents one subcompartment state in one condition. The ribbons connecting
strata show how bins flow between states. The width of each ribbon is
proportional to the number of 100 kb bins following that particular transition
path across WT --> T1 --> C1.

---

## Overview

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `12_subcompartment_alluvial.R` | For every genomic segment in the BED files, and for each whole translocation, counts subcompartment transitions and draws an alluvial diagram | Subcompartment bedGraph, translocation BED files | One JPEG per segment + one JPEG per whole translocation |

---

## Plots produced

### Part 1 — Per-segment plots (one per BED row)

| Segment | Conditions |
|---------|-----------|
| chr3 fragment – Der(17)t(3;17) | WT, T1, C1 |
| chr17 body – Der(17)t(3;17) | WT, T1, C1 |
| chr3 body – Der(3)t(3;17) | WT, T1, C1 |
| chr17 fragment – Der(3)t(3;17) | WT, T1, C1 |
| chr6 body – t(6;19) | WT, T1, C1 |
| chr19 fragment – t(6;19) | WT, T1, C1 |
| chr2 fragment – t(2;10) | **WT and C1 only** |
| chr10 body – t(2;10) | **WT and C1 only** |

### Part 2 — Whole-translocation plots (fragment + body combined)

| Translocation | Conditions |
|--------------|-----------|
| Der(17)t(3;17) | WT, T1, C1 |
| Der(3)t(3;17) | WT, T1, C1 |
| t(6;19) | WT, T1, C1 |
| t(2;10) | **WT and C1 only** |

t(2;10) only appears in C1, not in T1. Plotting T1 for t(2;10) coordinates
would show the unaffected genome at those positions, which would be
misleading, so T1 is excluded. This is controlled by the
`WT_C1_ONLY` setting in the CONFIG section.

Segments that appear in both the T1 and C1 BED files (e.g. the shared
t(3;17) and t(6;19) regions) are only plotted once.

---

## Requirements

**R ≥ 4.1**

```r
install.packages(c("tidyverse", "ggalluvial"))
```

---

## Input files

| File | Description | Source |
|------|-------------|--------|
| `GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph` | 100 kb subcompartment annotations for WT, T1, and C1 | GEO: GSE246947 |
| `verify_T1_translocations.bed` | Segment coordinates for T1 | This study |
| `verify_C1_translocations.bed` | Segment coordinates for C1 | This study |

### Translocation BED format

Tab-separated with one row per segment. Paired rows share a `transloc_id`.
Must include columns: `chrom`, `start`, `end`, `label`, `transloc_id`.

```
chrom   start       end         label                   transloc_id
chr3    0           58600000    T1_transloc_1           T1
chr17   22700001    83257441    T1_transloc_1_partner   T1
```

The `label` column values must match the keys in `LABEL_NAMES` in CONFIG.

### Subcompartment bedGraph format

Tab-separated with one row per 100 kb bin. Must include columns
`MCF10A_WT.state`, `MCF10A_T1.state`, and `MCF10A_C1.state` containing
subcompartment labels such as A0, A1, A2, A3, B0, B1, B2, B3.

---

## How to run

### Step 1 — Edit the CONFIG section

Open `12_subcompartment_alluvial.R` and edit the variables in the
`CONFIG SECTION` near the top of the script:

```r
# Subcompartment file
INPUT_FILE <- "/path/to/subcompartments.bedGraph"

# BED files for translocations, one for each condition
T1_BED <- "/path/to/verify_T1_translocations.bed"
C1_BED <- "/path/to/verify_C1_translocations.bed"

# Maps from "label" column in the BED files to readable segment names
LABEL_NAMES <- c(
  "T1_transloc_1"         = "chr3 fragment – Der(17)t(3;17)",
  "T1_transloc_1_partner" = "chr17 body – Der(17)t(3;17)",
  ...
)

# Maps from transloc_id to readable translocation name
T1_NAME_MAP <- c(
  "T1" = "Der(17)t(3;17)",
  "T2" = "Der(3)t(3;17)",
  "T3" = "t(6;19)"
)
C1_NAME_MAP <- c(
  "T1" = "t(2;10)",
  "T2" = "Der(17)t(3;17)",
  "T3" = "Der(3)t(3;17)",
  "T4" = "t(6;19)"
)

# Translocations to plot for WT and C1 only (skipping T1)
WT_C1_ONLY <- c("t(2;10)")

# Output folder
OUTPUT_DIR <- "/path/to/results/subcompartment_alluvial"
```

Things to check before running:
- **`LABEL_NAMES`**: maps each `label` value in the BED files to a
  human-readable segment name used in plot titles and filenames. Update
  if your BED files use different label values.
- **`T1_NAME_MAP` / `C1_NAME_MAP`**: map each `transloc_id` to a
  whole-translocation name (e.g. `"T1"` --> `"Der(17)t(3;17)"`). Must
  match the `transloc_id` column in your BED files.
- **`WT_C1_ONLY`**: lists translocation names that should only be plotted
  for WT and C1. Currently set to `c("t(2;10)")`.

### Step 2 — Run

```bash
Rscript 12_subcompartment_alluvial.R
```

Or from RStudio: open the script and click **Source**.

The script prints progress for each segment as it is processed.

---

## Output files

All files are saved to `OUTPUT_DIR`.

**Per-segment plots:**

| File | Segment |
|------|---------|
| `alluvial_chr3_fragment_Der_17_t_3_17.jpeg` | chr3 fragment |
| `alluvial_chr17_body_Der_17_t_3_17.jpeg` | chr17 body |
| `alluvial_chr3_body_Der_3_t_3_17.jpeg` | chr3 body |
| `alluvial_chr17_fragment_Der_3_t_3_17.jpeg` | chr17 fragment |
| `alluvial_chr6_body_t_6_19.jpeg` | chr6 body |
| `alluvial_chr19_fragment_t_6_19.jpeg` | chr19 fragment |
| `alluvial_chr2_fragment_t_2_10.jpeg` | chr2 fragment (WT + C1) |
| `alluvial_chr10_body_t_2_10.jpeg` | chr10 body (WT + C1) |

**Whole-translocation plots:**

| File | Translocation |
|------|--------------|
| `alluvial_Der_17_t_3_17_whole_translocation.jpeg` | Der(17)t(3;17) |
| `alluvial_Der_3_t_3_17_whole_translocation.jpeg` | Der(3)t(3;17) |
| `alluvial_t_6_19_whole_translocation.jpeg` | t(6;19) |
| `alluvial_t_2_10_whole_translocation.jpeg` | t(2;10) (WT + C1) |

Each plot shows:
- **Strata** coloured by subcompartment state, red shades for A-compartment
  (A0 = light pink --> A3 = dark red), blue shades for B-compartment
  (B0 = light blue --> B3 = dark navy)
- **Ribbons** showing how bins transition between states across conditions,
  scaled by the number of bins following each path
- **Reversed y-axis** so A-compartment (active) appears at the top
- **Two axes** (WT, C1) for t(2;10) plots; **three axes** (WT, T1, C1) for all others

---

## Design decisions

**Why two types of plot?**
Per-segment plots let you examine whether the fragment or the body region
drives any subcompartment switching independently. The whole-translocation
plot gives an overview of the entire event and is easier to include as a
single thesis figure.

**How are whole-translocation regions built?**
For each `transloc_id`, the two paired BED rows (fragment and body) are
grouped by chromosome and their coordinate ranges are merged. Bins from
both chromosome regions are then pooled before counting transitions.

**Why is t(2;10) plotted for WT and C1 only?**
The t(2;10) translocation does not occur in T1. Including T1 would show the
unaffected genome at those coordinates, which would be misleading. The
`WT_C1_ONLY` setting lists translocation names that should use a two-axis
(WT + C1) plot. The check uses `fixed()` string matching rather than a regex
pattern, because the parentheses in `"t(2;10)"` would otherwise be
interpreted as regex capture groups, causing the match to silently fail.

**Colour scheme**
Matches scripts 10 and 11: A-compartment in shades of red, B-compartment
in shades of blue. Any state not in the palette appears in gray.
