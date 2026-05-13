# Subcompartment Classification Analysis

Script for classifying genomic bins in verified translocated regions as **retained**,
**adopted**, or **other** relative to the wildtype chromatin subcompartment
state, and testing whether those fractions differ significantly from the
genome-wide background.

---

## Biological background

Chromatin subcompartments (A1, A2, A3, B1, B2, B3) reflect the transcriptional
and epigenetic state of genomic regions. A-compartment regions are generally
transcriptionally active, while B-compartment regions are generally inactive.

When a chromosomal translocation moves a segment of DNA onto a new chromosome,
it places that segment in a new chromatin environment. There are three possible
outcomes for each 100 kb bin in the translocated region:

- **Retained**: the bin keeps the same subcompartment state it had in the
  wildtype. The translocation did not change its chromatin identity.
- **Adopted**: the bin switches to match the dominant subcompartment state
  of its new chromosomal neighborhood. It has taken on the identity of its
  new environment.
- **Other**: the bin changed subcompartment, but not to the dominant state
  of the new neighborhood.

The "adopted" category is of particular biological interest: a high adoption
rate suggests the translocated segment is being reshaped by its new chromatin
environment, which may in turn explain changes in gene expression.

---

## Overview

This is a single self-contained script.

| Script | What it does | Input | Output |
|--------|-------------|-------|--------|
| `10_subcompartment_classification.py` | Classifies translocated bins as retained/adopted/other at subcompartment and A/B levels; runs binomial tests against a genome-wide control | Subcompartment bedGraph, translocation BED files, neighbor BED files | 6 plots (3 per condition) + `bins_annotation.csv` + `all_binomial_stats.csv` |

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
| `GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph` | 100 kb subcompartment annotations for WT, T1, and C1 | GEO: GSE246947 |
| `verify_T1_translocations.bed` | Genomic coordinates of verified translocated regions in T1 | This study |
| `verify_C1_translocations.bed` | Genomic coordinates of verified translocated regions in C1 | This study |
| `verify_T1_neighbor.bed` | neighboring regions on the receiving chromosome for T1 | This study |
| `verify_T1_neighbor_originchr.bed` | neighboring regions on the origin chromosome for T1 | This study |
| `verify_C1_neighbor.bed` | neighboring regions on the receiving chromosome for C1 | This study |
| `verify_C1_neighbor_originchr.bed` | neighboring regions on the origin chromosome for C1 | This study |

### Subcompartment bedGraph format

Tab-separated with one row per 100 kb bin. Must include columns named
`MCF10A_WT.state`, `MCF10A_T1.state`, and `MCF10A_C1.state` containing
subcompartment labels such as A1, A2, A3, B1, B2, B3.

### Translocation BED format

Tab-separated, must include at least:

```
chrom   start   end   label   transloc_id
```

`transloc_id` values must match the keys in `T1_NAME_MAP` and `C1_NAME_MAP`
in the CONFIG section.

---

## How to run

### Step 1: Edit the CONFIG section

Open `10_subcompartment_classification.py` and edit the variables in the
`CONFIG SECTION` near the top of the script:

```python
# Subcompartment annotations
SUBCOMPARTMENT_FILE = "/path/to/subcompartments.bedGraph"

# Translocation BED files
T1_TRANSLOC_BED = "/path/to/verify_T1_translocations.bed"
C1_TRANSLOC_BED = "/path/to/verify_C1_translocations.bed"

# Neighbor BED files
# - "new" neighbor: where the translocated piece ENDED UP
# - "old" neighbor: where the translocated piece CAME FROM
T1_NEW_NEIGHBOR = "/path/to/verify_T1_neighbor.bed"
T1_OLD_NEIGHBOR = "/path/to/verify_T1_neighbor_originchr.bed"
C1_NEW_NEIGHBOR = "/path/to/verify_C1_neighbor.bed"
C1_OLD_NEIGHBOR = "/path/to/verify_C1_neighbor_originchr.bed"

# Maps from transloc_id (in BED files) to readable names for plot labels
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

# Output folder
OUTPUT_FOLDER = "/path/to/results/subcompartment_classification"
```

The keys in `T1_NAME_MAP` and `C1_NAME_MAP` must match the `transloc_id`
column in your BED files. Check by opening a BED file and reading that column.

### Step 2: Run the script

```bash
python 10_subcompartment_classification.py
```

---

## Output files

### Plots

Six plots are produced in total, three per condition (T1 and C1):

| File | Description |
|------|-------------|
| `T1_subcompartment.png` | Fraction of bins retained/adopted/other at full subcompartment resolution (A1, A2, A3, B1, B2, B3) |
| `T1_collapsed.png` | Same at collapsed A/B level |
| `T1_direction.png` | Fraction of bins switching A->B, B->A, or retained |
| `C1_subcompartment.png` | Same plots for C1 |
| `C1_collapsed.png` | |
| `C1_direction.png` | |

Each plot includes a "Genome control" bar showing the background rate of
switching across all bins not in any translocated region.

### Data files

| File | Description |
|------|-------------|
| `bins_annotation.csv` | Per-bin classifications (retained/adopted/other/direction) for all translocated bins across T1 and C1 |
| `all_binomial_stats.csv` | Binomial test results for every translocation x category x analysis level combination |

---

## Design decisions

**What does "adopted" mean?**
A bin is classified as adopted if it changes subcompartment state to match
the most common (dominant) subcompartment in its new chromosomal neighborhood
(defined by the neighbor BED file for the receiving chromosome). This is the
most biologically meaningful category: it indicates the translocated segment
has taken on the chromatin identity of its new environment.

**Why is the dominant state of the neighborhood used as the adopted state?**
The translocated segment lands next to existing chromatin with an established
compartment identity. If the segment adopts the same compartment as the
majority of its new neighbors, it is most likely being shaped by local
chromatin propagation mechanisms (e.g. spreading of histone modifications),
rather than changing at random.

**Why a binomial test instead of Wilcoxon or Mann-Whitney U?**
The data here is fundamentally different from the distance analyzes in
scripts 05–09. There, each gene or model produced a continuous distance
value, and we asked whether two distributions of distances were different.
Here, each bin produces a categorical label ("retained", "adopted", or
"other"), and we are counting *how many* bins fall into each category. The
question becomes: out of *N* total bins, *K* were classified as adopted,
is *K* larger than we'd expect given the genome-wide background rate? That
is exactly what a binomial test is designed to answer. Wilcoxon and
Mann-Whitney U operate on continuous values and cannot be applied to these
counts directly.

**Why a binomial test against the genome control?**
The genome-wide control provides the background rate of subcompartment
switching, i.e. the fraction of all non-translocated bins that change state
between WT and the condition. This accounts for the fact that some switching
happens throughout the genome simply due to noise or cell-line differences
unrelated to the translocation. The binomial test then asks: is the fraction
of adopted (or retained, or other) bins in this translocation significantly
higher or lower than the background rate?

**Why is "adopted" not defined for control bins?**
Control bins have no defined new neighborhood (they are not translocated),
so there is no meaningful "adopted" state. Control bins are therefore
classified as only "retained" or "other". This means the genome control bar
in each plot has no "adopted" segment, which is expected.

**Collapsed A/B analysis**
The collapsed analysis reduces the six subcompartment labels to two (A and B)
for a simpler view of large-scale compartment switching. This is useful because
A1/A2/A3 differences are more subtle than the A/B boundary, and the collapsed
view is easier to interpret statistically.

