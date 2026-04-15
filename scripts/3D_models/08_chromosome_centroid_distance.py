
"""
08_chromosome_centroid_distance.py
====================================
Computes and compares the 3D distance between the centroids of two chromosomes
across wildtype (WT), T1, and C1 Chrom3D structural models.

What this script does
---------------------
Chromosomal translocations physically connect two chromosomes. One consequence
is that the two chromosomes involved may be pulled into closer spatial proximity
in the nucleus — they need to be near each other for the translocation to have
occurred, and may remain spatially associated afterwards.

This script tests whether the two translocated chromosomes (e.g. chr3 and
chr17) are significantly closer together in the translocated cell lines (T1,
C1) compared to the wildtype (WT), by:

1. Parsing all Chrom3D .cmm structural model files for each condition.
2. Computing the centroid (mean x, y, z position) of each chromosome across
   all its beads in each model.
3. Calculating the Euclidean distance between the two chromosome centroids
   in each model.
4. Plotting the distribution of distances across models as a box plot with
   individual model points overlaid.
5. Running pairwise Mann-Whitney U tests to assess whether distances differ
   significantly between conditions.

Why chromosome centroids?
--------------------------
Each chromosome is represented as a chain of genomic beads in Chrom3D. The
centroid — the mean 3D position of all beads on a chromosome — gives a single
summary coordinate for the chromosome's location in the nucleus. The distance
between two centroids is therefore a measure of how spatially separated the
two chromosomes are overall, averaged across their full lengths.

Why Mann-Whitney U instead of a t-test?
-----------------------------------------
The distribution of centroid distances across structural models is not assumed
to be normally distributed. The Mann-Whitney U test is a non-parametric
alternative that tests whether one condition tends to produce larger or smaller
distances than another, without assuming a particular distribution shape.

Usage
-----
    1. Edit the CONFIG section below:
       - Set CHR_A and CHR_B to the two chromosomes involved in the
         translocation (e.g. "chr3" and "chr17").
       - Point the cmm_dirs to your model directories for WT, T1, and C1.
       - Set output_dir to where you want plots saved.
    2. Run: python 08_chromosome_centroid_distance.py

Dependencies
------------
    pandas, numpy, matplotlib, seaborn, scipy
    Install with: pip install pandas numpy matplotlib seaborn scipy

"""

import os
import re
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


# =============================================================================
# CONFIG – edit all settings here before running
# =============================================================================
CONFIG = {
    # The two chromosomes to measure distance between.
    # Set these to the chromosomes involved in the translocation.
    "chr_A": "chr3",
    "chr_B": "chr17",

    # Directories containing .cmm structural model files for each condition
    "cmm_dirs": {
        "WT": "/path/to/cmm/WT",
        "T1": "/path/to/cmm/T1",
        "C1": "/path/to/cmm/C1",
    },

    # Where to save the output plot (created automatically)
    "output_dir": "/path/to/results/plots/chromosome_centroid_distance619",
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def parse_cmm(file_path: str) -> pd.DataFrame:
    """
    Parse a UCSF Chimera / Chrom3D .cmm marker file and return a DataFrame
    with columns: chr, x, y, z.

    Each <marker> line with a beadID of the form 'chrN:start-end' represents
    one genomic bead positioned at (x, y, z) in 3D space. Only the chromosome
    and coordinates are needed here; start/end positions are not used because
    we are computing chromosome-level centroids, not gene-level positions.
    <link> elements and all other attributes are ignored.
    """
    markers = []
    with open(file_path, "r") as fh:
        for line in fh:
            if "<marker" not in line:
                continue
            attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
            if "beadID" not in attrs:
                continue
            try:
                chrom_part, _ = attrs["beadID"].split(":")
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, x, y, z])
            except Exception:
                continue  # skip malformed lines silently

    return pd.DataFrame(markers, columns=["chr", "x", "y", "z"])


def compute_centroid(beads_df: pd.DataFrame, chrom: str) -> np.ndarray | None:
    """
    Compute the centroid (mean 3D position) of all beads on a given chromosome.

    Returns None if the chromosome has no beads in this model, which would
    cause this model to be skipped for that chromosome pair.
    """
    sub = beads_df[beads_df["chr"] == chrom]
    if sub.empty:
        return None
    return sub[["x", "y", "z"]].mean().values


# =============================================================================
# Main – compute centroid distances across all models and conditions
# =============================================================================

chr_A = CONFIG["chr_A"]
chr_B = CONFIG["chr_B"]

results = []

for cond, cmm_dir in CONFIG["cmm_dirs"].items():
    cmm_files = sorted(glob.glob(os.path.join(cmm_dir, "*.cmm")))
    print(f"\nProcessing {cond}: {len(cmm_files)} models found")

    for f in cmm_files:
        beads = parse_cmm(f)
        if beads.empty:
            continue

        centroid_A = compute_centroid(beads, chr_A)
        centroid_B = compute_centroid(beads, chr_B)

        # Skip this model if either chromosome has no beads
        if centroid_A is None or centroid_B is None:
            print(f"  Skipping {os.path.basename(f)}: "
                  f"missing beads for {chr_A} or {chr_B}")
            continue

        # Euclidean distance between the two chromosome centroids
        dist = np.linalg.norm(centroid_A - centroid_B)

        results.append({
            "Condition": cond,
            "Distance":  dist,
            "Model":     os.path.basename(f),
        })

df = pd.DataFrame(results)
print(f"\nTotal: {len(df)} model-condition pairs computed")
print(df.groupby("Condition")["Distance"].describe().round(3))

# =============================================================================
# Box plot with individual model points overlaid
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# Use a consistent condition order: WT first, then T1, C1
condition_order = [c for c in ["WT", "T1", "C1"] if c in df["Condition"].unique()]

fig, ax = plt.subplots(figsize=(6, 6))

sns.boxplot(
    data=df,
    x="Condition", y="Distance",
    order=condition_order,
    ax=ax,
)
# Overlay individual model points so the reader can see the raw distribution
sns.stripplot(
    data=df,
    x="Condition", y="Distance",
    order=condition_order,
    jitter=True,
    color="black", alpha=0.5, size=4,
    ax=ax,
)

ax.set_title(f"Distance between {chr_A} and {chr_B} centroids")
ax.set_ylabel("Centroid distance (model units)")
ax.set_xlabel("Condition")
fig.tight_layout()

plot_path = os.path.join(CONFIG["output_dir"], "chromosome_distance_boxplot.png")
fig.savefig(plot_path)
plt.close(fig)
print(f"\nPlot saved: {plot_path}")

# =============================================================================
# Pairwise Mann-Whitney U tests
# =============================================================================

def compare(df: pd.DataFrame, cond1: str, cond2: str) -> None:
    """
    Run a two-sided Mann-Whitney U test comparing centroid distances between
    two conditions and print the result.

    The Mann-Whitney U test is used rather than a t-test because the distance
    distributions across models are not assumed to be normally distributed.
    It tests whether one condition tends to produce systematically larger or
    smaller distances than the other.
    """
    d1 = df[df["Condition"] == cond1]["Distance"]
    d2 = df[df["Condition"] == cond2]["Distance"]

    if len(d1) == 0 or len(d2) == 0:
        print(f"  {cond1} vs {cond2}: insufficient data")
        return

    stat, p = mannwhitneyu(d1, d2, alternative="two-sided")
    median1  = d1.median()
    median2  = d2.median()
    direction = "closer" if median2 < median1 else "further"
    print(f"  {cond1} vs {cond2}: U={stat:.0f}, p={p:.4g}  "
          f"(median {cond1}={median1:.3f}, median {cond2}={median2:.3f} "
          f"— {cond2} {direction} than {cond1})")

print("\nStatistical tests (Mann-Whitney U, two-sided):")
compare(df, "WT", "T1")
compare(df, "WT", "C1")
compare(df, "T1", "C1")
