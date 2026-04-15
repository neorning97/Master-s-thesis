"""
09_nuclear_radial_positioning.py
==================================
Script for analysing the radial nuclear positioning of translocated genomic
segments across wildtype (WT), T1, and C1 Chrom3D structural models.

What this script does
---------------------
Chromosomal translocations move segments of DNA from one chromosome to
another. This repositioning may also change where those segments sit in 3D
nuclear space. Specifically, whether they move toward the nuclear interior
(associated with active A-compartments) or toward the nuclear periphery
(associated with inactive B-compartments).

This script tracks the radial position (distance to the nuclear centre) of
each translocated genomic segment across three conditions:
  - WT  : wildtype models, before any translocation
  - T1  : premalignant translocated cell line
  - C1  : malignant translocated cell line

For each Chrom3D structural model and each genomic segment, the script:
1. Finds all beads overlapping the segment's genomic coordinates.
2. Computes the Euclidean distance of each bead from the origin (0, 0, 0),
   which represents the nuclear centre in Chrom3D models.
3. Averages those distances across all beads in the segment, giving one
   radial distance value per segment per model.
4. Aggregates across models and produces four types of plot plus a
   statistics table.

Hypothesis
----------
Based on compartment analysis, two specific predictions are tested:

  - Der(17) [contains the chr3-derived segment]: expected to move toward
    the nuclear INTERIOR (smaller radial distance) from WT -> T1 -> C1,
    consistent with adoption of an active A-compartment environment.

  - Der(3) [contains the chr17-derived segment]: expected to move toward
    the nuclear PERIPHERY (larger radial distance) from WT -> T1 -> C1,
    consistent with adoption of an inactive B-compartment environment.

Output plots
------------
  1. trajectory_radial_distance.png
     Mean radial distance per segment across WT -> T1 -> C1 (one line per
     translocation segment with error bars showing ±SEM).

  2. trajectory_chr3_chr17_and_fragments.png
     Same as plot 1 but showing whole chr3 and chr17 alongside the small
     translocated fragments, to contextualise the segment-level shifts
     against the chromosome-wide radial positions.

Statistical test
----------------
Wilcoxon signed-rank test, paired across models. For each segment and each
condition pair (WT vs T1, WT vs C1, T1 vs C1), the per-model mean radial
distances are compared. Pairing is done by model index (model_0000,
model_0001, ...) so the same structural ensemble is compared across
conditions. Results are saved to radial_positioning_stats.tsv.

Usage
-----
    1. Edit the CONFIG section below to point to your files.
    2. Run: python 09_nuclear_radial_positioning.py

Dependencies
------------
    pandas, numpy, matplotlib, scipy
    Install with: pip install pandas numpy matplotlib scipy

"""

import os
import re
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon


# =============================================================================
# CONFIG – edit all paths and settings here before running
# =============================================================================
CONFIG = {
    # Directories containing Chrom3D .cmm structural model files
    "cmm_dirs": {
        "WT": "/path/to/cmm/WT",
        "T1": "/path/to/cmm/T1",
        "C1": "/path/to/cmm/C1",
    },

    # BED files defining translocated segments (one row per segment).
    # Required columns: chrom, start, end, transloc_id, label
    "transloc_beds": {
        "T1": "/path/to/verify_T1_translocations.bed",
        "C1": "/path/to/verify_C1_translocations.bed",
    },

    # Maps transloc_id values in each BED file to human-readable segment names.
    # Update these to match the transloc_id values in your BED files.
    "transloc_name_maps": {
        "T1": {
            "T1": "Der(17)t(3;17)",
            "T2": "Der(3)t(3;17)",
            "T3": "t(6;19)",
        },
        "C1": {
            "T1": "t(2;10)",
            "T2": "Der(17)t(3;17)",
            "T3": "Der(3)t(3;17)",
            "T4": "t(6;19)",
        },
    },

    # Colour scheme: one colour per named translocation segment
    "segment_colors": {
        "Der(3)t(3;17)":          "#e41a1c",   # red
        "Der(17)t(3;17)":         "#2171b5",   # blue
        "t(6;19)":                "#756bb1",   # purple
        "t(2;10)":                "#31a354",   # green
        "chr3":                   "#e41a1c",   # red
        "chr17":                  "#2171b5",   # blue
        "chr3 fragment in Der(17)":  "#fd8d3c",  # orange
        "chr17 fragment in Der(3)":  "#31a354",  # green
    },

    # Where to save all output files (created automatically)
    "output_dir": "/path/to/results/nuclear_positioning",
}

# Order in which conditions are displayed in all plots
CONDITIONS = ["WT", "T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def parse_cmm(file_path: str) -> pd.DataFrame:
    """
    Parse a UCSF Chimera / Chrom3D .cmm marker file and return a DataFrame
    with columns: chr, start, end, x, y, z.

    Each <marker> line with a beadID of the form 'chrN:start-end' represents
    one genomic bead positioned at (x, y, z) in 3D space. <link> elements
    and all other attributes are ignored.
    """
    markers = []
    with open(file_path) as fh:
        for line in fh:
            if "<marker" not in line:
                continue
            attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
            if "beadID" not in attrs:
                continue
            try:
                chrom_part, coords = attrs["beadID"].split(":")
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part
                start, end = map(int, coords.split("-"))
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, start, end, x, y, z])
            except Exception:
                continue  # skip malformed lines silently

    return pd.DataFrame(markers, columns=["chr", "start", "end", "x", "y", "z"])


def radial_distance(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    Compute Euclidean distance from the nuclear centre (origin = 0, 0, 0).

    In Chrom3D models the nucleus is modelled as a sphere centred on the
    origin, so this distance is a direct measure of radial nuclear position.
    The nuclear radius was set to 5 µm for this dataset.
    """
    return np.sqrt(x**2 + y**2 + z**2)


def get_beads_for_segment(
    beads_df: pd.DataFrame,
    chrom: str,
    seg_start: int,
    seg_end: int,
) -> pd.DataFrame:
    """
    Return all beads that overlap a given genomic segment [seg_start, seg_end)
    on the specified chromosome.

    A bead overlaps if any part of it falls within the segment, using the
    same half-open interval logic as BED format.
    """
    mask = (
        (beads_df["chr"] == chrom) &
        (beads_df["start"] < seg_end) &
        (beads_df["end"]   > seg_start)
    )
    return beads_df[mask].copy()


def build_segment_table(config: dict) -> pd.DataFrame:
    """
    Build a unified table of all genomic segments to track, combining:
      1. Translocated segments from the T1 and C1 BED files
      2. Small chr3 and chr17 fragments extracted from the t(3;17) BED rows
      3. Whole-chromosome entries for chr3 and chr17 (for context)

    Returns a DataFrame with columns:
      label, chrom, seg_start, seg_end
    """
    records = []

    # ── Translocation segments from BED files ────────────────────────────────
    for cond_key, bed_path in config["transloc_beds"].items():
        name_map = config["transloc_name_maps"][cond_key]
        trans = pd.read_csv(bed_path, sep="\t")
        trans[["start", "end"]] = trans[["start", "end"]].astype(int)
        trans["chrom"] = trans["chrom"].astype(str).apply(
            lambda x: x if x.startswith("chr") else "chr" + x
        )

        for tid, label in name_map.items():
            rows = trans[trans["transloc_id"] == tid]
            for _, seg in rows.iterrows():
                records.append({
                    "label":     label,
                    "chrom":     seg["chrom"],
                    "seg_start": seg["start"],
                    "seg_end":   seg["end"],
                })

            # Extract the small chr3 and chr17 fragments from t(3;17) BED rows
            if label == "Der(17)t(3;17)":
                frag = rows[rows["chrom"] == "chr3"]
                for _, seg in frag.iterrows():
                    records.append({
                        "label":     "chr3 fragment in Der(17)",
                        "chrom":     "chr3",
                        "seg_start": seg["start"],
                        "seg_end":   seg["end"],
                    })
            elif label == "Der(3)t(3;17)":
                frag = rows[rows["chrom"] == "chr17"]
                for _, seg in frag.iterrows():
                    records.append({
                        "label":     "chr17 fragment in Der(3)",
                        "chrom":     "chr17",
                        "seg_start": seg["start"],
                        "seg_end":   seg["end"],
                    })

    # ── Whole chromosomes (for context) ──────────────────────────────────────
    for chrom in ["chr3", "chr17"]:
        records.append({
            "label":     chrom,
            "chrom":     chrom,
            "seg_start": 0,
            "seg_end":   300_000_000,  # larger than any human chromosome
        })

    segments_df = pd.DataFrame(records).drop_duplicates(
        subset=["label", "chrom", "seg_start", "seg_end"]
    )
    return segments_df.reset_index(drop=True)


# =============================================================================
# Main – collect per-bead radial distances across all models and segments
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Build segment table ───────────────────────────────────────────────────────
print("Building segment table...")
segments_df = build_segment_table(CONFIG)
print(f"  {len(segments_df)} unique genomic segments to track:")
print(segments_df[["label", "chrom", "seg_start", "seg_end"]].to_string(index=False))

# ── Sanity check: print chromosome ranges in the first model per condition ────
for cond in ["T1", "C1"]:
    cmm_files = sorted(glob.glob(os.path.join(CONFIG["cmm_dirs"][cond], "*.cmm")))
    if cmm_files:
        test_beads = parse_cmm(cmm_files[0])
        print(f"\n{cond} — chromosomes in first CMM file:")
        print(sorted(test_beads["chr"].unique()))
        for ch in ["chr3", "chr17"]:
            sub = test_beads[test_beads["chr"] == ch]
            if not sub.empty:
                print(f"  {ch}: {sub['start'].min():,} – {sub['end'].max():,}")

# ── Main loop: collect per-bead radial distances ──────────────────────────────
print("\nProcessing CMM models...")
all_bead_records = []

for cond in CONDITIONS:
    cmm_dir   = CONFIG["cmm_dirs"][cond]
    cmm_files = sorted(glob.glob(os.path.join(cmm_dir, "*.cmm")))
    print(f"\n  {cond}: {len(cmm_files)} CMM files")

    if not cmm_files:
        print(f"  WARNING: no CMM files found for {cond}, skipping.")
        continue

    for model_idx, f in enumerate(cmm_files):
        # Use a zero-padded integer index as the model ID so that models pair
        # correctly across conditions regardless of filename differences.
        # This assumes CMM files are sorted consistently across conditions,
        # which sorted(glob(...)) guarantees if naming is consistent.
        model_id = f"model_{model_idx:04d}"
        beads    = parse_cmm(f)

        if beads.empty:
            print(f"    {model_id}: no beads parsed, skipping.")
            continue

        beads["radial_dist"] = radial_distance(
            beads["x"].values, beads["y"].values, beads["z"].values
        )

        for _, seg in segments_df.iterrows():
            seg_beads = get_beads_for_segment(
                beads, seg["chrom"], seg["seg_start"], seg["seg_end"]
            )
            if seg_beads.empty:
                continue

            for _, bead in seg_beads.iterrows():
                all_bead_records.append({
                    "condition":   cond,
                    "model":       model_id,
                    "label":       seg["label"],
                    "chrom":       seg["chrom"],
                    "seg_start":   seg["seg_start"],
                    "seg_end":     seg["seg_end"],
                    "bead_start":  bead["start"],
                    "bead_end":    bead["end"],
                    "radial_dist": bead["radial_dist"],
                })

df_beads = pd.DataFrame(all_bead_records)

if df_beads.empty:
    print("\nERROR: No bead data collected. Check CMM paths and segment coordinates.")
    raise SystemExit(1)

# Save raw per-bead data for reference
raw_out = os.path.join(CONFIG["output_dir"], "per_bead_radial_distances.csv")
df_beads.to_csv(raw_out, index=False)
print(f"\nSaved per-bead data: {raw_out}  ({len(df_beads)} rows)")

# ── Aggregate: mean radial distance per (condition, label, model) ─────────────
# This gives one value per model per condition per segment, which is what
# the Wilcoxon test uses — paired across models.
df_agg = (
    df_beads
    .groupby(["condition", "label", "model"], as_index=False)["radial_dist"]
    .mean()
    .rename(columns={"radial_dist": "mean_radial"})
)

# Aggregate across models for trajectory plots (mean ± SEM)
df_traj = (
    df_agg
    .groupby(["condition", "label"])["mean_radial"]
    .agg(mean="mean", sem="sem", n="count")
    .reset_index()
)
df_traj["condition"] = pd.Categorical(
    df_traj["condition"], categories=CONDITIONS, ordered=True
)
df_traj = df_traj.sort_values(["label", "condition"])

print("\nModel counts per condition (should be equal for pairing to work):")
print(df_agg.groupby("condition")["model"].nunique())


# =============================================================================
# Statistics: Wilcoxon signed-rank test, paired across models
# =============================================================================

print("\n--- Wilcoxon signed-rank tests (paired across models) ---")
stat_rows = []

for label in df_agg["label"].unique():
    sub = df_agg[df_agg["label"] == label]

    def get_vals(c):
        return sub[sub["condition"] == c].set_index("model")["mean_radial"]

    for c1_name, c2_name in [("WT", "T1"), ("WT", "C1"), ("T1", "C1")]:
        s1     = get_vals(c1_name)
        s2     = get_vals(c2_name)
        common = s1.index.intersection(s2.index)

        if len(common) < 5:
            print(f"  {label} | {c1_name} vs {c2_name}: "
                  f"too few paired models ({len(common)}), skipping.")
            stat_rows.append({
                "segment":       label,
                "comparison":    f"{c1_name}_vs_{c2_name}",
                "n_models":      len(common),
                "mean_diff":     np.nan,
                "median_diff":   np.nan,
                "wilcoxon_stat": np.nan,
                "p_value":       np.nan,
            })
            continue

        v1   = s1[common].values
        v2   = s2[common].values
        diff = v1 - v2

        try:
            stat, p = wilcoxon(v1, v2)
        except ValueError:
            stat, p = np.nan, np.nan

        # Positive diff = c1_name has larger radial distance (more peripheral)
        direction = "↑ periphery" if np.median(diff) < 0 else "↓ interior"
        print(f"  {label} | {c1_name} vs {c2_name}: "
              f"Δmedian={np.median(diff):+.4f} ({direction}), "
              f"W={stat:.1f}, p={p:.4g}, n={len(common)}")

        stat_rows.append({
            "segment":       label,
            "comparison":    f"{c1_name}_vs_{c2_name}",
            "n_models":      len(common),
            "mean_diff":     np.mean(diff),
            "median_diff":   np.median(diff),
            "wilcoxon_stat": stat,
            "p_value":       p,
        })

stats_df  = pd.DataFrame(stat_rows)
stats_out = os.path.join(CONFIG["output_dir"], "radial_positioning_stats.tsv")
stats_df.to_csv(stats_out, sep="\t", index=False)
print(f"\nSaved statistics table: {stats_out}")


# =============================================================================
# Plotting helper
# =============================================================================

def plot_trajectory(ax, df_traj_sub, labels, colors):
    """
    Draw a trajectory plot: one line per segment showing mean radial distance
    (± SEM) across WT → T1 → C1.

    Each point is the mean across all structural models for that condition.
    Error bars show ±1 standard error of the mean across models.
    """
    for label in labels:
        sub   = df_traj_sub[df_traj_sub["label"] == label].sort_values("condition")
        color = colors.get(label, "#636363")
        x     = [CONDITIONS.index(c) for c in sub["condition"]]
        ax.errorbar(
            x, sub["mean"].values, yerr=sub["sem"].values,
            marker="o", linewidth=2, markersize=7,
            capsize=4, color=color, label=label,
        )
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["WT", "T1\n(premalignant)", "C1\n(malignant)"], fontsize=11)
    ax.set_ylabel("Mean radial distance from nuclear centre (µm)", fontsize=10)
    ax.legend(title="Segment", bbox_to_anchor=(1.02, 1), loc="upper left")


# =============================================================================
# Plot 1 – Trajectory: translocation segments across WT → T1 → C1
# =============================================================================

print("\nGenerating trajectory plot...")

transloc_labels  = ["Der(3)t(3;17)", "Der(17)t(3;17)", "t(6;19)", "t(2;10)"]
df_traj_transloc = df_traj[df_traj["label"].isin(transloc_labels)]

fig, ax = plt.subplots(figsize=(7, 5))
plot_trajectory(ax, df_traj_transloc, transloc_labels, CONFIG["segment_colors"])
ax.set_title(
    "Nuclear radial positioning of translocated segments\nWT → T1 → C1",
    fontweight="bold",
)
plt.tight_layout()
plt.savefig(
    os.path.join(CONFIG["output_dir"], "trajectory_radial_distance.png"),
    dpi=300, bbox_inches="tight",
)
plt.close()
print("  Saved: trajectory_radial_distance.png")


# =============================================================================
# Plot 2 – Chr3/Chr17 whole chromosomes + small fragments combined
#           Provides context: how does the segment-level shift compare with
#           the overall radial position of the full chromosome?
# =============================================================================

print("Generating combined chr3/chr17 trajectory plot...")

all_chrom_labels = [
    "chr3", "chr17",
    "chr3 fragment in Der(17)",
    "chr17 fragment in Der(3)",
]
df_traj_chrom = df_traj[df_traj["label"].isin(all_chrom_labels)]

fig, ax = plt.subplots(figsize=(8, 5))
plot_trajectory(ax, df_traj_chrom, all_chrom_labels, CONFIG["segment_colors"])
ax.set_title(
    "Nuclear radial positioning of chr3 and chr17\nWT → T1 → C1",
    fontweight="bold",
)
plt.tight_layout()
plt.savefig(
    os.path.join(CONFIG["output_dir"], "trajectory_chr3_chr17_and_fragments.png"),
    dpi=300, bbox_inches="tight",
)
plt.close()
print("  Saved: trajectory_chr3_chr17_and_fragments.png")


print(f"\nAll plots saved to: {CONFIG['output_dir']}")
print(f"Statistics saved to: {stats_out}")
