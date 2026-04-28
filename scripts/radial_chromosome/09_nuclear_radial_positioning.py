"""
09_nuclear_radial_positioning.py
==================================
This script tracks where translocated DNA segments sit in the nucleus,
across three conditions: WT (normal), T1 (premalignant), C1 (malignant).

The big question:
- After a translocation, do the moved DNA segments end up CLOSER to the
  centre of the nucleus (the "active" zone where genes are turned on)
  or FURTHER from the centre (the "inactive" zone near the edge)?
- We want to see if the position changes from WT --> T1 --> C1.

The hypothesis we're testing:
- Der(17): contains a chunk of chr3 that ended up on chr17.
           We expect this to move TOWARD the centre (more active).
- Der(3):  contains a chunk of chr17 that ended up on chr3.
           We expect this to move AWAY from the centre (more inactive).

How we do it:
1. For each 3D model file (.cmm), read all the bead positions.
2. For each translocated segment (defined in BED files), find all the
   beads that overlap with that segment.
3. Calculate each bead's distance from the nuclear centre (0, 0, 0).
4. Average those distances per segment per model.
5. Make trajectory plots showing how the average changes across
   WT --> T1 --> C1.
6. Run Wilcoxon signed-rank tests to check if the changes are
   statistically significant.

About the BED files:
- A BED file is a simple tab-separated text file describing genomic
  regions. Each row has: chromosome, start, end, plus extra columns.
- We use them here to define which segments were translocated.

Edit the paths and settings below before running.

Required libraries: pandas, numpy, matplotlib, scipy
Install them with: pip install pandas numpy matplotlib scipy
"""

import os                         
import re                          
import glob                        
import numpy as np                 
import pandas as pd                
import matplotlib.pyplot as plt    
from scipy.stats import wilcoxon   

# =============================================================================
# CONFIG SECTION - Edit these paths and settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# Folders containing the 3D model files (.cmm) for each condition
# -----------------------------------------------------------------------------
WT_CMM_FOLDER = "/path/to/cmm/WT"
T1_CMM_FOLDER = "/path/to/cmm/T1"
C1_CMM_FOLDER = "/path/to/cmm/C1"

# Put them in a dictionary so we can loop through them later
cmm_folders = {
    "WT": WT_CMM_FOLDER,
    "T1": T1_CMM_FOLDER,
    "C1": C1_CMM_FOLDER,
}

# -----------------------------------------------------------------------------
# BED files describing the translocated segments
# -----------------------------------------------------------------------------
# A BED file row tells us "this segment from chromX is at position start-end"
# Required columns in these files: chrom, start, end, transloc_id, label
T1_BED_FILE = "/path/to/verify_T1_translocations.bed"
C1_BED_FILE = "/path/to/verify_C1_translocations.bed"

bed_files = {
    "T1": T1_BED_FILE,
    "C1": C1_BED_FILE,
}

# -----------------------------------------------------------------------------
# Maps from transloc_id (used in the BED files) to nice human-readable names
# -----------------------------------------------------------------------------
# In the BED files, each translocation has a short ID like "T1", "T2", etc.
# These maps turn those into readable names like "Der(17)t(3;17)".
# Update these to match your actual BED files.
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

# Group these together for looping later
name_maps = {
    "T1": T1_NAME_MAP,
    "C1": C1_NAME_MAP,
}

# -----------------------------------------------------------------------------
# Colours for the plots - one colour per segment label
# -----------------------------------------------------------------------------
# Hex codes like "#e41a1c" are colours: a way to specify a colour using
# 6 hex digits (red, green, blue). You can find tools online to pick them.
segment_colors = {
    "Der(3)t(3;17)":             "#e41a1c",   # red
    "Der(17)t(3;17)":            "#2171b5",   # blue
    "t(6;19)":                   "#756bb1",   # purple
    "t(2;10)":                   "#31a354",   # green
    "chr3":                      "#e41a1c",   # red
    "chr17":                     "#2171b5",   # blue
    "chr3 fragment in Der(17)":  "#fd8d3c",   # orange
    "chr17 fragment in Der(3)":  "#31a354",   # green
}

# -----------------------------------------------------------------------------
# Where to save all the output files (plots and stats table)
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/results/plots/nuclear_positioning"

# Order of conditions in the plots (left to right)
CONDITIONS = ["WT", "T1", "C1"]


# =============================================================================
# Helper function: Read a .cmm file and pull out the bead positions
# =============================================================================
# Same as in the previous scripts: parses the .cmm file and returns
# a DataFrame with chr, start, end, x, y, z columns.

def read_cmm_file(file_path):
    """
    Open a .cmm file and read out all the bead information.
    Returns a pandas DataFrame with columns: chr, start, end, x, y, z
    """
    # Empty list to fill with bead info
    bead_list = []

    # Open the file for reading, line by line
    with open(file_path) as my_file:
        for line in my_file:

            # Skip lines that don't have a marker
            if "<marker" not in line:
                continue

            # Find all attribute="value" pairs and turn them into a dict
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            # Skip if there's no beadID
            if "beadID" not in attrs:
                continue

            # Try to extract the values, skip the line if anything goes wrong
            try:
                # The beadID looks like "chr1:1000-2000"
                bead_id = attrs["beadID"]
                chrom_part, coords_part = bead_id.split(":")

                # Make sure chromosome name starts with "chr"
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part

                # Get start and end as integers
                start_str, end_str = coords_part.split("-")
                start = int(start_str)
                end = int(end_str)

                # Get x, y, z as floats (NaN if missing)
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))

                # Add to our list
                bead_list.append([chrom_part, start, end, x, y, z])

            except Exception:
                # Skip badly-formatted lines
                continue

    # Turn the list into a DataFrame
    beads_df = pd.DataFrame(
        bead_list,
        columns=["chr", "start", "end", "x", "y", "z"]
    )
    return beads_df


# =============================================================================
# Helper function: Compute distance from each bead to the nuclear centre
# =============================================================================
# In Chrom3D models, the nucleus is a sphere centred at (0, 0, 0).
# So the distance from a bead to the centre is just sqrt(x^2 + y^2 + z^2).
# This function takes three arrays (x values, y values, z values) and
# returns an array of distances.

def radial_distance(x, y, z):
    """
    Calculate the distance from each (x, y, z) point to the origin (0,0,0).
    """
    # ** is the exponent operator in Python: x**2 means "x squared"
    # numpy makes this work on whole arrays at once. Much faster than
    # using a Python for loop.
    distances = np.sqrt(x**2 + y**2 + z**2)
    return distances


# =============================================================================
# Helper function: Find all beads that overlap a genomic segment
# =============================================================================
# A "segment" is just a stretch of DNA: chromosome + start position + end position.
# A bead also has chromosome + start + end.
# A bead "overlaps" a segment if they share any DNA.
#
# Two intervals [a, b) and [c, d) overlap if a < d AND b > c.
# (This is the standard half-open interval logic used in BED files.)

def get_beads_for_segment(beads_df, chrom, seg_start, seg_end):
    """
    Return all beads on the given chromosome that overlap the segment
    from seg_start to seg_end.
    """
    # Build a boolean mask: True for beads that match all three conditions
    same_chrom = (beads_df["chr"] == chrom)
    starts_before_seg_ends = (beads_df["start"] < seg_end)
    ends_after_seg_starts = (beads_df["end"] > seg_start)

    # Combine with &: all three must be true
    mask = same_chrom & starts_before_seg_ends & ends_after_seg_starts

    # Filter and return a copy
    return beads_df[mask].copy()


# =============================================================================
# Helper function: Build the table of all segments we want to track
# =============================================================================
# We want to track three kinds of segments:
#   1. The translocated segments listed in the BED files (e.g. Der(17)t(3;17))
#   2. The small chr3/chr17 fragments inside the t(3;17) translocations
#      (so we can see where THAT specific bit of DNA went)
#   3. The whole chr3 and chr17 (for context: where does the whole
#      chromosome sit overall?)

def build_segment_table():
    """
    Build a DataFrame listing every segment we want to track.
    Columns: label, chrom, seg_start, seg_end
    """
    # We'll collect each segment as a small dict in this list,
    # then turn the list into a DataFrame at the end.
    records = []

    # -------------------------------------------------------------------------
    # Part 1 + 2: Translocations and their chr3/chr17 fragments from BED files
    # -------------------------------------------------------------------------
    # Loop through each condition's BED file
    for cond_key, bed_path in bed_files.items():

        # Get the name map for this condition
        name_map = name_maps[cond_key]

        # Read the BED file (tab-separated)
        trans = pd.read_csv(bed_path, sep="\t")

        # Make sure start and end are integers
        trans["start"] = trans["start"].astype(int)
        trans["end"] = trans["end"].astype(int)

        # Make sure chromosome names start with "chr"
        # Same trick as in earlier scripts
        trans["chrom"] = trans["chrom"].astype(str)
        fixed_chroms = []
        for c in trans["chrom"]:
            if c.startswith("chr"):
                fixed_chroms.append(c)
            else:
                fixed_chroms.append("chr" + c)
        trans["chrom"] = fixed_chroms

        # Loop through each transloc_id --> readable label pair in the map
        # .items() gives us both the key (tid) and value (label)
        for tid, label in name_map.items():

            # Get all the rows for this translocation ID
            rows = trans[trans["transloc_id"] == tid]

            # Add each row as a segment to our records list
            for _, seg in rows.iterrows():
                one_record = {
                    "label": label,
                    "chrom": seg["chrom"],
                    "seg_start": seg["start"],
                    "seg_end": seg["end"],
                }
                records.append(one_record)

            # ------------------------------------------------------------------
            # Special case: if this is one of the t(3;17) translocations,
            # also add the small chr3 or chr17 fragment as its own segment
            # so we can track it separately.
            # ------------------------------------------------------------------
            if label == "Der(17)t(3;17)":
                # Der(17) contains a chr3 fragment, find the chr3 rows
                frag_rows = rows[rows["chrom"] == "chr3"]
                for _, seg in frag_rows.iterrows():
                    fragment_record = {
                        "label": "chr3 fragment in Der(17)",
                        "chrom": "chr3",
                        "seg_start": seg["start"],
                        "seg_end": seg["end"],
                    }
                    records.append(fragment_record)

            elif label == "Der(3)t(3;17)":
                # Der(3) contains a chr17 fragment, find the chr17 rows
                frag_rows = rows[rows["chrom"] == "chr17"]
                for _, seg in frag_rows.iterrows():
                    fragment_record = {
                        "label": "chr17 fragment in Der(3)",
                        "chrom": "chr17",
                        "seg_start": seg["start"],
                        "seg_end": seg["end"],
                    }
                    records.append(fragment_record)

    # -------------------------------------------------------------------------
    # Part 3: Whole chromosomes (for context)
    # -------------------------------------------------------------------------
    # We add chr3 and chr17 as "segments" that span the whole chromosome.
    # We use 300,000,000 as the end because no human chromosome is bigger
    # than that, so this guarantees we'll catch every bead.
    # The underscore in 300_000_000 is just for readability (Python ignores it).
    for chrom in ["chr3", "chr17"]:
        whole_chrom_record = {
            "label": chrom,
            "chrom": chrom,
            "seg_start": 0,
            "seg_end": 300_000_000,
        }
        records.append(whole_chrom_record)

    # Turn the list of dictionaries into a DataFrame
    segments_df = pd.DataFrame(records)

    # Remove any duplicate rows
    # (a segment defined by the same label, chrom, start, end is the same segment)
    segments_df = segments_df.drop_duplicates(
        subset=["label", "chrom", "seg_start", "seg_end"]
    )

    # Reset the row numbers so they're 0, 1, 2, ... again
    segments_df = segments_df.reset_index(drop=True)

    return segments_df


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# -----------------------------------------------------------------------------
# Step 1: Build the segment table
# -----------------------------------------------------------------------------
print("Building segment table...")
segments_df = build_segment_table()

print("  " + str(len(segments_df)) + " unique genomic segments to track:")
# Print the table without the row numbers (index=False)
print(segments_df[["label", "chrom", "seg_start", "seg_end"]].to_string(index=False))

# -----------------------------------------------------------------------------
# Step 2: Sanity check - print the chromosome ranges in the first model
# -----------------------------------------------------------------------------
# This is just to help us verify that the BED file coordinates actually
# match what's in the .cmm files. If they don't, we won't find any beads.
print("")
print("Sanity check: chromosomes in the first .cmm file of T1 and C1...")

for cond in ["T1", "C1"]:
    folder = cmm_folders[cond]
    search_pattern = os.path.join(folder, "*.cmm")
    cmm_files = sorted(glob.glob(search_pattern))

    # Skip if there are no files
    if len(cmm_files) == 0:
        continue

    # Read just the first file
    first_file = cmm_files[0]
    test_beads = read_cmm_file(first_file)

    # Show what chromosomes are in it
    print("")
    print(cond + " — chromosomes in first CMM file:")
    print(sorted(test_beads["chr"].unique()))

    # For chr3 and chr17, show the start-end range
    for ch in ["chr3", "chr17"]:
        sub = test_beads[test_beads["chr"] == ch]
        if not sub.empty:
            min_start = sub["start"].min()
            max_end = sub["end"].max()
            # Format with commas for readability (e.g. 1,000,000)
            print("  " + ch + ": {:,} – {:,}".format(min_start, max_end))


# -----------------------------------------------------------------------------
# Step 3: Main loop - for each model, get the radial distance of each bead
# -----------------------------------------------------------------------------
# We'll go through each condition, each model, and each segment, and
# collect one row per bead-in-segment with its distance to the nuclear centre.

print("")
print("Processing CMM models...")

# This list will hold one record per bead-in-segment
all_bead_records = []

# Loop through each condition
for cond in CONDITIONS:

    # Get this condition's folder of model files
    cmm_folder = cmm_folders[cond]
    search_pattern = os.path.join(cmm_folder, "*.cmm")
    cmm_files = sorted(glob.glob(search_pattern))

    print("")
    print("  " + cond + ": " + str(len(cmm_files)) + " CMM files")

    # Skip if no files were found
    if len(cmm_files) == 0:
        print("  WARNING: no CMM files found for " + cond + ", skipping.")
        continue

    # Loop through each model file, with an index counter
    # enumerate gives us both the index (0, 1, 2...) and the file path
    for model_idx, f in enumerate(cmm_files):

        # Build a model_id like "model_0000", "model_0001", etc.
        # The {:04d} format means: integer with at least 4 digits,
        # padded with zeros (so 5 becomes "0005").
        # We use the same numbering across conditions so we can later
        # PAIR up the same model number across WT, T1, C1 (which we need
        # for the Wilcoxon signed-rank test).
        model_id = "model_{:04d}".format(model_idx)

        # Read the .cmm file
        beads = read_cmm_file(f)

        # Skip if no beads were parsed
        if beads.empty:
            print("    " + model_id + ": no beads parsed, skipping.")
            continue

        # ---------------------------------------------------------------------
        # Calculate the radial distance for every bead in this model
        # ---------------------------------------------------------------------
        # We add a new column "radial_dist" to the beads table.
        # We use .values to get the underlying numpy arrays, this is faster
        # than using the pandas Series directly.
        x_values = beads["x"].values
        y_values = beads["y"].values
        z_values = beads["z"].values
        beads["radial_dist"] = radial_distance(x_values, y_values, z_values)

        # ---------------------------------------------------------------------
        # For each segment, find the beads that overlap and record them
        # ---------------------------------------------------------------------
        for _, seg in segments_df.iterrows():

            # Get all beads that overlap this segment
            seg_beads = get_beads_for_segment(
                beads,
                seg["chrom"],
                seg["seg_start"],
                seg["seg_end"]
            )

            # Skip if no beads overlap this segment
            if seg_beads.empty:
                continue

            # Save one record per overlapping bead
            for _, bead in seg_beads.iterrows():
                one_record = {
                    "condition": cond,
                    "model": model_id,
                    "label": seg["label"],
                    "chrom": seg["chrom"],
                    "seg_start": seg["seg_start"],
                    "seg_end": seg["seg_end"],
                    "bead_start": bead["start"],
                    "bead_end": bead["end"],
                    "radial_dist": bead["radial_dist"],
                }
                all_bead_records.append(one_record)

# Turn our list of records into a DataFrame
df_beads = pd.DataFrame(all_bead_records)

# Stop if we got nothing, something is wrong with the inputs
if df_beads.empty:
    print("")
    print("ERROR: No bead data collected. "
          "Check CMM paths and segment coordinates.")
    # Exit the script with an error code (1 means "something went wrong")
    raise SystemExit(1)

# -----------------------------------------------------------------------------
# Step 4: Save the raw per-bead data
# -----------------------------------------------------------------------------
raw_out = os.path.join(OUTPUT_FOLDER, "per_bead_radial_distances.csv")
df_beads.to_csv(raw_out, index=False)
print("")
print("Saved per-bead data: " + raw_out + "  (" + str(len(df_beads)) + " rows)")


# -----------------------------------------------------------------------------
# Step 5: Average per (condition, segment, model)
# -----------------------------------------------------------------------------
# A segment can have several beads. We want one number per model per segment
# per condition: the average radial distance across all beads in the segment.
# We get this with groupby + mean.

# Group by condition, label (segment), and model. Take the mean of radial_dist.
df_agg = df_beads.groupby(
    ["condition", "label", "model"],
    as_index=False
)["radial_dist"].mean()

# Rename the column to make its meaning clearer
df_agg = df_agg.rename(columns={"radial_dist": "mean_radial"})


# -----------------------------------------------------------------------------
# Step 6: Aggregate across models for the trajectory plots
# -----------------------------------------------------------------------------
# For the plot, we want one point per (condition, segment), with:
#   - the mean across all models
#   - the SEM (standard error of the mean) for the error bars
#   - the count of how many models contributed

# Group by condition and segment, then calculate three stats on mean_radial
grouped = df_agg.groupby(["condition", "label"])["mean_radial"]
df_traj = grouped.agg(
    mean="mean",   # mean across models
    sem="sem",     # standard error of the mean
    n="count"      # how many models
).reset_index()

# Make "condition" a Categorical so it sorts in our preferred order
# (WT, T1, C1) instead of alphabetically.
df_traj["condition"] = pd.Categorical(
    df_traj["condition"],
    categories=CONDITIONS,
    ordered=True
)

# Sort by label, then condition (so each segment's row appears in WT-T1-C1 order)
df_traj = df_traj.sort_values(["label", "condition"])

# Print model counts to make sure we have equal numbers for pairing
print("")
print("Model counts per condition (should be equal for pairing to work):")
print(df_agg.groupby("condition")["model"].nunique())


# =============================================================================
# Step 7: Run Wilcoxon signed-rank tests
# =============================================================================
# For each segment and each pair of conditions (WT vs T1, WT vs C1, T1 vs C1),
# we test whether the radial distances are significantly different.
# This is a PAIRED test, so we have to match up models (model_0000 in WT
# is paired with model_0000 in T1, etc.).

print("")
print("--- Wilcoxon signed-rank tests (paired across models) ---")

# We'll collect one row per test in this list
stat_rows = []

# Get a list of unique segment labels
unique_labels = df_agg["label"].unique()

# The condition pairs we want to test
condition_pairs = [("WT", "T1"), ("WT", "C1"), ("T1", "C1")]

# Loop through each segment
for label in unique_labels:

    # Get only the rows for this segment
    sub = df_agg[df_agg["label"] == label]

    # Loop through each pair of conditions
    for c1_name, c2_name in condition_pairs:

        # Get the values for condition 1, indexed by model name
        # (so we can match up the same model across conditions)
        s1_rows = sub[sub["condition"] == c1_name]
        s1 = s1_rows.set_index("model")["mean_radial"]

        # Same for condition 2
        s2_rows = sub[sub["condition"] == c2_name]
        s2 = s2_rows.set_index("model")["mean_radial"]

        # Find the models that are in BOTH conditions (the paired ones)
        # .intersection gives us the model names that appear in both
        common = s1.index.intersection(s2.index)

        # If we have fewer than 5 paired models, we don't have enough data
        # to do a meaningful test, record NaN values and move on.
        if len(common) < 5:
            print("  " + label + " | " + c1_name + " vs " + c2_name + ": "
                  "too few paired models (" + str(len(common)) + "), skipping.")

            # Save a row with NaN values so the output table is complete
            no_test_row = {
                "segment": label,
                "comparison": c1_name + "_vs_" + c2_name,
                "n_models": len(common),
                "mean_diff": np.nan,
                "median_diff": np.nan,
                "wilcoxon_stat": np.nan,
                "p_value": np.nan,
            }
            stat_rows.append(no_test_row)
            continue

        # Get the actual values, in the same order, as numpy arrays
        v1 = s1[common].values
        v2 = s2[common].values

        # Calculate the differences (will use this for direction reporting)
        diff = v1 - v2

        # Run the Wilcoxon signed-rank test.
        # Wrap in try/except because it can fail if all differences are zero.
        try:
            stat, p = wilcoxon(v1, v2)
        except ValueError:
            stat, p = np.nan, np.nan

        # Figure out the direction of change.
        # If the median diff is negative, c1 was SMALLER than c2,
        # i.e. moving from c1 to c2 means moving FURTHER from centre (peripheral).
        # If positive, c2 is smaller, i.e. moving toward the interior.
        median_diff = np.median(diff)
        if median_diff < 0:
            direction = "↑ periphery"
        else:
            direction = "↓ interior"

        # Print the result.
        # {:+.4f} = sign-included, 4 decimal places
        # {:.1f}  = 1 decimal place
        # {:.4g}  = up to 4 significant digits
        print("  " + label + " | " + c1_name + " vs " + c2_name + ": "
              + "Δmedian={:+.4f}".format(median_diff)
              + " (" + direction + "), "
              + "W={:.1f}, ".format(stat)
              + "p={:.4g}, ".format(p)
              + "n=" + str(len(common)))

        # Save the result
        result_row = {
            "segment": label,
            "comparison": c1_name + "_vs_" + c2_name,
            "n_models": len(common),
            "mean_diff": np.mean(diff),
            "median_diff": median_diff,
            "wilcoxon_stat": stat,
            "p_value": p,
        }
        stat_rows.append(result_row)

# Save the stats table
stats_df = pd.DataFrame(stat_rows)
stats_out = os.path.join(OUTPUT_FOLDER, "radial_positioning_stats.tsv")
stats_df.to_csv(stats_out, sep="\t", index=False)
print("")
print("Saved statistics table: " + stats_out)


# =============================================================================
# Plotting helper function
# =============================================================================
# We use the same plot style for both plots (only the segment list differs),
# so we put the plotting code in a helper function.

def plot_trajectory(ax, df_traj_sub, labels, colors):
    """
    Draw a trajectory plot on the given axes (ax):
    one line per segment showing mean radial distance (with SEM error bars)
    across WT --> T1 --> C1.
    """
    # One line per segment label
    for label in labels:

        # Get the rows for this segment, in WT-->T1-->C1 order
        sub = df_traj_sub[df_traj_sub["label"] == label]
        sub = sub.sort_values("condition")

        # Get the colour for this segment.
        # .get(...) returns the colour if the label is in the dict,
        # otherwise it returns the default "#636363" (grey).
        color = colors.get(label, "#636363")

        # The x-axis position for each condition: 0, 1, 2 for WT, T1, C1
        # We use a list comprehension to get the index of each condition
        x = []
        for c in sub["condition"]:
            x.append(CONDITIONS.index(c))

        # Draw the line with error bars
        # marker="o": put a dot at each point
        # capsize=4: add little caps to the error bars
        # yerr=...: vertical error bars
        ax.errorbar(
            x,
            sub["mean"].values,
            yerr=sub["sem"].values,
            marker="o",
            linewidth=2,
            markersize=7,
            capsize=4,
            color=color,
            label=label,
        )

    # Set the x-axis tick positions and labels
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(
        ["WT", "T1\n(premalignant)", "C1\n(malignant)"],
        fontsize=11
    )

    # Y-axis label
    ax.set_ylabel(
        "Mean radial distance from nuclear centre (µm)",
        fontsize=10
    )

    # Add a legend to the right of the plot
    ax.legend(
        title="Segment",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )


# =============================================================================
# Plot 1: Trajectory of translocated segments across WT --> T1 --> C1
# =============================================================================

print("")
print("Generating trajectory plot...")

# Which segments to include in this plot
transloc_labels = ["Der(3)t(3;17)", "Der(17)t(3;17)", "t(6;19)", "t(2;10)"]

# Filter the trajectory data to just these segments
df_traj_transloc = df_traj[df_traj["label"].isin(transloc_labels)]

# Create the figure and plot
fig, ax = plt.subplots(figsize=(7, 5))
plot_trajectory(ax, df_traj_transloc, transloc_labels, segment_colors)

ax.set_title(
    "Nuclear radial positioning of translocated segments\nWT → T1 → C1",
    fontweight="bold"
)

# tight_layout adjusts spacing, bbox_inches="tight" prevents legend cutoff
plt.tight_layout()

plot_path_1 = os.path.join(OUTPUT_FOLDER, "trajectory_radial_distance.png")
plt.savefig(plot_path_1, dpi=300, bbox_inches="tight")
plt.close()

print("  Saved: trajectory_radial_distance.png")


# =============================================================================
# Plot 2: Whole chr3/chr17 with their fragments overlaid
# =============================================================================
# This puts the fragments next to the whole chromosomes for context.

print("Generating combined chr3/chr17 trajectory plot...")

all_chrom_labels = [
    "chr3",
    "chr17",
    "chr3 fragment in Der(17)",
    "chr17 fragment in Der(3)",
]

df_traj_chrom = df_traj[df_traj["label"].isin(all_chrom_labels)]

fig, ax = plt.subplots(figsize=(8, 5))
plot_trajectory(ax, df_traj_chrom, all_chrom_labels, segment_colors)

ax.set_title(
    "Nuclear radial positioning of chr3 and chr17\nWT → T1 → C1",
    fontweight="bold"
)

plt.tight_layout()

plot_path_2 = os.path.join(OUTPUT_FOLDER, "trajectory_chr3_chr17_and_fragments.png")
plt.savefig(plot_path_2, dpi=300, bbox_inches="tight")
plt.close()

print("  Saved: trajectory_chr3_chr17_and_fragments.png")


# =============================================================================
# Done!
# =============================================================================
print("")
print("All plots saved to: " + OUTPUT_FOLDER)
print("Statistics saved to: " + stats_out)
print("Done!")
