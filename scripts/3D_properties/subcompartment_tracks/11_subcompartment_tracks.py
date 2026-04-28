"""
11_subcompartment_tracks.py
============================
This script makes pictures showing the subcompartment "color pattern" along
each translocated chromosome.

Recap on subcompartments:
- Subcompartments are different "neighborhoods" along the genome.
- They're labeled A0, A1, A2, A3 (more active) and B0, B1, B2, B3 (less active).
- Each 100 kb chunk of the genome is assigned to one subcompartment.

What this script makes:
For each translocation, it makes a figure with 4 horizontal "tracks":
  Track 1: the full original chromosome A (where the fragment came from)
  Track 2: the full original chromosome B (where the fragment ended up)
  Track 3: the new derivative chromosome (fragment + body joined together,
           with a green line marking the breakpoint)
  Track 4: a legend showing what each color means

Each track is a long colored bar. The colors come from the subcompartment
of each 100 kb chunk:
- Reds = A subcompartments (active)
- Blues = B subcompartments (inactive)

A bit about translocations:
- Each translocation has TWO pieces: a "fragment" (smaller) and a
  "body" (larger).
- They get joined together to form a new "derivative" chromosome.
- We figure out which is which by comparing their lengths.

We make one figure per translocation per condition (WT, T1, C1).
For example, t(3;17) gets 3 figures (one per condition), but t(2;10)
only gets 2 (because it doesn't exist in T1).

Edit the file paths below before running.

Required libraries: pandas, matplotlib
Install them with: pip install pandas matplotlib
"""

import os                         
import pandas as pd                
import matplotlib.pyplot as plt    


# =============================================================================
# CONFIG SECTION - Edit these paths and settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# Subcompartment file (bedGraph)
# -----------------------------------------------------------------------------
# Tells us what subcompartment each 100 kb genomic bin belongs to,
# in WT, T1, and C1.
SUBCOMPARTMENT_FILE = "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph"

# -----------------------------------------------------------------------------
# BED files for each condition's translocations
# -----------------------------------------------------------------------------
# Each BED file row is one segment of a translocation.
# Two rows with the same transloc_id form one translocation
# (one fragment + one body).
T1_TRANSLOC_BED = "/path/to/verify_T1_translocations.bed"
C1_TRANSLOC_BED = "/path/to/verify_C1_translocations.bed"

# Group them in a dict for looping
transloc_beds = {
    "T1": T1_TRANSLOC_BED,
    "C1": C1_TRANSLOC_BED,
}

# -----------------------------------------------------------------------------
# Conditions to plot, with their column names in the subcompartment file
# -----------------------------------------------------------------------------
# We need this map because the column names in the file aren't just
# "WT", "T1", "C1", they have a longer prefix.
condition_columns = {
    "WT": "MCF10A_WT.state",
    "T1": "MCF10A_T1.state",
    "C1": "MCF10A_C1.state",
}

# -----------------------------------------------------------------------------
# Maps from transloc_id --> human-readable name
# -----------------------------------------------------------------------------
T1_NAME_MAP = {
    "T1": "Der(17)t(3;17)",   # chr3 fragment goes onto chr17 body
    "T2": "Der(3)t(3;17)",    # chr17 fragment goes onto chr3 body
    "T3": "t(6;19)",          # chr19 fragment goes onto chr6 body
}
C1_NAME_MAP = {
    "T1": "t(2;10)",          # chr2 fragment goes onto chr10 body
    "T2": "Der(17)t(3;17)",
    "T3": "Der(3)t(3;17)",
    "T4": "t(6;19)",
}

# Group them for looping
name_maps = {
    "T1": T1_NAME_MAP,
    "C1": C1_NAME_MAP,
}

# -----------------------------------------------------------------------------
# Which conditions to plot for each translocation
# -----------------------------------------------------------------------------
# Some translocations don't exist in some cell lines, so we don't try to
# plot those. For example, t(2;10) is only in C1, so we plot WT and C1
# (skipping T1).
transloc_conditions = {
    "t(2;10)":        ["WT", "C1"],
    "Der(17)t(3;17)": ["WT", "T1", "C1"],
    "Der(3)t(3;17)":  ["WT", "T1", "C1"],
    "t(6;19)":        ["WT", "T1", "C1"],
}

# -----------------------------------------------------------------------------
# Output folder for plots
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/results/plots/subcompartment_tracks"

# -----------------------------------------------------------------------------
# Colour scheme: shades of red for A, shades of blue for B
# -----------------------------------------------------------------------------
# Lower numbers = lighter shades, higher numbers = darker shades
subcomp_colors = {
    "A0": "#fcaeae",   # light red
    "A1": "#fb6a6a",
    "A2": "#ef3b2c",
    "A3": "#99000d",   # dark red
    "B0": "#deebf7",   # light blue
    "B1": "#9ecae1",
    "B2": "#4292c6",
    "B3": "#084594",   # dark blue
}


# =============================================================================
# Helper function: Load a BED file and pair up its rows
# =============================================================================
# Each translocation has 2 rows in the BED file (sharing a transloc_id).
# We figure out which one is the "fragment" (shorter) and which is the
# "body" (longer) by looking at how long each row is.

def load_translocation_pairs(bed_path, name_map):
    """
    Load a translocation BED file and pair up its rows.
    Returns a list of dicts, each describing one translocation
    with keys: name, transloc_id, fragment, body
    """
    # Read the BED file (tab-separated)
    bed = pd.read_csv(bed_path, sep="\t")

    # Make sure start and end are integers
    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["end"].astype(int)

    # Add a column that tells us how long each row is
    bed["span"] = bed["end"] - bed["start"]

    # We'll collect each translocation as a dict in this list
    pairs = []

    # Get the unique translocation IDs
    unique_tids = bed["transloc_id"].unique()

    # Loop through each translocation ID
    for tid in unique_tids:

        # Get the rows for this translocation ID
        rows = bed[bed["transloc_id"] == tid]

        # We expect exactly 2 rows per translocation.
        # If we don't get 2, print a warning and skip.
        if len(rows) != 2:
            print("  Warning: transloc_id " + tid + " has "
                  + str(len(rows)) + " rows, expected 2 — skipping.")
            continue

        # Sort by span (length) so the shorter row comes first
        rows_sorted = rows.sort_values("span")

        # The first row (shortest) is the fragment, the second is the body.
        # .iloc[0] and .iloc[1] get the rows by position.
        fragment = rows_sorted.iloc[0]
        body = rows_sorted.iloc[1]

        # Get the readable name from the map.
        # If tid isn't in the map, just use tid itself as a fallback.
        readable_name = name_map.get(tid, tid)

        # Build the dict for this translocation
        one_pair = {
            "name": readable_name,
            "transloc_id": tid,
            "fragment": fragment,
            "body": body,
        }

        pairs.append(one_pair)

    return pairs


# =============================================================================
# Helper function: Get bins on a chromosome (optionally within a range)
# =============================================================================
# Used to extract a slice of the subcompartment data for plotting.

def get_chrom_bins(sub_df, chrom, start=None, end=None):
    """
    Return all subcompartment bins on the given chromosome.
    If start and/or end are provided, also filter to that range.
    """
    # Start with a mask that just checks the chromosome
    mask = (sub_df["chrom"] == chrom)

    # If start was provided, also require bin start >= start
    if start is not None:
        mask = mask & (sub_df["start"] >= start)

    # If end was provided, also require bin end <= end
    if end is not None:
        mask = mask & (sub_df["end"] <= end)

    return sub_df[mask]


# =============================================================================
# Helper function: Draw one subcompartment track (one horizontal bar)
# =============================================================================
# A "track" is a long horizontal bar made of many small colored rectangles,
# one per 100 kb bin. Each rectangle's color shows the bin's subcompartment.
#
# A track can be made of one OR two segments. The derivative chromosome
# track is made of two segments (the fragment + the body).

def plot_track(ax, segments, condition_col, title, x_label_suffix=""):
    """
    Draw a horizontal subcompartment track on the given axes.

    Parameters:
    - ax:              the matplotlib axes to draw on
    - segments:        a list of (DataFrame, original_start) tuples.
                       Each tuple is one piece of genome to draw.
                       original_start is the genomic position where the
                       segment STARTS. We need this to translate genomic
                       coords into plot coords.
    - condition_col:   the column name with the subcompartment states
    - title:           the plot title
    - x_label_suffix:  extra text to add to the x-axis label
    """
    # cursor = current x position on the plot.
    # We start at 0 and move it forward as we draw segments.
    cursor = 0

    # Loop through each segment to draw
    for df_seg, original_start in segments:

        # Make a copy so we don't modify the original DataFrame
        df_seg = df_seg.copy()

        # ---------------------------------------------------------------------
        # Figure out where each bin should be drawn on the plot
        # ---------------------------------------------------------------------
        # The bin's "start" and "end" are in genomic coordinates.
        # But on the plot, we want the segment to START at our cursor.
        # So we shift each bin's coordinates by (cursor - original_start).
        # If cursor=0 and original_start=10000, shift = -10000:
        #   a bin from 10000-20000 becomes 0-10000 on the plot.
        shift = cursor - original_start
        df_seg["plot_start"] = df_seg["start"] + shift
        df_seg["plot_end"] = df_seg["end"] + shift

        # ---------------------------------------------------------------------
        # Draw each bin as a colored rectangle
        # ---------------------------------------------------------------------
        for _, row in df_seg.iterrows():

            # Look up the color for this bin's subcompartment
            # If the state isn't in our color dict, use grey ("#cccccc")
            state = row[condition_col]
            color = subcomp_colors.get(state, "#cccccc")

            # fill_between fills a region between two y values.
            # We give it the x values [plot_start, plot_end] and fill
            # from y=0 to y=1, making a colored rectangle.
            # edgecolor="none" prevents thin black lines around each bin.
            ax.fill_between(
                [row["plot_start"], row["plot_end"]],
                0,
                1,
                color=color,
                edgecolor="none"
            )

        # Move the cursor to the end of this segment so the next segment
        # (if any) starts where this one ended.
        seg_length = df_seg["end"].max() - original_start
        cursor = cursor + seg_length

    # -------------------------------------------------------------------------
    # Set up the plot's axes and title
    # -------------------------------------------------------------------------
    # X-axis: from 0 to total length (where the cursor ended up)
    ax.set_xlim(0, cursor)

    # Y-axis: from 0 to 1 (the height of the colored rectangles)
    ax.set_ylim(0, 1)

    # Hide the y-axis tick marks (they don't mean anything for tracks)
    ax.set_yticks([])

    # Add the title and x-axis label
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Position (bp)" + x_label_suffix)

    # -------------------------------------------------------------------------
    # If we drew TWO segments, mark the breakpoint with a green vertical line
    # -------------------------------------------------------------------------
    if len(segments) == 2:

        # The breakpoint is where segment 1 ends.
        # segments[0] is the first (fragment) segment.
        # segments[0][0] is the DataFrame, segments[0][1] is original_start.
        seg1_df, seg1_original_start = segments[0]
        seg1_length = seg1_df["end"].max() - seg1_original_start

        # axvline draws a vertical line at the given x position
        ax.axvline(
            x=seg1_length,
            color="green",
            linestyle="-",
            linewidth=5,
            label="Breakpoint"
        )

        # Show a small legend for the breakpoint
        ax.legend(fontsize=8)


# =============================================================================
# Helper function: Draw the legend track
# =============================================================================
# This makes the bottom track of each figure, a colored legend showing
# which color corresponds to which subcompartment.

def add_legend(ax):
    """
    Draw a legend showing all subcompartment colors.
    The axes itself is hidden, we just use it as a place for the legend.
    """
    # We "draw" empty lines (with no data) just so they appear in the legend
    # with the right color. ax.plot([], [], ...) makes an invisible line.
    for state, color in subcomp_colors.items():
        ax.plot([], [], color=color, label=state, linewidth=6)

    # Add the legend, positioned to the right of the plot
    # bbox_to_anchor=(1.05, 1) means: place at x=1.05, y=1 in axes units
    # (so just to the right of the axis, at the top)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    # Hide the axes frame and ticks
    ax.axis("off")


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# -----------------------------------------------------------------------------
# Step 1: Load the subcompartment data
# -----------------------------------------------------------------------------
print("Loading subcompartment file...")

sub = pd.read_csv(SUBCOMPARTMENT_FILE, sep="\t")
print("Loaded subcompartment data: " + str(len(sub)) + " bins")

# -----------------------------------------------------------------------------
# Step 2: Load all translocation pairs from both BED files
# -----------------------------------------------------------------------------
# We collect them into one big dictionary keyed by translocation NAME
# (e.g. "Der(17)t(3;17)"), so we don't process the same translocation
# twice if it appears in both T1 and C1 BED files.

print("")
print("Loading translocation pairs...")

all_pairs = {}  # key = readable name, value = pair dict

# Loop through each condition's BED file
for cond_key, bed_path in transloc_beds.items():

    # Get the name map for this condition
    this_name_map = name_maps[cond_key]

    # Load and pair the rows in this BED file
    pairs = load_translocation_pairs(bed_path, this_name_map)

    # Add each pair to our master dict, only if we haven't seen it yet
    for p in pairs:
        pair_name = p["name"]

        if pair_name not in all_pairs:
            all_pairs[pair_name] = p

            # Print a summary so we can verify what was loaded.
            # The {:,} format adds thousands separators (e.g. 22,700,000).
            frag = p["fragment"]
            body = p["body"]

            summary = (
                "  Loaded: " + pair_name + "  "
                + "fragment=" + frag["chrom"] + ":"
                + "{:,}".format(frag["start"]) + "–"
                + "{:,}".format(frag["end"]) + "  "
                + "body=" + body["chrom"] + ":"
                + "{:,}".format(body["start"]) + "–"
                + "{:,}".format(body["end"])
            )
            print(summary)

# -----------------------------------------------------------------------------
# Step 3: Make figures - one per translocation per condition
# -----------------------------------------------------------------------------
# Loop through each translocation we loaded
for trans_name, pair in all_pairs.items():

    # Pull out the fragment and body info for convenience
    frag = pair["fragment"]
    body = pair["body"]

    frag_chrom = frag["chrom"]
    body_chrom = body["chrom"]

    # -------------------------------------------------------------------------
    # Get the subcompartment data we need for this translocation
    # -------------------------------------------------------------------------
    # 1. The full original fragment chromosome (for context, track 1)
    frag_chrom_all = get_chrom_bins(sub, frag_chrom)

    # 2. The full original body chromosome (for context, track 2)
    body_chrom_all = get_chrom_bins(sub, body_chrom)

    # 3. Just the fragment piece (for the derivative chromosome, track 3)
    frag_bins = get_chrom_bins(sub, frag_chrom, frag["start"], frag["end"])

    # 4. Just the body piece (for the derivative chromosome, track 3)
    body_bins = get_chrom_bins(sub, body_chrom, body["start"], body["end"])

    # -------------------------------------------------------------------------
    # Figure out which conditions to plot for this translocation
    # -------------------------------------------------------------------------
    # If we have an entry for this translocation in our config, use it.
    # Otherwise, default to plotting all conditions.
    if trans_name in transloc_conditions:
        conds_to_plot = transloc_conditions[trans_name]
    else:
        # list() turns the dict's keys into a list
        conds_to_plot = list(condition_columns.keys())

    # -------------------------------------------------------------------------
    # Make one figure per condition
    # -------------------------------------------------------------------------
    for cond_name in conds_to_plot:

        # Look up the column name for this condition
        col = condition_columns.get(cond_name)

        # Skip if this condition isn't in our config
        if col is None:
            print("  Skipping " + cond_name + ": not found in conditions config.")
            continue

        print("")
        print("Plotting " + trans_name + " — " + cond_name)

        # ---------------------------------------------------------------------
        # Create the figure with 4 stacked subplots (axes)
        # ---------------------------------------------------------------------
        # plt.subplots(4, 1, ...) creates a figure with 4 rows and 1 column.
        # axes is a list of 4 Axes objects we can draw on.
        # figsize is the size in inches.
        fig, axes = plt.subplots(4, 1, figsize=(14, 9))

        # Add the overall figure title
        fig.suptitle(
            "Subcompartment tracks " + trans_name + " — " + cond_name,
            fontsize=13,
            fontweight="bold"
        )

        # ---------------------------------------------------------------------
        # Track 1 (axes[0]): the full original fragment chromosome
        # ---------------------------------------------------------------------
        # We pass a list with ONE segment (the whole chromosome, starting at 0)
        plot_track(
            ax=axes[0],
            segments=[(frag_chrom_all, 0)],
            condition_col=col,
            title=frag_chrom + " (original) — " + cond_name,
            x_label_suffix=" [" + frag_chrom + "]"
        )

        # ---------------------------------------------------------------------
        # Track 2 (axes[1]): the full original body chromosome
        # ---------------------------------------------------------------------
        plot_track(
            ax=axes[1],
            segments=[(body_chrom_all, 0)],
            condition_col=col,
            title=body_chrom + " (original) — " + cond_name,
            x_label_suffix=" [" + body_chrom + "]"
        )

        # ---------------------------------------------------------------------
        # Track 3 (axes[2]): the derivative chromosome (fragment + body)
        # ---------------------------------------------------------------------
        # This one is made of TWO segments joined together.
        # We pass each segment with its original genomic start position.

        # Convert positions to megabases (Mb) for the title text
        # We use // for integer division (drops the decimal)
        # and / for normal division (keeps decimals).
        # 1000000 is one million.
        frag_mb_start_int = frag["start"] // 1000000
        frag_mb_end = frag["end"] / 1000000
        body_mb_start = body["start"] / 1000000
        body_mb_end = body["end"] / 1000000

        # Build the title using string concatenation
        # {:.1f} = 1 decimal place
        track3_title = (
            "Derivative " + trans_name + " — " + cond_name + " | "
            + "[" + frag_chrom + ":"
            + str(frag_mb_start_int) + "–"
            + "{:.1f}".format(frag_mb_end) + "Mb | "
            + body_chrom + ":"
            + "{:.1f}".format(body_mb_start) + "–"
            + "{:.1f}".format(body_mb_end) + "Mb]"
        )

        plot_track(
            ax=axes[2],
            segments=[
                (frag_bins, frag["start"]),   # fragment first
                (body_bins, body["start"]),   # then body
            ],
            condition_col=col,
            title=track3_title
        )

        # ---------------------------------------------------------------------
        # Track 4 (axes[3]): the legend
        # ---------------------------------------------------------------------
        add_legend(axes[3])

        # ---------------------------------------------------------------------
        # Adjust spacing and save the figure
        # ---------------------------------------------------------------------
        plt.tight_layout()

        # Build a safe filename, we replace characters that don't play nice
        # in filenames (semicolons, parentheses)
        safe_name = trans_name.replace(";", "_")
        safe_name = safe_name.replace("(", "")
        safe_name = safe_name.replace(")", "")

        # Build the full output path
        filename = "subcompartment_" + safe_name + "_" + cond_name + ".png"
        out_path = os.path.join(OUTPUT_FOLDER, filename)

        # Save and close
        # bbox_inches="tight" trims any extra whitespace around the figure
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

        print("  Saved: " + out_path)

# Done!
print("")
print("Done!")
