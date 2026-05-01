"""
10_subcompartment_classification.py
=====================================
This script looks at WHERE in the nucleus each piece of translocated DNA
ends up, in terms of "subcompartments".

What are subcompartments?
- The nucleus has different "neighborhoods" called subcompartments.
- They're labeled A1, A2, A3, B1, B2, B3.
- A subcompartments are generally ACTIVE (genes turned on).
- B subcompartments are generally INACTIVE (genes turned off).
- Sometimes we group them into just A or B for simpler analysis.

The big question:
After a translocation, does the moved DNA:
- KEEP its original subcompartment (RETAINED)?
- ADOPT the subcompartment of its new neighborhood (ADOPTED)?
- Switch to something else entirely (OTHER)?

How we do it:
1. Load the subcompartment annotations, a file telling us what
   subcompartment each 100kb chunk of the genome belongs to.
2. For each translocated region:
   a. Find all the 100kb bins that overlap it.
   b. Find the dominant subcompartment of the new neighborhood.
   c. For each bin, classify it as retained, adopted, or other.
3. Do the same at the simpler A/B level.
4. Also track A-->B and B-->A switches specifically.
5. Build a "control" group from the rest of the genome (bins NOT in
   any translocation) to use as a background rate.
6. Run binomial tests: is each translocation switching MORE than the
   genome-wide background rate would predict?
7. Save plots and a stats table.

About the binomial test:
- A binomial test asks: "Out of N tries, I got K successes. Is that
  more (or less) than I'd expect by chance, given the background rate?"
- For us, "successes" = bins that adopted (or retained, or other).
- "Background rate" = the same fraction in the control bins.

Edit the file paths below before running.

Required libraries: pandas, numpy, matplotlib, scipy
Install them with: pip install pandas numpy matplotlib scipy
"""

import os                            
import numpy as np                   
import pandas as pd                  
import matplotlib.pyplot as plt      
from scipy.stats import binomtest    

# Increase font sizes for all plots
plt.rcParams.update({
    "font.size": 22,         # default text size
    "axes.titlesize": 30,    # plot title
    "axes.labelsize": 30,    # x and y axis labels
    "xtick.labelsize": 22,   # x tick labels
    "ytick.labelsize": 20,   # y tick labels
    "legend.fontsize": 20,   # legend text
    "legend.title_fontsize": 21  # legend title
})

# =============================================================================
# CONFIG SECTION - Edit these paths before running
# =============================================================================

# -----------------------------------------------------------------------------
# Subcompartment file (bedGraph)
# -----------------------------------------------------------------------------
# This file tells us what subcompartment each 100kb genomic bin belongs
# to in WT, T1, and C1. Expected columns include:
# chrom, start, end, MCF10A_WT.state, MCF10A_T1.state, MCF10A_C1.state
SUBCOMPARTMENT_FILE = "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph"

# -----------------------------------------------------------------------------
# BED files for each condition's translocated regions
# -----------------------------------------------------------------------------
T1_TRANSLOC_BED = "/path/to/verify_T1_translocations.bed"
C1_TRANSLOC_BED = "/path/to/verify_C1_translocations.bed"

# Group them in a dict for looping
transloc_beds = {
    "T1": T1_TRANSLOC_BED,
    "C1": C1_TRANSLOC_BED,
}

# -----------------------------------------------------------------------------
# BED files for "neighbor" regions
# -----------------------------------------------------------------------------
# - "new" neighbor: where the translocated piece ENDED UP
#   (so we can figure out what subcompartment to "adopt")
# - "old" neighbor: where the translocated piece CAME FROM
T1_NEW_NEIGHBOR = "/path/to/verify_T1_neighbor.bed"
T1_OLD_NEIGHBOR = "/path/to/verify_T1_neighbor_originchr.bed"
C1_NEW_NEIGHBOR = "/path/to/verify_C1_neighbor.bed"
C1_OLD_NEIGHBOR = "/path/to/verify_C1_neighbor_originchr.bed"

# Nested dict so we can do neighbor_beds[condition]["new"] etc.
neighbor_beds = {
    "T1": {"new": T1_NEW_NEIGHBOR, "old": T1_OLD_NEIGHBOR},
    "C1": {"new": C1_NEW_NEIGHBOR, "old": C1_OLD_NEIGHBOR},
}

# -----------------------------------------------------------------------------
# Maps from transloc_id (in the BED files) to readable names for the plots
# -----------------------------------------------------------------------------
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

# Group them for looping
name_maps = {
    "T1": T1_NAME_MAP,
    "C1": C1_NAME_MAP,
}

# -----------------------------------------------------------------------------
# Output folder for plots and tables
# -----------------------------------------------------------------------------
OUTPUT_FOLDER = "/path/to/plots/subcompartment_classification"

# -----------------------------------------------------------------------------
# A small number used to avoid divide-by-zero when calculating odds ratios
# -----------------------------------------------------------------------------
# If a category has 0 bins, we'd divide by zero.
# Adding a tiny epsilon (0.000001) avoids that.
EPSILON = 1e-6

# Conditions to process (we don't process WT here, it's just the reference)
CONDITIONS = ["T1", "C1"]


# =============================================================================
# Helper function: Find subcompartment bins overlapping a region
# =============================================================================
# Same overlap idea as in script 09:
# Two intervals [a, b) and [c, d) overlap if a < d AND b > c.

def get_overlapping_bins(sub_df, region_chrom, region_start, region_end):
    """
    Find all subcompartment bins overlapping a single region.
    Returns the rows of sub_df that overlap.
    """
    same_chrom = (sub_df["chrom"] == region_chrom)
    starts_before_end = (sub_df["start"] < region_end)
    ends_after_start = (sub_df["end"] > region_start)

    # Combine the three conditions with & (AND)
    mask = same_chrom & starts_before_end & ends_after_start

    return sub_df[mask]


# =============================================================================
# Helper function: Find the most common subcompartment in a set of bins
# =============================================================================
# This is what we use to figure out the "adopted" state, i.e. what
# subcompartment dominates the new neighborhood.

def find_dominant_state(bins_df, condition, collapsed=False):
    """
    Find the most common subcompartment state in a set of bins.
    If collapsed=True, use the A/B level instead of A1/A2/B1/B2/B3.
    Returns None if the bins DataFrame is empty.
    """
    # If there are no bins, we can't pick a dominant state
    if bins_df.empty:
        return None

    # Pick the right column name based on whether we want A/B or A1/A2/etc.
    if collapsed:
        column_name = "MCF10A_" + condition + ".state_collapsed"
    else:
        column_name = "MCF10A_" + condition + ".state"

    # value_counts() counts how many times each value appears
    # Then idxmax() returns the value with the highest count
    counts = bins_df[column_name].value_counts()
    most_common = counts.idxmax()

    return most_common


# =============================================================================
# Helper function: Classify one bin (retained / adopted / other)
# =============================================================================

def classify_one_bin(row, condition, dominant, collapsed=False):
    """
    Classify a single bin by comparing its WT state to its condition state.
    - "retained": same as in WT
    - "adopted":  changed to match the dominant state of new neighborhood
    - "other":    changed to something else

    If 'dominant' is None (used for the genome control), bins can only
    be "retained" or "other" (never "adopted").
    """
    # Pick the right column names depending on collapsed or not
    if collapsed:
        wt_col = "MCF10A_WT.state_collapsed"
        cond_col = "MCF10A_" + condition + ".state_collapsed"
    else:
        wt_col = "MCF10A_WT.state"
        cond_col = "MCF10A_" + condition + ".state"

    # Get the values for this bin
    wt_state = row[wt_col]
    cond_state = row[cond_col]

    # Classify
    if cond_state == wt_state:
        return "retained"
    elif dominant is not None and cond_state == dominant:
        return "adopted"
    else:
        return "other"


# =============================================================================
# Helper function: Get the direction of A/B switching
# =============================================================================

def get_change_direction(row, condition):
    """
    Look at the A/B status of one bin and classify it as:
    - "A_to_B"  : was A in WT, now B
    - "B_to_A"  : was B in WT, now A
    - "retained": no A/B change
    """
    wt_state = row["MCF10A_WT.state_collapsed"]
    cond_state = row["MCF10A_" + condition + ".state_collapsed"]

    if wt_state == "A" and cond_state == "B":
        return "A_to_B"
    elif wt_state == "B" and cond_state == "A":
        return "B_to_A"
    else:
        return "retained"


# =============================================================================
# Helper function: Get fractions of each category
# =============================================================================
# Used to get the genome-wide background rate (e.g. "20% of genome control
# bins are 'other'").

def get_category_fractions(df, column):
    """
    Calculate what fraction of rows fall into each category in a column.
    Returns a dict like {'retained': 0.7, 'other': 0.2, 'adopted': 0.1}.
    """
    n = len(df)

    # value_counts gives raw counts of each value
    counts = df[column].value_counts()

    # Build a dict of category --> fraction
    fractions = {}
    for category in counts.index:
        fractions[category] = counts[category] / n

    return fractions


# =============================================================================
# Helper function: Make a stacked bar plot
# =============================================================================
# Each bar represents one translocation, and the bar is divided into colored
# segments showing what fraction is "retained", "adopted", or "other".

def make_stacked_bar_plot(df, y_columns, filename, title,
                          colors=None, use_colormap=False):
    """
    Save a stacked bar plot.
    - df:           DataFrame with one row per translocation
    - y_columns:    Which columns (categories) to stack (e.g. retained/adopted/other)
    - filename:     Filename for the output PNG (just the name, not the path)
    - title:        Plot title
    - colors:       Optional list of colors (one per y_column)
    - use_colormap: If True, use matplotlib's "tab10" color scheme instead
    """
    # Make a new figure (8x5 inches)
    fig, ax = plt.subplots(figsize=(12, 8))

    # We're using pandas' built-in plotting (df.plot()).
    # It draws a bar plot on the axes we give it.
    # x="transloc_id" puts translocation names on the x-axis.
    # stacked=True stacks the categories on top of each other.
    if use_colormap:
        # Use a built-in matplotlib color scheme
        df.plot(
            x="transloc_id",
            y=y_columns,
            kind="bar",
            stacked=True,
            ax=ax,
            colormap="tab10"
        )
    else:
        # Use the colors we explicitly provided
        df.plot(
            x="transloc_id",
            y=y_columns,
            kind="bar",
            stacked=True,
            ax=ax,
            color=colors
        )

    # Add labels and title
    ax.set_ylabel("Fraction of bins")
    ax.set_xlabel("")
    ax.set_title(title)

    # Adjust spacing and save
    plt.tight_layout()
    full_path = os.path.join(OUTPUT_FOLDER, filename)
    fig.savefig(full_path, dpi=300)
    plt.close()


# =============================================================================
# Main code starts here
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# -----------------------------------------------------------------------------
# Step 1: Load the subcompartment annotations
# -----------------------------------------------------------------------------
print("Loading subcompartment annotations...")

sub = pd.read_csv(SUBCOMPARTMENT_FILE, sep="\t")

# Make sure start and end are integers
sub["start"] = sub["start"].astype(int)
sub["end"] = sub["end"].astype(int)

# -----------------------------------------------------------------------------
# Step 2: Add "collapsed" A/B columns
# -----------------------------------------------------------------------------
# The full state is something like "A1" or "B2".
# The collapsed state is just the FIRST letter: "A" or "B".
# We make collapsed columns for WT, T1, C1.
# .str[0] gets the first character of each string.

for cond in ["WT", "T1", "C1"]:
    full_col = "MCF10A_" + cond + ".state"
    collapsed_col = "MCF10A_" + cond + ".state_collapsed"
    sub[collapsed_col] = sub[full_col].str[0]


# -----------------------------------------------------------------------------
# Step 3: Process each condition
# -----------------------------------------------------------------------------
# We'll collect all bins from all translocations here for saving later
all_bins_collected = []

# We'll collect statistics from binomial tests here
stats_results = []

# Loop through each condition (T1 and C1)
for cond in CONDITIONS:

    print("")
    print("===== Processing " + cond + " =====")

    # -------------------------------------------------------------------------
    # Step 3a: Load the BED files for this condition
    # -------------------------------------------------------------------------
    trans = pd.read_csv(transloc_beds[cond], sep="\t")
    new_nb = pd.read_csv(neighbor_beds[cond]["new"], sep="\t")
    old_nb = pd.read_csv(neighbor_beds[cond]["old"], sep="\t")

    # Make sure start and end are integers in all three tables
    for one_df in (trans, new_nb, old_nb):
        one_df["start"] = one_df["start"].astype(int)
        one_df["end"] = one_df["end"].astype(int)

    # Get the name map for this condition (transloc_id --> readable label)
    trans_name_map = name_maps[cond]

    # -------------------------------------------------------------------------
    # Step 3b: Set up storage for this condition's results
    # -------------------------------------------------------------------------
    # These lists will hold one dict per translocation for plotting
    results_sub = []         # full subcompartment results
    results_collapsed = []   # collapsed A/B results
    results_direction = []   # A-->B / B-->A / retained results

    # We'll also keep a list of bins from THIS condition so we can build
    # the genome control later (everything NOT in a translocation)
    all_trans_bins_list = []

    # We need to keep the raw counts (not fractions) for the binomial test.
    # So we store them in dicts keyed by translocation ID.
    trans_counts_sub = {}
    trans_counts_collapsed = {}
    trans_counts_dir = {}

    # -------------------------------------------------------------------------
    # Step 3c: Loop through each translocation
    # -------------------------------------------------------------------------
    # Get the list of unique translocation IDs in this BED file
    unique_tids = trans["transloc_id"].unique()

    for tid in unique_tids:

        # Get the rows for this translocation
        trans_rows = trans[trans["transloc_id"] == tid]
        new_rows = new_nb[new_nb["transloc_id"] == tid]

        # ---------------------------------------------------------------------
        # Find all subcompartment bins overlapping this translocation
        # ---------------------------------------------------------------------
        # A translocation can have multiple regions. For each region, we
        # find the overlapping bins, and combine them all into one DataFrame.
        bins_pieces = []
        for _, r in trans_rows.iterrows():
            piece = get_overlapping_bins(sub, r["chrom"], r["start"], r["end"])
            bins_pieces.append(piece)

        # pd.concat sticks them all together into one DataFrame
        # ignore_index=True gives fresh row numbers
        trans_bins = pd.concat(bins_pieces, ignore_index=True)

        # ---------------------------------------------------------------------
        # Same for the new neighborhood
        # ---------------------------------------------------------------------
        new_bins_pieces = []
        for _, r in new_rows.iterrows():
            piece = get_overlapping_bins(sub, r["chrom"], r["start"], r["end"])
            new_bins_pieces.append(piece)
        new_bins = pd.concat(new_bins_pieces, ignore_index=True)

        # Drop duplicate bins (same bin can overlap multiple regions)
        trans_bins = trans_bins.drop_duplicates(subset=["chrom", "start", "end"])
        new_bins = new_bins.drop_duplicates(subset=["chrom", "start", "end"])

        # If no bins were found, skip this translocation
        if trans_bins.empty:
            print("  " + tid + ": no overlapping bins, skipping.")
            continue

        # ---------------------------------------------------------------------
        # Find the dominant subcompartment in the new neighborhood
        # ---------------------------------------------------------------------
        # This will be the "adopted" state for our classification
        new_dom = find_dominant_state(new_bins, cond, collapsed=False)
        new_dom_collapsed = find_dominant_state(new_bins, cond, collapsed=True)

        # ---------------------------------------------------------------------
        # Classify each bin
        # ---------------------------------------------------------------------
        # We use .apply() to run our classify function on each row.
        # axis=1 means "process each row" (axis=0 would be each column).
        # We use a lambda to pass extra arguments to our function.

        # First: full subcompartment level (A1 / A2 / A3 / B1 / B2 / B3)
        trans_bins["behavior"] = trans_bins.apply(
            lambda row: classify_one_bin(row, cond, new_dom, collapsed=False),
            axis=1
        )

        # Second: collapsed A/B level
        trans_bins["behavior_collapsed"] = trans_bins.apply(
            lambda row: classify_one_bin(row, cond, new_dom_collapsed, collapsed=True),
            axis=1
        )

        # Third: direction of switch (A_to_B / B_to_A / retained)
        trans_bins["change_direction"] = trans_bins.apply(
            lambda row: get_change_direction(row, cond),
            axis=1
        )

        # Add columns identifying which translocation and condition this is
        trans_bins["transloc_id"] = tid
        trans_bins["condition"] = cond

        # Save for the combined output later
        all_bins_collected.append(trans_bins)

        # Save just the position columns for the genome-control logic later
        position_cols = trans_bins[["chrom", "start", "end"]]
        all_trans_bins_list.append(position_cols)

        # ---------------------------------------------------------------------
        # Calculate counts and fractions for plotting and statistics
        # ---------------------------------------------------------------------
        n = len(trans_bins)

        # Raw counts (used for binomial tests later)
        counts_sub = trans_bins["behavior"].value_counts()
        counts_collapsed = trans_bins["behavior_collapsed"].value_counts()
        counts_dir = trans_bins["change_direction"].value_counts()

        trans_counts_sub[tid] = counts_sub
        trans_counts_collapsed[tid] = counts_collapsed
        trans_counts_dir[tid] = counts_dir

        # Fractions (used for plotting)
        # We build a dict starting with the transloc_id, then add each
        # category's fraction.
        # The ** operator unpacks a dict into keyword args (like spreading
        # in JavaScript). Here we use it to merge dicts.
        # (counts / n) gives a Series; .to_dict() turns it into a dict.

        fractions_sub_dict = (counts_sub / n).to_dict()
        row_for_sub = {"transloc_id": tid}
        row_for_sub.update(fractions_sub_dict)  # add all the categories
        results_sub.append(row_for_sub)

        fractions_collapsed_dict = (counts_collapsed / n).to_dict()
        row_for_collapsed = {"transloc_id": tid}
        row_for_collapsed.update(fractions_collapsed_dict)
        results_collapsed.append(row_for_collapsed)

        fractions_dir_dict = (counts_dir / n).to_dict()
        row_for_dir = {"transloc_id": tid}
        row_for_dir.update(fractions_dir_dict)
        results_direction.append(row_for_dir)

    # -------------------------------------------------------------------------
    # Step 3d: Build the genome-wide control group
    # -------------------------------------------------------------------------
    # The control = all bins NOT in any translocated region.
    # This gives us a background rate of switching.

    # Combine all translocation bin positions and remove duplicates
    all_trans_bins_df = pd.concat(all_trans_bins_list).drop_duplicates()

    # Use merge with how="left" and indicator=True to find rows in `sub`
    # that DON'T have a match in `all_trans_bins_df`.
    # The indicator column tells us if each row was in:
    #   "left_only" (only in sub), this is what we want for control
    #   "right_only" (only in all_trans_bins_df)
    #   "both" (in both)
    merged = sub.merge(
        all_trans_bins_df,
        on=["chrom", "start", "end"],
        how="left",
        indicator=True
    )

    # Keep only rows that are NOT in any translocation
    control_bins = merged[merged["_merge"] == "left_only"]

    # Drop the indicator column, we don't need it anymore
    control_bins = control_bins.drop(columns="_merge")

    # ---------------------------------------------------------------------
    # Classify the control bins
    # ---------------------------------------------------------------------
    # Note: dominant=None because the control bins don't have a defined
    # "new neighborhood". So they can only be "retained" or "other" -
    # never "adopted".

    control_bins["behavior"] = control_bins.apply(
        lambda row: classify_one_bin(row, cond, None, collapsed=False),
        axis=1
    )

    control_bins["behavior_collapsed"] = control_bins.apply(
        lambda row: classify_one_bin(row, cond, None, collapsed=True),
        axis=1
    )

    control_bins["change_direction"] = control_bins.apply(
        lambda row: get_change_direction(row, cond),
        axis=1
    )

    # Get the background fractions for the binomial test
    control_p_sub = get_category_fractions(control_bins, "behavior")
    control_p_collapsed = get_category_fractions(control_bins, "behavior_collapsed")
    control_p_dir = get_category_fractions(control_bins, "change_direction")

    # -------------------------------------------------------------------------
    # Step 3e: Run binomial tests for each translocation and category
    # -------------------------------------------------------------------------
    # For each translocation x category, ask: is the observed count
    # significantly different from what we'd expect given the background rate?

    # Loop through each translocation
    for tid in trans_counts_sub:

        # Get the total number of bins in this translocation
        n = int(trans_counts_sub[tid].sum())

        # We have three "levels" of analysis to run tests on.
        # For each level, we have:
        #   - the actual counts (counts_dict)
        #   - the background rates (control_p)
        # We loop through the three levels.
        levels_to_test = [
            ("subcompartment", trans_counts_sub[tid],       control_p_sub),
            ("collapsed",      trans_counts_collapsed[tid], control_p_collapsed),
            ("direction",      trans_counts_dir[tid],       control_p_dir),
        ]

        for level_name, counts_dict, control_p in levels_to_test:

            # Loop through each category in this level
            # (e.g. "retained", "adopted", "other" for the subcompartment level)
            for cat in counts_dict.index:

                # k = number of bins in this category for this translocation
                k = int(counts_dict.get(cat, 0))

                # p0 = background probability of this category (from control)
                p0 = control_p.get(cat, 0)

                # Run the binomial test
                # Question: out of n trials, getting k successes, when expected
                # probability is p0, is this surprising?
                test = binomtest(k, n, p=p0)

                # ---------------------------------------------------------
                # Calculate the log2 odds ratio
                # ---------------------------------------------------------
                # Odds = P / (1 - P), where P is a probability.
                # The odds ratio compares two odds.
                # log2 of the odds ratio is positive if the translocation
                # is enriched (more than background), negative if depleted.
                # We add EPSILON to avoid divide-by-zero if any value is 0.
                odds_trans = (k + EPSILON) / (n - k + EPSILON)
                odds_control = (p0 + EPSILON) / (1 - p0 + EPSILON)
                log2_or = np.log2(odds_trans / odds_control)

                # Get the readable name for this translocation
                # If tid isn't in the map, use tid itself as a fallback
                readable_name = trans_name_map.get(tid, tid)

                # Store all the info in a dict
                one_stat_row = {
                    "condition": cond,
                    "translocation": readable_name,
                    "level": level_name,
                    "category": cat,
                    "n_bins": n,
                    "k_bins": k,
                    "fraction": k / n,
                    "expected_fraction": p0,
                    "p_value": test.pvalue,
                    "log2_odds_ratio": log2_or,
                }
                stats_results.append(one_stat_row)

    # -------------------------------------------------------------------------
    # Step 3f: Add a "Genome control" bar to the plots
    # -------------------------------------------------------------------------
    # The control group provides a reference bar in the bar plots so the
    # reader can see if the translocations differ from the genome background.

    # value_counts(normalize=True) gives fractions instead of raw counts
    control_fractions_sub = control_bins["behavior"].value_counts(normalize=True).to_dict()
    control_row_sub = {"transloc_id": "Genome control"}
    control_row_sub.update(control_fractions_sub)
    results_sub.append(control_row_sub)

    control_fractions_collapsed = control_bins["behavior_collapsed"].value_counts(normalize=True).to_dict()
    control_row_collapsed = {"transloc_id": "Genome control"}
    control_row_collapsed.update(control_fractions_collapsed)
    results_collapsed.append(control_row_collapsed)

    control_fractions_dir = control_bins["change_direction"].value_counts(normalize=True).to_dict()
    control_row_dir = {"transloc_id": "Genome control"}
    control_row_dir.update(control_fractions_dir)
    results_direction.append(control_row_dir)

    # -------------------------------------------------------------------------
    # Step 3g: Convert results lists to DataFrames and apply readable names
    # -------------------------------------------------------------------------
    # We've been storing results as lists of dicts. Now convert to DataFrames
    # and replace transloc_ids with readable names.

    # Subcompartment level
    final_sub = pd.DataFrame(results_sub).fillna(0)
    # If the transloc_id is in our name map, replace it with the readable name
    # Otherwise keep it as-is. We use a list comprehension for clarity.
    new_names = []
    for tid in final_sub["transloc_id"]:
        new_names.append(trans_name_map.get(tid, tid))
    final_sub["transloc_id"] = new_names

    # Collapsed level
    final_collapsed = pd.DataFrame(results_collapsed).fillna(0)
    new_names = []
    for tid in final_collapsed["transloc_id"]:
        new_names.append(trans_name_map.get(tid, tid))
    final_collapsed["transloc_id"] = new_names

    # Direction level
    final_direction = pd.DataFrame(results_direction).fillna(0)
    new_names = []
    for tid in final_direction["transloc_id"]:
        new_names.append(trans_name_map.get(tid, tid))
    final_direction["transloc_id"] = new_names

    # -------------------------------------------------------------------------
    # Step 3h: Save the three plots for this condition
    # -------------------------------------------------------------------------

    # Plot 1: full subcompartments
    make_stacked_bar_plot(
        df=final_sub,
        y_columns=["retained", "adopted", "other"],
        filename=cond + "_subcompartment.png",
        title=cond + " vs WT: Subcompartments",
        use_colormap=True
    )

    # Plot 2: collapsed A/B
    make_stacked_bar_plot(
        df=final_collapsed,
        y_columns=["retained", "adopted", "other"],
        filename=cond + "_collapsed.png",
        title=cond + " vs WT: Collapsed A/B",
        use_colormap=True
    )

    # Plot 3: direction of switching
    # Using specific colors here: red for A-->B, blue for B-->A, grey for retained
    make_stacked_bar_plot(
        df=final_direction,
        y_columns=["A_to_B", "B_to_A", "retained"],
        filename=cond + "_direction.png",
        title=cond + " vs WT: A↔B direction",
        colors=["#FF0800", "#1A43BF", "#cccccc"]
    )

    print("  Plots saved for " + cond)


# =============================================================================
# Step 4: Save combined output files
# =============================================================================

# Combine all the per-bin tables from all conditions into one big DataFrame
all_bins_df = pd.concat(all_bins_collected)
bins_out = os.path.join(OUTPUT_FOLDER, "bins_annotation.csv")
all_bins_df.to_csv(bins_out, index=False)

# Save the binomial test results
stats_df = pd.DataFrame(stats_results)
stats_out = os.path.join(OUTPUT_FOLDER, "all_binomial_stats.csv")
stats_df.to_csv(stats_out, index=False)

print("")
print("  Bin annotations: " + bins_out)
print("  Statistics:      " + stats_out)
print("Done!")
