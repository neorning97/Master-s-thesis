"""
08_chromosome_centroid_distance.py
====================================
This script measures how far apart TWO whole chromosomes are in 3D,
and compares that distance across three conditions (WT, T1, C1).

The big question:
- A translocation physically joins two chromosomes together.
- Are those two chromosomes spatially closer to each other in
  translocated cells (T1, C1) than in normal cells (WT)?
- They had to be near each other for the translocation to happen,
  and they might still be close afterwards.

How we answer it:
1. For each .cmm 3D model file, read all the bead positions.
2. For each chromosome, calculate its "centroid", the average position
   of all its beads. This gives one (x, y, z) point that represents
   the whole chromosome.
3. Calculate the distance between the two chromosomes' centroids.
4. Repeat across all models for each condition (WT, T1, C1).
5. Make a box plot to compare the three conditions.
6. Run a Mann-Whitney U test to check if the differences are statistically
   significant. (We use this test instead of a t-test because the data
   isn't necessarily normally distributed.)

Edit the file paths and chromosome names below before running.

Required libraries: pandas, numpy, matplotlib, seaborn, scipy
Install them with: pip install pandas numpy matplotlib seaborn scipy
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
# CONFIG SECTION - Edit these settings before running the script
# =============================================================================

# The two chromosomes involved in the translocation.
# Change these to match your specific translocation.
CHR_A = "chr3"
CHR_B = "chr17"

# Folders containing the .cmm 3D model files for each condition
WT_FOLDER = "/path/to/cmm/WT"
T1_FOLDER = "/path/to/cmm/T1"
C1_FOLDER = "/path/to/cmm/C1"

# Where to save the output plot
OUTPUT_FOLDER = "/path/to/results/plots/chromosome_centroid_distance"

# Put the conditions and folders in a dictionary so we can loop through them.
# The keys (WT, T1, C1) are the condition names, and the values are the folders.
condition_folders = {
    "WT": WT_FOLDER,
    "T1": T1_FOLDER,
    "C1": C1_FOLDER,
}


# =============================================================================
# Helper function: Read a .cmm file and pull out the bead positions
# =============================================================================
# This is similar to scripts 05 and 06, but a bit simpler because here
# we only need the chromosome and 3D coordinates. We don't need the
# start/end positions (we're working at the chromosome level, not gene level).

def read_cmm_file(file_path):
    """
    Open a .cmm file and read the bead positions.
    Returns a pandas DataFrame with columns: chr, x, y, z
    """
    # Empty list to fill with bead info
    bead_list = []

    # Open and read the file line by line
    with open(file_path, "r") as my_file:
        for line in my_file:

            # Skip lines that don't have a marker
            if "<marker" not in line:
                continue

            # Find all the attribute="value" pairs in the line and turn
            # them into a dictionary
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            # If there's no beadID, we can't use this marker
            if "beadID" not in attrs:
                continue

            # Try to extract what we need (use try/except for safety)
            try:
                # The beadID looks like "chr1:1000-2000".
                # We split on ":" to get "chr1" and "1000-2000".
                # We don't actually need the "1000-2000" part here,
                # so we use _ to throw it away.
                bead_id = attrs["beadID"]
                chrom_part, _ = bead_id.split(":")

                # Make sure chromosome name starts with "chr"
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part

                # Get the 3D coordinates as floats (decimal numbers)
                # If they're missing, use NaN (not a number)
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))

                # Add to our list
                bead_list.append([chrom_part, x, y, z])

            except Exception:
                # Skip any badly-formatted lines and keep going
                continue

    # Turn the list into a DataFrame and return it
    beads_df = pd.DataFrame(bead_list, columns=["chr", "x", "y", "z"])
    return beads_df


# =============================================================================
# Helper function: Compute the centroid of one chromosome
# =============================================================================
# The "centroid" of a chromosome is the average position of all its beads.
# It's like asking "where is the middle of all these dots in 3D space?"
# We get this by averaging the x values, y values, and z values separately.

def compute_centroid(beads_df, chrom):
    """
    Calculate the centroid (mean 3D position) of all beads belonging
    to a particular chromosome.
    Returns a numpy array [x, y, z], or None if the chromosome
    has no beads in this model.
    """
    # Get only the beads on this chromosome
    sub = beads_df[beads_df["chr"] == chrom]

    # If there are no beads for this chromosome, return None
    # (the calling code will check for None and skip this model)
    if sub.empty:
        return None

    # Get just the x, y, z columns
    coords = sub[["x", "y", "z"]]

    # Calculate the mean of each column. .mean() gives us a Series
    # with three values (one for each axis).
    means = coords.mean()

    # Convert it to a numpy array so we can do math with it more easily
    centroid = means.values

    return centroid


# =============================================================================
# Main code - Step 1: Calculate centroid distances for each model
# =============================================================================
# We'll go through each condition's folder, read each model file,
# calculate the distance between the two chromosomes' centroids,
# and store the result.

# A list to collect our results, each entry will be a small dictionary
# with the condition, distance, and model filename.
results_list = []

# Loop through each condition (WT, T1, C1)
# .items() gives us both the key (condition name) and value (folder path)
for cond_name, cmm_folder in condition_folders.items():

    # Find all .cmm files in this condition's folder, sorted alphabetically
    search_pattern = os.path.join(cmm_folder, "*.cmm")
    cmm_files = sorted(glob.glob(search_pattern))

    print("")
    print("Processing " + cond_name + ": " + str(len(cmm_files)) + " models found")

    # Loop through each model file
    for f in cmm_files:

        # Get just the filename for nicer messages
        file_name = os.path.basename(f)

        # Read the .cmm file using our helper function
        beads = read_cmm_file(f)

        # Skip if the file had no beads
        if beads.empty:
            continue

        # Compute the centroid of each of the two chromosomes
        centroid_A = compute_centroid(beads, CHR_A)
        centroid_B = compute_centroid(beads, CHR_B)

        # If either chromosome wasn't found in this model, skip it
        # (Sometimes a chromosome might be missing from a model)
        if centroid_A is None or centroid_B is None:
            print("  Skipping " + file_name + ": missing beads for "
                  + CHR_A + " or " + CHR_B)
            continue

        # ---------------------------------------------------------------------
        # Calculate the distance between the two centroids
        # ---------------------------------------------------------------------
        # The distance between two points (x1, y1, z1) and (x2, y2, z2) is:
        #     sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
        # That's the same as the length (norm) of the vector that goes
        # from one point to the other, i.e. the difference between them.
        # np.linalg.norm does this calculation for us.
        difference = centroid_A - centroid_B
        dist = np.linalg.norm(difference)

        # Build a small dictionary with this result and add it to our list
        one_result = {
            "Condition": cond_name,
            "Distance": dist,
            "Model": file_name,
        }
        results_list.append(one_result)

# Turn our list of dictionaries into a DataFrame.
# Each dictionary becomes a row, and the keys become column names.
df = pd.DataFrame(results_list)

print("")
print("Total: " + str(len(df)) + " model-condition pairs computed")

# Print some summary statistics for each condition.
# .groupby("Condition") groups the rows by condition.
# Then .describe() gives count, mean, std, min, max, etc.
# .round(3) rounds to 3 decimal places so it's easier to read.
print(df.groupby("Condition")["Distance"].describe().round(3))


# =============================================================================
# Main code - Step 2: Make the box plot
# =============================================================================

# Make sure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# We want the conditions to appear in a specific order on the plot:
# WT first, then T1, then C1. But only include conditions that
# actually have data (in case one of them was empty).
# This is a list comprehension, it builds a list by checking each value.
preferred_order = ["WT", "T1", "C1"]
conditions_in_data = df["Condition"].unique()

condition_order = []
for c in preferred_order:
    if c in conditions_in_data:
        condition_order.append(c)

# -----------------------------------------------------------------------------
# Create the plot
# -----------------------------------------------------------------------------
# fig is the whole figure (the canvas), ax is the actual plot area.
# figsize=(6, 6) makes it 6x6 inches.
fig, ax = plt.subplots(figsize=(6, 6))

# Draw the box plot.
# A box plot shows the median, quartiles, and outliers for each group.
sns.boxplot(
    data=df,
    x="Condition",
    y="Distance",
    order=condition_order,
    ax=ax,
)

# Overlay the individual data points on top of the box plot.
# This lets the reader see the actual distribution, not just summary boxes.
# jitter=True scatters the points horizontally a bit so they don't overlap.
sns.stripplot(
    data=df,
    x="Condition",
    y="Distance",
    order=condition_order,
    jitter=True,
    color="black",
    alpha=0.5,    # transparency so points don't hide each other
    size=4,       # small dots
    ax=ax,
)

# Add labels and a title
ax.set_title("Distance between " + CHR_A + " and " + CHR_B + " centroids")
ax.set_ylabel("Centroid distance (model units)")
ax.set_xlabel("Condition")

# tight_layout prevents labels from getting cut off
fig.tight_layout()

# Build the output path and save the figure
plot_path = os.path.join(OUTPUT_FOLDER, "chromosome_distance_boxplot.png")
fig.savefig(plot_path)

# Close the figure to free memory
plt.close(fig)

print("")
print("Plot saved: " + plot_path)


# =============================================================================
# Main code - Step 3: Run the statistical tests
# =============================================================================
# We compare each pair of conditions to see if their distance distributions
# are significantly different.
# We use Mann-Whitney U because the data may not be normally distributed.
# (A t-test would assume a normal distribution.)

def compare_two_conditions(df, cond1, cond2):
    """
    Compare the distance distributions for two conditions using
    a Mann-Whitney U test, and print the result.
    """
    # Get the distance values for each condition
    distances_1 = df[df["Condition"] == cond1]["Distance"]
    distances_2 = df[df["Condition"] == cond2]["Distance"]

    # If either group has no data, we can't run the test
    if len(distances_1) == 0 or len(distances_2) == 0:
        print("  " + cond1 + " vs " + cond2 + ": insufficient data")
        return

    # Run the test.
    # alternative="two-sided" means we test for ANY difference (either bigger
    # OR smaller), we don't pick a direction in advance.
    # The test returns:
    #   stat: the U statistic
    #   p:    the p-value (smaller = more significant difference)
    stat, p = mannwhitneyu(distances_1, distances_2, alternative="two-sided")

    # Calculate the median for each group so we can describe the direction
    median_1 = distances_1.median()
    median_2 = distances_2.median()

    # Figure out whether cond2 has a smaller or larger median than cond1
    # (we can use a simple if/else for clarity)
    if median_2 < median_1:
        direction = "closer"
    else:
        direction = "further"

    # Print the results.
    # The format codes:
    #   {:.0f}: whole number, no decimals
    #   {:.4g}: up to 4 significant digits (e.g. 0.0123 or 1.23e-05)
    #   {:.3f}: 3 digits after the decimal point
    print("  " + cond1 + " vs " + cond2
          + ": U={:.0f}, p={:.4g}".format(stat, p)
          + "  (median " + cond1 + "={:.3f},".format(median_1)
          + " median " + cond2 + "={:.3f}".format(median_2)
          + " — " + cond2 + " " + direction + " than " + cond1 + ")")


# Now actually run the three pairwise tests
print("")
print("Statistical tests (Mann-Whitney U, two-sided):")
compare_two_conditions(df, "WT", "T1")
compare_two_conditions(df, "WT", "C1")
compare_two_conditions(df, "T1", "C1")

print("")
print("Done!")
