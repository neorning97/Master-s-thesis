"""
05_compute_tc_distance_to_center.py
=====================================
This script measures how far translocated genes are from the centre of the
nucleus, using 3D models of chromosomes (Chrom3D .cmm files).

The idea:
- The nucleus is modeled as a sphere with its centre at point (0, 0, 0).
- Each gene has a position in 3D space (x, y, z coordinates).
- The distance from (0, 0, 0) to the gene's position tells us how close
  the gene is to the nuclear centre.
- Genes near the centre are usually "active" (turned on).
- Genes near the edge are usually "inactive" (turned off).

We do this for each .cmm model file we have, then average the results.

Run script 01 first to make the gene expression files.
Then edit the file paths below and run this script.

Required libraries: pandas, numpy
Install them with: pip install pandas numpy
"""

import os       
import re       
import glob     
import numpy as np   
import pandas as pd  

# =============================================================================
# CONFIG SECTION - Edit these file paths before running the script
# =============================================================================
# I'm putting all my settings in one dictionary at the top so I can easily
# change them later without digging through the code.

# T1 settings - paths for the T1 condition
T1_GENE_FILE = "/path/to/T1_genes_expression.tsv"   # input from script 01
T1_CMM_FOLDER = "/path/to/cmm/T1"                   # folder with .cmm files
T1_OUTPUT = "/path/to/T1_distance_to_center.tsv"    # main output file
T1_OUTPUT_AVG = "/path/to/T1_distance_to_center_agg.tsv"  # averaged output

# C1 settings - paths for the C1 condition
C1_GENE_FILE = "/path/to/C1_genes_expression.tsv"
C1_CMM_FOLDER = "/path/to/cmm/C1"
C1_OUTPUT = "/path/to/C1_distance_to_center.tsv"
C1_OUTPUT_AVG = "/path/to/C1_distance_to_center_agg.tsv"

# Now I'll put all of this into a list so I can loop through both conditions
# Each item in the list is a dictionary with the settings for one condition
conditions_list = [
    {
        "name": "T1",
        "gene_file": T1_GENE_FILE,
        "cmm_folder": T1_CMM_FOLDER,
        "output": T1_OUTPUT,
        "output_avg": T1_OUTPUT_AVG,
    },
    {
        "name": "C1",
        "gene_file": C1_GENE_FILE,
        "cmm_folder": C1_CMM_FOLDER,
        "output": C1_OUTPUT,
        "output_avg": C1_OUTPUT_AVG,
    },
]


# =============================================================================
# Helper function: Read a .cmm file and pull out the bead positions
# =============================================================================
# A .cmm file is a text file with information about "beads". Each bead
# represents a small piece of a chromosome at a specific 3D position.
# We need to extract: which chromosome, where on the chromosome (start/end),
# and the 3D coordinates (x, y, z).

def read_cmm_file(file_path):
    """
    Open a .cmm file and read out all the bead information.
    Returns a pandas DataFrame (a table) with columns:
    chr, start, end, x, y, z
    """
    # Make an empty list to store the bead info as we find it
    bead_list = []

    # Open the file for reading
    # Using "with open" makes sure the file is closed automatically when done
    with open(file_path) as my_file:

        # Loop through the file one line at a time
        for line in my_file:

            # We only care about lines that contain "<marker"
            # If a line doesn't have this, skip to the next line
            if "<marker" not in line:
                continue

            # Each marker line has attributes like beadID="chr1:1000-2000" x="1.5" y="2.0" z="3.0"
            # Use a regular expression to find all attribute="value" pairs
            # The result is a list of pairs like [("beadID", "chr1:1000-2000"), ("x", "1.5"), ...]
            # Then we turn that list into a dictionary for easier access
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            # If this marker doesn't have a beadID, we can't use it --> skip it
            if "beadID" not in attrs:
                continue

            # Now try to extract the info we want
            # Use try/except in case something is malformed in the file
            try:
                # The beadID looks like "chr1:1000-2000"
                # Split on the colon to get "chr1" and "1000-2000"
                bead_id = attrs["beadID"]
                chrom_part, coords_part = bead_id.split(":")

                # Make sure the chromosome name starts with "chr"
                # (sometimes it might just be "1" instead of "chr1")
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part

                # Split "1000-2000" into start=1000 and end=2000
                # Then convert them from text to integers
                start_str, end_str = coords_part.split("-")
                start = int(start_str)
                end = int(end_str)

                # Get the x, y, z coordinates and convert them to numbers (floats)
                # If they're missing for some reason, use "nan"
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))

                # Add this bead's info to our list
                bead_list.append([chrom_part, start, end, x, y, z])

            except Exception:
                # If anything went wrong with this line, just skip it
                # and keep going with the rest of the file
                continue

    # Now turn our list of beads into a pandas DataFrame (a table)
    # Each inner list becomes a row, and we name the columns
    beads_df = pd.DataFrame(
        bead_list,
        columns=["chr", "start", "end", "x", "y", "z"]
    )

    return beads_df


# =============================================================================
# Helper function: Find the closest bead to each gene
# =============================================================================
# Each gene has a position on a chromosome (a start and end).
# Each bead also has a position on a chromosome.
# For each gene, we want to find the bead that is closest to it
# (on the same chromosome) and get that bead's 3D coordinates.

def find_bead_coords_for_genes(genes_df, beads_df):
    """
    For each gene in genes_df, find the closest bead on the same chromosome
    in beads_df, and return the (x, y, z) coordinates of that bead.

    Returns a numpy array with shape (number_of_genes, 3).
    Each row is the [x, y, z] of the bead closest to that gene.
    """
    n_genes = len(genes_df)

    # Make an empty array to hold the results
    # It has n_genes rows and 3 columns (for x, y, z)
    # We fill it with "nan" so any genes we can't match will stay as nan
    coords_array = np.full((n_genes, 3), np.nan)

    # Loop through each gene one at a time
    # enumerate gives us both the index (i) and the row data
    # iterrows() gives us each row of the DataFrame
    for i, (_, gene_row) in enumerate(genes_df.iterrows()):

        # Which chromosome is this gene on?
        gene_chrom = gene_row["chrom"]

        # Get only the beads that are on the same chromosome as this gene
        # This is like filtering: keep only rows where chr matches
        beads_on_same_chr = beads_df[beads_df["chr"] == gene_chrom]

        # If there are no beads on this chromosome, we can't match the gene
        # Just skip it (the row will stay as nan)
        if len(beads_on_same_chr) == 0:
            continue

        # Calculate the midpoint of the gene (average of start and end)
        gene_midpoint = (gene_row["start"] + gene_row["end"]) / 2

        # Calculate the midpoint of each bead
        bead_midpoints = (beads_on_same_chr["start"] + beads_on_same_chr["end"]) / 2

        # For each bead, find how far its midpoint is from the gene's midpoint
        # We use abs() to get the absolute distance (no negative numbers)
        distances = abs(bead_midpoints - gene_midpoint)

        # Find which bead has the smallest distance (the closest one)
        # idxmin() gives us the index/row label of the minimum value
        closest_bead_index = distances.idxmin()

        # Look up that bead's x, y, z coordinates
        closest_bead = beads_on_same_chr.loc[closest_bead_index]
        x = closest_bead["x"]
        y = closest_bead["y"]
        z = closest_bead["z"]

        # Store these coordinates in our results array
        coords_array[i] = [x, y, z]

    return coords_array


# =============================================================================
# Main code
# =============================================================================

# Loop through each condition (T1 and C1) and process them one by one
for condition in conditions_list:

    # Get the name and paths for this condition
    cond_name = condition["name"]
    gene_file = condition["gene_file"]
    cmm_folder = condition["cmm_folder"]
    output_file = condition["output"]
    output_avg_file = condition["output_avg"]

    print("")
    print("===== Processing " + cond_name + " =====")

    # -------------------------------------------------------------------------
    # Step 1: Load the gene expression file (made by script 01)
    # -------------------------------------------------------------------------
    # This file is tab-separated
    genes = pd.read_csv(gene_file, sep="\t")

    # Make sure all chromosome names start with "chr"
    # We do this by going through each value in the chrom column
    # and adding "chr" if it's missing
    # First, make sure the column is text (string)
    genes["chrom"] = genes["chrom"].astype(str)

    # Now go through each value and fix it if needed
    # We use a list comprehension here, it's like a one-line for loop
    fixed_chroms = []
    for c in genes["chrom"]:
        if c.startswith("chr"):
            fixed_chroms.append(c)
        else:
            fixed_chroms.append("chr" + c)
    genes["chrom"] = fixed_chroms

    # -------------------------------------------------------------------------
    # Step 2: Filter to keep only the genes we care about
    # -------------------------------------------------------------------------
    # We only want translocated genes ("inside_transloc") and controls ("control").
    # We don't want the "neighbor" genes because we're asking about
    # the translocated genes themselves moving in the nucleus.
    # The .isin() function checks if each value is in our list of allowed values.
    keep_these = ["inside_transloc", "control"]
    genes = genes[genes["region_type"].isin(keep_these)]

    # .copy() makes a fresh copy so pandas doesn't give us warnings later
    # .reset_index(drop=True) renumbers the rows from 0, 1, 2, ... again.
    # When we filtered above, pandas kept the OLD row numbers (so they have
    # gaps now, like 2, 5, 7, 12...). If we don't reset them, when we later
    # try to add the "distances" numpy array as a new column, pandas can
    # get confused about which gene each distance belongs to and give us
    # NaN values or wrong matches.
    genes = genes.copy().reset_index(drop=True)

    print("Loaded " + str(len(genes)) + " genes (inside_transloc + control)")

    # -------------------------------------------------------------------------
    # Step 3: Find all the .cmm model files in the folder
    # -------------------------------------------------------------------------
    # glob.glob finds files matching a pattern, we want all .cmm files
    # os.path.join builds the path correctly (handles / vs. \ for different OSes)
    # sorted() puts them in alphabetical order
    search_pattern = os.path.join(cmm_folder, "*.cmm")
    cmm_files = sorted(glob.glob(search_pattern))

    print("Found " + str(len(cmm_files)) + " CMM files to process")

    # -------------------------------------------------------------------------
    # Step 4: Process each .cmm file (each is a different 3D model)
    # -------------------------------------------------------------------------
    # We'll keep all the results from each file in this list,
    # then combine them at the end
    all_results_list = []

    # Loop through each .cmm file
    for f in cmm_files:

        # os.path.basename gets just the filename (not the full path)
        # for cleaner printing
        file_name = os.path.basename(f)
        print("  Processing model: " + file_name)

        # Read the .cmm file using our helper function
        beads = read_cmm_file(f)

        # If the file had no beads, skip it
        if beads.empty:
            print("    No markers found, skipping")
            continue

        # For each gene, find the closest bead and get its 3D coordinates
        coords = find_bead_coords_for_genes(genes, beads)

        # ---------------------------------------------------------------------
        # Calculate the distance from each bead to the centre (0, 0, 0)
        # ---------------------------------------------------------------------
        # The Euclidean distance from (x, y, z) to (0, 0, 0) is:
        #     sqrt(x^2 + y^2 + z^2)
        # numpy has a function np.linalg.norm that does this for us.
        # axis=1 means "calculate one distance for each row".
        distances = np.linalg.norm(coords, axis=1)

        # ---------------------------------------------------------------------
        # Build a results table for this model
        # ---------------------------------------------------------------------
        # We want to keep some info from the genes table plus our new distance.
        # First, copy the columns we want to keep.
        # We also reset_index here so the row numbers are 0, 1, 2, ... which
        # matches the positions in the "distances" numpy array below.
        keep_cols = ["gene_id", "gene_name", "chrom", "region_type", "transloc_id"]
        result_for_this_model = genes[keep_cols].copy().reset_index(drop=True)

        # Add the distance column
        result_for_this_model["distance_to_center"] = distances

        # Add a column that says which model this came from
        result_for_this_model["model"] = file_name

        # Add this table to our list of results
        all_results_list.append(result_for_this_model)

    # -------------------------------------------------------------------------
    # Step 5: Check if we got any results, and if not, skip to next condition
    # -------------------------------------------------------------------------
    if len(all_results_list) == 0:
        print("  No results — skipping.")
        continue

    # -------------------------------------------------------------------------
    # Step 6: Combine all the per-model results into one big table
    # -------------------------------------------------------------------------
    # pd.concat sticks all the DataFrames together, one on top of another
    # ignore_index=True gives us fresh row numbers (instead of duplicates)
    big_results_df = pd.concat(all_results_list, ignore_index=True)

    # -------------------------------------------------------------------------
    # Step 7: Save the per-model results to a TSV file
    # -------------------------------------------------------------------------
    # Make sure the output folder exists. exist_ok=True means
    # "don't crash if the folder already exists".
    output_folder = os.path.dirname(output_file)
    os.makedirs(output_folder, exist_ok=True)

    # Save as a tab-separated file (sep="\t")
    # index=False means don't write the row numbers as a column
    big_results_df.to_csv(output_file, sep="\t", index=False)
    print("  Saved per-model distances: " + output_file)

    # -------------------------------------------------------------------------
    # Step 8: Calculate average distance for each gene across all models
    # -------------------------------------------------------------------------
    # Each gene appears once per model, so we have many distance values per gene.
    # We want to average them, also get the standard deviation, and count
    # how many models we used.
    #
    # groupby groups rows together that have the same values in the listed columns.
    # Then we tell it which column to do calculations on (distance_to_center)
    # and what calculations to do (mean, std, count).
    group_cols = ["gene_name", "chrom", "region_type", "transloc_id"]

    grouped = big_results_df.groupby(group_cols, as_index=False)

    # Calculate stats on the distance_to_center column
    # The named arguments (mean_distance=...) become the new column names
    averaged_df = grouped["distance_to_center"].agg(
        mean_distance="mean",   # average distance across models
        std_distance="std",     # standard deviation (how spread out the values are)
        n_models="count"        # how many models we had data from
    )

    # Save the averaged table to a file
    averaged_df.to_csv(output_avg_file, sep="\t", index=False)
    print("  Saved aggregated distances: " + output_avg_file)

# When the loop is done, we've processed both conditions and we're finished
print("")
print("Done!")
