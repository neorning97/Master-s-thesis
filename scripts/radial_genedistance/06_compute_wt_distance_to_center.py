"""
06_compute_wt_distance_to_center.py
=====================================
This script measures how far genes are from the centre of the nucleus
in WILDTYPE (WT) cells. WT means "normal": cells that have NOT had a
translocation.

Why do we need this?
- Script 05 measured distances in translocated cells (T1 and C1).
- But to know if the translocation actually MOVED the genes, we need to
  compare with where those same genes were before the translocation.
- That's what this script does: it measures the same genes in WT cells.
- Then script 07 will compare the two and tell us if the genes moved.

Important difference from script 05:
- In WT cells, there's no rearrangement, so we use the NORMAL gene
  coordinates from a reference genome file (a GTF file).
- In script 05, we used the translocated coordinates from the gene
  expression file.

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
# I'm putting all my settings as variables at the top so I can easily
# change them later without digging through the code.

# Gene expression files (made by script 01), we use these to know
# which genes we care about
T1_GENE_FILE = "/path/to/T1_genes_expression.tsv"
C1_GENE_FILE = "/path/to/C1_genes_expression.tsv"

# GTF file = a standard genome annotation file that tells us where each
# gene is located in a NORMAL (non-translocated) genome.
# WT cells aren't rearranged, so we use these standard coordinates.
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# Folder containing the WT .cmm 3D model files
WT_CMM_FOLDER = "/path/to/cmm/WT"

# Where to save our output
OUTPUT_FILE = "/path/to/WT_distance_to_center_agg.tsv"


# =============================================================================
# Helper function: Read a .cmm file and pull out the bead positions
# =============================================================================
# A .cmm file is a text file with information about "beads". Each bead
# represents a small piece of a chromosome at a specific 3D position.
# This is the same function as in script 05.

def read_cmm_file(file_path):
    """
    Open a .cmm file and read out all the bead information.
    Returns a pandas DataFrame (a table) with columns:
    chr, start, end, x, y, z
    """
    # Make an empty list to store the bead info as we find it
    bead_list = []

    # Open the file for reading
    # "with open" makes sure the file is closed automatically when done
    with open(file_path) as my_file:

        # Loop through the file one line at a time
        for line in my_file:

            # We only care about lines that contain "<marker"
            # If a line doesn't have this, skip to the next line
            if "<marker" not in line:
                continue

            # Each marker line has attributes like:
            #   beadID="chr1:1000-2000" x="1.5" y="2.0" z="3.0"
            # Use a regular expression to find all attribute="value" pairs.
            # The result is a list of pairs like
            #   [("beadID", "chr1:1000-2000"), ("x", "1.5"), ...]
            # Then we turn that list into a dictionary for easier access.
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            # If this marker doesn't have a beadID, we can't use it --> skip
            if "beadID" not in attrs:
                continue

            # Now try to extract the info we want.
            # Use try/except in case something is malformed in the file.
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

                # Get the x, y, z coordinates and convert to numbers (floats)
                # If they're missing for some reason, use "nan" (not a number)
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
# Main code
# =============================================================================

print("Computing WT distance to nuclear centre...")

# -----------------------------------------------------------------------------
# Step 1: Load the gene lists from T1 and C1
# -----------------------------------------------------------------------------
# We want to study the same genes that were in script 05.
# Read both gene files using pandas.

print("")
print("Step 1: Loading gene lists...")

t1_genes = pd.read_csv(T1_GENE_FILE, sep="\t")
c1_genes = pd.read_csv(C1_GENE_FILE, sep="\t")

# We only want "inside_transloc" and "control" genes
# (same as in script 05, we don't include "neighbor" genes)
keep_these = ["inside_transloc", "control"]

# Filter both tables
t1_genes = t1_genes[t1_genes["region_type"].isin(keep_these)]
c1_genes = c1_genes[c1_genes["region_type"].isin(keep_these)]

# We need to combine the genes from both T1 and C1 into one list,
# but we don't want duplicates (same gene appearing twice).
# Also, for each gene, we want to remember if it's "inside_transloc"
# or "control".
#
# A dictionary works well here: the key is the gene name,
# and the value is the region type.
# If the same gene is in both T1 and C1, the dictionary will just
# keep one copy automatically.

gene_region = {}  # empty dictionary to start

# Go through T1 genes first
for index, row in t1_genes.iterrows():
    gene_name = row["gene_name"]
    region_type = row["region_type"]
    gene_region[gene_name] = region_type

# Then go through C1 genes
for index, row in c1_genes.iterrows():
    gene_name = row["gene_name"]
    region_type = row["region_type"]
    gene_region[gene_name] = region_type

# Get the list of unique gene names (the dictionary keys)
genes_of_interest = list(gene_region.keys())

print("Total genes for WT mapping: " + str(len(genes_of_interest)))

# -----------------------------------------------------------------------------
# Step 2: Load the GTF file to get the gene coordinates
# -----------------------------------------------------------------------------
# A GTF file lists where every gene is located in the genome.
# It's a tab-separated file with these columns (in this order):
#   1. chrom (which chromosome)
#   2. source (where the info came from, like "ensembl")
#   3. feature (what type of thing this row is, like "gene" or "exon")
#   4. start (where it starts on the chromosome)
#   5. end (where it ends on the chromosome)
#   6. score
#   7. strand (+ or -, which direction the gene reads)
#   8. frame
#   9. attribute (extra info like the gene name, in a weird format)

print("")
print("Step 2: Loading GTF file...")

# Names for the 9 columns of a GTF file
gtf_column_names = ["chrom", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attribute"]

# Read the GTF file. Some special arguments:
#   sep="\t": file is tab-separated
#   comment="#": lines starting with # are comments and should be skipped
#   header=None: the file doesn't have a header row of column names
#   names=...: so we provide our own column names
gtf_df = pd.read_csv(
    GTF_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=gtf_column_names
)

# We only want rows where the feature is "gene"
# (a GTF can also have rows for exons, transcripts, etc., but we don't need them)
gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()

# The "attribute" column has text that looks like:
#   gene_id "ENSG00000123"; gene_name "ABC1"; gene_biotype "protein_coding"
# We need to pull out just the gene name.
# We use a regex to find the part inside the quotes after gene_name
# str.extract grabs the part inside the parentheses ( ).
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')

# Now keep only the rows for genes we actually care about
# (this makes the table much smaller and faster to work with)
gtf_df = gtf_df[gtf_df["gene_name"].isin(genes_of_interest)].copy()

# Make sure all chromosome names start with "chr"
# Same trick as in script 05, using a list comprehension
gtf_df["chrom"] = gtf_df["chrom"].astype(str)

fixed_chroms = []
for c in gtf_df["chrom"]:
    if c.startswith("chr"):
        fixed_chroms.append(c)
    else:
        fixed_chroms.append("chr" + c)
gtf_df["chrom"] = fixed_chroms

# Make sure start and end are numbers (not text).
# pd.to_numeric converts text to numbers.
# errors="coerce" means: if something can't be converted, make it NaN
# instead of crashing.
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")

# Reset the row numbers since we did filtering above
# (this avoids any index alignment issues later)
gtf_df = gtf_df.reset_index(drop=True)

print("GTF entries found for these genes: " + str(len(gtf_df)))

# -----------------------------------------------------------------------------
# Step 3: Find all the WT .cmm model files
# -----------------------------------------------------------------------------
# glob.glob finds all files matching a pattern. We want all .cmm files.
# os.path.join builds the path correctly for any operating system.
# sorted() puts them in alphabetical order.

print("")
print("Step 3: Finding WT .cmm files...")

search_pattern = os.path.join(WT_CMM_FOLDER, "*.cmm")
cmm_files = sorted(glob.glob(search_pattern))

print("Found " + str(len(cmm_files)) + " WT CMM files to process")

# -----------------------------------------------------------------------------
# Step 4: Process each WT model file
# -----------------------------------------------------------------------------
# For each gene, we'll collect all the distance values it gets across
# all the models. We use a dictionary where:
#   key = gene name
#   value = list of distance values (one per model)
#
# When we encounter a gene for the first time, we need to start with
# an empty list. We could do that with a check each time, OR we can
# use defaultdict from collections, but to keep things simple I'll
# just check if the gene is already in the dictionary.

print("")
print("Step 4: Processing each WT model file...")

# Empty dictionary, will fill it up as we go
distances_per_gene = {}

# Loop through each .cmm file
for f in cmm_files:

    # Get just the filename (not the full path) for nicer printing
    file_name = os.path.basename(f)
    print("  Processing model: " + file_name)

    # Read the .cmm file using our helper function
    beads = read_cmm_file(f)

    # If the file had no beads, skip it
    if beads.empty:
        print("    No markers found, skipping")
        continue

    # -------------------------------------------------------------------------
    # For each gene, find the closest bead and calculate its distance
    # -------------------------------------------------------------------------
    # Loop through each gene in the GTF table
    for index, gene_row in gtf_df.iterrows():

        # Get info about this gene
        gene_name = gene_row["gene_name"]
        gene_chrom = gene_row["chrom"]
        gene_start = gene_row["start"]
        gene_end = gene_row["end"]

        # Calculate the gene's midpoint
        gene_midpoint = (gene_start + gene_end) / 2

        # Get only the beads on the same chromosome as this gene
        beads_on_same_chr = beads[beads["chr"] == gene_chrom]

        # If there are no beads on this chromosome, we can't match the gene
        # Skip it and move on
        if len(beads_on_same_chr) == 0:
            continue

        # Calculate the midpoint of each bead on this chromosome
        bead_midpoints = (beads_on_same_chr["start"] + beads_on_same_chr["end"]) / 2

        # For each bead, find how far its midpoint is from the gene's midpoint
        # abs() gives us the absolute distance (no negative numbers)
        midpoint_diffs = abs(bead_midpoints - gene_midpoint)

        # Find which bead has the smallest distance (the closest one)
        # idxmin() gives the row label of the minimum value
        closest_bead_index = midpoint_diffs.idxmin()

        # Look up that bead's x, y, z coordinates
        closest_bead = beads_on_same_chr.loc[closest_bead_index]
        x = closest_bead["x"]
        y = closest_bead["y"]
        z = closest_bead["z"]

        # ---------------------------------------------------------------------
        # Calculate the distance from this bead to the centre (0, 0, 0)
        # ---------------------------------------------------------------------
        # The Euclidean distance from (x, y, z) to (0, 0, 0) is:
        #     sqrt(x^2 + y^2 + z^2)
        # We can use np.linalg.norm to do this for us.
        # We pass the coordinates as a small array [x, y, z].
        coords = np.array([x, y, z])
        dist = np.linalg.norm(coords)

        # ---------------------------------------------------------------------
        # Add this distance to our dictionary for this gene
        # ---------------------------------------------------------------------
        # If we've never seen this gene before, start with an empty list
        if gene_name not in distances_per_gene:
            distances_per_gene[gene_name] = []

        # Add this distance to the gene's list
        distances_per_gene[gene_name].append(dist)

# -----------------------------------------------------------------------------
# Step 5: Calculate the average distance for each gene across all models
# -----------------------------------------------------------------------------
# Each gene now has a list of distances (one per model).
# We want to make a final table with:
#   - gene_name
#   - region_type ("inside_transloc" or "control")
#   - mean_WT (average distance across all models)
#   - std_WT (standard deviation: how spread out the values are)
#   - n_models (how many models we got data from)

print("")
print("Step 5: Calculating averages...")

# We'll build up a list of dictionaries, one per gene.
# Then turn that list into a DataFrame at the end.
results_list = []

# Loop through each gene in our distances dictionary
# .items() gives us both the key (gene_name) and value (list of distances)
for gene_name, distance_list in distances_per_gene.items():

    # Calculate stats using numpy
    mean_dist = np.mean(distance_list)   # average
    std_dist = np.std(distance_list)     # standard deviation
    n_models = len(distance_list)        # how many models contributed

    # Look up the region type for this gene
    # .get() is like ["..."] but lets us provide a default if the key is missing
    region = gene_region.get(gene_name, "unknown")

    # Make a small dictionary with this gene's results
    one_result = {
        "gene_name": gene_name,
        "region_type": region,
        "mean_WT": mean_dist,
        "std_WT": std_dist,
        "n_models": n_models,
    }

    # Add it to our results list
    results_list.append(one_result)

# Turn the list of dictionaries into a DataFrame
# Each dictionary becomes one row, and the keys become column names
wt_results_df = pd.DataFrame(results_list)

# -----------------------------------------------------------------------------
# Step 6: Save the results to a TSV file
# -----------------------------------------------------------------------------

# Make sure the output folder exists.
# exist_ok=True means "don't crash if it already exists"
output_folder = os.path.dirname(OUTPUT_FILE)
os.makedirs(output_folder, exist_ok=True)

# Save the table.
# sep="\t": tab separated
# index=False: don't write row numbers as a column
wt_results_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("")
print("Aggregated to " + str(len(wt_results_df)) + " unique genes")
print("Saved WT distances to: " + OUTPUT_FILE)
print("")
print("Done!")
