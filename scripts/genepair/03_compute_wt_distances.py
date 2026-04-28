"""
03_compute_wt_distances.py
===========================
This script calculates 3D distances between gene pairs in the WILDTYPE (WT)
cells. WT means "normal", cells that have NOT had a translocation.

Why we need this:
- Script 02 calculated distances in the translocated cells (T1, C1).
- This script does the same for normal (WT) cells.
- Comparing the two tells us: did the translocation MOVE genes closer
  together or push them apart?

Important:
- We use the EXACT SAME gene pairs as script 02. That way, the comparison
  is fair: we're measuring the same pairs in both cases.
- We use the standard reference gene coordinates from a GTF file (not
  rearranged), because WT cells aren't rearranged.

Run script 02 first to make the T1 and C1 distance files.
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
# CONFIG SECTION - Edit these paths before running
# =============================================================================

# -----------------------------------------------------------------------------
# Distance files from script 02 (used to get the gene pairs we want)
# -----------------------------------------------------------------------------
T1_DISTANCES_FILE = "/path/to/T1_distances.tsv"
C1_DISTANCES_FILE = "/path/to/C1_distances.tsv"

# -----------------------------------------------------------------------------
# GTF file with reference (non-rearranged) gene coordinates
# -----------------------------------------------------------------------------
GTF_FILE = "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# -----------------------------------------------------------------------------
# Folder with WT 3D model files
# -----------------------------------------------------------------------------
WT_CMM_FOLDER = "/path/to/cmm/WT"

# -----------------------------------------------------------------------------
# Where to save the output
# -----------------------------------------------------------------------------
OUTPUT_FILE = "/path/to/results/WT_distances_agg.tsv"


# =============================================================================
# Helper function: Add "chr" prefix to a chromosome name if it's missing
# =============================================================================

def fix_chrom_name(c):
    """If c starts with 'chr', return as-is. Otherwise prepend 'chr'."""
    c = str(c).strip()
    if c.startswith("chr"):
        return c
    else:
        return "chr" + c


# =============================================================================
# Helper function: Read a .cmm file
# =============================================================================
# Same as in script 02

def read_cmm_file(file_path):
    """
    Open a .cmm file and read the bead positions.
    Returns a DataFrame with columns: chr, start, end, x, y, z
    """
    bead_list = []

    with open(file_path, "r") as my_file:
        for line in my_file:

            # Skip lines without bead markers
            if "<marker" not in line:
                continue

            # Parse all attribute="value" pairs into a dict
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            if "beadID" not in attrs:
                continue

            try:
                # The beadID format is "chr1:1000-2000"
                bead_id = attrs["beadID"]
                chrom_part, coords_part = bead_id.split(":")

                # Make sure chromosome name starts with "chr"
                chrom_part = fix_chrom_name(chrom_part)

                # Get start and end as integers
                start_str, end_str = coords_part.split("-")
                start = int(start_str)
                end = int(end_str)

                # Get x, y, z coordinates as floats
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))

                bead_list.append([chrom_part, start, end, x, y, z])

            except Exception:
                # Skip badly-formatted lines
                continue

    beads_df = pd.DataFrame(
        bead_list,
        columns=["chr", "start", "end", "x", "y", "z"]
    )
    return beads_df


# =============================================================================
# Main code starts here
# =============================================================================

print("Computing WT distances...")


# -----------------------------------------------------------------------------
# Step 1: Load the distance files from script 02
# -----------------------------------------------------------------------------
# These files give us the gene PAIRS we need to recompute in WT.

df_t1 = pd.read_csv(T1_DISTANCES_FILE, sep="\t")
df_c1 = pd.read_csv(C1_DISTANCES_FILE, sep="\t")


# -----------------------------------------------------------------------------
# Step 2: Build a lookup of gene --> region type
# -----------------------------------------------------------------------------
# We need to remember whether each gene is "inside_transloc", "neighbor",
# or "control", to put back in the final output.
#
# A dictionary works well here:
# - Keys = gene names
# - Values = region type ("inside_transloc", "neighbor", or "control")

gene_region_type = {}

# Loop through both DataFrames (T1 and C1)
for df in [df_t1, df_c1]:

    # Process gene1 column.
    # zip pairs up the values from two columns row by row.
    for gene, region in zip(df["gene1"], df["gene1_region_type"]):
        gene_region_type[gene] = region

    # Process gene2 column
    for gene, region in zip(df["gene2"], df["gene2_region_type"]):
        gene_region_type[gene] = region


# Get the list of all unique gene names
all_genes = list(gene_region_type.keys())
print("Total unique genes across T1 + C1: " + str(len(all_genes)))


# -----------------------------------------------------------------------------
# Step 3: Build the set of unique gene PAIRS we need to measure
# -----------------------------------------------------------------------------
# We use a Python set to automatically deduplicate (the same pair might
# appear in both T1 and C1 distance files).

gene_pairs = set()

for df in [df_t1, df_c1]:
    # zip pairs up gene1 and gene2 row by row.
    # set.update adds each pair to the set.
    for g1, g2 in zip(df["gene1"], df["gene2"]):
        gene_pairs.add((g1, g2))

print("Total gene pairs to compute: " + str(len(gene_pairs)))


# -----------------------------------------------------------------------------
# Step 4: Load the GTF file to get gene coordinates
# -----------------------------------------------------------------------------
# WT cells aren't rearranged, so we use the standard reference gene
# coordinates from the GTF.

print("Loading GTF file...")

# GTF files have 9 columns, we provide names since GTF has no header
gtf_columns = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
]

# Read the GTF, skipping comment lines (starting with #)
gtf_df = pd.read_csv(
    GTF_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=gtf_columns
)

# Keep only "gene" rows (skip transcripts, exons, etc.)
gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()

# Pull out the gene_name from the attribute column.
# The attribute column looks like: gene_id "ENSG..."; gene_name "ABC1"; ...
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')

# Keep only the genes we actually need (saves memory)
gtf_df = gtf_df[gtf_df["gene_name"].isin(all_genes)].copy()

# Make sure chromosome names start with "chr"
# Using a list comprehension (one-line for loop)
gtf_df["chrom"] = [fix_chrom_name(c) for c in gtf_df["chrom"]]

# Make sure start/end are numbers
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")

# Reset row numbers (gone gappy after filtering)
gtf_df = gtf_df.reset_index(drop=True)


# -----------------------------------------------------------------------------
# Step 5: Find all the WT model files
# -----------------------------------------------------------------------------
search_pattern = os.path.join(WT_CMM_FOLDER, "*.cmm")
cmm_files = sorted(glob.glob(search_pattern))
print("Found " + str(len(cmm_files)) + " CMM files to process")


# -----------------------------------------------------------------------------
# Step 6: Process each model and collect distances
# -----------------------------------------------------------------------------
# For each gene pair, we want to get a list of distance values
# (one per model). We store these in a dictionary:
#   key = (gene1, gene2) tuple
#   value = list of distances across models

distances_dict = {}

# Also keep track of which chromosome each gene was mapped to.
# This is the same across all models (since the genes are at fixed
# positions), so we only need to remember it once.
chr_dict_cache = {}

# Loop through each WT model
for f in cmm_files:

    file_name = os.path.basename(f)
    print("")
    print("Processing model: " + file_name)

    # ---------------------------------------------------------------------
    # Step 6a: Read the .cmm file
    # ---------------------------------------------------------------------
    beads = read_cmm_file(f)

    if beads.empty:
        print("  No markers found, skipping")
        continue

    # ---------------------------------------------------------------------
    # Step 6b: Build a fast bead lookup
    # ---------------------------------------------------------------------
    # For each chromosome, we precompute the bead midpoints and 3D
    # coordinates so we can look them up quickly.
    bead_index = {}

    for chrom, bdf in beads.groupby("chr"):

        # Calculate the midpoint of each bead on this chromosome
        midpoints = ((bdf["start"] + bdf["end"]) / 2).to_numpy()

        # Get the 3D coordinates as a numpy array
        coordinates = bdf[["x", "y", "z"]].to_numpy()

        # Store in the dict
        bead_index[chrom] = {
            "mids": midpoints,
            "coords": coordinates,
        }

    # We don't need the raw beads DataFrame anymore
    del beads

    # ---------------------------------------------------------------------
    # Step 6c: Find the 3D coordinates for every gene in this model
    # ---------------------------------------------------------------------
    # We build two dictionaries:
    #   coords_dict: gene name --> [x, y, z]
    #   chr_dict:    gene name --> chromosome
    coords_dict = {}
    chr_dict = {}

    # Loop through each gene in the GTF table
    for _, row in gtf_df.iterrows():

        gene = row["gene_name"]
        chrom = row["chrom"]

        # Calculate the gene's midpoint
        mid = (row["start"] + row["end"]) / 2

        # Look up the beads for this chromosome.
        # .get() returns None if the chromosome isn't in the dict.
        entry = bead_index.get(chrom)
        if entry is None:
            continue

        # Find the bead with the closest midpoint to this gene's midpoint
        differences = np.abs(entry["mids"] - mid)
        nearest_idx = differences.argmin()

        # Save this gene's 3D coordinates and chromosome
        coords_dict[gene] = entry["coords"][nearest_idx]
        chr_dict[gene] = chrom

    # Save the chromosome mapping the first time around
    # (it's the same across all models, so we only need to do this once).
    # The dict will be empty {} until the first model finishes.
    if len(chr_dict_cache) == 0:
        chr_dict_cache = chr_dict

    # ---------------------------------------------------------------------
    # Step 6d: Calculate distances for all gene pairs in this model
    # ---------------------------------------------------------------------
    # Loop through each gene pair we need to measure
    for gene1, gene2 in gene_pairs:

        # Only proceed if BOTH genes were successfully mapped
        if gene1 not in coords_dict:
            continue
        if gene2 not in coords_dict:
            continue

        # Calculate the Euclidean distance between the two 3D positions
        # np.linalg.norm of a difference vector gives its length:
        # sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
        diff = coords_dict[gene1] - coords_dict[gene2]
        dist = np.linalg.norm(diff)

        # Store the distance.
        # If we've never seen this pair before, start with an empty list.
        pair_key = (gene1, gene2)
        if pair_key not in distances_dict:
            distances_dict[pair_key] = []

        # Add this distance to the list
        distances_dict[pair_key].append(dist)

    # Free memory before the next model
    del bead_index


# -----------------------------------------------------------------------------
# Step 7: Aggregate distances across all models
# -----------------------------------------------------------------------------
# For each gene pair, calculate:
# - mean distance across all models
# - standard deviation (how much it varies)
# - count of models

print("")
print("Aggregating distances across models...")

# We'll build a list of dicts (one per gene pair), then turn it into a DataFrame
records = []

# Loop through each gene pair in our distances dictionary
# .items() gives both the key (the pair) and value (list of distances)
for pair, dist_list in distances_dict.items():

    # Unpack the tuple
    gene1, gene2 = pair

    # Calculate stats
    mean_dist = np.mean(dist_list)
    std_dist = np.std(dist_list)
    n_models = len(dist_list)

    # Look up region types and chromosomes.
    # .get(key, default) returns the default if key is missing.
    g1_region = gene_region_type.get(gene1, "unknown")
    g2_region = gene_region_type.get(gene2, "unknown")
    g1_chrom = chr_dict_cache.get(gene1, "")
    g2_chrom = chr_dict_cache.get(gene2, "")

    # Build a small dict for this pair's results
    one_record = {
        "gene1": gene1,
        "gene1_region_type": g1_region,
        "chr1": g1_chrom,
        "gene2": gene2,
        "gene2_region_type": g2_region,
        "chr2": g2_chrom,
        "mean_distance": mean_dist,
        "std_distance": std_dist,
        "n_models": n_models,
    }

    records.append(one_record)

# Turn the list of dicts into a DataFrame
# Each dict becomes a row; the keys become column names
summary_df = pd.DataFrame(records)


# -----------------------------------------------------------------------------
# Step 8: Save the results
# -----------------------------------------------------------------------------
# Make sure the output folder exists
output_folder = os.path.dirname(OUTPUT_FILE)
os.makedirs(output_folder, exist_ok=True)

# Save the table as TSV
summary_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("")
print("Aggregated to " + str(len(summary_df)) + " unique gene pairs")
print("Saved to: " + OUTPUT_FILE)
print("Done!")
