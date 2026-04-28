"""
02_compute_tc_distances.py
===========================
This script calculates 3D distances between pairs of genes in the translocated
cell lines (T1 and C1).

The idea:
- After a translocation, a gene from chromosome A might end up physically
  next to a gene on chromosome B in the 3D space of the nucleus.
- We want to MEASURE these distances, to see how close translocated genes
  end up to their new neighbors.
- We do this for many different 3D models (each .cmm file is one possible
  3D model of the nucleus), then average the distances.

What this script does for each condition (T1, C1):
1. Loads the gene table from script 01.
2. For each .cmm file (3D model):
   a. Reads all the bead positions.
   b. Finds which bead each gene sits in.
   c. Calculates pairwise distances between:
      - Translocated genes <--> Their new neighbor genes
      - Control genes      <--> The same neighbor genes (for comparison)
3. Combines all models and averages the distances.

Why control vs neighbor pairs too?
- Translocated-neighbor pairs tell us "how close did the moved genes end up?"
- Control-neighbor pairs tell us "how close would RANDOM genes be to those
  same neighbors?"
- Comparing the two tells us if the translocated genes are SPECIFICALLY
  close to their new neighbors, or if it's just normal nuclear closeness.

About the SCENARIO setting:
- "inter" = only count pairs where the genes are on DIFFERENT chromosomes.
- "intra" = only count pairs where the genes are on the SAME chromosome.

Run script 01 first to make the gene expression files.
Then edit the file paths below and run this script.

Required libraries: pandas, numpy, scipy
Install them with: pip install pandas numpy scipy
"""

import os                                      
import re                                      
import glob                                    
import numpy as np                             
import pandas as pd                            
from scipy.spatial.distance import cdist       


# =============================================================================
# CONFIG SECTION - Edit these settings before running
# =============================================================================

# -----------------------------------------------------------------------------
# Which type of pairs to keep
# -----------------------------------------------------------------------------
# "inter" = only INTERchromosomal (different chromosomes)
# "intra" = only INTRAchromosomal (same chromosome)
SCENARIO = "inter"

# -----------------------------------------------------------------------------
# Input files for T1
# -----------------------------------------------------------------------------
T1_GENE_FILE = "/path/to/T1_genes_expression.tsv"
T1_NEIGHBOR_BED = "/path/to/T1_neighbor.bed"
T1_CMM_FOLDER = "/path/to/cmm/T1"
T1_OUTPUT = "/path/to/results/T1_distances.tsv"
T1_OUTPUT_AGG = "/path/to/results/T1_distances_agg.tsv"

# -----------------------------------------------------------------------------
# Input files for C1
# -----------------------------------------------------------------------------
C1_GENE_FILE = "/path/to/C1_genes_expression.tsv"
C1_NEIGHBOR_BED = "/path/to/C1_neighbor.bed"
C1_CMM_FOLDER = "/path/to/cmm/C1"
C1_OUTPUT = "/path/to/results/C1_distances.tsv"
C1_OUTPUT_AGG = "/path/to/results/C1_distances_agg.tsv"

# Group everything in a list so we can loop through both conditions
conditions_list = [
    {
        "name": "T1",
        "gene_file": T1_GENE_FILE,
        "neighbor_bed": T1_NEIGHBOR_BED,
        "cmm_folder": T1_CMM_FOLDER,
        "output": T1_OUTPUT,
        "output_agg": T1_OUTPUT_AGG,
    },
    {
        "name": "C1",
        "gene_file": C1_GENE_FILE,
        "neighbor_bed": C1_NEIGHBOR_BED,
        "cmm_folder": C1_CMM_FOLDER,
        "output": C1_OUTPUT,
        "output_agg": C1_OUTPUT_AGG,
    },
]


# =============================================================================
# Helper function: Add "chr" prefix to a chromosome name if it's missing
# =============================================================================
# Some files use "chr1", others just use "1". This helper normalizes them.

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

def read_cmm_file(file_path):
    """
    Open a .cmm file and read out all the bead information.
    Returns a DataFrame with columns: chr, start, end, x, y, z
    """
    bead_list = []

    with open(file_path, "r") as my_file:
        for line in my_file:

            # Only process lines containing bead markers
            if "<marker" not in line:
                continue

            # Parse all attribute="value" pairs into a dict
            found_attrs = re.findall(r'(\w+)="([^"]+)"', line)
            attrs = dict(found_attrs)

            if "beadID" not in attrs:
                continue

            try:
                # The beadID is in the format "chr1:1000-2000"
                bead_id = attrs["beadID"]
                chrom_part, coords_part = bead_id.split(":")

                # Make sure chromosome name starts with "chr"
                chrom_part = fix_chrom_name(chrom_part)

                # Get start and end as integers
                start_str, end_str = coords_part.split("-")
                start = int(start_str)
                end = int(end_str)

                # Get x, y, z coordinates as floats (NaN if missing)
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
# Helper function: For each gene, find the nearest bead and its 3D coordinates
# =============================================================================

def map_genes_to_beads(genes_df, bead_index, name_col):
    """
    For each gene in genes_df, find the nearest bead on the same chromosome
    and return its 3D coordinates.

    bead_index is a dict: chromosome --> {"mids": ..., "coords": ..., "chr": ...}
    name_col is the column with gene names (either "gene_name" or "gene_id")

    Returns three things:
    - coords: a 2D numpy array of shape (n_genes, 3) with x, y, z values
    - names:  a 1D array of gene names
    - chr_mapped: a 1D array of which chromosome each gene was mapped to
    """
    # If there are no genes, return empty arrays
    if genes_df.shape[0] == 0:
        empty_coords = np.zeros((0, 3))
        empty_names = np.array([], dtype=str)
        empty_chrs = np.array([], dtype=str)
        return empty_coords, empty_names, empty_chrs

    # Get gene names as a numpy array
    names = genes_df[name_col].to_numpy(str)

    # Calculate the midpoint of each gene
    # (e.g. for a gene at start=100, end=300, the midpoint is 200)
    mids = ((genes_df["start"] + genes_df["end"]) / 2).to_numpy()

    # Make empty arrays to fill with our results
    n_genes = len(names)
    coords = np.full((n_genes, 3), np.nan)            # n_genes rows, 3 columns
    chr_mapped = np.full(n_genes, "", dtype=object)   # n_genes empty strings

    # Loop through each gene
    # zip pairs up the chromosome and midpoint values
    chroms = genes_df["chrom"].values
    for i in range(n_genes):
        chrom = chroms[i]
        mid = mids[i]

        # Look up the beads on this chromosome
        # .get() returns None if the chromosome isn't in the dict
        entry = bead_index.get(chrom)
        if entry is None:
            continue   # leave coords[i] as NaN

        # Find the bead with the closest midpoint to this gene's midpoint
        # entry["mids"] is an array of all bead midpoints on this chromosome
        # np.abs(...).argmin() returns the index of the smallest absolute difference
        differences = np.abs(entry["mids"] - mid)
        nearest_idx = differences.argmin()

        # Save the bead's coordinates and chromosome
        coords[i] = entry["coords"][nearest_idx]
        chr_mapped[i] = entry["chr"]

    return coords, names, chr_mapped


# =============================================================================
# Helper function: Build a fast lookup table for beads
# =============================================================================
# We split the beads by chromosome and pre-compute the midpoints, so we can
# do fast lookups in map_genes_to_beads() above.

def build_bead_index(beads_df):
    """
    Convert a beads DataFrame into a dict for fast lookup.
    Returns dict: chromosome --> {"mids": numpy array, "coords": numpy array, "chr": str}
    """
    bead_index = {}

    # Group rows by chromosome.
    # .groupby returns (group_name, group_dataframe) pairs.
    for chrom, bdf in beads_df.groupby("chr"):

        # Calculate the midpoint of each bead on this chromosome
        midpoints = ((bdf["start"] + bdf["end"]) / 2).to_numpy()

        # Get the 3D coordinates as a numpy array
        coordinates = bdf[["x", "y", "z"]].to_numpy()

        # Store everything in the dict
        bead_index[chrom] = {
            "mids": midpoints,
            "coords": coordinates,
            "chr": chrom,
        }

    return bead_index


# =============================================================================
# Main code starts here
# =============================================================================

# Loop through each condition
for cond_info in conditions_list:

    # Pull out the values for this condition
    cond_name = cond_info["name"]
    gene_file = cond_info["gene_file"]
    neighbor_bed_file = cond_info["neighbor_bed"]
    cmm_folder = cond_info["cmm_folder"]
    output_file = cond_info["output"]
    output_agg_file = cond_info["output_agg"]

    print("")
    print("===== Processing " + cond_name + " =====")

    # -------------------------------------------------------------------------
    # Step 1: Load the gene expression and neighbor BED files
    # -------------------------------------------------------------------------
    df_genes = pd.read_csv(gene_file, sep="\t")
    df_neighbor_bed = pd.read_csv(neighbor_bed_file, sep="\t")

    print("Loaded " + str(len(df_genes)) + " genes, "
          + str(len(df_neighbor_bed)) + " neighbor regions")

    # -------------------------------------------------------------------------
    # Step 2: Pick the right gene name column
    # -------------------------------------------------------------------------
    # Use gene_name if available, otherwise fall back to gene_id
    if "gene_name" in df_genes.columns:
        name_col = "gene_name"
    else:
        name_col = "gene_id"

    # -------------------------------------------------------------------------
    # Step 3: Clean up the gene table
    # -------------------------------------------------------------------------
    # Make sure chromosomes start with "chr".
    # We use a list comprehension here (a one-line for loop):
    df_genes["chrom"] = [fix_chrom_name(c) for c in df_genes["chrom"]]

    # Convert start/end to numbers (NaN if invalid)
    df_genes["start"] = pd.to_numeric(df_genes["start"], errors="coerce")
    df_genes["end"] = pd.to_numeric(df_genes["end"], errors="coerce")

    # Clean region_type column: strip whitespace and lowercase
    df_genes["region_type"] = df_genes["region_type"].str.strip().str.lower()

    # -------------------------------------------------------------------------
    # Step 4: Clean up the neighbor BED file too
    # -------------------------------------------------------------------------
    # Note: neighbor BED uses "neighbor_chr" (the receiving chromosome)
    df_neighbor_bed["neighbor_chr"] = [
        fix_chrom_name(c) for c in df_neighbor_bed["neighbor_chr"]
    ]
    df_neighbor_bed["chrom"] = df_neighbor_bed["chrom"].astype(str).str.strip()

    # -------------------------------------------------------------------------
    # Step 5: Find all 3D model files
    # -------------------------------------------------------------------------
    search_pattern = os.path.join(cmm_folder, "*.cmm")
    cmm_files = sorted(glob.glob(search_pattern))
    print("Found " + str(len(cmm_files)) + " CMM files to process")

    # -------------------------------------------------------------------------
    # Step 6: Set up the output file
    # -------------------------------------------------------------------------
    # We'll write each model's results to disk one at a time, instead of
    # keeping them all in memory. This saves memory but means we need
    # to manage the file ourselves.

    # Make sure the output folder exists
    output_folder = os.path.dirname(output_file)
    os.makedirs(output_folder, exist_ok=True)

    # We'll only write the column headers for the FIRST model
    # (otherwise we'd have headers in the middle of the file).
    wrote_header = False

    # Track if we got any results at all
    any_results = False

    # Open the output file once for writing.
    # "w" means write mode (overwrite if file exists).
    out_fh = open(output_file, "w")

    # -------------------------------------------------------------------------
    # Step 7: Loop through each 3D model
    # -------------------------------------------------------------------------
    for f in cmm_files:

        file_name = os.path.basename(f)
        print("")
        print("Processing model: " + file_name)

        # ---------------------------------------------------------------------
        # Step 7a: Read the .cmm file
        # ---------------------------------------------------------------------
        beads = read_cmm_file(f)
        if beads.empty:
            print("  No markers found, skipping")
            continue

        # ---------------------------------------------------------------------
        # Step 7b: Build the fast bead lookup
        # ---------------------------------------------------------------------
        bead_index = build_bead_index(beads)

        # We don't need the raw beads DataFrame anymore --> free the memory
        del beads

        # ---------------------------------------------------------------------
        # Step 7c: Split genes by region type
        # ---------------------------------------------------------------------
        trans_genes = df_genes[df_genes["region_type"] == "inside_transloc"].copy()
        neighbor_genes = df_genes[df_genes["region_type"] == "neighbor"].copy()
        control_genes = df_genes[df_genes["region_type"] == "control"].copy()

        # ---------------------------------------------------------------------
        # Step 7d: Map each gene to its nearest bead in this model
        # ---------------------------------------------------------------------
        trans_coords, trans_names, trans_chr = map_genes_to_beads(
            trans_genes, bead_index, name_col
        )
        neigh_coords, neigh_names, neigh_chr = map_genes_to_beads(
            neighbor_genes, bead_index, name_col
        )
        ctrl_coords, ctrl_names, ctrl_chr = map_genes_to_beads(
            control_genes, bead_index, name_col
        )

        print("  Mapped " + str(len(trans_names)) + " trans, "
              + str(len(neigh_names)) + " neighbor, "
              + str(len(ctrl_names)) + " control genes")

        # We'll collect all the gene-pair results for this model here
        model_results = []

        # ---------------------------------------------------------------------
        # Step 7e: Compute pairwise distances for each translocation
        # ---------------------------------------------------------------------
        # For each translocation, we need to:
        # 1. Find its neighbor regions
        # 2. For each neighbor region, find which neighbor genes are in it
        # 3. Calculate distances from translocated genes to those neighbors
        # 4. Calculate distances from control genes to those neighbors

        unique_tids = trans_genes["transloc_id"].unique()

        for tid in unique_tids:

            # Get this translocation's neighbor regions
            neighbors_tid = df_neighbor_bed[df_neighbor_bed["transloc_id"] == tid]

            # Loop through each neighbor region
            for _, nb_row in neighbors_tid.iterrows():

                # The neighbor region's coordinates
                n_chr = nb_row["neighbor_chr"]
                nb_start = nb_row["start"]
                nb_end = nb_row["end"]

                # Find which neighbor genes are in this region.
                # The mask is a True/False Series, one entry per neighbor gene.
                same_chrom = (neighbor_genes["chrom"] == n_chr)
                starts_before_end = (neighbor_genes["start"] <= nb_end)
                ends_after_start = (neighbor_genes["end"] >= nb_start)
                neigh_mask = same_chrom & starts_before_end & ends_after_start

                # Skip if no neighbor genes in this region
                # np.any returns True if at least one element is True
                if not np.any(neigh_mask):
                    continue

                # Get the matching neighbor genes
                # reset_index() gives them clean row numbers 0, 1, 2, ...
                neigh_sub = neighbor_genes[neigh_mask].reset_index()

                # -------------------------------------------------------------
                # Translocated <--> Neighbor pairs
                # -------------------------------------------------------------
                # Find the translocated genes for THIS translocation
                trans_mask = (trans_genes["transloc_id"] == tid)

                # We'll save the distance matrix for use in the control step below.
                # Default to None in case there are no translocated genes.
                D_trans = None

                if np.any(trans_mask):

                    # Get the matching translocated genes
                    trans_sub = trans_genes[trans_mask].reset_index()

                    # ---------------------------------------------------------
                    # Calculate ALL pairwise distances at once
                    # ---------------------------------------------------------
                    # cdist computes a matrix of distances between all pairs
                    # of points in two sets. Much faster than a Python for loop.
                    # The result is a 2D matrix:
                    #   D[i, j] = distance from trans gene i to neighbor gene j

                    # First, we need to extract just the rows of trans_coords
                    # and neigh_coords that are part of this translocation.
                    # We use boolean indexing with .to_numpy() to convert the mask.
                    trans_coords_sub = trans_coords[trans_mask.to_numpy()]
                    neigh_coords_sub = neigh_coords[neigh_mask.to_numpy()]

                    D_trans = cdist(trans_coords_sub, neigh_coords_sub,
                                    metric="euclidean")

                    # ---------------------------------------------------------
                    # Build a results table for these pairs
                    # ---------------------------------------------------------
                    # We need to make a long-format table with one row per pair.
                    # np.meshgrid creates "all combinations" of indices.
                    # If we have 3 trans genes and 2 neighbor genes,
                    # we get 6 pairs (3 x 2 combinations).
                    t_idx, n_idx = np.meshgrid(
                        trans_sub.index, neigh_sub.index, indexing="ij"
                    )

                    # .ravel() flattens the 2D grids into 1D arrays
                    t_idx_flat = t_idx.ravel()
                    n_idx_flat = n_idx.ravel()

                    # Build the results DataFrame for these trans-neighbor pairs
                    pairs_df = pd.DataFrame({
                        "translocation": tid,
                        "gene1": trans_sub.loc[t_idx_flat, name_col].values,
                        "chr1": trans_sub.loc[t_idx_flat, "chrom"].values,
                        "gene1_region_type": "inside_transloc",
                        "gene1_chr_status": "translocated",
                        "gene2": neigh_sub.loc[n_idx_flat, name_col].values,
                        "chr2": neigh_sub.loc[n_idx_flat, "chrom"].values,
                        "gene2_region_type": "neighbor",
                        "gene2_chr_status": "non_translocated",
                        "distance": D_trans.ravel(),
                    })

                    model_results.append(pairs_df)

                # -------------------------------------------------------------
                # Control <--> Neighbor pairs
                # -------------------------------------------------------------
                # Same idea, but for control genes
                ctrl_mask = (control_genes["transloc_id"] == tid)

                if not np.any(ctrl_mask):
                    continue

                ctrl_sub = control_genes[ctrl_mask].reset_index()

                # Calculate all distances
                ctrl_coords_sub = ctrl_coords[ctrl_mask.to_numpy()]
                neigh_coords_sub2 = neigh_coords[neigh_mask.to_numpy()]
                D_ctrl = cdist(ctrl_coords_sub, neigh_coords_sub2,
                               metric="euclidean")

                # -------------------------------------------------------------
                # Match the number of control pairs to translocated pairs
                # -------------------------------------------------------------
                # If we have many more control pairs than translocated pairs,
                # we randomly sample down so the comparison is fair.
                if D_trans is not None:
                    n_trans_pairs = D_trans.size
                else:
                    n_trans_pairs = 0

                n_ctrl_pairs = D_ctrl.size

                # Build the index grids for control pairs
                ctrl_idx, neigh_idx = np.meshgrid(
                    ctrl_sub.index, neigh_sub.index, indexing="ij"
                )
                ctrl_idx_flat = ctrl_idx.ravel()
                neigh_idx_flat = neigh_idx.ravel()
                D_ctrl_flat = D_ctrl.ravel()

                # If there are more control pairs than trans pairs, subsample
                if n_ctrl_pairs > n_trans_pairs > 0:

                    # Use a reproducible random number generator.
                    # The seed comes from the translocation ID so the same
                    # tid always gives the same sample.
                    seed = abs(hash(tid)) % (2**32)
                    rng = np.random.default_rng(seed)

                    # Pick n_trans_pairs random indices, no duplicates
                    sample_idx = rng.choice(
                        len(ctrl_idx_flat),
                        size=n_trans_pairs,
                        replace=False
                    )

                    # Subset using the sampled indices
                    ctrl_idx_sample = ctrl_idx_flat[sample_idx]
                    neigh_idx_sample = neigh_idx_flat[sample_idx]
                    D_ctrl_sample = D_ctrl_flat[sample_idx]

                else:
                    # No subsampling needed --> keep all pairs
                    ctrl_idx_sample = ctrl_idx_flat
                    neigh_idx_sample = neigh_idx_flat
                    D_ctrl_sample = D_ctrl_flat

                # -------------------------------------------------------------
                # Build the control-neighbor results DataFrame
                # -------------------------------------------------------------
                pairs_df = pd.DataFrame({
                    "translocation": tid,
                    "gene1": ctrl_sub.loc[ctrl_idx_sample, name_col].values,
                    "chr1": ctrl_sub.loc[ctrl_idx_sample, "chrom"].values,
                    "gene1_region_type": "control",
                    "gene1_chr_status": "non_translocated",
                    "gene2": neigh_sub.loc[neigh_idx_sample, name_col].values,
                    "chr2": neigh_sub.loc[neigh_idx_sample, "chrom"].values,
                    "gene2_region_type": "neighbor",
                    "gene2_chr_status": "non_translocated",
                    "distance": D_ctrl_sample,
                })

                model_results.append(pairs_df)

        # ---------------------------------------------------------------------
        # Step 7f: If no results for this model, skip on
        # ---------------------------------------------------------------------
        if len(model_results) == 0:
            del bead_index
            continue

        # ---------------------------------------------------------------------
        # Step 7g: Combine all the per-translocation results for this model
        # ---------------------------------------------------------------------
        df_model = pd.concat(model_results, ignore_index=True)

        # ---------------------------------------------------------------------
        # Step 7h: Filter by inter/intra chromosomal scenario
        # ---------------------------------------------------------------------
        # Separate the trans-neighbor pairs from control-neighbor pairs
        # (we treat them differently in the filtering).

        # A pair is "trans-neighbor" if either gene1 or gene2 is inside_transloc
        # and the other is neighbor.
        is_trans_g1 = (df_model["gene1_region_type"] == "inside_transloc") & \
                      (df_model["gene2_region_type"] == "neighbor")
        is_trans_g2 = (df_model["gene1_region_type"] == "neighbor") & \
                      (df_model["gene2_region_type"] == "inside_transloc")
        is_trans_neigh = is_trans_g1 | is_trans_g2

        # A pair is "control-neighbor" if either gene1 or gene2 is control
        # and the other is neighbor.
        is_ctrl_g1 = (df_model["gene1_region_type"] == "control") & \
                     (df_model["gene2_region_type"] == "neighbor")
        is_ctrl_g2 = (df_model["gene1_region_type"] == "neighbor") & \
                     (df_model["gene2_region_type"] == "control")
        is_ctrl_neigh = is_ctrl_g1 | is_ctrl_g2

        df_trans_m = df_model[is_trans_neigh]
        df_ctrl_m = df_model[is_ctrl_neigh]

        # Apply the inter/intra filter to translocated-neighbor pairs.
        # For interchromosomal: keep only pairs on different chromosomes.
        # For intrachromosomal: keep only pairs on the same chromosome.
        if SCENARIO == "inter":
            df_trans_m = df_trans_m[df_trans_m["chr1"] != df_trans_m["chr2"]]
        elif SCENARIO == "intra":
            df_trans_m = df_trans_m[df_trans_m["chr1"] == df_trans_m["chr2"]]

        # Control-neighbor pairs are ALWAYS interchromosomal
        # (control genes are picked from chromosomes not involved in the translocation,
        # so they're on different chromosomes from the neighbors by definition).
        df_ctrl_m = df_ctrl_m[df_ctrl_m["chr1"] != df_ctrl_m["chr2"]]

        # Combine the filtered translocated and control results
        df_filtered = pd.concat([df_trans_m, df_ctrl_m], ignore_index=True)

        # ---------------------------------------------------------------------
        # Step 7i: Remove non-canonical chromosomes (scaffolds, patches, etc.)
        # ---------------------------------------------------------------------
        # We only want chromosomes 1-22, X, and Y.
        # Build a set of valid chromosome names (without "chr" prefix).
        # range(1, 23) gives us 1, 2, ..., 22, then we add X and Y.
        valid_numbers = set()
        for i in range(1, 23):
            valid_numbers.add(str(i))
        valid_chroms = valid_numbers | {"X", "Y"}    # | is set union

        # Apply the filter to both chromosome columns
        for col in ["chr1", "chr2"]:

            # Remove the "chr" prefix to compare with valid_chroms
            clean = df_filtered[col].str.replace("chr", "", case=False)

            # Keep only rows where the cleaned chromosome is in valid_chroms
            df_filtered = df_filtered[clean.isin(valid_chroms)]

            # Re-add the "chr" prefix to the kept rows
            clean = df_filtered[col].str.replace("chr", "", case=False)
            df_filtered[col] = "chr" + clean

        # Reset row numbers (gone gappy after filtering)
        df_filtered = df_filtered.reset_index(drop=True)

        # ---------------------------------------------------------------------
        # Step 7j: Write this model's results to disk
        # ---------------------------------------------------------------------
        # We append to the open file. Only write headers for the first model.
        df_filtered.to_csv(
            out_fh,
            sep="\t",
            index=False,
            header=(not wrote_header)
        )
        wrote_header = True
        any_results = True

        # Free memory before going to the next model
        del df_model, df_trans_m, df_ctrl_m, df_filtered, model_results, bead_index

    # Close the output file (we're done with this condition)
    out_fh.close()

    # -------------------------------------------------------------------------
    # Step 8: Aggregate across all models
    # -------------------------------------------------------------------------
    # If we got no results, skip aggregation
    if not any_results:
        print("No gene pairs found for this condition, skipping")
        continue

    print("")
    print("Saved filtered distances to: " + output_file)

    # Re-read everything we wrote, so we can aggregate
    df_all = pd.read_csv(output_file, sep="\t")
    print("Total raw gene pair distances: " + str(len(df_all)))

    # ---------------------------------------------------------------------
    # Group by gene pair and compute mean/std/count of distances
    # ---------------------------------------------------------------------
    # For each unique gene pair, we want to know:
    # - the mean distance across all models
    # - the standard deviation across models
    # - how many models contributed
    group_cols = [
        "gene1", "chr1", "gene1_region_type",
        "gene2", "chr2", "gene2_region_type"
    ]

    df_agg = df_all.groupby(group_cols, as_index=False)["distance"].agg(
        mean_distance="mean",
        std_distance="std",
        n_models="count"
    )

    print("Aggregated to " + str(len(df_agg)) + " unique gene pairs")

    # Make sure the output folder exists
    agg_folder = os.path.dirname(output_agg_file)
    os.makedirs(agg_folder, exist_ok=True)

    # Save the aggregated table
    df_agg.to_csv(output_agg_file, sep="\t", index=False)
    print("Saved aggregated distances to: " + output_agg_file)


# Done!
print("")
print("Done!")
