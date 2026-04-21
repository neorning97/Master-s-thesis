"""
02_compute_tc_distances.py
===========================
Script 2 of 4 – 3D distance calculation for translocated conditions (T1, C1).

For each condition, this script:
1. Reads the gene expression TSV from script 01.
2. For each 3D model (.cmm file), maps genes to their nearest bead.
3. Computes pairwise distances between translocated/control genes and neighbors.
4. Writes each model's distances to disk immediately (to save memory).
5. Aggregates distances across models and saves a summary TSV.

Usage:  Edit the paths below, then run: python 02_compute_tc_distances.py
Dependencies:  pip install pandas numpy scipy
"""

import os
import re
import glob

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


# ============================================================================================
# CONFIG – edit these paths before running, and change SCENARIO depending on what you analyze
# ============================================================================================

# "inter" = only interchromosomal pairs, "intra" = only intrachromosomal pairs
SCENARIO = "inter"

CONDITIONS = {
    "T1": {
        "gene_expression_file": "/path/to/T1_genes_expression.tsv",
        "neighbor_bed_file":    "/path/to/T1_neighbor.bed",
        "cmm_dir":              "/path/to/cmm/T1",
        "output_filtered":      "/path/to/T1_distances.tsv",
        "output_agg":           "/path/to/T1_distances_agg.tsv",
    },
    "C1": {
        "gene_expression_file": "/path/to/C1_genes_expression.tsv",
        "neighbor_bed_file":    "/path/to/C1_neighbor.bed",
        "cmm_dir":              "/path/to/cmm/C1",
        "output_filtered":      "/path/to/C1_distances.tsv",
        "output_agg":           "/path/to/C1_distances_agg.tsv",
    },
}


# =============================================================================
# Process each condition independently
# =============================================================================

for cond, paths in CONDITIONS.items():
    print(f"\n===== Processing {cond} =====")

    #-----------------------------------------------------------------
    # Load gene data and neighbor regions
    #-----------------------------------------------------------------
    df_genes = pd.read_csv(paths["gene_expression_file"], sep="\t")
    df_neighbor_bed = pd.read_csv(paths["neighbor_bed_file"], sep="\t")
    print(f"Loaded {len(df_genes)} genes, {len(df_neighbor_bed)} neighbor regions")

    # Use gene_name if available, otherwise fall back to gene_id
    name_col = "gene_name" if "gene_name" in df_genes.columns else "gene_id"

    #-----------------------------------------------------------------
    # Clean and standardize gene coordinates
    #-----------------------------------------------------------------
    # Makes sure chromosome notation is consistent (e.g., chr1, chrX, etc.)
    df_genes["chrom"] = df_genes["chrom"].astype(str).str.strip().apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )
    # Convert coordinates to numeric (invalid values --> NaN)
    df_genes["start"] = pd.to_numeric(df_genes["start"], errors="coerce")
    df_genes["end"] = pd.to_numeric(df_genes["end"], errors="coerce")

    # Normalize region_type (e.g., remove whitespace, lowercase)
    df_genes["region_type"] = df_genes["region_type"].str.strip().str.lower()

    # Clean neighbor BED chromosome columns
    df_neighbor_bed["neighbor_chr"] = df_neighbor_bed["neighbor_chr"].astype(str).str.strip().apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )
    df_neighbor_bed["chrom"] = df_neighbor_bed["chrom"].astype(str).str.strip()

    #-----------------------------------------------------------------
    # Find all 3D genome model files (.cmm)
    #-----------------------------------------------------------------
    cmm_files = sorted(glob.glob(os.path.join(paths["cmm_dir"], "*.cmm")))
    print(f"Found {len(cmm_files)} CMM files to process")

    # Prepare output file (write incrementally to save memory)
    os.makedirs(os.path.dirname(paths["output_filtered"]), exist_ok=True)
    wrote_header = False    # Track if header already written
    any_results = False     # Track if any data was generated

    # Open output file once and append model results as we go
    with open(paths["output_filtered"], "w") as out_fh:

        #-----------------------------------------------------------------
        # Loop over each 3D model
        #-----------------------------------------------------------------
        for f in cmm_files:
            print(f"\nProcessing model: {os.path.basename(f)}")

            #-----------------------------------------------------------------
            # Parse .cmm file (3D bead position)
            #-----------------------------------------------------------------
            markers = []
            with open(f, "r") as fh:
                for line in fh:

                    # Only process lines containing bead markers
                    if "<marker" not in line:
                        continue
                    # Extract attributes from XML-like line
                    attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
                    # Skip if no beadID
                    if "beadID" not in attrs:
                        continue
                    try:
                        # beadID format: chr:start-end
                        chrom_part, coords = attrs["beadID"].split(":")
                        # Ensure chromosome naming is consistent
                        if not chrom_part.startswith("chr"):
                            chrom_part = "chr" + chrom_part
                        # Extract genomic coordinates
                        start_i, end_i = map(int, coords.split("-"))
                        # Extract 3D coordinates (x, y, z)
                        x = float(attrs.get("x", "nan"))
                        y = float(attrs.get("y", "nan"))
                        z = float(attrs.get("z", "nan"))
                        markers.append([chrom_part, start_i, end_i, x, y, z])
                    except Exception:
                        # Skip malformed entries
                        continue

            # Convert to DataFrame
            beads = pd.DataFrame(markers, columns=["chr", "start", "end", "x", "y", "z"])
            if beads.empty:
                print("  No markers found, skipping")
                continue

            #-----------------------------------------------------------------
            # Build fast lookup: chromosome --> bead midpoints + coordinates
            #-----------------------------------------------------------------
            bead_index = {}
            for chrom, bdf in beads.groupby("chr"):
                bead_index[chrom] = {
                    "mids": ((bdf["start"] + bdf["end"]) / 2).to_numpy(),   # Genomic midpoints
                    "coords": bdf[["x", "y", "z"]].to_numpy(),              # 3D positions
                    "chr": chrom,
                }
            # Free memory
            del beads

            #-----------------------------------------------------------------
            # Split genes by region type
            #-----------------------------------------------------------------
            trans_genes = df_genes[df_genes["region_type"] == "inside_transloc"].copy()
            neighbor_genes = df_genes[df_genes["region_type"] == "neighbor"].copy()
            control_genes = df_genes[df_genes["region_type"] == "control"].copy()

            #-----------------------------------------------------------------
            # Map genes --> nearest 3D bead
            #-----------------------------------------------------------------
            def map_genes(genes_df):
                """
                Maps each gene to the nearest bead (by genomic midpoint).
                Returns:
                    coords: Nx3 array of 3D coordinates
                    names: gene names
                    chr_mapped: chromosome assigned
                """
                # Handles empty input
                if genes_df.shape[0] == 0:
                    return np.zeros((0, 3)), np.array([], dtype=str), np.array([], dtype=str)

                names = genes_df[name_col].to_numpy(str)
                # Compute gene midpoints
                mids = ((genes_df["start"] + genes_df["end"]) / 2).to_numpy()

                coords = np.full((len(names), 3), np.nan)
                chr_mapped = np.full(len(names), "", dtype=object)

                # Loop though genes
                for i, (chrom, mid) in enumerate(zip(genes_df["chrom"], mids)):
                    entry = bead_index.get(chrom)
                    if entry is None:
                        continue
                    # Find closest bead midpoint
                    nearest = np.abs(entry["mids"] - mid).argmin()

                    coords[i] = entry["coords"][nearest]
                    chr_mapped[i] = entry["chr"]

                return coords, names, chr_mapped
            
            # Apply mapping to each gene category
            trans_coords, trans_names, trans_chr = map_genes(trans_genes)
            neigh_coords, neigh_names, neigh_chr = map_genes(neighbor_genes)
            ctrl_coords, ctrl_names, ctrl_chr = map_genes(control_genes)

            print(f"  Mapped {len(trans_names)} trans, "
                  f"{len(neigh_names)} neighbor, "
                  f"{len(ctrl_names)} control genes")

            model_results = []

            #-----------------------------------------------------------------
            # Compute pairwise distances
            #-----------------------------------------------------------------
            for tid in trans_genes["transloc_id"].unique():
                # Get neighbor regions linked to this translocation
                neighbors_tid = df_neighbor_bed[df_neighbor_bed["transloc_id"] == tid]

                for _, nb_row in neighbors_tid.iterrows():
                    # Define neighbor region
                    n_chr = nb_row["neighbor_chr"]
                    nb_start = nb_row["start"]
                    nb_end = nb_row["end"]

                    # Find neighbor genes in this region
                    neigh_mask = (
                        (neighbor_genes["chrom"] == n_chr) &
                        (neighbor_genes["start"] <= nb_end) &
                        (neighbor_genes["end"] >= nb_start)
                    )
                    if not np.any(neigh_mask):
                        continue

                    neigh_sub = neighbor_genes[neigh_mask].reset_index()

                    #-----------------------------------------------------------------
                    # Translocated <--> Neighbor pairs
                    #-----------------------------------------------------------------
                    trans_mask = trans_genes["transloc_id"] == tid
                    D_trans = None
                    if np.any(trans_mask):
                        trans_sub = trans_genes[trans_mask].reset_index()

                        # cdist computes all pairwise distances without a huge intermediate array
                        D_trans = cdist(
                            trans_coords[trans_mask.to_numpy()],
                            neigh_coords[neigh_mask.to_numpy()],
                            metric="euclidean",
                        )

                        # Create all pair combinations
                        t_idx, n_idx = np.meshgrid(trans_sub.index, neigh_sub.index, indexing="ij")
                        
                        model_results.append(pd.DataFrame({
                            "translocation": tid,
                            "gene1": trans_sub.loc[t_idx.ravel(), name_col].values,
                            "chr1": trans_sub.loc[t_idx.ravel(), "chrom"].values,
                            "gene1_region_type": "inside_transloc",
                            "gene1_chr_status": "translocated",
                            "gene2": neigh_sub.loc[n_idx.ravel(), name_col].values,
                            "chr2": neigh_sub.loc[n_idx.ravel(), "chrom"].values,
                            "gene2_region_type": "neighbor",
                            "gene2_chr_status": "non_translocated",
                            "distance": D_trans.ravel(),
                        }))

                    #-----------------------------------------------------------------
                    # Control <--> Neighbor pairs
                    #-----------------------------------------------------------------
                    ctrl_mask = control_genes["transloc_id"] == tid
                    if not np.any(ctrl_mask):
                        continue

                    ctrl_sub = control_genes[ctrl_mask].reset_index()
                    D_ctrl = cdist(
                        ctrl_coords[ctrl_mask.to_numpy()],
                        neigh_coords[neigh_mask.to_numpy()],
                        metric="euclidean",
                    )

                    # Match number of control pairs to translocated pairs
                    n_trans_pairs = D_trans.size if D_trans is not None else 0
                    n_ctrl_pairs = D_ctrl.size

                    if n_ctrl_pairs > n_trans_pairs > 0:
                        ctrl_idx, neigh_idx = np.meshgrid(ctrl_sub.index, neigh_sub.index, indexing="ij")
                        ctrl_idx_flat = ctrl_idx.ravel()
                        neigh_idx_flat = neigh_idx.ravel()

                        rng = np.random.default_rng(abs(hash(tid)) % (2**32))
                        sample_idx = rng.choice(len(ctrl_idx_flat), size=n_trans_pairs, replace=False)

                        ctrl_idx_sample = ctrl_idx_flat[sample_idx]
                        neigh_idx_sample = neigh_idx_flat[sample_idx]
                        D_ctrl_sample = D_ctrl.ravel()[sample_idx]
                    else:
                        ctrl_idx_sample, neigh_idx_sample = np.meshgrid(
                            ctrl_sub.index, neigh_sub.index, indexing="ij"
                        )
                        ctrl_idx_sample = ctrl_idx_sample.ravel()
                        neigh_idx_sample = neigh_idx_sample.ravel()
                        D_ctrl_sample = D_ctrl.ravel()

                    model_results.append(pd.DataFrame({
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
                    }))

            if not model_results:
                del bead_index
                continue

            #-----------------------------------------------------------------
            # Combine and filter by scenario (inter- and intrachromosomal)
            #-----------------------------------------------------------------
            df_model = pd.concat(model_results, ignore_index=True)

            # Separate transloc-neighbor and control-neighbor pairs
            is_trans_neigh = (
                ((df_model["gene1_region_type"] == "inside_transloc") & (df_model["gene2_region_type"] == "neighbor")) |
                ((df_model["gene1_region_type"] == "neighbor") & (df_model["gene2_region_type"] == "inside_transloc"))
            )
            is_ctrl_neigh = (
                ((df_model["gene1_region_type"] == "control") & (df_model["gene2_region_type"] == "neighbor")) |
                ((df_model["gene1_region_type"] == "neighbor") & (df_model["gene2_region_type"] == "control"))
            )

            df_trans_m = df_model[is_trans_neigh]
            df_ctrl_m = df_model[is_ctrl_neigh]

            # Filter inter/intra chromosomal pairs
            if SCENARIO == "inter":
                df_trans_m = df_trans_m[df_trans_m["chr1"] != df_trans_m["chr2"]]
            elif SCENARIO == "intra":
                df_trans_m = df_trans_m[df_trans_m["chr1"] == df_trans_m["chr2"]]

            # Control-neighbor pairs are always interchromosomal
            df_ctrl_m = df_ctrl_m[df_ctrl_m["chr1"] != df_ctrl_m["chr2"]]

            df_filtered = pd.concat([df_trans_m, df_ctrl_m], ignore_index=True)

            # Remove non-canonical chromosomes (scaffolds, patches, etc.)
            valid_chroms = {str(i) for i in range(1, 23)} | {"X", "Y"}
            for col in ("chr1", "chr2"):
                clean = df_filtered[col].str.replace("chr", "", case=False)
                df_filtered = df_filtered[clean.isin(valid_chroms)]
                df_filtered[col] = "chr" + clean[clean.isin(valid_chroms)]
            df_filtered = df_filtered.reset_index(drop=True)

            # Write to disk and free memory
            df_filtered.to_csv(out_fh, sep="\t", index=False, header=not wrote_header)
            wrote_header = True
            any_results = True

            del df_model, df_trans_m, df_ctrl_m, df_filtered, model_results, bead_index

    #-----------------------------------------------------------------
    # Aggregate across models
    #-----------------------------------------------------------------
    if not any_results:
        print("No gene pairs found for this condition, skipping")
        continue

    print(f"\nSaved filtered distances to: {paths['output_filtered']}")

    # -- Aggregate across all models (mean, std, count per gene pair) --
    df_all = pd.read_csv(paths["output_filtered"], sep="\t")
    print(f"Total raw gene pair distances: {len(df_all)}")

    group_cols = ["gene1", "chr1", "gene1_region_type", "gene2", "chr2", "gene2_region_type"]
    df_agg = df_all.groupby(group_cols, as_index=False)["distance"].agg(
        mean_distance="mean",
        std_distance="std",
        n_models="count",
    )
    print(f"Aggregated to {len(df_agg)} unique gene pairs")

    os.makedirs(os.path.dirname(paths["output_agg"]), exist_ok=True)
    df_agg.to_csv(paths["output_agg"], sep="\t", index=False)
    print(f"Saved aggregated distances to: {paths['output_agg']}")

print("\nDone!")
