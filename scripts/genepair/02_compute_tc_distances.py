"""
02_compute_tc_distances.py
===========================
Script 2 of 4 for gene pair distance calculations – 3D distance calculation for translocated conditions (T1, C1).

What this script does
---------------------
For each translocated condition (T1 and C1), this script:

1. Reads the gene expression TSV produced by script 01.
2. For each 3D structural model (.cmm file), maps every gene to its nearest
   genomic bead in the model by midpoint position.
3. Computes pairwise Euclidean distances between:
      - Translocated genes  <-->  neighbor genes  (the comparison of interest)
      - Control genes       <-->  neighbor genes  (the baseline)
4. Writes each model's distances to disk immediately to avoid holding all
   models in memory at once (important on machines with limited RAM).
5. After all models are processed, aggregates distances across models
   (mean, std, count) and saves a summary TSV.

Why stream to disk?
-------------------
With ~4000 translocated genes x ~900 neighbor genes, the distance matrix
per model is large. Accumulating all 10 models in memory before writing
would exceed available RAM on most laptops. Instead, each model's results
are written as soon as they are computed and then freed from memory.

Why use cdist instead of numpy broadcasting?
--------------------------------------------
The intuitive numpy approach (A[:,None,:] - B[None,:,:]) silently creates a
temporary (N, M, 3) array in memory before computing the norm. With thousands
of genes this intermediate array alone can exhaust RAM. scipy's cdist computes
the same result without that intermediate, keeping peak memory proportional
to the output matrix (N, M) rather than (N, M, 3).

Why pre-build a bead index?
----------------------------
For each gene we need to find the nearest bead on the same chromosome. Doing
this by filtering a DataFrame on every gene is slow and creates many temporary
objects. Instead, we convert the bead coordinates to numpy arrays once per
model (build_bead_index), then all lookups are simple array operations.

Output of this script is used by scripts 03 and 04.

Usage
-----
    1. Run script 01 first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 02_compute_tc_distances.py

Dependencies
------------
    pandas, numpy, scipy
    Install with: pip install pandas numpy scipy

"""

import os
import re
import glob

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    # Per-condition inputs
    "conditions": {
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
    },

    # "inter" = only interchromosomal pairs (translocated chr ≠ neighbor chr)
    # "intra" = only intrachromosomal pairs (same chromosome)
    "scenario": "inter",
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def parse_cmm(file_path: str) -> pd.DataFrame:
    """
    Parse a UCSF Chimera .cmm marker file and return a DataFrame of beads
    with columns: chr, start, end, x, y, z.

    Each <marker> line with a beadID of the form 'chrN:start-end' represents
    one genomic bead in the 3D structural model, positioned at (x, y, z).
    """
    markers = []
    with open(file_path, "r") as fh:
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
                start_i, end_i = map(int, coords.split("-"))
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, start_i, end_i, x, y, z])
            except Exception:
                continue  # skip malformed lines silently

    return pd.DataFrame(markers, columns=["chr", "start", "end", "x", "y", "z"])


def build_bead_index(beads_by_chr: dict) -> dict:
    """
    Pre-convert each chromosome's bead DataFrame into numpy arrays for fast
    nearest-bead lookup.

    For each gene we need to find the closest bead by genomic midpoint.
    Pre-extracting midpoints and coordinates as numpy arrays (once per model)
    avoids slow and memory-heavy DataFrame filtering inside the gene loop.

    Returns a dict: chromosome → {'mids': array, 'coords': array, 'chr': str}
    """
    index = {}
    for chrom, bdf in beads_by_chr.items():
        index[chrom] = {
            "mids":   ((bdf["start"] + bdf["end"]) / 2).to_numpy(),
            "coords": bdf[["x", "y", "z"]].to_numpy(),
            "chr":    chrom,
        }
    return index


def map_genes_to_beads(
    genes_df: pd.DataFrame,
    bead_index: dict,
    name_col: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    For each gene, find the nearest bead on the same chromosome (by midpoint)
    and record its 3D coordinates.

    Returns
    -------
    coords    : (N, 3) float array – x, y, z of the matched bead (NaN if none)
    names     : (N,)  str array   – gene names or IDs
    chr_mapped: (N,)  str array   – chromosome of the matched bead
    """
    if genes_df.shape[0] == 0:
        return np.zeros((0, 3)), np.array([], dtype=str), np.array([], dtype=str)

    names      = genes_df[name_col].to_numpy(str)
    mids       = ((genes_df["start"] + genes_df["end"]) / 2).to_numpy()
    coords     = np.full((len(names), 3), np.nan)
    chr_mapped = np.full(len(names), "", dtype=object)

    for i, (chrom, mid) in enumerate(zip(genes_df["chrom"], mids)):
        entry = bead_index.get(chrom)
        if entry is None:
            continue
        nearest       = np.abs(entry["mids"] - mid).argmin()
        coords[i]     = entry["coords"][nearest]
        chr_mapped[i] = entry["chr"]

    return coords, names, chr_mapped


def pairwise_distances(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Compute the Euclidean distance between every pair of 3D points in A and B.
    Returns an (|A|, |B|) matrix where entry [i, j] = distance(A[i], B[j]).

    Uses scipy's cdist rather than numpy broadcasting (A[:,None,:] - B[None,:,:])
    because cdist avoids creating a large (N, M, 3) intermediate array, which
    would spike memory usage when N and M are both in the thousands.
    """
    if A.size == 0 or B.size == 0:
        return np.zeros((A.shape[0], B.shape[0]))
    return cdist(A, B, metric="euclidean")


# =============================================================================
# Main loop over conditions
# =============================================================================

for cond, paths in CONFIG["conditions"].items():
    print(f"\n===== Processing {cond} =====")

    # Load gene expression and neighbor BED files
    df_genes        = pd.read_csv(paths["gene_expression_file"], sep="\t")
    df_neighbor_bed = pd.read_csv(paths["neighbor_bed_file"], sep="\t")
    print(f"Loaded {len(df_genes)} genes, {len(df_neighbor_bed)} neighbor regions")

    # Use gene_name if available, otherwise fall back to gene_id
    name_col = "gene_name" if "gene_name" in df_genes.columns else "gene_id"

    # Ensure consistent chromosome notation and numeric coordinates
    df_genes["chrom"] = df_genes["chrom"].astype(str).str.strip().apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )
    df_genes["start"]       = pd.to_numeric(df_genes["start"], errors="coerce")
    df_genes["end"]         = pd.to_numeric(df_genes["end"],   errors="coerce")
    df_genes["region_type"] = df_genes["region_type"].str.strip().str.lower()

    df_neighbor_bed["neighbor_chr"] = df_neighbor_bed["neighbor_chr"].astype(str).str.strip().apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )
    df_neighbor_bed["chrom"] = df_neighbor_bed["chrom"].astype(str).str.strip()

    # Find all .cmm structural model files for this condition
    cmm_files = sorted(glob.glob(os.path.join(paths["cmm_dir"], "*.cmm")))
    print(f"Found {len(cmm_files)} CMM files to process")

    # Open the output file once and stream each model's results directly to
    # disk so we never hold all models in memory simultaneously
    os.makedirs(os.path.dirname(paths["output_filtered"]), exist_ok=True)
    wrote_header = False
    any_results  = False

    with open(paths["output_filtered"], "w") as out_fh:
        for f in cmm_files:
            print(f"\nProcessing model: {os.path.basename(f)}")
            beads = parse_cmm(f)
            if beads.empty:
                print("  No markers found, skipping")
                continue

            # Build numpy bead index and immediately free the raw DataFrame
            beads_by_chr = {c: b for c, b in beads.groupby("chr")}
            bead_index   = build_bead_index(beads_by_chr)
            del beads, beads_by_chr

            # Separate genes by region type
            trans_genes    = df_genes[df_genes["region_type"] == "inside_transloc"].copy()
            neighbor_genes = df_genes[df_genes["region_type"] == "neighbor"].copy()
            control_genes  = df_genes[df_genes["region_type"] == "control"].copy()

            # Map each gene group to its nearest bead in this model
            trans_coords, trans_names, trans_chr = map_genes_to_beads(trans_genes,    bead_index, name_col)
            neigh_coords, neigh_names, neigh_chr = map_genes_to_beads(neighbor_genes,  bead_index, name_col)
            ctrl_coords,  ctrl_names,  ctrl_chr  = map_genes_to_beads(control_genes,   bead_index, name_col)

            print(f"  Mapped {len(trans_names)} trans genes, "
                  f"{len(neigh_names)} neighbor genes, "
                  f"{len(ctrl_names)} control genes")

            model_results = []

            # ── Per-translocation distance calculation ────────────────────────
            for tid in trans_genes["transloc_id"].unique():
                neighbors_tid = df_neighbor_bed[df_neighbor_bed["transloc_id"] == tid]

                for _, nb_row in neighbors_tid.iterrows():
                    n_chr    = nb_row["neighbor_chr"]
                    nb_start = nb_row["start"]
                    nb_end   = nb_row["end"]

                    # Select neighbor genes that fall within this neighbor region
                    neigh_mask = (
                        (neighbor_genes["chrom"] == n_chr) &
                        (neighbor_genes["start"] <= nb_end) &
                        (neighbor_genes["end"]   >= nb_start)
                    )
                    if not np.any(neigh_mask):
                        continue

                    neigh_sub = neighbor_genes[neigh_mask].reset_index()

                    # ── Translocated genes ↔ neighbors ──────────────────────
                    trans_mask = trans_genes["transloc_id"] == tid
                    D_trans    = None
                    if np.any(trans_mask):
                        trans_sub = trans_genes[trans_mask].reset_index()
                        D_trans   = pairwise_distances(
                            trans_coords[trans_mask.to_numpy()],
                            neigh_coords[neigh_mask.to_numpy()],
                        )
                        t_idx, n_idx = np.meshgrid(trans_sub.index, neigh_sub.index, indexing="ij")
                        model_results.append(pd.DataFrame({
                            "translocation":     tid,
                            "gene1":             trans_sub.loc[t_idx.ravel(), name_col].values,
                            "chr1":              trans_sub.loc[t_idx.ravel(), "chrom"].values,
                            "gene1_region_type": "inside_transloc",
                            "gene1_chr_status":  "translocated",
                            "gene2":             neigh_sub.loc[n_idx.ravel(), name_col].values,
                            "chr2":              neigh_sub.loc[n_idx.ravel(), "chrom"].values,
                            "gene2_region_type": "neighbor",
                            "gene2_chr_status":  "non_translocated",
                            "distance":          D_trans.ravel(),
                        }))

                    # ── Control genes ↔ neighbors (sub-sampled to match trans count) ──
                    ctrl_mask = control_genes["transloc_id"] == tid
                    if not np.any(ctrl_mask):
                        continue

                    ctrl_sub  = control_genes[ctrl_mask].reset_index()
                    D_ctrl    = pairwise_distances(
                        ctrl_coords[ctrl_mask.to_numpy()],
                        neigh_coords[neigh_mask.to_numpy()],
                    )
                    n_trans_pairs = D_trans.size if D_trans is not None else 0
                    n_ctrl_pairs  = D_ctrl.size

                    if n_ctrl_pairs > n_trans_pairs > 0:
                        # Sub-sample control pairs so the group sizes match.
                        # The seed is derived from the translocation ID so
                        # results are reproducible across runs.
                        ctrl_idx, neigh_idx = np.meshgrid(ctrl_sub.index, neigh_sub.index, indexing="ij")
                        ctrl_idx_flat  = ctrl_idx.ravel()
                        neigh_idx_flat = neigh_idx.ravel()

                        rng        = np.random.default_rng(abs(hash(tid)) % (2**32))
                        sample_idx = rng.choice(len(ctrl_idx_flat), size=n_trans_pairs, replace=False)

                        ctrl_idx_sample  = ctrl_idx_flat[sample_idx]
                        neigh_idx_sample = neigh_idx_flat[sample_idx]
                        D_ctrl_sample    = D_ctrl.ravel()[sample_idx]
                    else:
                        ctrl_idx_sample, neigh_idx_sample = np.meshgrid(
                            ctrl_sub.index, neigh_sub.index, indexing="ij"
                        )
                        ctrl_idx_sample  = ctrl_idx_sample.ravel()
                        neigh_idx_sample = neigh_idx_sample.ravel()
                        D_ctrl_sample    = D_ctrl.ravel()

                    model_results.append(pd.DataFrame({
                        "translocation":     tid,
                        "gene1":             ctrl_sub.loc[ctrl_idx_sample, name_col].values,
                        "chr1":              ctrl_sub.loc[ctrl_idx_sample, "chrom"].values,
                        "gene1_region_type": "control",
                        "gene1_chr_status":  "non_translocated",
                        "gene2":             neigh_sub.loc[neigh_idx_sample, name_col].values,
                        "chr2":              neigh_sub.loc[neigh_idx_sample, "chrom"].values,
                        "gene2_region_type": "neighbor",
                        "gene2_chr_status":  "non_translocated",
                        "distance":          D_ctrl_sample,
                    }))

            if not model_results:
                del bead_index
                continue

            # ── Scenario filtering ────────────────────────────────────────────
            df_model = pd.concat(model_results, ignore_index=True)

            # Separate transloc-neighbor pairs from control-neighbor pairs
            df_trans_m = df_model[
                ((df_model["gene1_region_type"] == "inside_transloc") & (df_model["gene2_region_type"] == "neighbor")) |
                ((df_model["gene1_region_type"] == "neighbor") & (df_model["gene2_region_type"] == "inside_transloc"))
            ]
            df_ctrl_m = df_model[
                ((df_model["gene1_region_type"] == "control") & (df_model["gene2_region_type"] == "neighbor")) |
                ((df_model["gene1_region_type"] == "neighbor") & (df_model["gene2_region_type"] == "control"))
            ]

            if CONFIG["scenario"] == "inter":
                df_trans_m = df_trans_m[df_trans_m["chr1"] != df_trans_m["chr2"]]
            elif CONFIG["scenario"] == "intra":
                df_trans_m = df_trans_m[df_trans_m["chr1"] == df_trans_m["chr2"]]

            # Control–neighbor pairs are always interchromosomal by design
            df_ctrl_m = df_ctrl_m[df_ctrl_m["chr1"] != df_ctrl_m["chr2"]]

            df_filtered_model = pd.concat([df_trans_m, df_ctrl_m], ignore_index=True)

            # Remove non-canonical chromosomes (scaffolds, patches, etc.)
            for col in ("chr1", "chr2"):
                clean = df_filtered_model[col].str.replace("chr", "", case=False)
                valid = {str(i) for i in range(1, 23)} | {"X", "Y"}
                df_filtered_model = df_filtered_model[clean.isin(valid)]
                df_filtered_model[col] = "chr" + clean[clean.isin(valid)]
            df_filtered_model = df_filtered_model.reset_index(drop=True)

            # Write this model's results to disk and free memory
            df_filtered_model.to_csv(out_fh, sep="\t", index=False, header=not wrote_header)
            wrote_header = True
            any_results  = True

            del df_model, df_trans_m, df_ctrl_m, df_filtered_model, model_results, bead_index

    if not any_results:
        print("No gene pairs found for this condition, skipping")
        continue

    print(f"\nSaved filtered distances to: {paths['output_filtered']}")

    # ── Aggregate across all models ───────────────────────────────────────────
    # Read the streamed file back once to compute mean/std/count per gene pair
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
