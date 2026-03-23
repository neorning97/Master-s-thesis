import pandas as pd
import numpy as np
import os
import glob
import re

# -------------------------------
# User paths for each condition
# -------------------------------
conditions = {
    "T1": {
        "cmm_dir": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/T1",
        "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_genes_expression_subplot.tsv",
    #    "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_genes_expression_moved.tsv",
        "neighbor_bed_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/T1_neighbor1.bed",
    #    "neighbor_bed_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/T1_neighbor_originchr1.bed",
        "output_filtered": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_subplot.tsv",
        "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_agg_subplot.tsv"
    #    "output_filtered": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_moved.tsv",
    #    "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_agg_moved.tsv"
    },
    "C1": {
        "cmm_dir": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/C1",
        "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_genes_expression_subplot.tsv",
    #    "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_genes_expression_moved.tsv",
        "neighbor_bed_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/C1_neighbor1.bed",
    #    "neighbor_bed_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/C1_neighbor_originchr1.bed",
        "output_filtered": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_subplot.tsv",
        "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_agg_subplot.tsv"
    #    "output_filtered": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_moved.tsv",
    #    "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_agg_moved.tsv"
    }
}

# -------------------------------
# Scenario parameter
# -------------------------------
# "inter" = interchromosomal (chr1 != chr2)
# "intra" = intrachromosomal (chr1 == chr2)
scenario = "inter"

# -------------------------------
# Functions
# -------------------------------
def parse_cmm(file_path):
    markers = []
    with open(file_path,"r") as fh:
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
                start_s, end_s = coords.split("-")
                start_i, end_i = int(start_s), int(end_s)
                x = float(attrs.get("x","nan"))
                y = float(attrs.get("y","nan"))
                z = float(attrs.get("z","nan"))
                markers.append([chrom_part, start_i, end_i, x, y, z])
            except:
                continue
    return pd.DataFrame(markers, columns=["chr","start","end","x","y","z"])

def map_genes_to_beads(genes_df, beads_by_chr, name_col):
    if genes_df.shape[0]==0:
        return np.zeros((0,3)), np.array([],dtype=str), np.array([],dtype=str)
    
    names = genes_df[name_col].to_numpy(str)
    mids = ((genes_df["start"]+genes_df["end"])/2).to_numpy()
    coords = np.full((len(names),3), np.nan)
    chr_mapped = np.full(len(names), "", dtype=object)

    for i,(chrom,mid) in enumerate(zip(genes_df["chrom"], mids)):
        beads_chr = beads_by_chr.get(chrom)
        if beads_chr is None or beads_chr.empty:
            continue
        bead_mids = ((beads_chr["start"]+beads_chr["end"])/2).to_numpy()
        idx = np.abs(bead_mids - mid).argmin()
        coords[i] = beads_chr.iloc[idx][["x","y","z"]].to_numpy()
        chr_mapped[i] = beads_chr.iloc[idx]["chr"]

    return coords, names, chr_mapped

def pairwise_distances(A,B):
    if A.size==0 or B.size==0:
        return np.zeros((A.shape[0], B.shape[0]))
    return np.linalg.norm(A[:,None,:]-B[None,:,:], axis=2)

# -------------------------------
# Main loop over conditions
# -------------------------------
for cond, paths in conditions.items():
    print(f"\n===== Processing {cond} =====")

    # Load gene expression and neighbor BED
    df_genes = pd.read_csv(paths["gene_expression_file"], sep="\t")
    df_neighbor_bed = pd.read_csv(paths["neighbor_bed_file"], sep="\t")
    print(f"Loaded {len(df_genes)} genes, {len(df_neighbor_bed)} neighbor regions")

    # Decide which column to use for gene names
    name_col = "gene_name" if "gene_name" in df_genes.columns else "gene_id"

    # Clean columns
    df_genes["chrom"] = df_genes["chrom"].astype(str).str.strip().apply(lambda x: x if x.startswith("chr") else "chr"+x)
    df_genes["start"] = pd.to_numeric(df_genes["start"], errors="coerce")
    df_genes["end"] = pd.to_numeric(df_genes["end"], errors="coerce")
    df_genes["region_type"] = df_genes["region_type"].str.strip().str.lower()

    df_neighbor_bed["neighbor_chr"] = df_neighbor_bed["neighbor_chr"].astype(str).str.strip().apply(lambda x: x if x.startswith("chr") else "chr"+x)
    df_neighbor_bed["chrom"] = df_neighbor_bed["chrom"].astype(str).str.strip()

    # Load CMM files
    cmm_files = sorted(glob.glob(os.path.join(paths["cmm_dir"], "*.cmm")))
    print(f"Found {len(cmm_files)} CMM files to process")

    all_results = []

    for f in cmm_files:
        print(f"\nProcessing model: {os.path.basename(f)}")
        beads = parse_cmm(f)
        if beads.empty:
            print("  No markers found, skipping")
            continue

        beads_by_chr = {c: b for c, b in beads.groupby("chr")}

        # Separate gene types
        trans_genes = df_genes[df_genes["region_type"]=="inside_transloc"].copy()
        neighbor_genes = df_genes[df_genes["region_type"]=="neighbor"].copy()
        control_genes = df_genes[df_genes["region_type"]=="control"].copy()

        # Map genes to 3D coordinates
        trans_coords, trans_names, trans_chr = map_genes_to_beads(trans_genes, beads_by_chr, name_col)
        neigh_coords, neigh_names, neigh_chr = map_genes_to_beads(neighbor_genes, beads_by_chr, name_col)
        ctrl_coords, ctrl_names, ctrl_chr = map_genes_to_beads(control_genes, beads_by_chr, name_col)

        print(f"  Mapped {len(trans_names)} trans genes, {len(neigh_names)} neighbor genes, {len(ctrl_names)} control genes")

        # Loop over translocations
        for tid in trans_genes["transloc_id"].unique():
            neighbors_tid = df_neighbor_bed[df_neighbor_bed["transloc_id"]==tid]

            for _, nb_row in neighbors_tid.iterrows():
                n_chr = nb_row["neighbor_chr"]
                nb_start = nb_row["start"]
                nb_end = nb_row["end"]

                neigh_mask = (
                    (neighbor_genes["chrom"] == n_chr) &
                    (neighbor_genes["start"] <= nb_end) &
                    (neighbor_genes["end"] >= nb_start)
                )
                if not np.any(neigh_mask):
                    continue

                trans_mask = trans_genes["transloc_id"] == tid
                if np.any(trans_mask):
                    D_trans = pairwise_distances(trans_coords[trans_mask], neigh_coords[neigh_mask])

                    trans_sub = trans_genes[trans_mask].reset_index()
                    neigh_sub = neighbor_genes[neigh_mask].reset_index()

                    t_idx, n_idx = np.meshgrid(trans_sub.index, neigh_sub.index, indexing='ij')
                    df_tmp = pd.DataFrame({
                        "translocation": tid,
                        "gene1": trans_sub.loc[t_idx.ravel(), name_col].values,
                        "chr1": trans_sub.loc[t_idx.ravel(), "chrom"].values,
                        "gene1_region_type": "inside_transloc",
                        "gene1_chr_status": "translocated",
                        "gene2": neigh_sub.loc[n_idx.ravel(), name_col].values,
                        "chr2": neigh_sub.loc[n_idx.ravel(), "chrom"].values,
                        "gene2_region_type": "neighbor",
                        "gene2_chr_status": "non_translocated",
                        "distance": D_trans.ravel()
                    })
                    all_results.append(df_tmp)

                ctrl_mask = control_genes["transloc_id"] == tid
                if not np.any(ctrl_mask):
                    continue

                ctrl_sub = control_genes[ctrl_mask].reset_index()
                neigh_sub = neighbor_genes[neigh_mask].reset_index()
                D_ctrl = pairwise_distances(ctrl_coords[ctrl_mask], neigh_coords[neigh_mask])

                n_trans_pairs = D_trans.size if np.any(trans_mask) else 0
                n_ctrl_pairs = D_ctrl.size

                if n_ctrl_pairs > n_trans_pairs > 0:
                    ctrl_idx, neigh_idx = np.meshgrid(ctrl_sub.index, neigh_sub.index, indexing='ij')
                    ctrl_idx_flat = ctrl_idx.ravel()
                    neigh_idx_flat = neigh_idx.ravel()

                    rng = np.random.default_rng(abs(hash(tid)) % (2**32))
                    sample_idx = rng.choice(len(ctrl_idx_flat), size=n_trans_pairs, replace=False)

                    ctrl_idx_sample = ctrl_idx_flat[sample_idx]
                    neigh_idx_sample = neigh_idx_flat[sample_idx]

                    D_ctrl_sample = D_ctrl.ravel()[sample_idx]
                else:
                    ctrl_idx_sample, neigh_idx_sample = np.meshgrid(ctrl_sub.index, neigh_sub.index, indexing='ij')
                    ctrl_idx_sample = ctrl_idx_sample.ravel()
                    neigh_idx_sample = neigh_idx_sample.ravel()
                    D_ctrl_sample = D_ctrl.ravel()

                df_tmp = pd.DataFrame({
                    "translocation": tid,
                    "gene1": ctrl_sub.loc[ctrl_idx_sample, name_col].values,
                    "chr1": ctrl_sub.loc[ctrl_idx_sample, "chrom"].values,
                    "gene1_region_type": "control",
                    "gene1_chr_status": "non_translocated",
                    "gene2": neigh_sub.loc[neigh_idx_sample, name_col].values,
                    "chr2": neigh_sub.loc[neigh_idx_sample, "chrom"].values,
                    "gene2_region_type": "neighbor",
                    "gene2_chr_status": "non_translocated",
                    "distance": D_ctrl_sample
                })
                all_results.append(df_tmp)

    # Combine all results
    if len(all_results) == 0:
        print("No gene pairs found for this condition, skipping")
        continue
    df = pd.concat(all_results, ignore_index=True)
    print(f"Computed {len(df)} raw gene pair distances")

    # -------------------------------
    # Separate trans-neighbor and control-neighbor
    # -------------------------------
    df_trans = df[((df["gene1_region_type"]=="inside_transloc") & (df["gene2_region_type"]=="neighbor")) |
                  ((df["gene1_region_type"]=="neighbor") & (df["gene2_region_type"]=="inside_transloc"))]

    df_ctrl = df[((df["gene1_region_type"]=="control") & (df["gene2_region_type"]=="neighbor")) |
                 ((df["gene1_region_type"]=="neighbor") & (df["gene2_region_type"]=="control"))]

    print("Unique chr1-chr2-region_type combinations BEFORE scenario filtering:")
    print(df_trans[["chr1","chr2","gene1_region_type","gene2_region_type"]].drop_duplicates())

    # Apply scenario filter to trans-neighbor pairs
    if scenario == "inter":
        df_trans = df_trans[df_trans["chr1"] != df_trans["chr2"]]
    elif scenario == "intra":
        df_trans = df_trans[df_trans["chr1"] == df_trans["chr2"]]

    # Control-neighbor always interchromosomal
    df_ctrl = df_ctrl[df_ctrl["chr1"] != df_ctrl["chr2"]]

    # Combine back
    df_filtered = pd.concat([df_trans, df_ctrl], ignore_index=True)
    print(f"{len(df_filtered)} rows remain after filtering")

    print("Sample of remaining gene pairs:")
    print(df_filtered.head())

    # Valid chromosomes
    df_filtered['chr1_clean'] = df_filtered['chr1'].str.replace('chr','', case=False)
    df_filtered['chr2_clean'] = df_filtered['chr2'].str.replace('chr','', case=False)
    valid_chroms = [str(i) for i in range(1,23)] + ["X", "Y"]
    df_filtered = df_filtered[df_filtered['chr1_clean'].isin(valid_chroms) & df_filtered['chr2_clean'].isin(valid_chroms)]
    df_filtered['chr1'] = "chr" + df_filtered['chr1_clean']
    df_filtered['chr2'] = "chr" + df_filtered['chr2_clean']
    df_filtered.drop(columns=['chr1_clean','chr2_clean'], inplace=True)

    # Aggregate distances
    group_cols = ["gene1", "chr1", "gene1_region_type", "gene2", "chr2", "gene2_region_type"]
    df_agg = df_filtered.groupby(group_cols, as_index=False)["distance"].agg(
        mean_distance="mean",
        std_distance="std",
        n_models="count"
    )
    print(f"Aggregated to {len(df_agg)} unique gene pairs")

    # Save outputs
    os.makedirs(os.path.dirname(paths["output_filtered"]), exist_ok=True)
    df_filtered.to_csv(paths["output_filtered"], sep="\t", index=False)
    print(f"Saved filtered distances to: {paths['output_filtered']}")

    os.makedirs(os.path.dirname(paths["output_agg"]), exist_ok=True)
    df_agg.to_csv(paths["output_agg"], sep="\t", index=False)
    print(f"Saved aggregated distances to: {paths['output_agg']}")
