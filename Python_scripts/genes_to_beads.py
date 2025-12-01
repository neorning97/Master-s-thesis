# Script that maps the genes in the neighbor and translocation regions to beads in the .cmm files,
# and then calculates the pair-wise spatial distance between the genes in the neighbor and translocation regions.
# The output files contain the spatial distance between these gene pairs.

{"id":"52341","variant":"standard","title":"Calculate distances between translocated and neighbor genes"}
import pandas as pd
import numpy as np
import os
import glob
import re

# -------------------------------
# User paths
# -------------------------------
cmm_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/T1"
gene_expression_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.tsv"
neighbor_bed_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/T1_neighbor.bed"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_translocation_neighbor_distances.tsv"

# -------------------------------
# Load files
# -------------------------------
df_genes = pd.read_csv(gene_expression_file, sep="\t")
df_neighbor_bed = pd.read_csv(neighbor_bed_file, sep="\t")
print(f"Loaded gene expression file with {len(df_genes)} rows and {len(df_genes.columns)} columns")
print(f"Loaded neighbor BED file with {len(df_neighbor_bed)} rows")

# Clean columns
df_genes["chrom"] = df_genes["chrom"].astype(str).str.strip()
df_genes["chrom"] = df_genes["chrom"].apply(lambda x: x if x.startswith("chr") else "chr"+x)
df_genes["start"] = pd.to_numeric(df_genes["start"], errors="coerce")
df_genes["end"] = pd.to_numeric(df_genes["end"], errors="coerce")

df_neighbor_bed["neighbor_chr"] = df_neighbor_bed["neighbor_chr"].astype(str).str.strip()
df_neighbor_bed["chrom"] = df_neighbor_bed["chrom"].astype(str).str.strip()
df_neighbor_bed["neighbor_chr"] = df_neighbor_bed["neighbor_chr"].apply(lambda x: x if x.startswith("chr") else "chr"+x)

# Normalize region_type
df_genes["region_type"] = df_genes["region_type"].str.strip().str.lower()

# -------------------------------
# Parse CMM files
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

# -------------------------------
# Map genes to nearest bead
# -------------------------------
def map_genes_to_beads(genes_df, beads_df):
    if genes_df.shape[0]==0:
        return np.zeros((0,3)), np.array([],dtype=str), np.array([],dtype=str)
    names = genes_df["gene_name"].to_numpy(str)
    mids = ((genes_df["start"]+genes_df["end"])/2).to_numpy()
    coords = np.full((len(names),3), np.nan)
    chr_mapped = np.full(len(names), "", dtype=object)
    for i,(chrom,mid) in enumerate(zip(genes_df["chrom"], mids)):
        beads_chr = beads_df[beads_df["chr"]==chrom]
        if beads_chr.empty:
            continue
        bead_mids = ((beads_chr["start"]+beads_chr["end"])/2).to_numpy()
        idx = np.abs(bead_mids - mid).argmin()
        coords[i] = beads_chr.loc[beads_chr.index[idx], ["x","y","z"]].to_numpy()
        chr_mapped[i] = beads_chr.iloc[idx]["chr"]
    return coords, names, chr_mapped

# -------------------------------
# Compute distances
# -------------------------------
def pairwise_distances(A,B):
    if A.size==0 or B.size==0:
        return np.zeros((A.shape[0], B.shape[0]))
    return np.linalg.norm(A[:,None,:]-B[None,:,:], axis=2)

# -------------------------------
# Loop over models
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
print(f"Found {len(cmm_files)} CMM files to process")
all_results = []

for f in cmm_files:
    print(f"\nProcessing model: {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No markers found in this model, skipping")
        continue
    
    # Translocation genes
    trans_genes = df_genes[df_genes["region_type"]=="inside_transloc"]
    
    # Neighbor genes
    neighbor_genes = df_genes[df_genes["region_type"]=="neighbor"].copy()
    print("Number of neighbor genes:", len(neighbor_genes))

    
    # Map beads
    trans_coords, trans_names, trans_chr = map_genes_to_beads(trans_genes, beads)
    neigh_coords, neigh_names, neigh_chr = map_genes_to_beads(neighbor_genes, beads)

    print(f"  Mapped {len(trans_names)} trans genes and {len(neigh_names)} neighbor genes to beads")
    
    # Only keep neighbor genes on recipient chromosome of corresponding translocation
    for tid in trans_genes["transloc_id"].unique():
        print(f"    Processing translocation {tid}")
        trans_mask = trans_genes["transloc_id"]==tid
        # neighbor regions for this translocation
        neighbors_tid = df_neighbor_bed[df_neighbor_bed["transloc_id"]==tid]
        
        for _, nb_row in neighbors_tid.iterrows():
            n_chr = nb_row["neighbor_chr"]
            nb_start = nb_row["start"]
            nb_end = nb_row["end"]
            # Only neighbor genes in this BED region
            neigh_mask = (
                (neighbor_genes["chrom"] == n_chr) &
                (neighbor_genes["end"] >= nb_start) &   # gene ends after region starts
                (neighbor_genes["start"] <= nb_end)     # gene starts before region ends
            )

            print(f"    Neighbor genes retained in region {n_chr}:{nb_start}-{nb_end}: {neigh_mask.sum()}")

            # Compute distances
            if np.any(trans_mask) and np.any(neigh_mask):
                D = pairwise_distances(trans_coords[trans_mask], neigh_coords[neigh_mask])
                t_names = trans_names[trans_mask]
                n_names = neigh_names[neigh_mask]
                for i,tname in enumerate(t_names):
                    for j,nname in enumerate(n_names):
                        all_results.append({
                            "translocation": tid,
                            "gene1": tname,
                            "chr1": trans_genes.loc[trans_mask,"chrom"].iloc[i],
                            "gene1_region_type": "inside_transloc",
                            "gene1_chr_status": "translocated",
                            "gene2": nname,
                            "chr2": neighbor_genes.loc[neigh_mask,"chrom"].iloc[j],
                            "gene2_region_type": "neighbor",
                            "gene2_chr_status": "non_translocated",
                            "distance": D[i,j]
                        })

# -------------------------------
# Save raw per-model results
# -------------------------------
summary_df = pd.DataFrame(all_results)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
summary_df.to_csv(output_file, sep="\t", index=False)
print("Distance calculations saved to:", output_file)

# -------------------------------
# Save aggregate distances across models
# -------------------------------
if len(all_results) > 0:
    results_df = pd.DataFrame(all_results)
    agg_df = results_df.groupby(
        ["translocation", "gene1", "chr1", "gene2", "chr2", "gene1_region_type", "gene2_region_type",
         "gene1_chr_status", "gene2_chr_status"]
    ).agg(
        mean_distance=("distance", "mean"),
        std_distance=("distance", "std"),
        n_models=("distance", "count")
    ).reset_index()
    
    agg_output_file = output_file.replace(".tsv", "_aggregated.tsv")
    agg_df.to_csv(agg_output_file, sep="\t", index=False)
    print("Aggregated distance calculations saved to:", agg_output_file)
else:
    print("No gene pairs found; aggregation skipped.")
