# Script that calculates the spatial distance between the genes in the translocated regions and the center
# of the nucleus (0,0,0).

import pandas as pd
import numpy as np
import os
import glob
import re

# -------------------------------
# User paths
# -------------------------------
cmm_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/T1"  # or T1, WT
gene_expression_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.tsv"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/T1_transloc_gene_center_distances.tsv"

# -------------------------------
# Load gene file and filter
# -------------------------------
df_genes = pd.read_csv(gene_expression_file, sep="\t")
df_transloc_genes = df_genes[df_genes["region_type"]=="inside_transloc"].copy()
print(f"Number of translocated genes: {len(df_transloc_genes)}")

# -------------------------------
# Parse CMM files to get 3D coordinates
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
                    chrom_part = "chr"+chrom_part
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
    names = genes_df["gene_name"].to_numpy(str)
    mids = ((genes_df["start"] + genes_df["end"])/2).to_numpy()
    coords = np.full((len(names),3), np.nan)
    for i,(chrom,mid) in enumerate(zip(genes_df["chrom"], mids)):
        beads_chr = beads_df[beads_df["chr"]==chrom]
        if beads_chr.empty:
            continue
        bead_mids = ((beads_chr["start"] + beads_chr["end"])/2).to_numpy()
        idx = np.abs(bead_mids - mid).argmin()
        coords[i] = beads_chr.iloc[idx][["x","y","z"]].to_numpy()
    return coords, names

# -------------------------------
# Process all CMM files
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
all_results = []

for f in cmm_files:
    print(f"Processing {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No beads found, skipping")
        continue
    trans_coords, trans_names = map_genes_to_beads(df_transloc_genes, beads)
    # Distance to nuclear center (0,0,0)
    distances = np.linalg.norm(trans_coords, axis=1)
    for name, d in zip(trans_names, distances):
        all_results.append({
            "gene_name": name,
            "distance_to_center": d
        })

# -------------------------------
# Save results
# -------------------------------
result_df = pd.DataFrame(all_results)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
result_df.to_csv(output_file, sep="\t", index=False)
print(f"Distances saved to {output_file}")

# Optional: aggregate across models
agg_df = result_df.groupby("gene_name").agg(
    mean_distance=("distance_to_center","mean"),
    std_distance=("distance_to_center","std"),
    n_models=("distance_to_center","count")
).reset_index()
agg_output_file = output_file.replace(".tsv","_aggregated.tsv")
agg_df.to_csv(agg_output_file, sep="\t", index=False)
print(f"Aggregated distances saved to {agg_output_file}")
