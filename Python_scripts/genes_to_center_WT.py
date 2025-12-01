# Script that calculates the spatial distance between genes, found in the translocated regions of T1 and C1, in WT,
# and the center of the nucleus (0,0,0).

import pandas as pd
import numpy as np
import os
import glob
import re

# -------------------------------
# User paths
# -------------------------------
cmm_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/WT"  # WT CMM files
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
t1_gene_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.tsv"
c1_gene_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/C1_genes_expression.tsv"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_T1C1_genes_center_distances.tsv"

# -------------------------------
# Load T1 and C1 translocated genes
# -------------------------------
df_t1 = pd.read_csv(t1_gene_file, sep="\t")
df_c1 = pd.read_csv(c1_gene_file, sep="\t")
df_t1_trans = df_t1[df_t1["region_type"]=="inside_transloc"]
df_c1_trans = df_c1[df_c1["region_type"]=="inside_transloc"]

# Combine genes and ensure uniqueness
genes_of_interest = pd.concat([df_t1_trans["gene_name"], df_c1_trans["gene_name"]]).unique()
print(f"Number of unique T1/C1 translocated genes: {len(genes_of_interest)}")

# -------------------------------
# Load GTF and extract coordinates
# -------------------------------
gtf_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=gtf_cols)
gtf = gtf[gtf["feature"]=="gene"].copy()
gtf["gene_name"] = gtf["attributes"].str.extract('gene_name "([^"]+)"')
gtf_filtered = gtf[gtf["gene_name"].isin(genes_of_interest)].copy()
gtf_filtered["chrom"] = gtf_filtered["chrom"].astype(str).apply(lambda x: x if x.startswith("chr") else "chr"+x)
print(f"Number of genes found in GTF: {len(gtf_filtered)}")

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
        coords[i] = beads_df.iloc[idx][["x","y","z"]].to_numpy()
    return coords, names

# -------------------------------
# Process all WT CMM files
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
all_results = []

for f in cmm_files:
    print(f"Processing {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No beads found, skipping")
        continue
    gene_coords, gene_names = map_genes_to_beads(gtf_filtered, beads)
    distances = np.linalg.norm(gene_coords, axis=1)
    for name, d in zip(gene_names, distances):
        all_results.append({
            "gene_name": name,
            "distance_to_center": d
        })

# -------------------------------
# Save per-model results
# -------------------------------
result_df = pd.DataFrame(all_results)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
result_df.to_csv(output_file, sep="\t", index=False)
print(f"Distances saved to {output_file}")

# -------------------------------
# Aggregate across models
# -------------------------------
agg_df = result_df.groupby("gene_name").agg(
    mean_distance=("distance_to_center","mean"),
    std_distance=("distance_to_center","std"),
    n_models=("distance_to_center","count")
).reset_index()
agg_output_file = output_file.replace(".tsv","_aggregated.tsv")
agg_df.to_csv(agg_output_file, sep="\t", index=False)
print(f"Aggregated distances saved to {agg_output_file}")
