# Script that calculates the spatial distance between genes, in the neighbor and translocation regions
# found in T1 and C1, for WT.

import pandas as pd
import numpy as np
import glob
import os
import re
from collections import defaultdict

# -------------------------------
# User paths
# -------------------------------
cmm_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/WT"
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/T1_translocation_neighbor_distances_aggregated_filtered.tsv"
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/C1_translocation_neighbor_distances_aggregated_filtered.tsv"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_translocation_neighbor_distances_aggregated.tsv"

# -------------------------------
# Load T1 and C1 gene pairs
# -------------------------------
df_t1 = pd.read_csv(t1_file, sep="\t")
df_c1 = pd.read_csv(c1_file, sep="\t")

# Get unique gene pairs (as tuples)
t1_pairs = set(zip(df_t1["gene1"], df_t1["gene2"]))
c1_pairs = set(zip(df_c1["gene1"], df_c1["gene2"]))
gene_pairs = t1_pairs.union(c1_pairs)

# Get all unique genes
all_genes = set(df_t1["gene1"]).union(df_t1["gene2"]).union(df_c1["gene1"]).union(df_c1["gene2"])

# -------------------------------
# Parse GTF to get WT gene coordinates
# -------------------------------
gtf_cols = ["chrom","source","feature","start","end","score","strand","frame","attribute"]
gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=gtf_cols)

gtf_df = gtf_df[gtf_df["feature"]=="gene"].copy()
gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')
gtf_df = gtf_df[gtf_df["gene_name"].isin(all_genes)].copy()

gtf_df["chrom"] = gtf_df["chrom"].astype(str).apply(lambda x: x if x.startswith("chr") else "chr"+x)
gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")

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

    coords_dict = {name: coord for name, coord in zip(names, coords)}
    chr_dict = {name: chrom for name, chrom in zip(names, chr_mapped)}
    return coords_dict, chr_dict

# -------------------------------
# Process all WT CMM files and store distances
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
print(f"Found {len(cmm_files)} CMM files to process")

# Use defaultdict(list) to accumulate distances for each gene pair
distances_dict = defaultdict(list)
chr_dict_cache = {}

for f in cmm_files:
    print(f"\nProcessing model: {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No markers found in this model, skipping")
        continue

    coords_dict, chr_dict = map_genes_to_beads(gtf_df, beads)
    print(f"Mapped {len(coords_dict)} genes to beads in WT model")

    # Cache chromosome info (all models assumed to have same mapping)
    if not chr_dict_cache:
        chr_dict_cache = chr_dict

    for gene1, gene2 in gene_pairs:
        if gene1 in coords_dict and gene2 in coords_dict:
            distance = np.linalg.norm(coords_dict[gene1] - coords_dict[gene2])
            distances_dict[(gene1, gene2)].append(distance)

# -------------------------------
# Compute mean and std distance for each pair
# -------------------------------
all_results = []
for (gene1, gene2), dist_list in distances_dict.items():
    all_results.append({
        "gene1": gene1,
        "chr1": chr_dict_cache.get(gene1, ""),
        "gene2": gene2,
        "chr2": chr_dict_cache.get(gene2, ""),
        "mean_distance": np.mean(dist_list),
        "std_distance": np.std(dist_list),
        "n_models": len(dist_list)
    })

# -------------------------------
# Save results
# -------------------------------
summary_df = pd.DataFrame(all_results)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
summary_df.to_csv(output_file, sep="\t", index=False)
print("WT distance calculations saved to:", output_file)
