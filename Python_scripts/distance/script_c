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
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_subplot.tsv"
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_subplot.tsv"
#t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_moved.tsv"
#c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_moved.tsv"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/WT_translocation_distances_agg_subplot.tsv"
#output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/WT_translocation_distances_agg_moved.tsv"


# -------------------------------
# Load T1 and C1 gene data
# -------------------------------
df_t1 = pd.read_csv(t1_file, sep="\t")
df_c1 = pd.read_csv(c1_file, sep="\t")

# -------------------------------
# Collect all unique genes and their region types
# -------------------------------
gene_region_type = {}

for df in [df_t1, df_c1]:
    for gene, region in zip(df["gene1"], df["gene1_region_type"]):
        gene_region_type[gene] = region
    for gene, region in zip(df["gene2"], df["gene2_region_type"]):
        gene_region_type[gene] = region

all_genes = list(gene_region_type.keys())
print(f"Total unique genes: {len(all_genes)}")

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
        coords[i] = beads_chr.iloc[idx][["x","y","z"]].to_numpy()
        chr_mapped[i] = beads_chr.iloc[idx]["chr"]

    coords_dict = {name: coord for name, coord in zip(names, coords)}
    chr_dict = {name: chrom for name, chrom in zip(names, chr_mapped)}
    return coords_dict, chr_dict

# -------------------------------
# Build gene pairs to compute distances
# -------------------------------
gene_pairs = set()
for df in [df_t1, df_c1]:
    # inside_transloc ↔ neighbor
    pairs = df[((df["gene1_region_type"]=="inside_transloc") & (df["gene2_region_type"]=="neighbor")) |
               ((df["gene1_region_type"]=="neighbor") & (df["gene2_region_type"]=="inside_transloc"))]
    gene_pairs.update(zip(pairs["gene1"], pairs["gene2"]))
    # control ↔ neighbor
    pairs = df[((df["gene1_region_type"]=="control") & (df["gene2_region_type"]=="neighbor")) |
               ((df["gene1_region_type"]=="neighbor") & (df["gene2_region_type"]=="control"))]
    gene_pairs.update(zip(pairs["gene1"], pairs["gene2"]))

print(f"Total gene pairs to compute: {len(gene_pairs)}")

# -------------------------------
# Process all WT CMM files and store distances
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
print(f"Found {len(cmm_files)} CMM files to process")

distances_dict = defaultdict(list)
chr_dict_cache = {}

for f in cmm_files:
    print(f"\nProcessing model: {os.path.basename(f)}")
    beads = parse_cmm(f)
    if beads.empty:
        print("  No markers found, skipping")
        continue

    coords_dict, chr_dict = map_genes_to_beads(gtf_df, beads)
    if not chr_dict_cache:
        chr_dict_cache = chr_dict

    for gene1, gene2 in gene_pairs:
        if gene1 in coords_dict and gene2 in coords_dict:
            distance = np.linalg.norm(coords_dict[gene1] - coords_dict[gene2])
            distances_dict[(gene1, gene2)].append(distance)

# -------------------------------
# Compute mean, std, and assign region types
# -------------------------------
all_results = []
for (gene1, gene2), dist_list in distances_dict.items():
    all_results.append({
        "gene1": gene1,
        "gene1_region_type": gene_region_type.get(gene1, "unknown"),
        "chr1": chr_dict_cache.get(gene1, ""),
        "gene2": gene2,
        "gene2_region_type": gene_region_type.get(gene2, "unknown"),
        "chr2": chr_dict_cache.get(gene2, ""),
        "mean_distance": np.mean(dist_list),
        "std_distance": np.std(dist_list),
        "n_models": len(dist_list)
    })

# -------------------------------
# Save WT distances with region type
# -------------------------------
summary_df = pd.DataFrame(all_results)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
summary_df.to_csv(output_file, sep="\t", index=False)
print("WT distance calculations saved to:", output_file)
