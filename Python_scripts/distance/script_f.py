import pandas as pd
import numpy as np
import os
import glob
import re
from collections import defaultdict

# -------------------------------
# Paths
# -------------------------------
cmm_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/WT"
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

#t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_genes_expression.tsv"
#c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_genes_expression.tsv"
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_genes_expression_moved.tsv"
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_genes_expression_moved.tsv"

#output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_WT_distance_to_center_agg.tsv"
output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_WT_distance_to_center_agg_moved.tsv"

# -------------------------------
# Load gene sets and region types
# -------------------------------
df_t1 = pd.read_csv(t1_file, sep="\t")
df_c1 = pd.read_csv(c1_file, sep="\t")

gene_region = {}
for df in [df_t1, df_c1]:
    for _, row in df.iterrows():
        gene_region[row["gene_name"]] = row["region_type"]

genes_of_interest = list(gene_region.keys())
print(f"Total genes for WT mapping: {len(genes_of_interest)}")

# -------------------------------
# Parse GTF
# -------------------------------
gtf_cols = ["chrom","source","feature","start","end","score","strand","frame","attribute"]
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=gtf_cols)

gtf = gtf[gtf["feature"] == "gene"].copy()
gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "([^"]+)"')
gtf = gtf[gtf["gene_name"].isin(genes_of_interest)].copy()

gtf["chrom"] = gtf["chrom"].astype(str)
gtf["chrom"] = gtf["chrom"].apply(lambda x: x if x.startswith("chr") else "chr"+x)
gtf["start"] = pd.to_numeric(gtf["start"])
gtf["end"] = pd.to_numeric(gtf["end"])

# -------------------------------
# CMM parsing
# -------------------------------
def parse_cmm(file_path):
    markers = []
    with open(file_path) as fh:
        for line in fh:
            if "<marker" not in line:
                continue
            attrs = dict(re.findall(r'(\w+)="([^"]+)"', line))
            if "beadID" not in attrs:
                continue
            try:
                chrom, coords = attrs["beadID"].split(":")
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                start, end = map(int, coords.split("-"))
                x = float(attrs.get("x","nan"))
                y = float(attrs.get("y","nan"))
                z = float(attrs.get("z","nan"))
                markers.append([chrom,start,end,x,y,z])
            except:
                continue
    return pd.DataFrame(markers, columns=["chr","start","end","x","y","z"])

# -------------------------------
# Map genes to beads
# -------------------------------
def map_genes_to_beads(genes_df, beads_df):
    coords = {}
    for _, row in genes_df.iterrows():
        beads_chr = beads_df[beads_df["chr"] == row["chrom"]]
        if beads_chr.empty:
            continue
        gene_mid = (row["start"] + row["end"]) / 2
        bead_mid = (beads_chr["start"] + beads_chr["end"]) / 2
        idx = np.abs(bead_mid - gene_mid).argmin()
        coords[row["gene_name"]] = beads_chr.iloc[idx][["x","y","z"]].to_numpy()
    return coords

# -------------------------------
# Main loop
# -------------------------------
cmm_files = sorted(glob.glob(os.path.join(cmm_dir,"*.cmm")))
distances = defaultdict(list)

for f in cmm_files:
    print("Processing", os.path.basename(f))
    beads = parse_cmm(f)
    if beads.empty:
        continue

    coords = map_genes_to_beads(gtf, beads)

    for gene, xyz in coords.items():
        distances[gene].append(np.linalg.norm(xyz))

# -------------------------------
# Aggregate
# -------------------------------
rows = []
for gene, dists in distances.items():
    rows.append({
        "gene_name": gene,
        "region_type": gene_region.get(gene, "unknown"),
        "mean_WT": np.mean(dists),
        "std_WT": np.std(dists),
        "n_models": len(dists)
    })

wt_df = pd.DataFrame(rows)
wt_df.to_csv(output_file, sep="\t", index=False)

print("WT distance-to-center saved to:", output_file)
