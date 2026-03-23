import pandas as pd
import numpy as np
import os
import glob
import re

# -------------------------------
# Paths per condition
# -------------------------------
conditions = {
    "T1": {
        "cmm_dir": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/T1",
        "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_genes_expression.tsv",
        "output": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_distance_to_center.tsv",
        "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_distance_to_center_agg.tsv",
        #"gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_genes_expression_moved.tsv",
        #"output": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_distance_to_center_moved.tsv",
        #"output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_distance_to_center_agg_moved.tsv",
    },
    "C1": {
        "cmm_dir": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/changed_colors/C1",
        "gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_genes_expression.tsv",
        "output": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_distance_to_center.tsv",
        "output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_distance_to_center_agg.tsv",
        #"gene_expression_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_genes_expression_moved.tsv",
        #"output": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_distance_to_center_moved.tsv",
        #"output_agg": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_distance_to_center_agg_moved.tsv",
    }
}

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
                chrom_part, coords = attrs["beadID"].split(":")
                if not chrom_part.startswith("chr"):
                    chrom_part = "chr" + chrom_part
                start, end = map(int, coords.split("-"))
                x = float(attrs.get("x", "nan"))
                y = float(attrs.get("y", "nan"))
                z = float(attrs.get("z", "nan"))
                markers.append([chrom_part, start, end, x, y, z])
            except:
                continue
    return pd.DataFrame(markers, columns=["chr","start","end","x","y","z"])

# -------------------------------
# Map genes to beads
# -------------------------------
def map_genes_to_beads(genes_df, beads_df):
    n = len(genes_df)
    coords = np.full((n, 3), np.nan)

    for i, (_, row) in enumerate(genes_df.iterrows()):
        beads_chr = beads_df[beads_df["chr"] == row["chrom"]]
        if beads_chr.empty:
            continue

        gene_mid = (row["start"] + row["end"]) / 2
        bead_mid = (beads_chr["start"] + beads_chr["end"]) / 2

        idx = np.abs(bead_mid - gene_mid).to_numpy().argmin()
        coords[i] = beads_chr.iloc[idx][["x","y","z"]].to_numpy()

    return coords

# -------------------------------
# Main loop
# -------------------------------
for cond, paths in conditions.items():
    print(f"\n===== {cond} =====")

    genes = pd.read_csv(paths["gene_expression_file"], sep="\t")
    genes["chrom"] = genes["chrom"].astype(str)
    genes["chrom"] = genes["chrom"].apply(lambda x: x if x.startswith("chr") else "chr"+x)

    # Only translocated + control genes
    genes = genes[genes["region_type"].isin(["inside_transloc","control"])].copy()

    cmm_files = sorted(glob.glob(os.path.join(paths["cmm_dir"], "*.cmm")))
    all_results = []

    for f in cmm_files:
        print("Processing", os.path.basename(f))
        beads = parse_cmm(f)
        if beads.empty:
            continue

        coords = map_genes_to_beads(genes, beads)

        # Distance to nuclear center
        dists = np.linalg.norm(coords, axis=1)

        tmp = genes[[
            "gene_id","gene_name","chrom",
            "region_type","transloc_id"
        ]].copy()

        tmp["distance_to_center"] = dists
        tmp["model"] = os.path.basename(f)
        all_results.append(tmp)

    df = pd.concat(all_results, ignore_index=True)

    os.makedirs(os.path.dirname(paths["output"]), exist_ok=True)
    df.to_csv(paths["output"], sep="\t", index=False)

    # Aggregate across models
    df_agg = (
        df.groupby(
            ["gene_name","chrom","region_type","transloc_id"],
            as_index=False
        )["distance_to_center"]
        .agg(mean_distance="mean", std_distance="std", n_models="count")
    )

    df_agg.to_csv(paths["output_agg"], sep="\t", index=False)
    print("Saved:", paths["output_agg"])
