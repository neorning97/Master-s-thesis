import pandas as pd
import bioframe
import os
import re

# -------------------------------
# Define conditions
# -------------------------------
conditions = {
    "T1": {
        "transloc_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/subplot_T1_transloc.bed",
        "neighbor_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/T1_neighbor1.bed",
    #    "neighbor_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/T1_neighbor_originchr1.bed",
        "expr_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/GSE246689_gene_tpm.tsv",
        "output_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_genes_expression_subplot.tsv"
    #    "output_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_genes_expression_moved.tsv"
    },
    "C1": {
        "transloc_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/subplot_C1_transloc.bed",
        "neighbor_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/C1_neighbor1.bed",
    #    "neighbor_bed": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/C1_neighbor_originchr1.bed",
        "expr_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/GSE246689_gene_tpm.tsv",
        "output_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_genes_expression_subplot.tsv"
    #    "output_file": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_genes_expression_moved.tsv"
    }
}

gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"

# -------------------------------
# Function to process one condition
# -------------------------------
def process_condition(transloc_bed, neighbor_bed, gtf_file, expr_file, output_file):

    # --- Load BED files ---
    transloc_df = pd.read_csv(transloc_bed, sep="\t")
    neighbor_df = pd.read_csv(neighbor_bed, sep="\t")

    for df in (transloc_df, neighbor_df):
        df["chrom"] = df["chrom"].astype(str).str.strip()
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)

    # --- Load GTF ---
    gtf_df = bioframe.read_table(gtf_file, schema="gtf")
    gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()
    gtf_df["start"] -= 1
    gtf_df["chrom"] = gtf_df["chrom"].astype(str).str.strip()

    gtf_df["gene_id"] = gtf_df["attributes"].str.extract('gene_id "([^"]+)"')
    gtf_df["gene_name"] = gtf_df["attributes"].str.extract('gene_name "([^"]+)"')
    gtf_df["gene_type"] = gtf_df["attributes"].str.extract('gene_type "([^"]+)"')

    gtf_df["gene_id"] = gtf_df["gene_id"].str.split(".").str[0]

    # --- Load expression ---
    expr_df = pd.read_csv(expr_file, sep="\t")
    if "gene" in expr_df.columns:
        expr_df.rename(columns={"gene": "gene_id"}, inplace=True)

    expr_df["gene_id"] = expr_df["gene_id"].astype(str).str.split(".").str[0]

    # IMPORTANT: remove gene_name from expression to avoid gene_name_x/y
    expr_df = expr_df.drop(columns=[c for c in ["gene_name"] if c in expr_df.columns])

    # --- FIX: Ensure replicate columns are separate ---
    # Example pattern: MCF10AC1_REP1MCF10AC1_REP3 -> split into two
    new_cols = []
    for c in expr_df.columns:
        if re.match(r'(MCF10\w+_REP\d+)(MCF10\w+_REP\d+)', c):
            # split concatenated columns
            parts = re.findall(r'(MCF10\w+_REP\d+)', c)
            new_cols.extend(parts)
        else:
            new_cols.append(c)

    # If the number of new columns differs from current, adjust
    if len(new_cols) != len(expr_df.columns):
        expr_df.columns = new_cols

    # --- Find inside_transloc and neighbor genes ---
    all_genes = []

    for _, row in transloc_df.iterrows():
        genes = gtf_df[
            (gtf_df["chrom"] == row["chrom"]) &
            (gtf_df["start"] < row["end"]) &
            (gtf_df["end"] > row["start"])
        ].copy()

        if not genes.empty:
            genes["region_type"] = "inside_transloc"
            genes["translocation_label"] = row["label"]
            genes["transloc_id"] = row["transloc_id"]
            all_genes.append(genes)


    for _, row in neighbor_df.iterrows():
        genes = gtf_df[
            (gtf_df["chrom"] == row["neighbor_chr"]) &
            (gtf_df["start"] < row["end"]) &
            (gtf_df["end"] > row["start"])
        ].copy()

        if not genes.empty:
            genes["region_type"] = "neighbor"
            genes["translocation_label"] = row["label"]
            genes["transloc_id"] = row["transloc_id"]
            all_genes.append(genes)

    if not all_genes:
        print("No genes found")
        return

    all_genes_df = pd.concat(all_genes, ignore_index=True)

    # --- Control genes ---
    trans_chroms = set(transloc_df["chrom"])
    control_pool = gtf_df[~gtf_df["chrom"].isin(trans_chroms)].copy()

    control_genes = []
    for tid in all_genes_df["transloc_id"].unique():
        n = all_genes_df.query(
            "transloc_id == @tid and region_type == 'inside_transloc'"
        ).shape[0]
        if n == 0:
            continue

        sampled = control_pool.sample(
            n=n, random_state=abs(hash(tid)) % (2**32)
        ).copy()

        sampled["region_type"] = "control"
        sampled["transloc_id"] = tid
        sampled["translocation_label"] = "control"
        control_genes.append(sampled)

    if control_genes:
        all_genes_df = pd.concat([all_genes_df] + control_genes, ignore_index=True)

    # --- Merge with expression ---
    merged_df = all_genes_df.merge(expr_df, on="gene_id", how="left")

    # --- Column order to match old file ---
    base_cols = [
        "chrom", "start", "end",
        "gene_id", "gene_name", "gene_type",
        "region_type", "translocation_label", "transloc_id",
        "source", "feature", "score", "strand", "frame"
    ]

    expr_cols = [c for c in merged_df.columns if c.startswith("MCF10A")]

    merged_df = merged_df[base_cols + expr_cols]

    # --- Save ---
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"Saved {output_file}")

# -------------------------------
# Run all conditions
# -------------------------------
for cond, paths in conditions.items():
    print(f"\nProcessing {cond}")
    process_condition(
        gtf_file=gtf_file,
        **paths
    )

