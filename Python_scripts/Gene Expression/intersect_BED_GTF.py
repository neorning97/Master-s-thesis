# CODE THAT IDENTIFIES GENES LOCATED NEAR GENOMIC TRANSLOCATION BREAKPOINTS (NEIGHBOR), GENES IN TRANSLOCATIONS, AND COMBINES THEM WITH GENE EXPRESSION DATA.
#   - GENES OUTSIDE THE TRANSLOCATED SEGMENT (NEIGHBOR) AND GENES INSIDE THE TRANSLOCATED SEGMENT (TRANSLOC) ARE CONSIDERED
#   - NEIGHBORHOODS RESPECT THE BREAKPOINT SIDE (STARTBP VS ENDBP)
#   - NEIGHBOR and TRANSLOCATE GENES ARE MERGED WITH GENE EXPRESSION DATA (TPM) FOR DOWNSTREAM ANALYSIS

import pandas as pd
import bioframe

# -------------------------------
# File paths
# -------------------------------
transloc_bed_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/T1_translocations.bed"
neighbor_bed_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/T1_neighbor.bed"
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
expr_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/GSE246689_gene_tpm.tsv"

output_tsv = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.tsv"
output_parquet = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.parquet"

# -------------------------------
# 1. Load BED files
# -------------------------------
transloc_df = pd.read_csv(transloc_bed_file, sep="\t")
neighbor_df = pd.read_csv(neighbor_bed_file, sep="\t")

# Clean formatting
for df in [transloc_df, neighbor_df]:
    df['chrom'] = df['chrom'].astype(str).str.strip()
    df['start'] = pd.to_numeric(df['start'], errors='coerce').astype(int)
    df['end'] = pd.to_numeric(df['end'], errors='coerce').astype(int)

print(f"Loaded {len(transloc_df)} translocated regions and {len(neighbor_df)} neighbor regions.")

# -------------------------------
# 2. Load GTF file
# -------------------------------
gtf_df = bioframe.read_table(gtf_file, schema="gtf")
gtf_df = gtf_df[gtf_df['feature'] == 'gene'].copy()

# Convert GTF to BED-like 0-based
gtf_df['start'] = gtf_df['start'] - 1
gtf_df['chrom'] = gtf_df['chrom'].astype(str).str.strip()

# Extract attributes
gtf_df['gene_id'] = gtf_df['attributes'].str.extract('gene_id "([^"]+)"')
gtf_df['gene_name'] = gtf_df['attributes'].str.extract('gene_name "([^"]+)"')
gtf_df['gene_type'] = gtf_df['attributes'].str.extract('gene_type "([^"]+)"')

print(f"Loaded {len(gtf_df)} genes from GTF.")

# -------------------------------
# 3. Load expression data
# -------------------------------
expr_df = pd.read_csv(expr_file, sep="\t")

# Ensure gene_id compatibility
if "gene_id" not in expr_df.columns:
    # Common GEO formats use 'gene' or 'GeneName'
    if "gene" in expr_df.columns:
        expr_df.rename(columns={"gene": "gene_id"}, inplace=True)
    elif "GeneName" in expr_df.columns:
        expr_df.rename(columns={"GeneName": "gene_id"}, inplace=True)

print(f"Loaded expression table with {expr_df.shape[0]} genes and {expr_df.shape[1]-1} samples.")

# Remove version suffix (ENSG00...1.12 → ENSG001...)

expr_df["gene_id"] = expr_df["gene_id"].astype(str).str.split(".").str[0].str.strip()
gtf_df["gene_id"] = gtf_df["gene_id"].astype(str).str.split(".").str[0].str.strip()

# -------------------------------
# 4. Find genes overlapping regions
# -------------------------------
all_genes_list = []

# --- Translocated regions ---
for _, row in transloc_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    label = row['label']
    transloc_id = row['transloc_id']

    inside_genes = gtf_df[
        (gtf_df['chrom'] == chrom) &
        (gtf_df['start'] < end) &
        (gtf_df['end'] > start)
    ].copy()

    if not inside_genes.empty:
        inside_genes['region_type'] = 'inside_transloc'
        inside_genes['translocation_label'] = label
        inside_genes['transloc_id'] = transloc_id
        all_genes_list.append(inside_genes)

# --- Neighbor regions ---
for _, row in neighbor_df.iterrows():
    chrom = row['neighbor_chr']
    start = row['start']
    end = row['end']
    label = row['label']
    transloc_id = row['transloc_id']  # FIXED: Preserve ID

    neighbor_genes = gtf_df[
        (gtf_df['chrom'] == chrom) &
        (gtf_df['start'] < end) &
        (gtf_df['end'] > start)
    ].copy()

    if not neighbor_genes.empty:
        neighbor_genes['region_type'] = 'neighbor'
        neighbor_genes['translocation_label'] = label
        neighbor_genes['transloc_id'] = transloc_id
        all_genes_list.append(neighbor_genes)

# Combine all genes
if all_genes_list:
    all_genes_df = pd.concat(all_genes_list, ignore_index=True)
else:
    all_genes_df = pd.DataFrame()
    print("No genes found in translocated or neighbor regions!")

# -------------------------------
# 5. Merge with expression data
# -------------------------------
if not all_genes_df.empty:

    # Drop 'gene_name' in expression table to avoid _x/_y after merge
    if 'gene_name' in expr_df.columns:
        expr_df = expr_df.drop(columns=['gene_name'])

    merged_df = all_genes_df.merge(expr_df, on="gene_id", how="left")

    # Remove unnecessary columns
    drop_cols = ['attributes', 'chrom_', 'start_', 'end_']
    merged_df = merged_df.drop(columns=[c for c in drop_cols if c in merged_df.columns])

    # Column ordering
    cols = [
        'chrom', 'start', 'end',
        'gene_id', 'gene_name', 'gene_type',
        'region_type', 'translocation_label', 'transloc_id'
    ]
    cols += [c for c in merged_df.columns if c not in cols]

    merged_df = merged_df[[c for c in cols if c in merged_df.columns]]

    # Save
    merged_df.to_csv(output_tsv, sep="\t", index=False)

    try:
        merged_df.to_parquet(output_parquet, index=False)
        parquet_saved = True
    except ImportError:
        parquet_saved = False

    print(f"Genes in translocated and neighbor regions saved to TSV: {output_tsv}")
    if parquet_saved:
        print(f"Parquet saved: {output_parquet}")

else:
    print("No genes to merge or save.")
