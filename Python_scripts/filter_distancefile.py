# Script that makes sure that the gene pairs in the file are not on the same chromosome
# --> that way we know that the gene pairs are only between neighbor and translocation regions.

import pandas as pd

# List of input files

input_files = [
    "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_translocation_neighbor_distances.tsv",
    "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_translocation_neighbor_distances_aggregated.tsv"
]

for input_file in input_files:
    output_file = input_file.replace(".tsv", "_filtered.tsv")
    
    # Load the data
    df = pd.read_csv(input_file, sep="\t")
    print(f"Loaded {len(df)} rows from {input_file}")
    
    # Keep only rows where chr1 != chr2
    df_filtered = df[df["chr1"] != df["chr2"]].copy()
    print(f"{len(df_filtered)} rows remain after filtering out same-chromosome pairs")
    
    # Save filtered data
    df_filtered.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered distances saved to {output_file}\n")
