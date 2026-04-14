"""
01_extract_gene_expression.py
==============================
Script 1 of 4 for gene pair distance calculations – Gene expression extraction.

What this script does
---------------------
For each translocated condition (T1 and C1), this script:

1. Reads the BED files that define which genomic regions are translocated
   and which regions are their neighbors on the receiving chromosome.
2. Overlaps those regions with a GTF gene annotation to find all genes
   that fall inside each region.
3. Samples a matched set of control genes from chromosomes that are NOT
   involved in the translocation. The control group is the same size as
   the translocated gene group, making statistical comparisons fair.
4. Attaches RNA-seq TPM expression values to every gene.
5. Saves the result as a TSV file, one per condition.

Output of this script is used by script 02.

Usage
-----
    1. Edit the CONFIG section below to point to your files.
    2. Run: python 01_extract_gene_expression.py

Dependencies
------------
    pandas, bioframe
    Install with: pip install pandas bioframe

"""

import os
import pandas as pd
import bioframe


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    # GTF annotation file (Ensembl protein-coding genes, GRCh38)
    "gtf_file":  "/path/to/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf",

    # RNA-seq TPM expression table (GEO accession GSE246689)
    "expr_file": "/path/to/GSE246689_gene_tpm.tsv",

    # Per-condition BED files
    # transloc_bed : regions that were translocated
    # neighbor_bed : neighboring regions on the receiving chromosome
    "conditions": {
        "T1": {
            "transloc_bed": "/path/to/T1_translocations.bed",
            "neighbor_bed": "/path/to/T1_neighbor",
        },
        "C1": {
            "transloc_bed": "/path/to/C1_translocations.bed",
            "neighbor_bed": "/path/to/C1_neighbor.bed",
        },
    },

    # Where to save the output TSV files (created automatically)
    "results_dir": "/path/to/results",
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def load_bed(path: str) -> pd.DataFrame:
    """Read a BED file and ensure correct column types."""
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str).str.strip()
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype(int)
    return df


def load_gtf_genes(gtf_file: str) -> pd.DataFrame:
    """
    Load a GTF annotation file and return one row per gene with columns:
    chrom, start (0-based), end, gene_id, gene_name, gene_type.

    GTF files use 1-based inclusive coordinates. We convert the start
    position to 0-based (BED-style) so coordinates are consistent.
    """
    gtf = bioframe.read_table(gtf_file, schema="gtf")
    gtf = gtf[gtf["feature"] == "gene"].copy()

    # Convert 1-based GTF start to 0-based
    gtf["start"] -= 1
    gtf["chrom"]  = gtf["chrom"].astype(str).str.strip()

    # Parse the free-text attributes column to extract gene metadata
    gtf["gene_id"]   = gtf["attributes"].str.extract(r'gene_id "([^"]+)"')
    gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"')
    gtf["gene_type"] = gtf["attributes"].str.extract(r'gene_type "([^"]+)"')

    # Strip Ensembl version suffixes (e.g. ENSG00000001234.5 → ENSG00000001234)
    # so gene IDs match between the GTF and expression table
    gtf["gene_id"] = gtf["gene_id"].str.split(".").str[0]

    return gtf


def load_expression(expr_file: str) -> pd.DataFrame:
    """
    Load the RNA-seq TPM expression table and standardise the gene ID column
    so it can be merged with the GTF data.
    """
    df = pd.read_csv(expr_file, sep="\t")

    # Rename 'gene' → 'gene_id' if necessary
    if "gene" in df.columns:
        df.rename(columns={"gene": "gene_id"}, inplace=True)

    # Strip version suffixes to match GTF gene IDs
    df["gene_id"] = df["gene_id"].astype(str).str.split(".").str[0]

    # Drop gene_name column from expression table to avoid duplicate columns
    # when merging with the GTF (which already has gene_name)
    df.drop(columns=[c for c in ["gene_name"] if c in df.columns], inplace=True)

    return df


def genes_in_region(gtf: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Return all genes in the GTF that overlap the genomic region [start, end)
    on the given chromosome. A gene overlaps if any part of it falls within
    the region (not just the midpoint).
    """
    return gtf[
        (gtf["chrom"] == chrom) &
        (gtf["start"]  < end)   &
        (gtf["end"]    > start)
    ].copy()


def sample_control_genes(
    gtf: pd.DataFrame,
    excluded_chroms: set,
    n: int,
    seed: int,
) -> pd.DataFrame:
    """
    Draw n random genes from chromosomes NOT involved in the translocation.

    Why controls? We need a baseline to compare against. By sampling the same
    number of genes as the translocated group, from unrelated chromosomes, we
    can test whether any distance effects are specific to the translocation or
    just a general property of the genome.

    The seed is derived deterministically from the translocation ID so results
    are reproducible across runs.
    """
    pool = gtf[~gtf["chrom"].isin(excluded_chroms)]
    return pool.sample(n=n, random_state=seed).copy()


# =============================================================================
# Main function
# =============================================================================

def extract_gene_expression(
    transloc_bed: str,
    neighbor_bed: str,
    gtf_file: str,
    expr_file: str,
    output_file: str,
) -> None:
    """
    For one condition, find all genes in translocated and neighboring regions,
    add matched control genes, attach TPM expression values, and save to TSV.

    Parameters
    ----------
    transloc_bed  : BED file of translocated genomic regions.
    neighbor_bed  : BED file of neighboring regions on the receiving chr.
    gtf_file      : GTF annotation file.
    expr_file     : RNA-seq TPM expression table.
    output_file   : Where to save the output TSV.
    """
    transloc_df = load_bed(transloc_bed)
    neighbor_df = load_bed(neighbor_bed)
    gtf_df      = load_gtf_genes(gtf_file)
    expr_df     = load_expression(expr_file)

    all_genes = []

    # ── Genes inside translocated regions ────────────────────────────────────
    for _, row in transloc_df.iterrows():
        hits = genes_in_region(gtf_df, row["chrom"], row["start"], row["end"])
        if not hits.empty:
            hits["region_type"]         = "inside_transloc"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"]         = row["transloc_id"]
            all_genes.append(hits)

    # ── Genes inside neighboring regions ────────────────────────────────────
    # The neighbor BED uses 'neighbor_chr' instead of 'chrom' because the
    # neighbor sits on the receiving chromosome, not the origin chromosome
    for _, row in neighbor_df.iterrows():
        hits = genes_in_region(gtf_df, row["neighbor_chr"], row["start"], row["end"])
        if not hits.empty:
            hits["region_type"]         = "neighbor"
            hits["translocation_label"] = row["label"]
            hits["transloc_id"]         = row["transloc_id"]
            all_genes.append(hits)

    if not all_genes:
        print(f"  No genes found – skipping {output_file}")
        return

    all_genes_df = pd.concat(all_genes, ignore_index=True)

    # ── Control genes (matched by count, from non-translocated chromosomes) ──
    trans_chroms = set(transloc_df["chrom"])
    control_rows = []

    for tid in all_genes_df["transloc_id"].unique():
        n_trans = (
            all_genes_df
            .query("transloc_id == @tid and region_type == 'inside_transloc'")
            .shape[0]
        )
        if n_trans == 0:
            continue

        # Seed derived from translocation ID for reproducibility
        seed    = abs(hash(tid)) % (2**32)
        sampled = sample_control_genes(gtf_df, trans_chroms, n=n_trans, seed=seed)
        sampled["region_type"]         = "control"
        sampled["transloc_id"]         = tid
        sampled["translocation_label"] = "control"
        control_rows.append(sampled)

    if control_rows:
        all_genes_df = pd.concat([all_genes_df] + control_rows, ignore_index=True)

    # ── Attach TPM expression values ──────────────────────────────────────────
    merged = all_genes_df.merge(expr_df, on="gene_id", how="left")

    # Reorder: gene metadata first, then all MCF10A expression columns
    base_cols = [
        "chrom", "start", "end",
        "gene_id", "gene_name", "gene_type",
        "region_type", "translocation_label", "transloc_id",
        "source", "feature", "score", "strand", "frame",
    ]
    expr_cols = [c for c in merged.columns if c.startswith("MCF10A")]
    merged    = merged[base_cols + expr_cols]

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    merged.to_csv(output_file, sep="\t", index=False)
    print(f"  Saved: {output_file}")


# =============================================================================
# Run
# =============================================================================

if __name__ == "__main__":
    results_dir = CONFIG["results_dir"]
    os.makedirs(results_dir, exist_ok=True)

    for cond, paths in CONFIG["conditions"].items():
        print(f"\nProcessing condition: {cond}")
        extract_gene_expression(
            transloc_bed=paths["transloc_bed"],
            neighbor_bed=paths["neighbor_bed"],
            gtf_file=CONFIG["gtf_file"],
            expr_file=CONFIG["expr_file"],
            output_file=os.path.join(results_dir, f"{cond}_genes_expression.tsv"),
        )

    print("\nDone. Output files:")
    for cond in CONFIG["conditions"]:
        print(f"  {results_dir}/{cond}_genes_expression.tsv")
