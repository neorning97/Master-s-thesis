#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
from typing import List, Optional, Set, Tuple

import hictkpy
import pandas as pd


def existing_file(arg: str) -> pathlib.Path:
    if (path := pathlib.Path(arg)).is_file():
        return path

    raise argparse.ArgumentTypeError(f'Not an existing file: "{arg}"')


def positive_int(arg) -> int:
    if (n := int(arg)) > 0:
        return n

    raise argparse.ArgumentTypeError("Not a positive int")


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "sample",
        type=str,
        help="Sample name.",
    )
    cli.add_argument(
        "hic-file",
        type=existing_file,
        help="Path to a Hi-C matrix in .hic or .[m]cool format.",
    )
    cli.add_argument(
        "--resolution",
        type=positive_int,
        help="Matrix resolution (only required when processing multi-resoltion files).",
    )
    cli.add_argument(
        "--domains",
        type=existing_file,
        help="Path to a BEDPE file with the list of domains to be processed.",
    )
    cli.add_argument(
        "--interaction-type",
        type=str,
        choices={"cis", "trans"},
        help="Interaction type.",
    )

    return cli


def import_chrom_pairs_from_bedpe(
    path: Optional[pathlib.Path],
    interaction_type: Optional[str],
) -> Optional[Set[Tuple[str, str]]]:
    if path is None:
        return None

    # We parse all 6 BEDPE columns even though we only need the 1st and 3rd to catch malformed files with too few columns
    cols = ["chrom1", "start", "end1", "chrom2", "start2", "end2"]
    df = pd.read_table(path, names=cols, usecols=list(range(len(cols))))[["chrom1", "chrom2"]]

    if interaction_type == "cis":
        df = df[df["chrom1"] == df["chrom2"]]
    elif interaction_type == "trans":
        df = df[df["chrom1"] != df["chrom2"]]

    df = df.astype(str).drop_duplicates()

    return set((c1, c2) for c1, c2 in df.itertuples(index=False))


def chrom_pair_has_interactions(f: hictkpy.File, chrom1: str, chrom2: str) -> bool:
    return any(True for _ in f.fetch(chrom1, chrom2))


def filter_chrom_pairs(
    f: hictkpy.File,
    chrom_pairs: List[Tuple[str, str]],
    bedpe: Optional[pathlib.Path],
    interaction_type: Optional[str],
) -> List[Tuple[str, str]]:
    """
    Return the intersection between the given chromosome pairs and the chromosome pairs imported from the BEDPE file.
    Skip chromosome pairs without interactions.
    """
    unique_pairs = set(chrom_pairs)
    if bedpe is not None:
        unique_pairs &= import_chrom_pairs_from_bedpe(bedpe, interaction_type)

    pairs = []
    for chrom1, chrom2 in unique_pairs:
        if chrom_pair_has_interactions(f, chrom1, chrom2):
            pairs.append((chrom1, chrom2))

    return list(sorted(pairs))


def make_chrom_pairs_cis(chroms: List[str]) -> List[Tuple[str, str]]:
    return [(chrom, chrom) for chrom in chroms]


def make_chrom_pairs_trans(chroms: List[str]) -> List[Tuple[str, str]]:
    pairs = []
    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i + 1 :]:
            pairs.append((chrom1, chrom2))

    return pairs


def import_chrom_pairs_from_matrix(
    hf: hictkpy.File,
    interaction_type: Optional[str],
) -> List[Tuple[str, str]]:
    """
    Import chromosome pairs from a Hi-C matrix.
    """
    chroms = list(hf.chromosomes().keys())

    if interaction_type == "trans":
        return make_chrom_pairs_trans(chroms)

    if interaction_type == "cis":
        return make_chrom_pairs_cis(chroms)

    assert interaction_type is None
    return make_chrom_pairs_cis(chroms) + make_chrom_pairs_trans(chroms)


def main():
    args = vars(make_cli().parse_args())

    sample = args["sample"]
    interaction_type = args["interaction_type"]

    hf = hictkpy.File(args["hic-file"], args["resolution"])

    chrom_pairs = import_chrom_pairs_from_matrix(
        hf,
        interaction_type,
    )

    chrom_pairs = filter_chrom_pairs(
        hf,
        chrom_pairs,
        bedpe=args["domains"],
        interaction_type=interaction_type,
    )

    for chrom1, chrom2 in chrom_pairs:
        print(f"{sample}\t{chrom1}\t{chrom2}")


if __name__ == "__main__":
    main()
