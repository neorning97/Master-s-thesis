#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import string

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "gtrack",
        type=existing_file,
        help="Path to a .gtrack file with the list of beads for a single chromosome copy.",
    )
    cli.add_argument(
        "ploidy",
        type=int,
        help="Target number of chromosome copies.",
    )

    return cli


def rename_edges(edges: str, suffix: str) -> str:
    if edges == ".":
        return edges

    new_edges = []
    for e in edges.split(";"):
        chrom, _, coords = e.partition(":")
        new_edges.append(f"{chrom}{suffix}:{coords}")

    return ";".join(new_edges)


def main():
    args = vars(make_cli().parse_args())

    beads = pd.read_table(
        args["gtrack"], comment="#", names=["chrom", "start", "end", "tid", "radius", "periphery", "edges"]
    )

    print("##gtrack version: 1.0")
    print("##track type: linked segments")
    print("###seqid\tstart\tend\tid\tradius\tperiphery\tedges")

    ploidy = args["ploidy"]
    charset = string.ascii_uppercase
    for _, row in beads.iterrows():
        for i in range(ploidy):
            suffix = f"_{charset[i]}"

            row_ = row.copy()

            row_["chrom"] += suffix
            row_["tid"] = row["tid"].replace(row["chrom"], row_["chrom"])
            row_["edges"] = rename_edges(row["edges"], suffix)

            print("\t".join((str(x) for x in row_)))


if __name__ == "__main__":
    main()
