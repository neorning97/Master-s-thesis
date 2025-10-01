#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
import sys

import pandas as pd
from pandas.api.types import is_integer_dtype


def existing_file(arg: str) -> pathlib.Path:
    if (path := pathlib.Path(arg)).is_file():
        return path

    raise argparse.ArgumentTypeError(f'Not an existing file: "{arg}"')


def make_cli():
    cli = argparse.ArgumentParser(
        "Preprocess domains in BED or BEDPE format.",
    )

    cli.add_argument(
        "file",
        type=existing_file,
        help="Path to a file in BED or BEDPE format.",
    )

    return cli


def validate_cols(df: pd.DataFrame):
    cols = [col for col in df.columns if not col.startswith("chrom")]
    invalid_cols = []
    for col in cols:
        dtype = df[col].dtype
        if not is_integer_dtype(dtype):
            invalid_cols.append(f"{col}: dtype={dtype}")

    if len(invalid_cols) != 0:
        invalid_cols = "\n - ".join(invalid_cols)
        raise RuntimeError(f"found the following invalid column(s):\n - {invalid_cols}")

    if len(cols) == 2:
        num_invalid_intervals = (df["start"] >= df["end"]).sum()
    else:
        num_invalid_intervals1 = (df["start1"] >= df["end1"]).sum()
        num_invalid_intervals2 = (df["start2"] >= df["end2"]).sum()

        num_invalid_intervals = (num_invalid_intervals1 | num_invalid_intervals2).sum()

    if num_invalid_intervals != 0:
        raise RuntimeError(
            f"found {num_invalid_intervals} (i.e., intervals whose start position is greater or equal to the end position)."
        )


def try_import_bedpe(path: pathlib.Path) -> pd.DataFrame:
    cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]

    df = pd.read_table(path, names=cols, usecols=list(range(len(cols))))

    df[["chrom1", "chrom2"]] = df[["chrom1", "chrom2"]].astype(str)

    validate_cols(df)

    return df


def try_import_bed3(path: pathlib.Path) -> pd.DataFrame:
    cols = ["chrom", "start", "end"]

    df = pd.read_table(path, names=cols, usecols=list(range(len(cols))))

    df["chrom"] = df["chrom"].astype(str)

    validate_cols(df)

    return df


def import_table(path: pathlib.Path) -> pd.DataFrame:
    try:
        return try_import_bedpe(path)
    except:  # noqa
        logging.warning('failed to parse "%s" as BEDPE. Falling back to BED', path)

    return try_import_bed3(path)


def setup_logger(level: str = "INFO"):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = make_cli().parse_args()

    df = import_table(args.file)

    domain_format = "BEDPE" if len(df.columns) == 6 else "BED"

    logging.info('successfully parsed %d %s domain(s) from file "%s"!', len(df), domain_format, args.file)

    try:
        df.to_csv(sys.stdout, sep="\t", index=False, header=False)
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    setup_logger()
    main()
