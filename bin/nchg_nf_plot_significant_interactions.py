#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import logging
import multiprocessing as mp
import pathlib
from typing import Any, Dict, List, Optional, Tuple

import h5py
import hictkpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow
from matplotlib.colors import LogNorm
from numpy.typing import NDArray


def existing_file(arg: str) -> pathlib.Path:
    if (path := pathlib.Path(arg)).is_file():
        return path

    raise argparse.ArgumentTypeError(f'Not an existing file: "{arg}"')


def probability(arg) -> float:
    if 0 <= (n := float(arg)) <= 1:
        return n

    raise argparse.ArgumentTypeError("Not a valid probability")


def positive_int(arg) -> int:
    if (n := int(arg)) > 0:
        return n

    raise argparse.ArgumentTypeError("Not a positive int")


def make_cli():
    cli = argparse.ArgumentParser(
        "Plot the significant interactions for the given .parquet file or chromosome chromosome pair."
    )

    cli.add_argument(
        "hic-file",
        type=pathlib.Path,
        help="Path to a Hi-C matrix in .hic, .mcool, or .cool format",
    )
    cli.add_argument(
        "sig-interactions",
        type=existing_file,
        help="Path to the .parquet or TSV file produced by either NCHG filter or NCHG view.",
    )
    cli.add_argument(
        "output-path",
        type=pathlib.Path,
        help="Path prefix where to store the output plot. Should be a valid file path without extension.",
    )
    cli.add_argument(
        "--plot-format",
        type=str,
        default="png",
        help="Format to use for plotting. Can be any file extension recognized by matplotlib. Ignored when --chrom1 or --chrom2 have been provided.",
    )
    cli.add_argument(
        "--chrom1",
        type=str,
        help="Name of the first chromosome. When not provided, all chromosomes found in the .parquet file will be processed.",
    )
    cli.add_argument(
        "--chrom2",
        type=str,
        help="Name of the second chromosome. When not provided, all chromosomes found in the .parquet file will be processed.",
    )
    cli.add_argument(
        "--expected-values",
        type=pathlib.Path,
        default=None,
        help="Path to the expected values file in HDF5 format produced by NCHG expected.",
    )
    cli.add_argument(
        "--resolution",
        type=positive_int,
        help="Hi-C matrix resolution.\n" "Required when input matrix is in .hic or .mcool format.",
    )
    cli.add_argument(
        "--pvalue",
        type=probability,
        default=1.0,
        help="Pvalue threshold used to classify significant interactions.",
    )
    cli.add_argument(
        "--log-ratio",
        type=float,
        default=-np.inf,
        help="Log-ratio threshold used to classify significant interactions.",
    )
    cli.add_argument(
        "--min-value",
        type=float,
        default=None,
        help="Lower bound for the color scale of the log-ratio heatmap.",
    )
    cli.add_argument(
        "--max-value",
        type=float,
        default=None,
        help="Upper bound for the color scale of the log-ratio heatmap.",
    )
    cli.add_argument(
        "--nproc",
        type=positive_int,
        default=1,
        help="Maximum number of parallel processes to use.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing file(s).",
    )

    return cli


def preprocess_data(df: pd.DataFrame) -> pd.DataFrame:
    df["chrom1"] = df["chrom1"].astype("string")
    df["chrom2"] = df["chrom2"].astype("string")

    chroms = pd.Series(df["chrom1"].unique().tolist() + df["chrom2"].unique().tolist()).unique()

    df.sort_values(["chrom1", "start1", "chrom2", "start2"], inplace=True)
    df["chrom1"] = pd.Categorical(df["chrom1"], categories=chroms)
    df["chrom2"] = pd.Categorical(df["chrom2"], categories=chroms)

    df.reset_index(drop=True, inplace=True)

    return df


def import_data(path: pathlib.Path, chrom1: Optional[str], chrom2: Optional[str]) -> pd.DataFrame:
    logging.info("importing data from %s...", path)
    try:
        df = pd.read_parquet(path)
        logging.info("read %d records .parquet file %s", len(df), path)
    except pyarrow.lib.ArrowInvalid:
        df = pd.read_table(path)
        logging.info("read %d records from text file %s", len(df), path)

    df = preprocess_data(df)
    if chrom1 is not None:
        assert chrom2 is not None
        logging.info("dropping interactions not belonging to %s:%s", chrom1, chrom2)
        df = df[(df["chrom1"] == chrom1) & (df["chrom2"] == chrom2)]

    return df


def import_expected_values(
    path: pathlib.Path,
) -> Dict[Tuple[str, str], Tuple[NDArray, NDArray]]:
    evs = {}
    logging.info("reading expected values from %s...", path)
    with h5py.File(path) as h5:
        chrom1 = h5["bin-masks/chrom1"][:]
        chrom2 = h5["bin-masks/chrom2"][:]
        offsets1 = h5["bin-masks/offsets1"][:]
        offsets2 = h5["bin-masks/offsets2"][:]

        values1 = h5["bin-masks/values1"][:]
        values2 = h5["bin-masks/values2"][:]

        for i, (chrom1, chrom2) in enumerate(zip(chrom1, chrom2)):
            i0 = offsets1[i]
            i1 = offsets1[i + 1]

            j0 = offsets2[i]
            j1 = offsets2[i + 1]
            evs[(chrom1.decode("utf-8"), chrom2.decode("utf-8"))] = (
                values1[i0:i1],
                values2[j0:j1],
            )

    return evs


def fetch_hic_matrix(
    path: pathlib.Path,
    resolution: Optional[int],
    chrom1: str,
    chrom2: str,
    expected_values,
) -> Tuple[int, NDArray]:
    logging.info("fetching %s:%s interactions from %s matrix at %d resolution...", chrom1, chrom2, path, resolution)
    f = hictkpy.File(path, resolution)
    m = f.fetch(chrom1, chrom2).to_numpy().astype(float)

    if expected_values is not None:
        mask1, mask2 = expected_values[(chrom1, chrom2)]
        m[mask1] = np.nan
        m[:, mask2] = np.nan

    return f.resolution(), m


def df_to_matrix(df: pd.DataFrame, bin_size: int, shape: Tuple[int, int]) -> NDArray:
    logging.info("populating matrix of significant interactions...")
    m = np.full(shape, np.nan)

    cols = ["start1", "end1", "start2", "end2", "log_ratio"]
    for _, (start1, end1, start2, end2, log_ratio) in df[cols].iterrows():
        i0 = int(start1 // bin_size)
        i1 = int(end1 // bin_size)
        j0 = int(start2 // bin_size)
        j1 = int(end2 // bin_size)

        if np.isfinite(log_ratio):
            m[i0:i1, j0:j1] = np.maximum(np.nan_to_num(m[i0:i1, j0:j1], nan=-np.inf), log_ratio)

    return m


def generate_chrom_pairs(df: pd.DataFrame) -> List[Tuple[str, str]]:
    return (
        df[["chrom1", "chrom2"]]
        .drop_duplicates()
        .apply(lambda row: (row["chrom1"], row["chrom2"]), axis="columns")
        .tolist()
    )


def process_chromosome_pair(
    chrom1: str,
    chrom2: str,
    df: pd.DataFrame,
    matrix_file: pathlib.Path,
    output_path: pathlib.Path,
    resolution: Optional[int],
    pvalue: float,
    log_ratio: float,
    evs: Optional[NDArray[float]],
    min_value: float,
    max_value: float,
):
    logging.info("processing %s:%s interactions...", chrom1, chrom2)
    if "pvalue_corrected" in df:
        pval_col = "pvalue_corrected"
    else:
        pval_col = "pvalue"

    resolution, obs_matrix = fetch_hic_matrix(matrix_file, resolution, chrom1, chrom2, evs)

    if obs_matrix.sum() == 0:
        logging.warn("file has no interactions for %s:%s matrix: no plots will be generated!", chrom1, chrom2)
        return

    logging.info("using %s column for filtering by p-value", pval_col)

    pval_significant = df[pval_col] < pvalue
    log_ratio_significant = (df["log_ratio"] >= log_ratio) & (~df["log_ratio"].isna())
    num_interactions = len(df)

    logging.info("filtering interactions using p-value=%g and log_ratio=%.2f...", pvalue, log_ratio)
    df = df[pval_significant & log_ratio_significant]
    logging.info("dropped %d/%d interactions", num_interactions - len(df), len(df))

    sig_matrix = df_to_matrix(df, resolution, obs_matrix.shape)  # noqa

    logging.info("generating plots for %s:%s...", chrom1, chrom2)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(2 * 6.4, 6.4))

    img0 = ax0.imshow(obs_matrix, norm=LogNorm(), interpolation="nearest")
    img1 = ax1.imshow(
        sig_matrix,
        vmin=min_value,
        vmax=max_value,
        interpolation="nearest",
        cmap="Reds",
    )

    ax0.set(xlabel=chrom2, ylabel=chrom1, title="Observed matrix")
    ax1.set(
        xlabel=chrom2,
        ylabel=chrom1,
        title="Significant interactions (Log-ratio)",
    )

    plt.colorbar(img0, ax=ax0)
    plt.colorbar(img1, ax=ax1)

    plt.tight_layout()

    logging.info("writing plots to %s...", output_path)
    fig.savefig(output_path, dpi=300)

    plt.close(fig)


def generate_placeholder_plot(path: pathlib.Path):
    fig, ax = plt.subplots()
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        "unable to find any significant interaction!",
        ha="center",
        va="center",
        fontsize=40,
        wrap=True,
    )

    logging.info("writing placeholder plot to %s...", path)
    fig.savefig(path, dpi=300)

    plt.close(fig)


def process_many(args: Dict[str, Any], df: pd.DataFrame, evs: Optional[NDArray[float]]):
    output_path = args["output-path"]
    extension = args["plot_format"]

    with mp.Pool(min(mp.cpu_count(), args["nproc"])) as pool:
        tasks = []
        for chrom1, chrom2 in generate_chrom_pairs(df):
            dff = df[(df["chrom1"] == chrom1) & (df["chrom2"] == chrom2)]
            dest = pathlib.Path(f"{output_path}.{chrom1}.{chrom2}.{extension}")
            if not args["force"] and dest.exists():
                raise RuntimeError(f"Refusing to overwrite existing file {output_path}. Pass --force to overwrite.")
            tasks.append(
                pool.apply_async(
                    process_chromosome_pair,
                    (
                        chrom1,
                        chrom2,
                        dff,
                        args["hic-file"],
                        dest,
                        args["resolution"],
                        args["pvalue"],
                        args["log_ratio"],
                        evs,
                        args["min_value"],
                        args["max_value"],
                    ),
                )
            )

        for t in tasks:
            t.get()


def main():
    args = vars(make_cli().parse_args())

    min_value = args["min_value"]
    max_value = args["max_value"]

    if min_value is not None and max_value is not None and min_value > max_value:
        raise RuntimeError("--min-value must be smaller than --max-value")

    chrom1 = args["chrom1"]
    chrom2 = args["chrom2"]
    if chrom1 is not None and chrom2 is None:
        chrom2 = chrom1

    df = import_data(args["sig-interactions"], chrom1, chrom2)

    evs = None
    if args["expected_values"] is not None:
        evs = import_expected_values(args["expected_values"])

    output_path = args["output-path"]
    extension = args["plot_format"]

    if chrom1 is not None:
        if chrom2 is None:
            chrom2 = chrom1

        output_path = pathlib.Path(f"{output_path}.{extension}")
        if not args["force"] and output_path.exists():
            raise RuntimeError(f"Refusing to overwrite existing file {output_path}. Pass --force to overwrite.")

        process_chromosome_pair(
            chrom1,
            chrom2,
            df,
            args["hic-file"],
            args["output-path"],
            args["resolution"],
            args["pvalue"],
            args["log_ratio"],
            evs,
            args["min_value"],
            args["max_value"],
        )
        return

    if len(df) == 0:
        output_path = pathlib.Path(f"{output_path}.placeholder.{extension}")
        if not args["force"] and output_path.exists():
            raise RuntimeError(f"Refusing to overwrite existing file {output_path}. Pass --force to overwrite.")
        generate_placeholder_plot(output_path)
        return

    process_many(args, df, evs)


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
