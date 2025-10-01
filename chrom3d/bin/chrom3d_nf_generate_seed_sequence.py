#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import hashlib
import pathlib
import random
from typing import List


def make_cli() -> argparse.ArgumentParser:
    def positive_int(s) -> float:
        if (x := int(s)) > 0:
            return x

        raise RuntimeError("Not a positive int")

    cli = argparse.ArgumentParser(
        "Use one or more file as source of entropy to generate a sequence of seeds. Seeds are printed to stdout."
    )
    cli.add_argument(
        "files",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more files to use as source of entropy",
    )

    cli.add_argument(
        "--number-of-seeds",
        type=positive_int,
        default=50,
        help="Number of seeds to generate.",
    )
    cli.add_argument(
        "--sep",
        type=str,
        default="\n",
        help="Seed separator.",
    )
    cli.add_argument(
        "--lower-bound",
        type=int,
        default=0,
        help="Lower bound for the uniform distribution from which seeds are sampled.",
    )
    cli.add_argument(
        "--upper-bound",
        type=int,
        default=(2**32) - 1,
        help="Upper bound for the uniform distribution from which seeds are sampled.",
    )

    return cli


def hash_files(files: List[pathlib.Path], chunk_size: int = 64 * 1024 * 1024) -> str:
    hasher = hashlib.sha512()

    for f in set(files):
        with open(f, "rb") as fp:
            while data := fp.read(chunk_size):
                hasher.update(data)

    return hasher.hexdigest()


def generate_seeds(lb: int, ub: int, num_seeds: int) -> List[int]:
    """
    Generate a sequence of seeds ensuring that there are no duplicate seeds, and that the resulting
    sequence has a deterministic (but random-looking) order
    """
    seeds = set()
    while len(seeds) != num_seeds:
        seeds.add(random.randint(lb, ub))

    seeds = list(sorted(seeds))
    random.shuffle(seeds)

    return seeds


def print_seeds(seeds: List[int], sep: str):
    padding = len(str(len(seeds)))

    # prints something like 012    1234567890
    seeds_str = (f"{i:0{padding}d}\t{n}" for i, n in enumerate(seeds))
    print(sep.join(seeds_str))


def main():
    args = vars(make_cli().parse_args())

    digest = hash_files(args["files"])

    random.seed(digest, version=2)

    seeds = generate_seeds(
        lb=args["lower_bound"],
        ub=args["upper_bound"],
        num_seeds=args["number_of_seeds"],
    )
    print_seeds(seeds, args["sep"])


if __name__ == "__main__":
    main()
