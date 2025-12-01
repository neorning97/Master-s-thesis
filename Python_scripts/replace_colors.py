# Script that changes the color of the beads to specific colors, so the same chromosomes have the same color on their beads across all models.

#!/usr/bin/env python3


import pandas as pd
import argparse
import pathlib
import re
import sys


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to a TSV with the color mappings.",
    )

    return cli


def read_lines():
    """
    Read file line by line from stdin and return lines as a list of strings
    """
    return list(sys.stdin.readlines())


def read_table(path):
    """
    Read a TSV like:
    chrom   red   green   blue
    chr1   1.0   1.0   1.0
    chr2   1.0   0.0   0.0
    """
    return pd.read_table(path)


def replace_color(lines, chrom, red, green, blue):
    # Use REGEX to match something like r="0.0822875" g="0.0589692" b="0.831684"
    pattern = re.compile(r'([rgb]=)"[\d.]+"\s+([rgb]=)"[\d.]+"\s+([rgb]=)"[\d.]+"')

    chrom_pattern = re.compile(rf'\b{chrom}([_-]\w+)?\b')

    for line_number, line in enumerate(lines):
        if not chrom_pattern.search(line):
            # Line does not refer to our chromosome: do nothing
            continue

        # replace colors and update list of lines
        new_line = pattern.sub(repl=rf'\1"{red}" \2"{green}" \3"{blue}"', string=line)
        lines[line_number] = new_line

    return lines


def write_lines(lines, path):
    with pathlib.Path(path).open("wt") as f:
        f.write("\n".join(lines))


def main():
    args = make_cli().parse_args()

    lines = read_lines()
    color_mappings = read_table(args.tsv)

    # Loop over rows from the table with color mappings
    for chrom, red, green, blue in color_mappings.itertuples(index=False):
        lines = replace_color(lines, chrom, red, green, blue)

    # Print lines with the new colors to stdout
    for line in lines:
        print(line, end="")


if __name__ == "__main__":
    main()
