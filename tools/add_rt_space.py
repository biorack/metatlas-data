#!/usr/bin/env python3
"""
Add an 'rt_space' column with value 'HF_Aug2019' to all input TSV files in place.
"""

import argparse
import csv
import sys
from pathlib import Path


RT_SPACE_VALUE = "HF_Aug2019"
RT_SPACE_COLUMN = "rt_space"


def add_rt_space(input_path):
    input_path = Path(input_path)

    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            print(f"Error: Could not read headers from {input_path}", file=sys.stderr)
            return

        input_columns = list(reader.fieldnames)
        rows = list(reader)

    if RT_SPACE_COLUMN in input_columns:
        print(f"Warning: '{RT_SPACE_COLUMN}' already exists in {input_path}, overwriting values.")
    else:
        input_columns.append(RT_SPACE_COLUMN)

    with open(input_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=input_columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in rows:
            row[RT_SPACE_COLUMN] = RT_SPACE_VALUE
            writer.writerow(row)

    print(f"Updated: {input_path}")


def main():
    parser = argparse.ArgumentParser(
        description=f"Add an '{RT_SPACE_COLUMN}' column with value '{RT_SPACE_VALUE}' to TSV files in place."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="One or more input TSV files to update in place.",
    )
    args = parser.parse_args()

    for input_file in args.input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"Error: File not found: {input_path}", file=sys.stderr)
            continue
        add_rt_space(input_path)


if __name__ == "__main__":
    main()