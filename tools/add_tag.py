#!/usr/bin/env python3
"""
Add a 'tags' column with empty values and rename 'standard_source' to 'source'
in all input TSV files in place.
"""

import argparse
import csv
import sys
from pathlib import Path


TAG_COLUMN = "tags"
TAG_VALUE = ""
RENAME = {"standard_source": "source"}


def update_tsv(input_path):
    input_path = Path(input_path)

    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            print(f"Error: Could not read headers from {input_path}", file=sys.stderr)
            return

        input_columns = list(reader.fieldnames)
        rows = list(reader)

    output_columns = [RENAME.get(col, col) for col in input_columns]

    if TAG_COLUMN in output_columns:
        print(f"Warning: '{TAG_COLUMN}' already exists in {input_path}, overwriting values.")
    else:
        output_columns.append(TAG_COLUMN)

    with open(input_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in rows:
            # Apply renames
            for old, new in RENAME.items():
                if old in row:
                    row[new] = row.pop(old)

            row[TAG_COLUMN] = TAG_VALUE
            writer.writerow(row)

    print(f"Updated: {input_path}")


def main():
    parser = argparse.ArgumentParser(
        description=f"Add a '{TAG_COLUMN}' column and rename 'standard_source' to 'source' in TSV files in place."
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
        update_tsv(input_path)


if __name__ == "__main__":
    main()