#!/usr/bin/env python3
"""
Reorder TSV columns to match the desired output format.
Handles various input column name variations across different file formats.
"""

import argparse
import csv
import os
import sys
from pathlib import Path


DESIRED_COLUMNS = [
    "label",
    "inchi_key",
    "adduct",
    "mz",
    "rt_peak",
    "rt_min",
    "rt_max",
    "mz_tolerance",
    "polarity",
    "identification_notes",
    "classes",
    "pathways",
    "standard_source",
]

RENAME_MAP = {
    "detected_polarity": "polarity",
    "file_name": "standard_source",
    "compound_classes": "classes",
    "compound_pathways": "pathways",
}

DROP_COLUMNS = {
    "old_label",
    "permanent_index",
    "inchi",
    "rt_units",
    "mz_tolerance_units",
    "mono_isotopic_molecular_weight",
    "pubchem_compound_id",
    "synonyms",
    "ms1_notes",
    "ms2_notes",
    "file_type",
    "rt_overlap",
    "intensity",
    "metatlas_score",
    "metatlas_num_matches",
    "nist_score",
    "nist_num_matches",
    "mona_score",
    "mona_num_matches",
    "experiment",
    "Peak shape",
    "RERUN?",
    "library",
    "verified",
    "name",
    "valid_msms",
    "msms_score",
    "code",
    "code_weight_with_msms",
    "code_rank_with_msms",
    "code_weight_without_msms",
    "code_rank_without_msms",
    "not_blank",
    "confidence_category",
    "top_two_by_intensity",
    "top_two_by_rank",
    "in_metatlas",
}


def normalize_columns(input_columns):
    """
    Returns:
      - effective_rename: dict of {old_col -> new_col} for columns being kept
      - effective_drop: set of column names to drop
    """
    col_set = set(input_columns)
    effective_rename = {}
    effective_drop = set()

    has_label = "label" in col_set

    for col in input_columns:
        if col in DROP_COLUMNS:
            effective_drop.add(col)
            continue

        if "-summary_" in col or col.startswith("mz_centroid_ms1") or col.startswith("i_msms") or col.startswith("mz_msms") or col.startswith("rt_msms"):
            effective_drop.add(col)
            continue

        if col == "compound_name":
            if not has_label:
                effective_rename[col] = "label"
            else:
                effective_drop.add(col)
            continue

        if col in RENAME_MAP:
            effective_rename[col] = RENAME_MAP[col]

    return effective_rename, effective_drop


def reorder_tsv(input_path, output_path=None):
    input_path = Path(input_path)

    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        input_columns = reader.fieldnames

        if input_columns is None:
            print(f"Error: Could not read headers from {input_path}", file=sys.stderr)
            sys.exit(1)

        input_columns = [c for c in input_columns if c.strip()]

        effective_rename, effective_drop = normalize_columns(input_columns)

        # Map new column names -> original column names for lookup during row writing
        available = {}
        for col in input_columns:
            if col in effective_drop:
                continue
            new_name = effective_rename.get(col, col)
            available[new_name] = col

        # Always output all DESIRED_COLUMNS (blank if missing), then any extras
        desired_set = set(DESIRED_COLUMNS)
        output_columns = list(DESIRED_COLUMNS)

        for col in input_columns:
            if col in effective_drop:
                continue
            new_name = effective_rename.get(col, col)
            if new_name not in desired_set and new_name not in output_columns:
                output_columns.append(new_name)

        if output_path:
            out_file = open(output_path, "w", newline="", encoding="utf-8")
        else:
            out_file = sys.stdout

        try:
            writer = csv.DictWriter(
                out_file,
                fieldnames=output_columns,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in reader:
                new_row = {}
                for col in input_columns:
                    if col in effective_drop:
                        continue
                    new_name = effective_rename.get(col, col)
                    new_row[new_name] = row.get(col, "")

                # Fill any missing desired columns with empty string
                for col in DESIRED_COLUMNS:
                    if col not in new_row:
                        new_row[col] = ""

                writer.writerow(new_row)

        finally:
            if output_path:
                out_file.close()

    if output_path:
        print(f"Written: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Reorder TSV columns to a standard format."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="One or more input TSV files to reorder.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help=(
            "Directory to write output files. If not specified and only one input "
            "file is given, output goes to stdout. If multiple files are given, "
            "defaults to writing alongside each input file with '_reordered' suffix."
        ),
    )
    parser.add_argument(
        "--suffix",
        default="_reordered",
        help="Suffix to append to output filenames (default: '_reordered').",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite input file in place (ignores --output-dir and --suffix).",
    )
    args = parser.parse_args()

    multiple = len(args.input_files) > 1

    for input_file in args.input_files:
        input_path = Path(input_file)

        if not input_path.exists():
            print(f"Error: File not found: {input_path}", file=sys.stderr)
            continue

        if args.overwrite:
            output_path = str(input_path)
        elif multiple or args.output_dir:
            if args.output_dir:
                out_dir = Path(args.output_dir)
                out_dir.mkdir(parents=True, exist_ok=True)
                output_path = str(out_dir / (input_path.stem + args.suffix + input_path.suffix))
            else:
                output_path = str(input_path.with_name(input_path.stem + args.suffix + input_path.suffix))
        else:
            output_path = None  # stdout

        reorder_tsv(input_path, output_path)


if __name__ == "__main__":
    main()