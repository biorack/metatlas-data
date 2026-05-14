#!/usr/bin/env python3
"""
Merge identification_notes from a CSV into TSV files based on matching inchi_key.
If the TSV row already has identification_notes, the new notes are appended.
"""

import argparse
import csv
import sys
from pathlib import Path


def load_notes(notes_csv_path):
    """
    Load the input CSV of inchi_key -> identification_notes.
    Returns a dict of {inchi_key: identification_notes}.
    """
    notes = {}
    with open(notes_csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        inchi_col, notes_col = header[0], header[1]
        for row in reader:
            if len(row) < 2:
                continue
            inchi_key = row[0].strip()
            note = row[1].strip()
            if inchi_key:
                notes[inchi_key] = note
    return notes


def merge_notes_into_tsv(input_path, output_path, notes):
    """
    Read a TSV, merge identification_notes from the notes dict by inchi_key,
    and write to output_path (or stdout if None).
    """
    input_path = Path(input_path)

    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        input_columns = reader.fieldnames

        if input_columns is None:
            print(f"Error: Could not read headers from {input_path}", file=sys.stderr)
            return

        if "inchi_key" not in input_columns:
            print(f"Warning: No inchi_key column in {input_path}, skipping.", file=sys.stderr)
            return

        rows = list(reader)

    # Add identification_notes column if not already present
    output_columns = list(input_columns)
    if "identification_notes" not in output_columns:
        output_columns.append("identification_notes")

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

        for row in rows:
            inchi_key = row.get("inchi_key", "").strip()
            new_note = notes.get(inchi_key, "")

            if new_note:
                existing = row.get("identification_notes", "").strip()
                if existing:
                    row["identification_notes"] = f"{existing}; {new_note}"
                else:
                    row["identification_notes"] = new_note

            writer.writerow(row)

    finally:
        if output_path:
            out_file.close()

    if output_path:
        print(f"Written: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge identification_notes from a CSV into TSV files by inchi_key."
    )
    parser.add_argument(
        "notes_csv",
        help="Input CSV file with columns: inchi_key, identification_notes.",
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="One or more input TSV files to update.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help=(
            "Directory to write output files. If not specified and only one input "
            "file is given, output goes to stdout. If multiple files are given, "
            "defaults to writing alongside each input file with '_updated' suffix."
        ),
    )
    parser.add_argument(
        "--suffix",
        default="_updated",
        help="Suffix to append to output filenames (default: '_updated').",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite input files in place (ignores --output-dir and --suffix).",
    )
    args = parser.parse_args()

    notes = load_notes(args.notes_csv)
    print(f"Loaded {len(notes)} inchi_key entries from {args.notes_csv}")

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

        merge_notes_into_tsv(input_path, output_path, notes)


if __name__ == "__main__":
    main()