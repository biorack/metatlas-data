import sys
import argparse
import pandas as pd
from pathlib import Path

FILTER_RULES = {
    "HILIC": [
        ("rt_peak", 0.8, 17.5),
    ],
    "C18": [
        ("rt_peak", 0.5, 10.0),
    ],
}

def detect_type(path: Path) -> str:
    for key in FILTER_RULES:
        if key in path.name:
            return key
    return None

def filter_tsv(input_path: Path, output_path: Path) -> None:
    ftype = detect_type(input_path)
    if ftype is None:
        print(f"[ERROR] Could not determine filter type for '{input_path.name}'. "
              f"Filename must contain one of: {list(FILTER_RULES)}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_path, sep="\t", header=0)
    original_len = len(df)

    for col_name, min_val, max_val in FILTER_RULES[ftype]:
        if col_name not in df.columns:
            print(f"[ERROR] Column '{col_name}' not found in {input_path.name}. "
                  f"Available columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)
        col = pd.to_numeric(df[col_name], errors="coerce")
        if min_val is not None:
            df = df[col >= min_val]
        if max_val is not None:
            col = pd.to_numeric(df[col_name], errors="coerce")
            df = df[col <= max_val]

    df.to_csv(output_path, sep="\t", index=False)
    print(f"[{ftype}] {input_path.name}: {original_len} -> {len(df)} rows => {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Filter TSV files by rt_peak column bounds.")
    parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True,
                        help="Input TSV file(s)")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("."),
                        help="Output directory (default: current directory)")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    for input_path in args.inputs:
        output_path = args.outdir / input_path.name
        filter_tsv(input_path, output_path)

if __name__ == "__main__":
    main()
