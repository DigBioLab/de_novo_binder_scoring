#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def main():
    ap = argparse.ArgumentParser(
        description="Merge run.csv with one or more metric CSVs on binder_id"
    )
    ap.add_argument("--run-csv", required=True, help="Path to run.csv (base file)")
    ap.add_argument(
        "--metric-csvs", nargs="+", required=True,
        help="Paths to metric CSVs (will be merged into run.csv by binder_id)"
    )
    ap.add_argument(
        "--out-dir", required=True,
        help="Directory to save merged CSV (creates merged_run.csv inside)"
    )
    args = ap.parse_args()

    # --- Load run.csv
    if not os.path.isfile(args.run_csv):
        sys.exit(f"ERROR: run.csv not found: {args.run_csv}")
    run_df = pd.read_csv(args.run_csv, dtype=str)
    if "binder_id" not in run_df.columns:
        sys.exit("ERROR: run.csv must contain a 'binder_id' column")

    merged_df = run_df.copy()
    seen_cols = set(merged_df.columns)

    # --- Process metric CSVs
    for csv_path in args.metric_csvs:
        if not os.path.isfile(csv_path):
            print(f"[warn] metric CSV not found: {csv_path}", file=sys.stderr)
            continue

        df = pd.read_csv(csv_path, dtype=str)
        if "binder_id" not in df.columns:
            print(f"[warn] skipping {csv_path}: missing binder_id column", file=sys.stderr)
            continue

        # Check overlapping columns
        overlap = set(df.columns) & seen_cols - {"binder_id"}
        if overlap:
            print(f"[warn] overlapping columns in {csv_path}: {', '.join(sorted(overlap))}", file=sys.stderr)

        # Check coverage of binder_ids
        missing = set(run_df["binder_id"]) - set(df["binder_id"])
        if missing:
            print(f"[warn] {csv_path} missing {len(missing)} binder_ids present in run.csv", file=sys.stderr)

        # Merge
        merged_df = merged_df.merge(df, on="binder_id", how="left")
        seen_cols |= set(df.columns)

    # --- Save output
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "merged_run.csv")
    merged_df.to_csv(out_path, index=False)
    print(f"[info] merged CSV written to {out_path}")

if __name__ == "__main__":
    main()
