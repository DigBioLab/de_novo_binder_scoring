#!/usr/bin/env python3
import os
import csv
import argparse
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces


def calculate_scores(native_path, model_path, chain_map):
    """
    Run DockQ: compare model to native (input) and return total score.
    """
    try:
        native_structure = load_PDB(native_path)
        model_structure  = load_PDB(model_path)
    except Exception as e:
        print(f"[dockq] load error:\n  native: {native_path}\n  model: {model_path}\n  {e}")
        return None

    try:
        _, total_score = run_on_all_native_interfaces(
            model_structure, native_structure, chain_map=chain_map
        )
        return total_score
    except Exception as e:
        print(f"[dockq] run error on {model_path}: {e}")
        return None


def analyze_against_input(input_folder: str, models: dict, output_csv: str, chain_map=None):
    """
    For each binder in input_folder (reference), compute DockQ against
    any provided model folders.

    models = { "af2": "/path/to/af2/pdbs", "cf": "/path/to/colab/pdbs" }
    """
    if chain_map is None:
        chain_map = {"A": "A", "B": "B"}

    binder_ids = sorted({
        os.path.splitext(f)[0]
        for f in os.listdir(input_folder)
        if f.lower().endswith(".pdb")
    })

    if not binder_ids:
        raise RuntimeError(f"No PDBs found in input folder: {input_folder}")

    rows = []
    for bid in binder_ids:
        native_path = os.path.join(input_folder, bid + ".pdb")
        if not os.path.isfile(native_path):
            continue

        row = {"binder_id": bid}
        for prefix, folder in models.items():
            model_pdb = os.path.join(folder, bid + ".pdb")
            colname = f"{prefix}_dockQ"

            if os.path.isfile(model_pdb):
                score = calculate_scores(native_path, model_pdb, chain_map)
                row[colname] = score
            else:
                row[colname] = None
                print(f"[miss] {prefix}: no pdb for {bid} at {model_pdb}")

        rows.append(row)

    # union of all column names
    fieldnames = ["binder_id"] + [f"{p}_dockQ" for p in models.keys()]

    with open(output_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"[ok] CSV written: {output_csv}")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute DockQ against required input folder, "
                    "with any number of model folders + prefixes."
    )
    ap.add_argument("--input_pdbs", required=True,
                    help="Path to REQUIRED reference/input PDBs (binder_id.pdb).")
    ap.add_argument("--model", action="append", required=True,
                    help="Model definition in the form prefix:/path/to/pdbs. "
                         "Can be given multiple times.")
    ap.add_argument("--out-csv", required=True, help="Destination CSV path.")
    ap.add_argument("--mapA", default="A", help="Map input chain A to this chain id in models")
    ap.add_argument("--mapB", default="B", help="Map input chain B to this chain id in models")
    return ap.parse_args()


def main():
    args = parse_args()

    models = {}
    for m in args.model:
        if ":" not in m:
            raise ValueError(f"--model must be prefix:/path form, got {m}")
        prefix, folder = m.split(":", 1)
        folder = os.path.abspath(folder)
        if not os.path.isdir(folder):
            raise FileNotFoundError(f"Model folder not found: {folder}")
        models[prefix] = folder

    chain_map = {"A": args.mapA, "B": args.mapB}
    analyze_against_input(args.input_pdbs, models, args.output_csv, chain_map)


if __name__ == "__main__":
    main()
