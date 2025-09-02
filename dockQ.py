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


def parse_folder_kv(value):
    """
    Parse --folder flags of the form name:path
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError("Folder must be given as name:path")
    name, path = value.split(":", 1)
    name = name.strip()
    path = os.path.abspath(os.path.expandvars(os.path.expanduser(path.strip())))
    if not name:
        raise argparse.ArgumentTypeError("Folder name is empty")
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Folder path does not exist: {path}")
    return name, path

def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute DockQ against required input folder, "
                    "with any number of model folders + prefixes."
    )
    ap.add_argument(
        "--input-pdbs",
        required=True,
        help="Path to REQUIRED reference/input PDBs (binder_id.pdb)."
    )
    ap.add_argument(
        "--folder",
        action="append",
        required=True,
        type=parse_folder_kv,
        help='Repeatable. Format "name:path". Example: --folder af3:/path/to/AF3/pdbs'
    )
    ap.add_argument("--out-csv", required=True, help="Destination CSV path.")
    ap.add_argument("--mapA", default="A", help="Map input chain A to this chain id in models")
    ap.add_argument("--mapB", default="B", help="Map input chain B to this chain id in models")
    return ap.parse_args()

def main():
    args = parse_args()

    # Convert list of (name, path) to dict
    models = dict(args.folder)

    chain_map = {"A": args.mapA, "B": args.mapB}
    analyze_against_input(args.input_pdbs, models, args.out_csv, chain_map)

if __name__ == "__main__":
    main()
