#!/usr/bin/env python3
import os
import argparse
import csv
import glob
import tempfile
from itertools import combinations
from pathlib import Path
import pandas as pd

# PyRosetta
import pyrosetta as pr
from pyrosetta import rosetta
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, OrResidueSelector
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover

# ---------------------------
# Helpers

def get_binder_id(pdb_path: str) -> str:
    return os.path.splitext(os.path.basename(pdb_path))[0]

def _chains_list(chains_str: str):
    return [c.strip() for c in chains_str.split(",") if c.strip()]

def _subpose_for_chains(pose: pr.Pose, chains_str: str) -> pr.Pose:
    sel = OrResidueSelector()
    for ch in _chains_list(chains_str):
        sel.add_residue_selector(ChainSelector(ch))
    subset = sel.apply(pose)
    res_idxs = get_residues_from_subset(subset)
    sub = pr.Pose()
    pose_from_pose(sub, pose, res_idxs)
    return sub

def align_by_chain(reference_pdb: str, mobile_pdb: str,
                   ref_chain: str, mob_chain: str,
                   out_path: str) -> None:
    """
    Superimpose mobile_pdb onto reference_pdb using a single chain as guide,
    write aligned mobile copy to out_path (original files remain untouched).
    """
    ref_pose = pr.pose_from_pdb(reference_pdb)
    mob_pose = pr.pose_from_pdb(mobile_pdb)

    mover = AlignChainMover()
    mover.pose(ref_pose)

    ref_chain_num = rosetta.core.pose.get_chain_id_from_chain(ref_chain, ref_pose)
    mob_chain_num = rosetta.core.pose.get_chain_id_from_chain(mob_chain, mob_pose)

    mover.source_chain(mob_chain_num)
    mover.target_chain(ref_chain_num)
    mover.apply(mob_pose)

    mob_pose.dump_pdb(out_path)

def rmsd_between(pdb_ref: str, pdb_mob: str, chains_str: str, bb_only: bool = False) -> float:
    """
    Load two PDBs, extract subposes for the given chains (e.g. "A,B"),
    compute RMSD with RMSDMetric.
    """
    ref = pr.pose_from_pdb(pdb_ref)
    mob = pr.pose_from_pdb(pdb_mob)

    ref_sub = _subpose_for_chains(ref, chains_str)
    mob_sub = _subpose_for_chains(mob, chains_str)

    m = RMSDMetric()
    m.set_comparison_pose(ref_sub)
    if bb_only:
        m.set_bb_only(True)
    return float(m.calculate(mob_sub))

# ---------------------------
# CLI parsing

def parse_folder_kv(value):
    """
    Parse --folder flags of the form name:path
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError("Folder must be given as name:path")
    name, path = value.split(":", 1)
    name = name.strip()
    path = os.path.expandvars(os.path.expanduser(path.strip()))
    if not name:
        raise argparse.ArgumentTypeError("Folder name is empty")
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Folder path does not exist: {path}")
    return name, path

# ---------------------------
# Main
def main():
    ap = argparse.ArgumentParser(description="Pairwise RMSDs across model folders (AlignChainMover + RMSDMetric).")
    ap.add_argument(
        "--folder",
        action="append",
        required=True,
        type=parse_folder_kv,
        help='Repeatable. Format "name:path". Example: --folder af3:/path/to/AF3/pdbs'
    )
    ap.add_argument("--out-csv", required=True, help="Output CSV path for RMSD table")
    ap.add_argument("--run-csv", help="Optional run CSV to merge results into")
    ap.add_argument("--binder-chain", default="A", help="Binder chain letter (default: A)")
    ap.add_argument("--target-chain", default="B", help="Target chain letter (default: B)")
    ap.add_argument("--bb-only", action="store_true", help="Backbone-only RMSD")
    ap.add_argument("--backup", action="store_true", help="Write run.csv.bak before overwrite")
    ap.add_argument("--verbose", action="store_true", help="Print more info while processing")
    args = ap.parse_args()

    pr.init('-ignore_unrecognized_res -ignore_zero_occupancy -mute all -corrections::beta_nov16 true')

    model_folders = dict(args.folder)  # {model_name: path}

    # binder_id -> {model: pdb_path}
    binder_map = {}
    for model, folder in model_folders.items():
        for pdb_path in glob.glob(os.path.join(folder, "*.pdb")):
            binder_id = get_binder_id(pdb_path)
            binder_map.setdefault(binder_id, {})[model] = pdb_path

    # Prepare output rows/columns
    rows = []
    all_cols = set(["binder_id"])
    bind_ch = args.binder_chain
    tgt_ch = args.target_chain
    tmpdir = tempfile.mkdtemp(prefix="rmsd_align_tmp_")

    try:
        for binder_id, model_to_pdb in sorted(binder_map.items()):
            if len(model_to_pdb) < 2:
                continue

            row = {"binder_id": binder_id}
            for (m1, p1), (m2, p2) in combinations(sorted(model_to_pdb.items()), 2):
                # Complex RMSD
                try:
                    tmp_p2_A = os.path.join(tmpdir, f"{binder_id}__{m2}__Aalign.pdb")
                    align_by_chain(p1, p2, ref_chain=bind_ch, mob_chain=bind_ch, out_path=tmp_p2_A)
                    r_complex = rmsd_between(p1, tmp_p2_A, chains_str=f"{bind_ch},{tgt_ch}", bb_only=args.bb_only)
                    row[f"RMSD_complex_{m1}_{m2}"] = round(r_complex, 3)
                    all_cols.add(f"RMSD_complex_{m1}_{m2}")
                except Exception as e:
                    row[f"RMSD_complex_{m1}_{m2}"] = f"ERR:{e}"
                    all_cols.add(f"RMSD_complex_{m1}_{m2}")

                # Binder chain RMSD after aligning by target
                try:
                    tmp_p2_B = os.path.join(tmpdir, f"{binder_id}__{m2}__Balign.pdb")
                    align_by_chain(p1, p2, ref_chain=tgt_ch, mob_chain=tgt_ch, out_path=tmp_p2_B)
                    r_a_after_b = rmsd_between(p1, tmp_p2_B, chains_str=bind_ch, bb_only=args.bb_only)
                    row[f"RMSD_chA_aft_chB_align_{m1}_{m2}"] = round(r_a_after_b, 3)
                    all_cols.add(f"RMSD_chA_aft_chB_align_{m1}_{m2}")
                except Exception as e:
                    row[f"RMSD_chA_aft_chB_align_{m1}_{m2}"] = f"ERR:{e}"
                    all_cols.add(f"RMSD_chA_aft_chB_align_{m1}_{m2}")

            rows.append(row)

    finally:
        try:
            for f in glob.glob(os.path.join(tmpdir, "*.pdb")):
                os.remove(f)
            os.rmdir(tmpdir)
        except Exception:
            pass

    # Write standalone output CSV
    headers = ["binder_id"] + sorted([c for c in all_cols if c != "binder_id"])
    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    if args.verbose:
        print(f"Written RMSD metrics to {args.out_csv} with {len(rows)} rows and {len(headers)} columns")

    # Merge into run CSV if provided
    if args.run_csv:
        run_csv_path = Path(args.run_csv)
        run_df = pd.read_csv(run_csv_path, dtype=str)
        if "binder_id" not in run_df.columns:
            raise SystemExit("--run-csv must contain a 'binder_id' column")

        res_df = pd.DataFrame(rows).fillna("")
        if args.backup:
            bak = run_csv_path.with_suffix(run_csv_path.suffix + ".bak")
            run_df.to_csv(bak, index=False)
            if args.verbose:
                print(f"Backed up {run_csv_path} to {bak}")

        merged_df = run_df.merge(res_df, on="binder_id", how="left")
        merged_df.to_csv(run_csv_path, index=False)
        if args.verbose:
            print(f"Updated {run_csv_path} with merged RMSD columns")

if __name__ == "__main__":
    main()


