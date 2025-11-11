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

from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector, ResiduePropertySelector
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, OrResidueSelector, NotResidueSelector,
    ResiduePropertySelector, AndResidueSelector
)
from pyrosetta.rosetta.core.chemical import ResidueProperty

# ---------------------------
# Helpers

def get_binder_id(pdb_path: str) -> str:
    return os.path.splitext(os.path.basename(pdb_path))[0]

def _chains_list(chains_str: str):
    return [c.strip() for c in chains_str.split(",") if c.strip()]


def _subpose_for_chains(pose: pr.Pose, chains_str: str, ignore_ions: bool = False) -> pr.Pose:
    # Base selector: all requested chains
    sel = OrResidueSelector()
    for ch in _chains_list(chains_str):
        sel.add_residue_selector(ChainSelector(ch))

    if ignore_ions:
        ion_sel = ResiduePropertySelector(ResidueProperty.METAL)
        not_ions = NotResidueSelector(ion_sel)
        # Combine: chain selector AND NOT ions
        sel = AndResidueSelector(sel, not_ions)

    subset = sel.apply(pose)  # still a vector1_bool
    res_idxs = get_residues_from_subset(subset)
    sub = pr.Pose()
    pose_from_pose(sub, pose, res_idxs)
    return sub

def _pose_without_ions(pose: pr.Pose) -> pr.Pose:
    """Return a copy of pose without metal ions (e.g. CA, MG, ZN)."""
    from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector, ResiduePropertySelector
    from pyrosetta.rosetta.core.chemical import ResidueProperty
    ion_sel = ResiduePropertySelector(ResidueProperty.METAL)
    not_ions = NotResidueSelector(ion_sel)
    subset = not_ions.apply(pose)
    res_idxs = get_residues_from_subset(subset)
    clean = pr.Pose()
    pose_from_pose(clean, pose, res_idxs)
    return clean


def align_by_chain(reference_pdb: str, mobile_pdb: str,
                   ref_chain: str, mob_chain: str,
                   out_path: str,
                   ignore_ions: bool = False) -> None:
    """
    Superimpose mobile_pdb onto reference_pdb using a single chain as guide,
    optionally removing ions first.
    """
    ref_pose = pr.pose_from_pdb(reference_pdb)
    mob_pose = pr.pose_from_pdb(mobile_pdb)
 

    if ignore_ions:
        ref_pose = _pose_without_ions(ref_pose)
        mob_pose = _pose_without_ions(mob_pose)

    # print(mob_pose)
    # print(ref_pose)

    mover = AlignChainMover()
    mover.pose(ref_pose)

    ref_chain_num = rosetta.core.pose.get_chain_id_from_chain(ref_chain, ref_pose)
    mob_chain_num = rosetta.core.pose.get_chain_id_from_chain(mob_chain, mob_pose)



    mover.source_chain(mob_chain_num)
    mover.target_chain(ref_chain_num)
    mover.apply(mob_pose)

    mob_pose.dump_pdb(out_path)

def rmsd_between(pdb_ref: str, pdb_mob: str, chains_str: str,
                 bb_only: bool = False, ignore_ions: bool = False) -> float:
    ref = pr.pose_from_pdb(pdb_ref)
    mob = pr.pose_from_pdb(pdb_mob)

    ref_sub = _subpose_for_chains(ref, chains_str, ignore_ions=ignore_ions)
    mob_sub = _subpose_for_chains(mob, chains_str, ignore_ions=ignore_ions)

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
    ap.add_argument("--update-runcsv",action="store_true",help="If set, update the run CSV (default: only write standalone out-csv)")

    args = ap.parse_args()

    pr.init('-ignore_unrecognized_res -ignore_zero_occupancy -mute all -corrections::beta_nov16 true')

    model_folders = dict(args.folder)  # {model_name: path}
    # Load ions_in_target info if run_csv is given
    ions_map = {}
    if args.run_csv:
        run_df = pd.read_csv(Path(args.run_csv), dtype=str).fillna("")
        if "binder_id" not in run_df.columns:
            raise SystemExit("--run-csv must contain a 'binder_id' column")
        if "ions_in_target" in run_df.columns:
            ions_map = dict(zip(run_df["binder_id"], run_df["ions_in_target"]))
            print(f"Loaded ions_in_target info for {len(ions_map)} binders from {args.run_csv}")

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
                    ignore_ions = bool(ions_map.get(binder_id))  # True if non-empty
                    print(f"RMSD for {binder_id} (ignore_ions={ignore_ions}): {ions_map.get(binder_id)}")
                    align_by_chain(p1, p2, ref_chain=bind_ch, mob_chain=bind_ch, out_path=tmp_p2_A, ignore_ions=ignore_ions)

                    ignore_ions = bool(ions_map.get(binder_id))  # True if non-empty string
                    #print(f"RMSD for {binder_id} (ignore_ions={ignore_ions}): {ions_map.get(binder_id)}")
                    r_complex = rmsd_between(
                        p1, tmp_p2_A,
                        chains_str=f"{bind_ch},{tgt_ch}",
                        bb_only=args.bb_only,
                        ignore_ions=ignore_ions
                    )
                    row[f"RMSD_complex_{m1}_{m2}"] = round(r_complex, 3)
                    if args.verbose:
                        print(f"RMSD complex for {binder_id} between {m1} and {m2}: {r_complex}")
                    all_cols.add(f"RMSD_complex_{m1}_{m2}")
                except Exception as e:
                    # If alignment or RMSD fails, record error but only alpha numerical value
                    e=str(e).replace(",", ";").replace("\n", " ").replace("\r", " ").replace("\t", " ")
                    if args.verbose:
                        print(f"Error computing RMSD for {binder_id} between {m1} and {m2}: {e}")
                    #print(f"RMSD for {binder_id} (ignore_ions={ignore_ions}): {ions_map.get(binder_id)}")
                    row[f"RMSD_complex_{m1}_{m2}"] = f"ERR:{e}"
                    all_cols.add(f"RMSD_complex_{m1}_{m2}")

                # Binder chain RMSD after aligning by target
                try:
                    tmp_p2_B = os.path.join(tmpdir, f"{binder_id}__{m2}__Balign.pdb")
                    ignore_ions = bool(ions_map.get(binder_id))  # True if non-empty
                    align_by_chain(p1, p2, ref_chain=tgt_ch, mob_chain=tgt_ch, out_path=tmp_p2_B, ignore_ions=ignore_ions)
                    r_a_after_b = rmsd_between(p1, tmp_p2_B, chains_str=bind_ch, bb_only=args.bb_only)
                    row[f"RMSD_chA_aft_chB_align_{m1}_{m2}"] = round(r_a_after_b, 3)
                    if args.verbose:
                        print(f"RMSD for {binder_id} between {m1} and {m2}: {r_a_after_b}")
                    all_cols.add(f"RMSD_chA_aft_chB_align_{m1}_{m2}")
                except Exception as e:
                    # If alignment or RMSD fails, record error but only alpha numerical values
                    if args.verbose:
                        print(f"Error computing RMSD for {binder_id} between {m1} and {m2}: {e}")
                    e=str(e).replace(",", ";").replace("\n", " ").replace("\r", " ").replace("\t", " ")
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
    
    if args.run_csv and not args.update_runcsv and args.verbose:
        print("Note: --run-csv provided but --update-runcsv not set, skipping merge")

    # Merge into run CSV only if explicitly requested
    if args.update_runcsv and args.run_csv:
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
    else:
        ions_map = {}

if __name__ == "__main__":
    main()


