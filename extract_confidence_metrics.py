#!/usr/bin/env python3
import os
import re
import csv
import glob
import json
import argparse
import shutil
import statistics
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, PDBParser

# =========================
# Utilities
# =========================

def _ensure_dir(d: str):
    os.makedirs(d, exist_ok=True)
    
def _make_atom_id_map(pose):
    amap = rosetta.core.id.AtomID_Map(rosetta.core.id.AtomID(), False)
    amap.resize(pose.size())
    return amap

def _merge_into_run(run_csv: str, df_new: pd.DataFrame):
    """
    Merge df_new into run.csv by 'binder_id'. Adds/overwrites columns present in df_new.
    """
    if not os.path.isfile(run_csv):
        raise FileNotFoundError(f"run.csv not found: {run_csv}")
    df_run = pd.read_csv(run_csv)
    if "binder_id" not in df_run.columns:
        raise ValueError("run.csv must contain 'binder_id'.")

    df_run = df_run.set_index("binder_id")
    df_new = df_new.set_index("binder_id")

    # only keep rows that exist in run.csv; warn about extras
    extras = [bid for bid in df_new.index if bid not in df_run.index]
    if extras:
        print(f"[merge] WARNING: {len(extras)} binder_id(s) not in run.csv; skipping first few: {extras[:5]}")

    df_new = df_new[df_new.index.isin(df_run.index)]
    for col in df_new.columns:
        df_run[col] = df_new[col]

    df_run.reset_index().to_csv(run_csv, index=False)
    print(f"[merge] Updated {run_csv} (+{len(df_new.columns)} columns)")

# ---------- PDB post-processing (inline from meta_analysis.py) ----------
                
def combine_subchains(folder_path):
    """
    Processes all PDB files in the given folder. For each file,
    it changes every chain with an ID other than "A" to "B" and
    saves the modified structure with the same file name.
    
    Parameters:
        folder_path (str): Path to the folder containing PDB files.
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    print(f"found argument combine subchains with folder {folder_path}")

    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(".pdb"):
            file_path = os.path.join(folder_path, filename)

            
            # Parse the structure from the PDB file
            structure = parser.get_structure(filename, file_path)
            
            # Iterate over all models and chains, renaming non-"A" chains to "B"

            for model in structure:
                for chain in model:
                    if chain.get_id() not in ["A","B"]:
                  
                        chain.id = "B"
            
            # Save the modified structure back to the same file (overwrites original)
            io.set_structure(structure)
            io.save(file_path)


def _renumber_structure(structure, starting_resid=1):
    """
    Renumber canonical residues (res.id[0] == ' ') in each chain, starting at starting_resid.
    HETATM/solvent keep original numbering.
    """
    for model in structure:
        for chain in model:
            next_id = starting_resid
            for res in chain:
                hetflag, old_id, icode = res.id  # tuple like (' ', 5, ' ')
                if hetflag == ' ':
                    res.id = (hetflag, next_id, ' ')
                    next_id += 1

def renumber_pdb_file(file_path, starting_resid=1):
    """
    Renumbers the residues of the PDB file at the given file_path, starting from starting_resid,
    and overwrites the original file with the updated content.

    Parameters
    ----------
    file_path : str
        Path to the input PDB file.
    starting_resid : int, optional
        The starting residue number (default is 1).

    Raises
    ------
    ValueError
        If the file is not found/readable or if the new residue number exceeds 9999.
    """
    if not os.path.isfile(file_path):
        raise ValueError(f"ERROR!! File not found or not readable: '{file_path}'")

    def pad_line(line):
        """Ensures the line is at least 80 characters long."""
        size_of_line = len(line)
        if size_of_line < 80:
            padding = 80 - size_of_line + 1
            line = line.rstrip('\n') + ' ' * padding + '\n'
        return line[:81]  # Exactly 80 characters plus the newline

    def run(lines, starting_resid):
        """Generator that yields renumbered PDB lines while handling TER cases."""
        prev_resid = None
        resid = starting_resid - 1  # adjust for first residue increment
        records = ('ATOM', 'HETATM', 'ANISOU')
        terminate_at_ter = False  # Flag to stop processing after a TER without more residues

        for i, line in enumerate(lines):
            line = pad_line(line)

            if terminate_at_ter:
                yield "END"

                break  # Stop processing after a standalone TER

            if line.startswith('MODEL'):
                resid = starting_resid - 1
                prev_resid = None
                yield line

            elif line.startswith('TER'):
                #print(line)
                # Look ahead to check if there are more residues
                next_lines = lines[i + 1:]  # Remaining lines
                has_more_residues = any(l.startswith(records) for l in next_lines)
                if has_more_residues:
                    yield line[:4] + " "*76 +"\n"
                if not has_more_residues:
                    yield line[:4] + " "*76 +"\n" +"END"
                

                if not has_more_residues:
                    terminate_at_ter = True  # Stop processing after this
                    #print("last res")

            elif line.startswith(records):
                line_resuid = line[17:27]
                if line_resuid != prev_resid:
                    prev_resid = line_resuid
                    resid += 1
                    if resid > 9999:
                        raise ValueError('Cannot set residue number above 9999.')
                yield line[:22] + str(resid).rjust(4) + line[26:]

            else:
                yield line

    # Read file, process, and overwrite
    with open(file_path, 'r') as f:
        lines = f.readlines()

    with open(file_path, 'w') as f:
        f.writelines(run(lines, starting_resid))


def reres(folder, starting_resid=1):
    """
    Applies residue renumbering to all PDB files in the given folder.
    """
    if not os.path.isdir(folder):
        raise ValueError(f"Provided folder does not exist or is not a directory: '{folder}'")
    for file_name in os.listdir(folder):
        if file_name.lower().endswith('.pdb'):
            file_path = os.path.join(folder, file_name)
            try:
                renumber_pdb_file(file_path, starting_resid)
            except Exception as e:
                print(f"[reres] Error processing '{file_path}': {e}")

def combine_subchains(folder_path: str):
    """
    For each PDB in folder: change every chain with an ID other than 'A' or 'B' to 'B'.
    Saves in-place (overwrites), using Biopython parsing (robust to formatting).
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    print(f"[subchains] combining non-A/B → B in {folder_path}")
    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(".pdb"):
            continue
        file_path = os.path.join(folder_path, filename)
        try:
            structure = parser.get_structure(filename, file_path)
            for model in structure:
                for chain in model:
                    if chain.id not in ("A", "B"):
                        chain.id = "B"
            io.set_structure(structure)
            io.save(file_path)
        except Exception as e:
            print(f"[subchains] WARN {filename}: {e}")

def postprocess_pdb_dir(pdb_dir: str, combine: bool):
    """
    New pipeline:
      (optional) combine_subchains → continuous residue renumber across chains → clean TERs → renumber atom serials
    """
    
    if not os.path.isdir(pdb_dir):
        return
    if combine:
        # collapse non-A/B to B, twice (mirrors your shell, helps remove internal TERs after merges)
        reres(pdb_dir)
        combine_subchains(pdb_dir)
        combine_subchains(pdb_dir)
        reres(pdb_dir)


# =========================
# ColabFold metrics
# =========================

def extract_colab_metrics(run_csv: str, colab_dir: str) -> pd.DataFrame:
    """
    Expects ColabFold outputs in:
      {colab_dir}/ptm_output  (json,pdb)

    Produces:
      - mean/std over actifptm/ptm/iptm across hits per binder_id
      - model_0 raw values
      - copies best rank_001 PDB to {colab_dir}/pdbs/{binder_id}.pdb
      - postprocesses PDBs: renumber + combine_subchains + renumber
    """
    ptm_folder    = os.path.join(colab_dir, "ptm_output")
    best_pdbs_dir = os.path.join(colab_dir, "pdbs")
    _ensure_dir(best_pdbs_dir)

    df = pd.read_csv(run_csv)
    if "binder_id" not in df.columns:
        raise ValueError("run.csv must contain 'binder_id'.")

    fieldnames = [
        "binder_id",
        "colab_actifptm_avg", "colab_actifptm_std", "colab_actifptm_model_0",
        "colab_ptm_avg",      "colab_ptm_std",      "colab_ptm_model_0",
        "colab_iptm_avg",     "colab_iptm_std",     "colab_iptm_model_0",
    ]
    updated = []

    for bid in df["binder_id"]:
        if not os.path.isdir(ptm_folder):
            print(f"[colab] Missing ptm_output dir: {ptm_folder}")
            break

        hits = [
            f for f in os.listdir(ptm_folder)
            if f.startswith(bid + "_") and f.endswith(".json")
               and "predicted_aligned_error" not in f
        ]

        all_actif, all_ptm, all_iptm = [], [], []
        best_iptm = -np.inf
        best_hit  = None
        model0_hit = None

        for fname in hits:
            path = os.path.join(ptm_folder, fname)
            try:
                data = json.load(open(path))
            except Exception as e:
                print(f"[colab] WARN cannot read {path}: {e}")
                continue

            a = data.get("actifptm")
            p = data.get("ptm")
            i = data.get("iptm")

            if a is not None: all_actif.append(a)
            if p is not None: all_ptm.append(p)
            if i is not None:
                all_iptm.append(i)
                if i > best_iptm:
                    best_iptm = i
                    best_hit  = fname

            if "rank_001_" in fname:
                model0_hit = fname

        def mean_std(vals):
            if not vals: return (None, None)
            if len(vals) == 1: return (float(vals[0]), 0.0)
            return (float(statistics.mean(vals)), float(statistics.stdev(vals)))

        mean_a, std_a = mean_std(all_actif)
        mean_p, std_p = mean_std(all_ptm)
        mean_i, std_i = mean_std(all_iptm)

        m0_a = m0_p = m0_i = None
        if model0_hit:
            try:
                data0 = json.load(open(os.path.join(ptm_folder, model0_hit)))
                m0_a = data0.get("actifptm")
                m0_p = data0.get("ptm")
                m0_i = data0.get("iptm")
                # copy rank_001 pdb as representative
                pdb_candidates = [
                    fname for fname in os.listdir(ptm_folder)
                    if fname.endswith(".pdb") and bid in fname and "rank_001" in fname
                ]
                if pdb_candidates:
                    src = os.path.join(ptm_folder, pdb_candidates[0])
                    dst = os.path.join(best_pdbs_dir, f"{bid}.pdb")
                    shutil.copy(src, dst)
            except Exception as e:
                print(f"[colab] WARN reading model_0 for {bid}: {e}")

        updated.append({
            "binder_id": bid,
            "colab_actifptm_avg":     mean_a,
            "colab_actifptm_std":     std_a,
            "colab_actifptm_model_0": m0_a,
            "colab_ptm_avg":          mean_p,
            "colab_ptm_std":          std_p,
            "colab_ptm_model_0":      m0_p,
            "colab_iptm_avg":         mean_i,
            "colab_iptm_std":         std_i,
            "colab_iptm_model_0":     m0_i,
        })

    # postprocess PDBs
    postprocess_pdb_dir(best_pdbs_dir, combine=True)

    df_out = pd.DataFrame(updated)
    # also drop a helper CSV
    out_csv = os.path.join(colab_dir, "colab_metrics.csv")
    if not df_out.empty:
        df_out.to_csv(out_csv, index=False)
        print(f"[colab] wrote {out_csv}")
    return df_out

# =========================
# Boltz1 metrics
# =========================

def convert_cif_to_pdb(cif_path, pdb_path, model_id=None):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(model_id or "model", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

def extract_boltz1_metrics(predictions_folder: str, boltz_dir: str) -> pd.DataFrame:
    """
    predictions_folder should look like:
      Boltz/outputs/<binder_id>/predictions/*_model_*.json (and *_model_0.cif)
    """
    rows = []
    pdb_dir = os.path.join(boltz_dir, "pdbs")
    _ensure_dir(pdb_dir)

    metrics = ["iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]
    
    predictions_folder = os.path.join(predictions_folder, "predictions")

    if not os.path.isdir(predictions_folder):
        print(f"[boltz1] missing predictions: {predictions_folder}")
        return pd.DataFrame(columns=["binder_id"])

    for binder_id in sorted(os.listdir(predictions_folder)):
        subdir = os.path.join(predictions_folder, binder_id)
        if not os.path.isdir(subdir):
            continue

        all_vals = {m: [] for m in metrics}
        model0_vals = {m: None for m in metrics}

        for fname in os.listdir(subdir):
            if not fname.endswith(".json"):
                continue
            if not re.match(r".*_model_\d+\.json$", fname):
                continue

            path = os.path.join(subdir, fname)
            try:
                data = json.load(open(path))
            except Exception as e:
                print(f"[boltz1] WARN parse {path}: {e}")
                continue

            for m in metrics:
                v = data.get(m)
                if v is not None:
                    all_vals[m].append(v)

            if fname.endswith("_model_0.json"):
                for m in metrics:
                    model0_vals[m] = data.get(m)

        def mean_std(lst):
            if not lst:
                return None, None
            arr = np.array(lst, dtype=float)
            return float(arr.mean()), float(arr.std(ddof=0))

        row = {"binder_id": binder_id}
        for m in metrics:
            avg, std = mean_std(all_vals[m])
            row[f"boltz1_{m}_avg"]     = avg
            row[f"boltz1_{m}_std"]     = std
            row[f"boltz1_{m}_model_0"] = model0_vals[m]
        rows.append(row)

        cif_path = os.path.join(subdir, f"{binder_id}_model_0.cif")
        if os.path.isfile(cif_path):
            out_pdb = os.path.join(pdb_dir, f"{binder_id}.pdb")
            try:
                convert_cif_to_pdb(cif_path, out_pdb, binder_id)
            except Exception as e:
                print(f"[boltz1] WARN CIF→PDB {binder_id}: {e}")

    # postprocess PDBs (combine chains too)
    postprocess_pdb_dir(pdb_dir, combine=True)

    df_out = pd.DataFrame(rows)
    out_csv = os.path.join(boltz_dir, "boltz1_metrics.csv")
    if not df_out.empty:
        df_out.to_csv(out_csv, index=False)
        print(f"[boltz1] wrote {out_csv}")
    return df_out

# =========================
# AF3 metrics
# =========================

def summarize_af3(af3_outputs_root: str, run_csv: str, pdb_out_dir: str) -> pd.DataFrame:
    """
    AF3/outputs structure:
      AF3/outputs/<binder_id_lower>/... summary_confidences.json + seed*/...json + *.cif

    Returns a DataFrame with:
      af3_iptm_model_0, af3_ptm_model_0, af3_iptm_avg/std, af3_ptm_avg/std
    """
    df_input = pd.read_csv(run_csv, dtype=str)
    if 'binder_id' not in df_input.columns:
        raise ValueError("run.csv must contain 'binder_id'")
    mapping = {bid.lower(): bid for bid in df_input['binder_id']}

    _ensure_dir(pdb_out_dir)
    rows = []

    if not os.path.isdir(af3_outputs_root):
        print(f"[af3] missing outputs dir: {af3_outputs_root}")
        return pd.DataFrame(columns=["binder_id"])

    for folder in sorted(glob.glob(os.path.join(af3_outputs_root, "*"))):
        if not os.path.isdir(folder):
            continue
        folder_name = os.path.basename(folder)  # usually lower-case id
        true_id = mapping.get(folder_name)
        if true_id is None:
            # allow also exact name
            true_id = mapping.get(folder_name.lower())
        if true_id is None:
            print(f"[af3] skip '{folder_name}' (not in run.csv)")
            continue

        top_summ = os.path.join(folder, f"{folder_name}_summary_confidences.json")
        if not os.path.isfile(top_summ):
            print(f"[af3] missing {top_summ}")
            continue
        data0 = json.load(open(top_summ))
        iptm0 = data0.get("iptm", np.nan)
        ptm0  = data0.get("ptm",  np.nan)

        seed_iptms, seed_ptms = [], []
        for seed_dir in glob.glob(os.path.join(folder, "seed*")):
            for sf in glob.glob(os.path.join(seed_dir, "*summary_confidences.json")):
                d = json.load(open(sf))
                if "iptm" in d: seed_iptms.append(d["iptm"])
                if "ptm"  in d: seed_ptms.append(d["ptm"])

        def stats(arr):
            if not arr: return (np.nan, np.nan)
            avg = round(float(np.mean(arr)), 2)
            std = round(float(np.std(arr, ddof=1)), 2) if len(arr) > 1 else 0.0
            return (avg, std)

        iptm_avg, iptm_std = stats(seed_iptms)
        ptm_avg,  ptm_std  = stats(seed_ptms)

        # convert first CIF to PDB
        cif_list = glob.glob(os.path.join(folder, "*.cif"))
        if cif_list:
            cif_path = cif_list[0]
            pdb_path = os.path.join(pdb_out_dir, f"{true_id}.pdb")
            try:
                convert_cif_to_pdb(cif_path, pdb_path, true_id)
            except Exception as e:
                print(f"[af3] WARN CIF→PDB {true_id}: {e}")

        rows.append({
            "binder_id":        true_id,
            "af3_iptm_model_0": iptm0,
            "af3_ptm_model_0":  ptm0,
            "af3_iptm_avg":     iptm_avg,
            "af3_iptm_std":     iptm_std,
            "af3_ptm_avg":      ptm_avg,
            "af3_ptm_std":      ptm_std
        })

    # postprocess PDBs (combine chains too)
    postprocess_pdb_dir(pdb_out_dir, combine=True)

    df_out = pd.DataFrame(rows).sort_values("binder_id")
    out_csv = os.path.join(os.path.dirname(pdb_out_dir), "af3_metrics.csv")
    if not df_out.empty:
        df_out.to_csv(out_csv, index=False)
        print(f"[af3] wrote {out_csv}")
    return df_out

# =========================
# AF2 metrics
# =========================

def _clean_stem(name: str) -> str:
    """Remove _af2pred and _relaxed from a filename stem."""
    stem = os.path.splitext(os.path.basename(name))[0]
    stem = re.sub(r"_af2pred$", "", stem)   # strip if at end
    stem = re.sub(r"_relaxed$", "", stem)   # strip if at end
    stem = re.sub(r"_relaxed", "", stem)    # also remove if in middle
    return stem

def extract_af2_metrics(af2_dir: str) -> pd.DataFrame:
    """
    Expects:
      AF2/outputs/out.sc    (Rosetta-style whitespace table; last field is decoy tag)
      AF2/outputs/*.pdb     (may end with _af2pred or _relaxed)

    Produces:
      - AF2/AF2_metrics.csv
      - renames/moves PDBs to AF2/pdbs/{binder_id}.pdb (strip _af2pred / _relaxed)
      - postprocesses PDBs: renumber (no chain combining)

    Returns:
      DataFrame keyed by binder_id with af2_* columns taken from out.sc.
    """
    outputs = os.path.join(af2_dir, "outputs")
    if not os.path.isdir(outputs):
        print(f"[af2] missing outputs dir: {outputs}")
        return pd.DataFrame(columns=["binder_id"])

    # Convert out.sc to CSV
    sc_path = os.path.join(af2_dir, "out.sc")
    df_metrics = pd.DataFrame(columns=["binder_id"])
    if os.path.isfile(sc_path):
        rows = []
        with open(sc_path, "r") as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if lines:
            header = re.split(r"\s+", re.sub(r"^SCORE:\s*", "", lines[0].strip()))
            for ln in lines[1:]:
                ln = re.sub(r"^SCORE:\s*", "", ln.strip())
                parts = re.split(r"\s+", ln)
                if len(parts) != len(header):
                    continue
                row = dict(zip(header, parts))
                tag = row.get(header[-1], "")
                bid = _clean_stem(tag)
                row_out = {"binder_id": bid}
                for k, v in row.items():
                    if k == header[-1]:
                        continue
                    try:
                        row_out[f"af2_{k}"] = float(v)
                    except Exception:
                        row_out[f"af2_{k}"] = v
                rows.append(row_out)
        if rows:
            df_metrics = pd.DataFrame(rows).drop_duplicates("binder_id")
    else:
        print(f"[af2] WARNING: out.sc not found: {sc_path}")

    # Move/rename PDBs
    pdbs_out = os.path.join(af2_dir, "pdbs")
    _ensure_dir(pdbs_out)

    for fp in glob.glob(os.path.join(outputs, "*.pdb")):
        stem = _clean_stem(fp)
        dst = os.path.join(pdbs_out, f"{stem}.pdb")
        try:
            shutil.move(fp, dst)
        except Exception:
            try:
                shutil.copy(fp, dst)
            except Exception as e:
                print(f"[af2] WARN moving {fp}: {e}")

    # postprocess PDBs (renumber only)
    postprocess_pdb_dir(pdbs_out, combine=True)

    out_csv = os.path.join(af2_dir, "AF2_metrics.csv")
    if not df_metrics.empty:
        df_metrics.to_csv(out_csv, index=False)
        print(f"[af2] wrote {out_csv}")
    return df_metrics

# =========================
# CLI Orchestrator
# =========================

def parse_args():
    p = argparse.ArgumentParser(
        description="Combine model metrics into run.csv; choose models with --models"
    )
    p.add_argument("--run_csv", required=True, help="Path to run.csv (will be updated in-place)")
    p.add_argument("--output_dir", required=True, help="Pipeline output base directory")
    p.add_argument("--models", default="colab,boltz1,af3,af2",
                   help="Comma-separated list among: colab,boltz1,af3,af2")
    # Optional overrides (defaults derive from output_dir)
    p.add_argument("--colab_dir", default=None, help="Override ColabFold base dir (default: <output_dir>/ColabFold)")
    p.add_argument("--boltz_dir", default=None, help="Override Boltz1 base dir (default: <output_dir>/Boltz)")
    p.add_argument("--af3_dir", default=None, help="Override AF3 base dir (default: <output_dir>/AF3)")
    p.add_argument("--af2_dir", default=None, help="Override AF2 base dir (default: <output_dir>/AF2)")
    return p.parse_args()

def main():
    args = parse_args()
    outdir = os.path.abspath(args.output_dir)
    models = [m.strip().lower() for m in args.models.split(",") if m.strip()]

    # Resolve default model dirs
    colab_dir = args.colab_dir or os.path.join(outdir, "ColabFold")
    boltz_dir = args.boltz_dir or os.path.join(outdir, "Boltz")
    af3_dir   = args.af3_dir   or os.path.join(outdir, "AF3")
    af2_dir   = args.af2_dir   or os.path.join(outdir, "AF2")

    print(f"[info] run_csv:   {args.run_csv}")
    print(f"[info] output_dir:{outdir}")
    print(f"[info] models:    {models}")

    # COLABFOLD
    if "colab" in models:
        df_colab = extract_colab_metrics(args.run_csv, colab_dir)
        if not df_colab.empty:
            _merge_into_run(args.run_csv, df_colab)

    # BOLTZ1
    if "boltz1" in models:
        predictions_folder = os.path.join(boltz_dir, "boltz_results_input_folder")
        df_boltz = extract_boltz1_metrics(predictions_folder, boltz_dir)
        if not df_boltz.empty:
            _merge_into_run(args.run_csv, df_boltz)

    # AF3
    if "af3" in models:
        af3_outputs_root = os.path.join(af3_dir, "outputs")
        af3_pdbs_dir     = os.path.join(af3_dir, "pdbs")
        df_af3 = summarize_af3(af3_outputs_root, args.run_csv, af3_pdbs_dir)
        if not df_af3.empty:
            _merge_into_run(args.run_csv, df_af3)

    # AF2
    if "af2" in models:
        df_af2 = extract_af2_metrics(af2_dir)
        if not df_af2.empty:
            _merge_into_run(args.run_csv, df_af2)

    print("[done] all requested models processed.")

if __name__ == "__main__":
    main()
