#!/usr/bin/env python3
import os
import re
import json
import glob
import shutil
import argparse
from typing import Dict, List, Tuple
from collections import Counter

import pandas as pd

# -------------------------
# Helpers
# -------------------------

def _read_csv(csv_file: str) -> pd.DataFrame:
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV not found: {csv_file}")
    return pd.read_csv(csv_file)

def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def _safe_str(x) -> str:
    return "" if pd.isna(x) or x is None else str(x)


def _chains_from_row(row: pd.Series) -> List[str]:
    """
    Deduce chain letters to process for this row:
    - binder_chain first (if present)
    - then chains listed in target_chains (if JSON-parsable)
    - else, any chains implied by 'msa_path_<X>' or 'target_subchain_<X>_seq' columns.
    """
    seen = []
    bc = str(row.get("binder_chain", "")).strip()
    if bc:
        seen.append(bc)

    tchains = []
    tc = row.get("target_chains", "")
    try:
        tchains = [c.strip() for c in json.loads(tc) if isinstance(c, str) and c.strip()]
    except Exception:
        # Fall back: infer from columns
        for c in row.index:
            m1 = re.match(r"msa_path_([A-Za-z])$", c)
            m2 = re.match(r"target_subchain_([A-Za-z])_seq$", c)
            if m1:
                tchains.append(m1.group(1))
            if m2:
                tchains.append(m2.group(1))

    # preserve order, de-dup, keep binder first
    for c in tchains:
        if c not in seen:
            seen.append(c)
    return seen

# -------------------------
# ColabFold inputs (renamed from reconstruct_complex_a3m)
# -------------------------

def write_colabfold_inputs(csv_file: str, outpath: str) -> None:
    """
    Builds complex-level A3Ms by stitching per-chain MSAs (or raw sequences when no_msa).
    - Dynamically handles chains (A/B/C/D/…).
    - If an .a3m path is missing, tries sibling FASTA, else falls back to the raw sequence.
    - Pads each chain block with '-' to align into the concatenated complex.
    """
    _ensure_dir(outpath)
    df = _read_csv(csv_file)
    csv_dir = os.path.dirname(os.path.abspath(csv_file))

    def _resolve_msa_path(row, chain: str) -> str:
        """Return an existing MSA-like file path for this chain, preferring .a3m, else FASTA fallback."""
        raw = _safe_str(row.get(f"msa_path_{chain}"))
        if not raw or raw.lower() == "no_msa":
            return ""
        msa_path = raw if os.path.isabs(raw) else os.path.normpath(os.path.join(csv_dir, raw))
        if os.path.isfile(msa_path):
            return msa_path
        # Fallback: .../unique_msa/msa/<name>.a3m  ->  .../unique_msa/<name>.fasta
        base = os.path.basename(msa_path)
        name, _ = os.path.splitext(base)
        parent = os.path.dirname(msa_path)
        grand = os.path.dirname(parent)
        fasta_guess = os.path.join(grand, f"{name}.fasta")
        return fasta_guess if os.path.isfile(fasta_guess) else ""

    for _, row in df.iterrows():
        binder_id = row["binder_id"]
        chain_lengths: List[int] = []
        chain_blocks: List[str] = []
        csv_chains: List[str] = []

        chain_counter = 101
        chains = _chains_from_row(row)
        if not chains:
            print(f"[colabfold] {binder_id}: no chains detected; skipping.")
            continue

        binder_chain = _safe_str(row.get("binder_chain"))

        for chain in chains:
            # choose sequence column
            if chain == binder_chain:
                seq_col = f"{binder_chain}_seq"
            else:
                tcol = f"target_subchain_{chain}_seq"
                if tcol in df.columns and _safe_str(row.get(tcol)):
                    seq_col = tcol
                elif f"{chain}_seq" in df.columns:
                    seq_col = f"{chain}_seq"
                else:
                    print(f"[colabfold] {binder_id}: no sequence column for chain {chain}; skipping chain.")
                    continue

            seq = _safe_str(row.get(seq_col))
            if not seq:
                print(f"[colabfold] {binder_id}: empty sequence for chain {chain}; skipping chain.")
                continue

            csv_chains.append(seq)
            chain_lengths.append(len(seq))

            msa_file = _resolve_msa_path(row, chain)
            block = ""
            if (not msa_file) or (os.path.splitext(msa_file)[1].lower() not in [".a3m", ".fasta", ".fa", ".faa"]):
                # No usable file -> single-sequence block
                block = f">{chain_counter}\n{seq}\n"
            else:
                with open(msa_file, "r") as f:
                    content = f.read()
                # Replace only the first header to set query id; drop comment lines
                lines = content.splitlines()
                replaced = False
                new_lines = []
                for ln in lines:
                    if not replaced and ln.startswith(">101"):
                        new_lines.append(f">{chain_counter}")
                        replaced = True
                    elif ln.startswith("#"):
                        continue
                    else:
                        new_lines.append(ln)
                if not replaced:
                    # set the first header if none matched >101
                    for i, ln in enumerate(new_lines):
                        if ln.startswith(">"):
                            new_lines[i] = f">{chain_counter}"
                            replaced = True
                            break
                    if not replaced:
                        new_lines = [f">{chain_counter}", seq]
                block = "\n".join(new_lines) + ("\n" if not new_lines or not new_lines[-1].endswith("\n") else "")

            chain_blocks.append(block.strip())
            chain_counter += 1

        if not chain_blocks:
            print(f"[colabfold] {binder_id}: no valid chain blocks; skipping.")
            continue

        # Header: "#<len1>,<len2>,...\t1,1,1,..."
        header_line = "#" + ",".join(str(l) for l in chain_lengths) + "\t" + ",".join(["1"] * len(chain_lengths))
        # Guide
        chain_numbers = [str(101 + i) for i in range(len(chain_lengths))]
        first_msa = ">101" + ("\t" + "\t".join(chain_numbers[1:]) if len(chain_numbers) > 1 else "") + "\n" + "".join(csv_chains)
        # Pad blocks
        padded_chain_blocks: List[str] = []
        total = len(chain_lengths)
        for i, block in enumerate(chain_blocks):
            prefix = "-" * sum(chain_lengths[:i])
            suffix = "-" * sum(chain_lengths[i+1:]) if i < total - 1 else ""
            new_lines = []
            for ln in block.splitlines():
                if ln.startswith(">"):
                    new_lines.append(ln)
                else:
                    new_lines.append(prefix + ln + suffix)
            padded_chain_blocks.append("\n".join(new_lines))

        out_txt = header_line + "\n" + first_msa + "\n" + "\n".join(padded_chain_blocks) + "\n"
        out_file = os.path.join(outpath, f"{binder_id}.a3m")
        with open(out_file, "w") as fh:
            fh.write(out_txt)

# -------------------------
# Boltz1 FASTA
# -------------------------

def write_boltz1_fastas(inpath: str, outpath: str):
    """
    Writes one FASTA per binder with entries:
      >{CHAIN}|protein|{msa_path or 'empty'}
      <sequence>
    Chains include binder_chain and any target subchains present.
    """
    _ensure_dir(outpath)
    df = _read_csv(inpath)

    for _, row in df.iterrows():
        binder_id = row["binder_id"]
        binder_chain = _safe_str(row.get("binder_chain"))
        chains = set(_chains_from_row(row))
        if binder_chain:
            chains.add(binder_chain)
        chains = sorted(chains)

        fasta_lines: List[str] = []
        for chain in chains:
            msa_path_col = f"msa_path_{chain}"
            # sequence selection
            if chain == binder_chain:
                seq = _safe_str(row.get(f"{binder_chain}_seq")) or _safe_str(row.get("A_seq"))
            else:
                seq = _safe_str(row.get(f"target_subchain_{chain}_seq")) or _safe_str(row.get(f"{chain}_seq"))
            if not seq:
                continue

            msa_path = _safe_str(row.get(msa_path_col))
            header = f">{chain}|protein|{'empty' if (not msa_path or msa_path.lower()=='no_msa') else msa_path}"
            fasta_lines.append(header)
            fasta_lines.append(seq)

        if not fasta_lines:
            print(f"[boltz] {binder_id}: no sequences, skipping.")
            continue

        out_file = os.path.join(outpath, f"{binder_id}.fasta")
        with open(out_file, "w") as fh:
            fh.write("\n".join(fasta_lines) + "\n")

# -------------------------
# AF3 JSON
# -------------------------

def get_predicted_binder_ids(outpath: str) -> Tuple[set, str]:
    """Find existing JSONs in outpath and return (ids, move_dir) to stash old ones.
       _previous_predictions is created one level above outpath (e.g. in AF3/)."""
    _ensure_dir(outpath)
    parent = os.path.dirname(outpath.rstrip(os.sep))  # go up one level
    prev_dir = os.path.join(parent, "_previous_predictions")
    _ensure_dir(prev_dir)

    ids = set()
    for p in glob.glob(os.path.join(outpath, "*.json")):
        stem = os.path.splitext(os.path.basename(p))[0]
        ids.add(stem)
    return ids, prev_dir

def write_af3_json(inpath: str, outpath: str):
    """
    Writes one JSON per binder for AF3 with:
      - binder as first sequence (binder_chain)
      - subsequent target subchains present in the CSV
      - MSA fields from msa_path_* if present, else empty fields.
      - optional ions from "ions_in_structure" column appended into "sequences"
        as objects of the form: {"ion": {"ion": "MG", "count": 2}}
    """
    df = _read_csv(inpath)

    ALLOWED_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}

    for _, row in df.iterrows():
        binder_id = row["binder_id"]

        binder_chain = _safe_str(row.get("binder_chain")) or "A"
        binder_seq = _safe_str(row.get(f"{binder_chain}_seq")) or _safe_str(row.get("A_seq"))
        if not binder_seq:
            print(f"[af3] {binder_id}: missing binder sequence; skipping.")
            continue

        sequences = [{
            "protein": {
                "id": binder_chain,
                "sequence": binder_seq,
                "templates": []
            }
        }]

        # add target subchains
        chains = _chains_from_row(row)
        for ch in chains:
            if ch == binder_chain:
                continue
            seq = _safe_str(row.get(f"target_subchain_{ch}_seq"))
            if not seq:
                continue
            sequences.append({
                "protein": {
                    "id": ch,
                    "sequence": seq,
                    "templates": []
                }
            })

         # determine last chain ID used
        all_chain_ids = [binder_chain] + [ch for ch in chains if ch != binder_chain]
        last_chain = max(all_chain_ids) if all_chain_ids else "A"
        start_ord = ord(last_chain) + 1

        # parse ions_in_target and append as ligands
        ions_raw = _safe_str(row.get("ions_in_target"))
        if ions_raw:
            try:
                ions_list = json.loads(ions_raw) if isinstance(ions_raw, str) else []
                if not isinstance(ions_list, list):
                    raise ValueError("ions_in_target is not a list")
            except Exception:
                print(f"[af3] {binder_id}: warning: could not parse ions_in_target: {ions_raw!r}; ignoring.")
                ions_list = []

            if ions_list:
                counts = Counter([str(x).upper().strip() for x in ions_list if x is not None])
                lig_id_ord = start_ord
                for ion_name, cnt in counts.items():
                    if ion_name in ALLOWED_IONS:
                        for _ in range(cnt):
                            lig_id = chr(lig_id_ord)
                            lig_id_ord += 1
                            sequences.append({
                                "ligand": {
                                    "id": lig_id,
                                    "ccdCodes": [ion_name]
                                }
                            })
                            print(f"[af3] {binder_id}: added ligand {ion_name} as ID {lig_id}")
                    else:
                        print(f"[af3] {binder_id}: warning: invalid/unsupported ion '{ion_name}' (ignored)")

        # attach MSA hints
        for seq in sequences:
            # ion entries don't have protein.id, handle accordingly
            if "protein" in seq:
                cid = seq["protein"]["id"]
                msa_col = f"msa_path_{cid}"
                msa_path = _safe_str(row.get(msa_col))
                if (not msa_path) or (msa_path.lower() == "no_msa"):
                    seq["protein"]["unpairedMsa"] = ""
                    seq["protein"]["pairedMsa"] = ""
                else:
                    seq["protein"]["unpairedMsaPath"] = msa_path
                    seq["protein"]["pairedMsa"] = ""

        json_data = {
            "name": binder_id,
            "sequences": sequences,
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 2
        }

        _ensure_dir(outpath)
        out_file = os.path.join(outpath, f"{binder_id}.json")
        with open(out_file, "w") as fh:
            json.dump(json_data, fh, indent=4)
        print(f"[af3] wrote: {out_file}")

# -------------------------
# CLI
# -------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate inputs for selected models: AF3, Boltz, ColabFold."
    )
    parser.add_argument(
        "--run-csv",
        required=True,
        dest="run_csv",
        help="Path to the CSV produced earlier (with *_seq and msa_path_* columns)."
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        dest="out_dir",
        help="Base output directory where model subfolders will be created."
    )
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["af3", "boltz", "colabfold"],
        help="Which models to generate inputs for. Default: all (af3 boltz colabfold)."
    )
    args = parser.parse_args()

    csv_file = args.run_csv
    out_base = args.out_dir
    models = set(args.models) if args.models else {"af3", "boltz", "colabfold"}

    print(f"=== Generate input files for selected models: {sorted(models)} ===")

    if "af3" in models:
        af3_dir = os.path.join(out_base, "AF3", "input_folder")
        _ensure_dir(af3_dir)
        write_af3_json(csv_file, af3_dir)

    if "boltz" in models:
        boltz_dir = os.path.join(out_base, "Boltz", "input_folder")
        _ensure_dir(boltz_dir)
        write_boltz1_fastas(csv_file, boltz_dir)

    if "colabfold" in models:
        cf_dir = os.path.join(out_base, "ColabFold", "input_folder")
        _ensure_dir(cf_dir)
        write_colabfold_inputs(csv_file, cf_dir)

if __name__ == "__main__":
    main()
