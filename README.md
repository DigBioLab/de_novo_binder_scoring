# de_novo_binder_scoring

This repository contains the scripts and analysis described in:  
[**Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders**](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1)

---

## Installation

Clone and set up the environment:

```bash
git clone https://github.com/DigBioLab/de_novo_binder_scoring.git
cd de_novo_binder_scoring

conda env create -f environment.yml
conda activate binder_scoring_env
chmod +x ./functions/DAlphaBall.gcc
```

**Note 1**: This repo uses [PyRosetta](https://www.pyrosetta.org/downloads) for two scripts (`compute_rosetta_metrics.py` and `rmsd.py`), which requires a license for commercial use.

**Note 2**: The used structure predictions tools (AF2 initial guess, ColabFold, Boltz and AF3) require seperate installations - these blocks are marked in the `example_run.sh` script.

---

## Usage


## Usage

### 1. Process inputs

Convert input PDBs into standardized inputs (`run.csv`, cleaned PDBs, and MSA FASTAs):

```bash
python process_inputs.py \
  --input_pdbs ./example_input/input_pdbs \
  --output_dir ./example_output
```

* Binder is expected as **chain A** (`A:no_msa` by default).
* Non-A chains are merged into **B** in PDB for downstream analysis (also works if target chains of input are already merged).
* Unique target sequences will get target IDs (`target_1`, `target_2`, …).
* Outputs: `run.csv`, `Binder_seq.fasta`, and `unique_msa/`.

#### Sanitation of names

All names are automatically sanitized: lowercase letters, underscores, and numbers are allowed; all non-alphanumerics are replaced with `_`.
Be careful that your input files do not sanitize to the same name (e.g., `abc.pdb`, `AbC.pdb`, and `abC.pdb` all become `abc.pdb`).

#### Input Modes

You can overwrite columns in `run.csv` for customization using the `--mode` flag:

```bash
python process_inputs.py --mode {pdb_only, seq_only_csv, hybrid}
```

* `pdb_only` (default): sequences and other columns are automatically inferred from input PDBs.
* `seq_only_csv`: sequences and columns are taken only from a CSV. Useful if you do not have input PDBs.

Example:

```bash
python process_inputs.py \
  --mode seq_only_csv \
  --input_csv ./example_input/input_sequence_only.csv \
  --output_dir ./example_output_seq
```

* `hybrid`: sequences are extracted from PDBs by default, but can be overwritten with CSV-specified sequences. To overwrite, the CSV must specify:

  * `target_chains` – all target chains (including those not overwritten)
  * `target_subchain_X_seq` – one column per chain to overwrite (e.g., `target_subchain_D_seq`)

Example:

```bash
python process_inputs.py \
  --mode hybrid \
  --input_pdbs ./example_input/input_pdbs \
  --input_csv ./example_input/input_overwrite.csv \
  --output_dir ./example_output_overwrite
```

* **Logging behavior**: overwritten sequences are stored in `pdb_extracted_trg_subch_{X}_not_used` to preserve original PDB info.

---

### 2. Generate MSAs

Generate MSAs using ColabFold (MMseqs2). Requires a separate ColabFold installation:

```bash
colabfold_batch ./example_outputs/unique_msa ./example_outputs/unique_msa/msa --msa-only
```

---

### 3. Prepare model inputs

Generate inputs for structure prediction models (AF3, Boltz, ColabFold):

```bash
python generate_model_inputs.py \
  --run-csv ./example_outputs/run.csv \
  --out-dir ./example_outputs
```

---

### 4. Relax input structures & compute Rosetta metrics

```bash
python compute_rosetta_metrics.py \
  --run-csv ./example_outputs/run.csv \
  --out-csv ./example_outputs/input_rosetta_metrics.csv \
  --folder input:./example_outputs/input_pdbs
```

---

### 5. AF2 initial guess

Run AF2 prediction on relaxed PDBs:

```bash
predict.py \
  -pdbdir ./example_outputs/input_pdbs/relaxed_pdbs \
  -scorefilename out.sc \
  -outsilent af2.silent
```

---

### 6. Run ColabFold

```bash
colabfold_batch ./example_outputs/ColabFold/input_folder ./example_outputs/ColabFold/ptm_output \
  --calc-extra-ptm --num-recycle 3 --num-models 3
```

Optional: remove generated PNGs:

```bash
find ./example_outputs/ColabFold/ptm_output -type f -name "*.png" -exec rm -f {} \;
```

---

### 7. Run Boltz

```bash
boltz predict ./example_outputs/Boltz/input_folder \
  --recycling_steps 10 \
  --diffusion_samples 3 \
  --write_full_pae \
  --out_dir ./example_outputs/Boltz
```

---

### 8. Run AF3

```bash
python run_alphafold.py \
  --input_dir=./example_outputs/AF3/input_folder \
  --model_dir=/path/to/alphafold3_weights \
  --db_dir=/path/to/alphafold3_database \
  --run_data_pipeline=False \
  --num_diffusion_samples=3 \
  --output_dir=./example_outputs/AF3/outputs
```

---

### 9. Extract confidence metrics

```bash
python extract_confidence_metrics.py \
  --run-csv ./example_outputs/run.csv \
  --out-dir ./example_outputs
```

---

### 10. Compute ipSAE and interface confidence metrics

```bash
python run_ipsae_batch.py \
  --run-csv ./example_outputs/run.csv \
  --out-csv ./example_outputs/ipsae_and_ipae.csv \
  --af3-dir ./example_outputs/AF3 \
  --boltz1-dir ./example_outputs/Boltz \
  --colab-dir ./example_outputs/ColabFold \
  --ipsae-script-path ./ipsae_w_ipae.py \
  --pae-cutoff 10 \
  --dist-cutoff 10 \
  --backup
```

---

### 11. Compute DockQ

```bash
python dockQ.py \
  --run-csv ./example_outputs/run.csv \
  --input-pdbs ./example_outputs/input_pdbs/ \
  --folder af3:./example_outputs/AF3/pdbs/ \
  --folder af2:./example_outputs/AF2/pdbs/ \
  --folder boltz:./example_outputs/Boltz/pdbs/ \
  --folder colab:./example_outputs/ColabFold/pdbs/ \
  --out-csv ./example_outputs/dockQ.csv \
  --backup \
  --verbose
```

---

### 12. Compute Rosetta metrics for model PDBs

```bash
python compute_rosetta_metrics.py \
  --run-csv ./example_outputs/run.csv \
  --out-csv ./example_outputs/rosetta_metrics.csv \
  --folder af3:./example_outputs/AF3/pdbs/ \
  --folder af2:./example_outputs/AF2/pdbs/ \
  --folder boltz1:./example_outputs/Boltz/pdbs/ \
  --folder colab:./example_outputs/ColabFold/pdbs/ \
  --folder input:./example_outputs/input_pdbs/
```

---

### 13. Compute RMSDs

```bash
python rmsd.py \
  --folder input:./example_outputs/input_pdbs/ \
  --folder af3:./example_outputs/AF3/pdbs/ \
  --folder af2:./example_outputs/AF2/pdbs/ \
  --folder boltz1:./example_outputs/Boltz/pdbs/ \
  --folder colab:./example_outputs/ColabFold/pdbs/ \
  --out-csv ./example_outputs/rmsd.csv
```

---

### 14. Compute PyMOL metrics

Requires the open-source PyMOL installation:

```bash
OUTPUT_DIR="$(pwd)/outputs"
PYMOL_DIR=$OUTPUT_DIR/pymol_files
mkdir -p "${PYMOL_DIR}"

# Create JSON file listing PDB directories
echo '{"input": "'$OUTPUT_DIR/input_pdbs'", "af2": "'$OUTPUT_DIR/AF2/pdbs'", "colab": "'$OUTPUT_DIR/ColabFold/pdbs'", "boltz1": "'$OUTPUT_DIR/Boltz/pdbs'", "af3": "'$OUTPUT_DIR/AF3/pdbs'"}' > "${PYMOL_DIR}/pdb_dirs.json"

# Run PyMOL analysis script
cd $OUTPUT_DIR
python -m pymol -c -d "run ../pymol_metrics.py"
```

---

### Full workflow example

See **`example_run.sh`** for a complete pipeline example including environment loading/unloading.

---



## Citation

If you use this code, please cite:  

**Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders**. *bioRxiv* (2025).  
DOI: [10.1101/2025.08.14.670059v1](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1)

---

### Additional citations

If you use any of the following tools or methods, please also cite:  

- **ColabFold (MSA generation and/or structure prediction)**  
  [10.1038/s41592-022-01488-1](https://www.nature.com/articles/s41592-022-01488-1)  

- **AF2 initial guess**  
  [https://doi.org/10.1038/s41467-023-38328-5](https://doi.org/10.1038/s41467-023-38328-5)  

- **Boltz-1**  
  [https://www.biorxiv.org/content/10.1101/2024.11.19.624167v4](https://www.biorxiv.org/content/10.1101/2024.11.19.624167v4)  

- **Boltz-2**  
  [https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1](https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1)  

- **AF3**  
  [https://doi.org/10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w)  

- **ipSAE**  
  [https://doi.org/10.1101/2025.02.10.637595](https://doi.org/10.1101/2025.02.10.637595)  

- **DockQ**  
  [https://doi.org/10.1093/bioinformatics/btae586](https://doi.org/10.1093/bioinformatics/btae586)  



