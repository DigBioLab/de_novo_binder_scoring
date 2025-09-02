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

⚠️ **Note**: This repo requires [PyRosetta](https://www.pyrosetta.org/downloads), which requires a license for commercial use.

⚠️ **Note**: This repo requires enviroments and licenses to run the structure predictions tools - these blocks are marked in the `example_run.sh` script.

---

## Usage

### 1. Process inputs

Convert PDBs into standardized inputs (`run.csv`, cleaned PDBs, and MSA FASTAs):

```bash
python process_inputs.py --input_pdbs ./example_input/input_pdbs   --output_dir ./outputs
```

- Binder is expected as **chain A** (`A:no_msa` by default).  
- Non-A chains are merged into **B** (or split into virtual B, C, D … if chain breaks are detected).  
- Unique target sequences are deduplicated (`target_1`, `target_2`, …).  
- Outputs: `run.csv`, `Binder_seq.fasta`, and `unique_msa/`.

#### Input Modes

You can control how sequences are taken into the run CSV with the `--mode` flag:

```bash
python process_inputs.py --mode {hybrid,pdb_only,seq_only_csv}
```

-----

`pdb_only` (default)

Sequences are extracted **only from the PDB files**.
The run CSV will be updated to reflect what is found in the structures.
Any sequences specified in the CSV are ignored.

-----

`seq_only_csv`

Sequences are taken **only from the CSV**.
Useful if you want full control over sequences (e.g., curated or designed sequences).

-----

`hybrid` (possibility of overwriting PDB extracted sequence of target chains)

In this mode, by default, sequences are extracted from the PDBs.
However one can also ***overwrite the PDB extracted seqeunces***, as well as any other columns. To overwrite the PDB extracted sequences with custom sequence, the input CSV must specify both:

  * `target_chains` please specify all target chains (also ones that should not be overwritten)
  * `target_subchain_X_seq` (one column per subchain that one wants to overwrite, e.g., `target_subchain_D_seq`)

Then the CSV-provided sequence will overwrite the extracted version for that chain.

**Additional behavior in hybrid mode:**

  * **Logging extracted sequences**
    When a CSV-provided sequence overwrites a PDB-extracted one, the extracted sequence is still stored in a dedicated column: `pdb_extracted_trg_subch_{X}_not_used`.
    This ensures transparency: you can always compare the provided sequence against what was found in the structure.

<!-- end list -->

---

### 2. Generate MSAs

We used the ColabFold server (MMseqs2) for MSAs. Example (local ColabFold):

```bash
colabfold_batch ./outputs/unique_msa ./outputs/unique_msa/msa --msa-only
```

---

### 3. Prepare model inputs

Generate input files for structure prediction models (AF3, Boltz-1/2, ColabFold):

```bash
python generate_model_inputs.py   --csv_file ./outputs/run.csv   --out_dir ./outputs   --models af3 boltz colabfold
```

The individual structure prediction tools need to be installed independenly, an example how they can be run can be found in the example_run.sh

For AF2 initial guess the pdbs have the correct format to use them directly as an input as specified here https://github.com/nrbennet/dl_binder_design.

---

### 4. Process outputs

#### Extract model metrics
```bash
python extract_confidence_metrics.py
--run_csv ./outputs/run.csv \
--output_dir ./outputs \
--colab_dir ./outputs/ColabFold \
--boltz_dir ./outputs/Boltz \
--af3_dir ./outputs/AF3 \
--af2_dir ./outputs/AF2
```

#### Compute ipSAE scores ([IPSAE](https://github.com/DunbrackLab/IPSAE)):
```bash
python run_ipsae_batch.py
--run-csv ./outputs/run.csv \
--af3-dir ./outputs/AF3 \
--boltz-dir ./outputs/Boltz \
--colab-dir ./outputs/ColabFold
```

#### Compute Rosetta metrics ([BindCraft-inspired](https://github.com/martinpacesa/BindCraft)):
```bash
python compute_rosetta_metrics.py
--run_csv ./outputs/run.csv \
--out_csv ./outputs/rosetta_metrics.csv \
--input input=./outputs/input_pdbs \
--input af3=./outputs/AF3/pdbs \
--input boltz=./outputs/Boltz/pdbs \
--input colab=./outputs/ColabFold/pdbs \
--input af2=./outputs/AF2/pdbs
```

#### Compute RMSD (requires at least two sets of PDBs):
```bash
python rmsd.py \ 
--folder input:./outputs/input_pdbs \
--folder af3:./outputs/AF3/pdbs \
--folder af2:./outputs/AF2/pdbs \
--folder boltz:./outputs/Boltz/pdbs \
--folder colab:./outputs/ColabFold/pdbs \
--out-csv ./outputs/rmsd.csv
```


#### Compute DockQ between input and predicted structure (requires at least two sets of PDBs):
```bash
python dockQ.py \
  --input_pdbs ./outputs/input_pdbs/ \
  --model af3:./outputs/AF3/pdbs/ \
  --model af2:./outputs/AF2/pdbs/ \
  --model boltz1:./outputs/Boltz/pdbs \
  --model colab:./outputs/ColabFold/pdbs \
  --output_csv ./outputs/dockQ.csv
```


#### Compute pymol metrics:
This will require separate installation of the opensource pymol package: https://github.com/schrodinger/pymol-open-source

Using the pymol environment navigate to de_novo_binder_scoring dir and then execute the following.
```bash
OUTPUT_DIR="$(pwd)/outputs"
PYMOL_DIR=$OUTPUT_DIR/pymol_files
mkdir -p "${PYMOL_DIR}"

# Create a JSON file that lists the directories for Pymol analysis.
echo '{"input": "'$OUTPUT_DIR/input_pdbs'", "af2": "'$OUTPUT_DIR/AF2/pdbs'", "colab": "'$OUTPUT_DIR/ColabFold/pdbs'", "boltz1": "'$OUTPUT_DIR/Boltz/pdbs'", "af3": "'$OUTPUT_DIR/AF3/pdbs'"}' > "${PYMOL_DIR}/pdb_dirs.json"

cd $OUTPUT_DIR
# Run Pymol in command-line mode to execute the analysis script.
python -m pymol -c -d "run ../pymol_metrics.py"
```
---
## Example run

A full workflow is provided in **`example_run.sh`**.

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



