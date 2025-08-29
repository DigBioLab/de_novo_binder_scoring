# de_novo_binder_scoring

This repo contains the scripts and analysis described in:  
[**Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders**](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1)

---

## Installation

```bash
git clone https://github.com/DigBioLab/de_novo_binder_scoring.git
cd de_novo_binder_scoring

conda env create -f environment.yml
conda activate binder_scoring_env
chmod +x ./functions/DAlphaBall.gcc
```

> **Note:** PyRosetta is optional but supported. It requires registration/licensing at [pyrosetta.org](https://www.pyrosetta.org/downloads).

---

## Usage

Process a folder of PDBs into standardized inputs (`run.csv`, MSAs, cleaned PDBs):

```bash
python process_inputs.py   --input_pdbs ./example_input/input_pdbs   --output_dir ./outputs
```

**Logic:**
- Binder is expected on **chain A** (`A:no_msa` by default).  
- Targets = chain B or merged non-A chains; chain breaks are split into virtual segments (B, C, D …).  
- Unique target sequences are grouped (`target_1`, `target_2`, …).  
- PDBs are renumbered in-place; MSAs must be generated separately (e.g. ColabFold server or local MMseqs2).

---

## Citation

If you use this code, please cite:  

**Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders**. *bioRxiv* (2025).  
DOI: [10.1101/2025.08.14.670059v1](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1)

