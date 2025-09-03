# ==============================================================================
# Variables & Directory Setup
# ==============================================================================

# >>>  Adjust ALL of the following paths to your environment <<<
# ==============================================================================
# Path to CONDA Install
CONDA_PATH="/path/to/conda/bin/activate"

# Output directory (must be full path)
OUTPUT_DIR="/path/to/outputs"

# Input directory containing PDB files
INPUT_PDBS="./example_input/input_pdbs"

# Path to this repository 
SCRIPT_DIR="/path/to/de_novo_binder_scoring"
# ==============================================================================

# Log directory setup
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "${LOG_DIR}"
OVERALL_START_TIME=$(date +%s)

# Count number of input PDB files 
count=$(ls -1 "${INPUT_PDBS}"/*.pdb | wc -l)
echo "Running ${count} structures..." > "${LOG_DIR}/log.txt"

# ==============================================================================
# 1. Pre-process input CSV and PDB Files
# ==============================================================================
# Activate the Python environment and run the meta_analysis script to convert PDB files to a CSV.
cd $SCRIPT_DIR

source "$CONDA_PATH"
conda activate binder_scoring_env

python process_inputs.py \
  --input_pdbs "${INPUT_PDBS}" \
  --output_dir "${OUTPUT_DIR}" 

# ==============================================================================
# 2. Relax input structures and compute Rosetta metrics
# ==============================================================================
START_TIME=$(date +%s)
echo -e "\nRelaxing input PDBs and computing Rosetta metrics" >> "${LOG_DIR}/log.txt"

python compute_rosetta_metrics.py \
  --run-csv "${OUTPUT_DIR}/run.csv" \
  --out-csv "${OUTPUT_DIR}/input_rosetta_metrics.csv" \
  --folder input:"${OUTPUT_DIR}/input_pdbs" \

conda deactivate
conda deactivate
END_TIME=$(date +%s)
echo "Structures relaxed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# ==============================================================================
# 3. Run Alphafold2 initial guess 
# >>> REQUIRES an environment to run af2_initial_guess <<<
# ==============================================================================

START_TIME=$(date +%s)
echo -e "\nRunning AF2 initial guess" >> "${LOG_DIR}/log.txt"


# >>> load af2 initial guess environment <<<

mkdir -p "${OUTPUT_DIR}/AF2" && cd "${OUTPUT_DIR}/AF2"

# Run AF2 prediction
predict.py -pdbdir "${OUTPUT_DIR}/input_pdbs/relaxed_pdbs" \
    -scorefilename "out.sc" \
    -outsilent "af2.silent"

END_TIME=$(date +%s)
echo "AF2 initial guess completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# >>> unload af2 initial guess environment <<<

# ==============================================================================
# 4A. Generate MSA files  
# >>> REQUIRES an environment to run ColabFold <<<
# ==============================================================================

cd $SCRIPT_DIR

# >>> load ColabFold environment <<<

colabfold_batch "${OUTPUT_DIR}/unique_msa" "${OUTPUT_DIR}/unique_msa/msa"  --msa-only

# >>> unload ColabFold environment <<<
#==============================================================================
# 4B. Generate model inputs
# ==============================================================================

source "$CONDA_PATH"
conda activate binder_scoring_env

python generate_model_inputs.py \
  --run-csv "${OUTPUT_DIR}/run.csv" \
  --out-dir  "${OUTPUT_DIR}"

conda deactivate
conda deactivate

# ==============================================================================
# 5. Run ColabFold 
# >>> REQUIRES an environment to run ColabFold <<<
# ==============================================================================
START_TIME=$(date +%s)
echo -e "\nRunning ColabFold" >> "${LOG_DIR}/log.txt"

# >>> load ColabFold environment <<<

colabfold_batch "${OUTPUT_DIR}/ColabFold/input_folder" "${OUTPUT_DIR}/ColabFold/ptm_output" --calc-extra-ptm --num-recycle 3  --num-models 3
find "${OUTPUT_DIR}/ColabFold/ptm_output" -type f -name "*.png" -exec rm -f {} \;

# >>> unload ColabFold environment <<<

END_TIME=$(date +%s)
echo "ColabFold completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# ==============================================================================
# 6. Run Boltz  
# >>> REQUIRES an environment to run Boltz <<<
# ==============================================================================
START_TIME=$(date +%s)
echo -e "\nRunning Boltz" >> "${LOG_DIR}/log.txt"

# >>> load Boltz environment <<<

boltz predict "${OUTPUT_DIR}/Boltz/input_folder" \
    --recycling_steps 10 \
    --diffusion_samples 3 \
    --write_full_pae \
    --out_dir "${OUTPUT_DIR}/Boltz"

# >>> unload Boltz environment <<<

END_TIME=$(date +%s)
echo "Boltz completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# ==============================================================================
# 7. Run AF3 
# >>> REQUIRES an environment to run AF3 <<<
# ==============================================================================
START_TIME=$(date +%s)
echo -e "\nRunning AF3" >> "${LOG_DIR}/log.txt"

# >>> load AF3 environment <<<

python run_alphafold.py \
--input_dir="${OUTPUT_DIR}/AF3/input_folder" \
--model_dir=/path/to/alphafold3_weights \
--db_dir=/path/to/alphafold3_database \
--run_data_pipeline=False \
--num_diffusion_samples=3 \
--output_dir="${OUTPUT_DIR}/AF3/outputs"

#  >>> unload AF3 environment <<<

END_TIME=$(date +%s)
echo "AF3 completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# ==============================================================================
# 8. Extracting confidence metrics and model_0 PDBs
# ==============================================================================
cd $SCRIPT_DIR
START_TIME=$(date +%s)
echo -e "\nExtracting confidence metrics" >> "${LOG_DIR}/log.txt"


source "$CONDA_PATH"
conda activate binder_scoring_env

python extract_confidence_metrics.py \
  --run-csv  "${OUTPUT_DIR}/run.csv" \
  --out-dir "${OUTPUT_DIR}"


# ==============================================================================
# 9. compute ipSAE and other interface confidence metrics
# ==============================================================================
echo -e "\nComputing ipSAE and other interface confidence metrics" >> "${LOG_DIR}/log.txt"
python run_ipsae_batch.py \
  --run-csv "${OUTPUT_DIR}/run.csv" \
  --out-csv "${OUTPUT_DIR}/ipsae_and_ipae.csv" \
  --af3-dir "${OUTPUT_DIR}/AF3" \
  --boltz1-dir "${OUTPUT_DIR}/Boltz" \
  --colab-dir "${OUTPUT_DIR}/ColabFold" \
  --ipsae-script-path ./ipsae_w_ipae.py \
  --pae-cutoff 10 --dist-cutoff 10 \

# ==============================================================================
# 10. Computing DockQ
# ==============================================================================
echo -e "\nComputing dockQ" >> "${LOG_DIR}/log.txt"

python dockQ.py \
  --input-pdbs "${OUTPUT_DIR}/input_pdbs/" \
  --folder af3:"${OUTPUT_DIR}/AF3/pdbs/" \
  --folder af2:"${OUTPUT_DIR}/AF2/pdbs/" \
  --folder boltz1:"${OUTPUT_DIR}/Boltz/pdbs" \
  --folder colab:"${OUTPUT_DIR}/ColabFold/pdbs" \
  --out-csv "${OUTPUT_DIR}/dockQ.csv"


# ==============================================================================
# 11. Compute Rosetta metrics
# ==============================================================================

START_TIME=$(date +%s)
echo -e "\nRelaxing model PDBs and computing Rosetta metrics" >> "${LOG_DIR}/log.txt"

python compute_rosetta_metrics.py \
  --run-csv "${OUTPUT_DIR}/run.csv" \
  --out-csv "${OUTPUT_DIR}/rosetta_metrics.csv" \
  --folder af3:"${OUTPUT_DIR}/AF3/pdbs" \
  --folder boltz1:"${OUTPUT_DIR}/Boltz/pdbs" \
  --folder colab:"${OUTPUT_DIR}/ColabFold/pdbs" \
  --folder af2:"${OUTPUT_DIR}/AF2/pdbs" \

END_TIME=$(date +%s)
echo "Structures relaxed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"
# ==============================================================================
# 12. Compute RMSDs
# ==============================================================================
echo -e "\nComputing RMSDs" >> "${LOG_DIR}/log.txt"

python rmsd.py \
  --folder input:"${OUTPUT_DIR}/input_pdbs/" \
  --folder af3:"${OUTPUT_DIR}/AF3/pdbs/" \
  --folder af2:"${OUTPUT_DIR}/AF2/pdbs/" \
  --folder boltz1:"${OUTPUT_DIR}/Boltz/pdbs" \
  --folder colab:"${OUTPUT_DIR}/ColabFold/pdbs" \
  --out-csv "${OUTPUT_DIR}/rmsd.csv"
conda deactivate
conda deactivate

# ==============================================================================
# 13. Pymol Metrics: Interface Analysis & Hydrogen Bonds
# >>> REQUIRES an environment to run Pymol <<<
# ==============================================================================
cd ${OUTPUT_DIR}
START_TIME=$(date +%s)
PYMOL_DIR="${OUTPUT_DIR}/pymol_files"
mkdir -p "${PYMOL_DIR}"

# >>> load Pymol environment <<<

# Create a JSON file that lists the directories for Pymol analysis.
echo '{"input": "'${OUTPUT_DIR}/input_pdbs'", "af2": "'${OUTPUT_DIR}/AF2/pdbs'", "colab": "'${OUTPUT_DIR}/ColabFold/pdbs'", "boltz1": "'${OUTPUT_DIR}/Boltz/pdbs'", "af3": "'${OUTPUT_DIR}/AF3/pdbs'"}' > "${PYMOL_DIR}/pdb_dirs.json"
# Run Pymol in command-line mode to execute the analysis script.
python -m pymol -c -d "run ${SCRIPT_DIR}/pymol_metrics.py"

# >>> unload Pymol environment <<<

END_TIME=$(date +%s)
echo "Pymol metrics calculated in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# ==============================================================================
# 14. Overall Execution Time Logging
# ==============================================================================
END_TIME=$(date +%s)
RUN_TIME=$(( END_TIME - OVERALL_START_TIME ))
echo -e "\nOverall pipeline execution time: ${RUN_TIME} seconds" >> "${LOG_DIR}/log.txt"
TIME_PER_STRUC=$(echo "scale=2; $RUN_TIME / $count" | bc)
echo "${TIME_PER_STRUC} seconds per design"  >> "${LOG_DIR}/log.txt"