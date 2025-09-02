#!/bin/sh
### General options
### –- specify queue --
#BSUB -q cabgpu
### -- set the job Name --
#BSUB -J example_job
### -- ask for number of cores (default: 1) --
#BSUB -n 16
### -- specify that the cores must be on the same host --
#BSUB -R "span[hosts=1]"
### -- Select the resources: 1 gpu in exclusive process mode --
#BSUB -gpu "num=1:mode=exclusive_process"
### -- set walltime limit: hh:mm --  maximum 24 hours for GPU-queues right now
#BSUB -W 24:00
# request 5GB of system-memory
#BSUB -R "rusage[mem=10GB]"
### -- set the email address --
##BSUB -u maxove@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o %J.out
#BSUB -e %J.err
# -- end of LSF options --

# ==============================================================================
# Variables & Directory Setup
# ==============================================================================

# Python virtual environment activation script
CONDA_ENV="/dtu/projects/RFdiffusion/scripts_collection/binder_desing_metrics/miniconda3/bin/activate"
AF3_WEIGHTS="/dtu/projects/RFdiffusion/AF3/3.0.0/model_parameter"
OUTPUT_DIR="./example_output_overwrite"
INPUT_PDBS="./example_input_overwrite/input_pdbs"
INPUT_CSV="./example_input/input_overwrite_example.csv"
SCRIPT_DIR="/work3/anryg/de_novo_binder_scoring"


LOG_DIR="${OUTPUT_DIR}/log"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

OVERALL_START_TIME=$(date +%s)

count=$(tail -n +1 "${INPUT_CSV}" | wc -l)
echo "Running ${count} structures..." > "${LOG_DIR}/log.txt"

# ==============================================================================
# 1. Pre-process input CSV and PDB Files
# ==============================================================================
# Activate the Python environment and run the meta_analysis script to convert PDB files to a CSV.
cd $SCRIPT_DIR

source "$CONDA_ENV"
conda activate binder_scoring_env

python process_inputs.py \
  --input_pdbs "${INPUT_PDBS}" \
  --output_dir "${OUTPUT_DIR}" \
  --mode "hybrid" \
  --input_csv "${INPUT_CSV}"

# # ==============================================================================
# # 2. Relax input structures and compute Rosetta metrics
# # ==============================================================================
# START_TIME=$(date +%s)
# echo -e "\nRelaxing input PDBs and computing Rosetta metrics" >> "${LOG_DIR}/log.txt"

# python compute_rosetta_metrics.py \
#   --run_csv "${OUTPUT_DIR}/run.csv" \
#   --out_csv "${OUTPUT_DIR}/input_rosetta_metrics.csv" \
#   --input input="${OUTPUT_DIR}/input_pdbs" \

# conda deactivate
# conda deactivate
# END_TIME=$(date +%s)
# echo "Structures relaxed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# # ==============================================================================
# # 3. Run Alphafold2 initial guess
# # ==============================================================================
# START_TIME=$(date +%s)
# echo -e "\nRunning AF2 initial guess" >> "${LOG_DIR}/log.txt"
# source /dtu/projects/RFdiffusion/setup.sh
# module load af2_binder_design
# mkdir -p "${OUTPUT_DIR}/AF2" && cd "${OUTPUT_DIR}/AF2"

# # Run AF2 prediction: generate a score file and a silent file from the input silent file.
# predict.py -pdbdir "${OUTPUT_DIR}/input_pdbs/relaxed_pdbs" \
#     -scorefilename "out.sc" \
#     -outsilent "af2.silent"

# module purge

# END_TIME=$(date +%s)
# echo "AF2 initial guess completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# # ==============================================================================
# # 4. Generate MSA files and model inputs
# # ==============================================================================
# cd $SCRIPT_DIR
# source /dtu/projects/RFdiffusion/setup.sh
# module load colabfold/1.5.6
# colabfold_batch "${OUTPUT_DIR}/unique_msa" "${OUTPUT_DIR}/unique_msa/msa"  --msa-only
# module purge


# source "$CONDA_ENV"
# conda activate binder_scoring_env

# python generate_model_inputs.py \
#   --csv_file "${OUTPUT_DIR}/run.csv" \
#   --out_dir  "${OUTPUT_DIR}"

# conda deactivate
# conda deactivate

# # ==============================================================================
# # 5. Run ColabFold
# # ==============================================================================
# START_TIME=$(date +%s)
# echo -e "\nRunning ColabFold" >> "${LOG_DIR}/log.txt"
# source /dtu/projects/RFdiffusion/setup.sh
# module load colabfold/1.5.6
# colabfold_batch "${OUTPUT_DIR}/ColabFold/input_folder" "${OUTPUT_DIR}/ColabFold/ptm_output" --calc-extra-ptm --num-recycle 3  --num-models 3
# find "${OUTPUT_DIR}/ColabFold/ptm_output" -type f -name "*.png" -exec rm -f {} \;
# module purge

# END_TIME=$(date +%s)
# echo "ColabFold completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# # ==============================================================================
# # 6. Run Boltz
# # ==============================================================================
# START_TIME=$(date +%s)
# echo -e "\nRunning Boltz" >> "${LOG_DIR}/log.txt"
# source /dtu/projects/RFdiffusion/boltz-1/1.0.0/boltx_env/bin/activate
# boltz predict "${OUTPUT_DIR}/Boltz/input_folder" \
#     --cache /dtu/projects/RFdiffusion/boltz-1/1.0.0 \
#     --recycling_steps 10 \
#     --diffusion_samples 3 \
#     --write_full_pae \
#     --out_dir "${OUTPUT_DIR}/Boltz"
# deactivate
# END_TIME=$(date +%s)
# echo "Boltz completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# # # ==============================================================================
# # # 7. Run AF3
# # # ==============================================================================
# # START_TIME=$(date +%s)
# # echo -e "\nRunning AF3" >> "${LOG_DIR}/log.txt"

# # source /dtu/projects/RFdiffusion/AF3/3.0.1/miniforge/bin/activate
# # source activate alphafold3.0.1_env
# # cd /dtu/projects/RFdiffusion/AF3/3.0.1/alphafold3

# # export PATH=/dtu/projects/RFdiffusion/AF3/3.0.0/hmmer/bin:$PATH
# # #export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
# # python run_alphafold.py \
# # --input_dir="${OUTPUT_DIR}/AF3/input_folder" \
# # --model_dir="${AF3_WEIGHTS}" \
# # --db_dir=/dtu/datasets1/alphafold3 \
# # --run_data_pipeline=False \
# # --num_diffusion_samples=3 \
# # --output_dir="${OUTPUT_DIR}/AF3/outputs"

# # #--flash_attention_implementation=xla \
# # conda deactivate
# # conda deactivate

# # END_TIME=$(date +%s)
# # echo "AF3 completed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"
# # cd $SCRIPT_DIR

# # ==============================================================================
# # 7. Extracting confidence metrics and model_0 PDBs
# # ==============================================================================
# cd $SCRIPT_DIR
# START_TIME=$(date +%s)
# echo -e "\nExtracting confidence metrics" >> "${LOG_DIR}/log.txt"


# source "$CONDA_ENV"
# conda activate binder_scoring_env

# python extract_confidence_metrics.py \
#   --run_csv  "${OUTPUT_DIR}/run.csv" \
#   --output_dir "${OUTPUT_DIR}"


# # ==============================================================================
# # 8. compute ipSAE and other interface confidence metrics
# # ==============================================================================
# echo -e "\nComputing ipSAE and other interface confidence metrics" >> "${LOG_DIR}/log.txt"
# python run_ipsae_batch.py \
#   --run-csv "${OUTPUT_DIR}/run.csv" \
#   --af3-dir "${OUTPUT_DIR}/AF3" \
#   --boltz1-dir "${OUTPUT_DIR}/Boltz" \
#   --colab-dir "${OUTPUT_DIR}/ColabFold" \
#   --ipsae-script-path ./ipsae_w_ipae.py \
#   --pae-cutoff 10 --dist-cutoff 10 \

# # ==============================================================================
# # 9. Computing DockQ
# # ==============================================================================
# echo -e "\nComputing dockQ" >> "${LOG_DIR}/log.txt"

# python dockQ.py \
#   --input_pdbs "${OUTPUT_DIR}/input_pdbs/" \
#   --model af3:"${OUTPUT_DIR}/AF3/pdbs/" \
#   --model af2:"${OUTPUT_DIR}/AF2/pdbs/" \
#   --model boltz1:"${OUTPUT_DIR}/Boltz/pdbs" \
#   --model colab:"${OUTPUT_DIR}/ColabFold/pdbs" \
#   --output_csv "${OUTPUT_DIR}/dockQ.csv"

# # ==============================================================================
# # 9. Computing Camsol
# # ==============================================================================
# echo -e "\nComputing Camsol" >> "${LOG_DIR}/log.txt"

# # Run CamSol to calculate intrinsic solubility scores.
# python camsol_intrinsic_linear.py "${OUTPUT_DIR}/Binder_seq.fasta" -out "${OUTPUT_DIR}/camsol_scores.txt"

# # ==============================================================================
# # 10. Compute Rosetta metrics
# # ==============================================================================
# source "$CONDA_ENV"
# conda activate binder_scoring_env
# START_TIME=$(date +%s)
# echo -e "\nRelaxing model PDBs and computing Rosetta metrics" >> "${LOG_DIR}/log.txt"

# python compute_rosetta_metrics.py \
#   --run_csv "${OUTPUT_DIR}/run.csv" \
#   --out_csv "${OUTPUT_DIR}/rosetta_metrics.csv" \
#   --input af3="${OUTPUT_DIR}/AF3/pdbs" \
#   --input boltz1="${OUTPUT_DIR}/Boltz/pdbs" \
#   --input colab="${OUTPUT_DIR}/ColabFold/pdbs" \
#   --input af2="${OUTPUT_DIR}/AF2/pdbs" \

# END_TIME=$(date +%s)
# echo "Structures relaxed in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"
# # ==============================================================================
# # 10. Compute RMSDs
# # ==============================================================================
# echo -e "\nComputing RMSDs" >> "${LOG_DIR}/log.txt"

# python rmsd.py \
#   --folder input:"${OUTPUT_DIR}/input_pdbs/" \
#   --folder af3:"${OUTPUT_DIR}/AF3/pdbs/" \
#   --folder af2:"${OUTPUT_DIR}/AF2/pdbs/" \
#   --folder boltz1:"${OUTPUT_DIR}/Boltz/pdbs" \
#   --folder colab:"${OUTPUT_DIR}/ColabFold/pdbs" \
#   --out-csv "${OUTPUT_DIR}/rmsd.csv"
# conda deactivate

# # ==============================================================================
# # 12. Pymol Metrics: Interface Analysis & Hydrogen Bonds
# # ==============================================================================
# cd ${OUTPUT_DIR}
# START_TIME=$(date +%s)
# PYMOL_DIR="${OUTPUT_DIR}/pymol_files"
# mkdir -p "${PYMOL_DIR}"
# #module load pymol
# source /dtu/projects/RFdiffusion/scripts_collection/binder_desing_metrics/pymol_venv/bin/activate
# /bin/activate
# # Create a JSON file that lists the directories for Pymol analysis.
# echo '{"input": "'${OUTPUT_DIR}/input_pdbs'", "af2": "'${OUTPUT_DIR}/AF2/pdbs'", "colab": "'${OUTPUT_DIR}/ColabFold/pdbs'", "boltz1": "'${OUTPUT_DIR}/Boltz/pdbs'", "af3": "'${OUTPUT_DIR}/AF3/pdbs'"}' > "${PYMOL_DIR}/pdb_dirs.json"
# # Run Pymol in command-line mode to execute the analysis script.
# python -m pymol -c -d "run ${SCRIPT_DIR}/pymol_metrics.py"
# module purge
# END_TIME=$(date +%s)
# echo "Pymol metrics calculated in $((END_TIME - START_TIME)) seconds" >> "${LOG_DIR}/log.txt"

# # ==============================================================================
# # 14. Overall Execution Time Logging
# # ==============================================================================
# END_TIME=$(date +%s)
# RUN_TIME=$(( END_TIME - OVERALL_START_TIME ))
# echo -e "\nOverall pipeline execution time: ${RUN_TIME} seconds" >> "${LOG_DIR}/log.txt"
# TIME_PER_STRUC=$(echo "scale=2; $RUN_TIME / $count" | bc)
# echo "${TIME_PER_STRUC} seconds per design"  >> "${LOG_DIR}/log.txt"