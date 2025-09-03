# Reproducibility Analyses

This directory contains the notebooks and scripts to reproduce the paper’s results.

## Contents
- **`metrics_analysis.ipynb`** – analysis of **individual** and **interaction** features.  
  **Input:** `./data/prepared_training_dataset.csv`
- **`logistic_regression_analysis.ipynb`** – logistic regression with **greedy feature selection** and an **XGBoost** baseline.  
  **Training scripts:** `./scripts/`  
  **Outputs:** `./data/lr_greedy_selection_outputs/`, `./data/xgboost_outputs/`

**Usage:** open the notebooks (in the order above) and run all cells; ensure inputs exist at the listed paths.
