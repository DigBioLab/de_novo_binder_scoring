# Reproducibility Analyses

This directory contains the notebooks and scripts required to reproduce the analyses and results presented in the paper.

## Contents

### Notebooks

- **`metrics_analysis.ipynb`**  
  Analysis of **individual** and **interaction** features.  

  **Input:**  
  - `./data/prepared_training_dataset.csv`

---

- **`logistic_regression_analysis.ipynb`**  
  Logistic regression with **greedy feature selection** and an **XGBoost** baseline.  

  **Training scripts:**  
  - `./scripts/`

  **Outputs:**  
  - `./data/lr_greedy_selection_outputs/`  
  - `./data/xgboost_outputs/`

---

- **`RMSD_agreement_features.ipynb`**  
  Analysis of RMSD agreement between structure predictors and their relationship to binding data.

---

- **`affinity_analysis.ipynb`**  
  Affinity analysis comparing predicted metrics and affinity-related feature behavior across datasets and model outputs.

  **Inputs:**  
  - `./data/all_affinity_data_with_predicted_metrics.csv`

---

- **`per_design_method_analysis.ipynb`**  
  Per-design-method analysis evaluating performance and predicted metrics across different design pipelines and generation methods.

  **Inputs:**  
  - original_data_a_bindcraft_a_boltzgen_prepared_column_names_updated.csv (on zenodo see below)

---

## Data

### Included locally

The following file should be placed in the `./data/` directory:

- `all_affinity_data_with_predicted_metrics.csv`
- `prepared_training_dataset.csv`
- `ranked_features_with_top50_interactions.csv`

---

### External dataset

The extended dataset with bindcraft and boltzgen is available on Zenodo:

- `original_data_a_bindcraft_a_boltzgen_prepared_column_names_updated.csv`  
  Zenodo DOI / Record: https://zenodo.org/records/15722219

Download the file and place it in the appropriate `./data/` directory before running the notebooks.

---
