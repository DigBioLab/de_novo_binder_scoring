import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from sklearn.metrics import average_precision_score, precision_recall_curve
import pickle
from pathlib import Path
from sklearn.model_selection import LeaveOneGroupOut, GridSearchCV

# --- Config ---
DATA_PATH = '../data/prepared_training_dataset.csv'
OUTPUT_DIR = '../data/xgboost_outputs'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# --- Load data ---
df = pd.read_csv(DATA_PATH)
df['binder'] = df['binder'].astype(int)
targets = df['target_id'].unique().tolist()

GENERAL_FEATURES = ['camsol_score', 'A_length', 'B_length', 'spatialPPI_poi']
MODEL_TAGS = ['all','af3', 'af2', 'boltz1', 'colab']


all_num_cols = df.select_dtypes(include=np.number).columns.tolist()
INPUT_FEATURES = [
    c for c in all_num_cols
    if ('input' in c.lower()) and not any(tag in c.lower() for tag in MODEL_TAGS)
]


logo = LeaveOneGroupOut()

# --- Metrics ---
def compute_precision_metrics(y_true, y_scores):
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    f1_scores = 2 * precision * recall / (precision + recall + 1e-8)
    return {
        'prec_at_0.2': np.max(precision[recall >= 0.2]) if np.any(recall >= 0.2) else np.nan,
        'prec_at_0.6': np.max(precision[recall >= 0.6]) if np.any(recall >= 0.6) else np.nan,
        'max_f1': np.max(f1_scores)
    }

# --- Nested CV with inner hyperparameter tuning and outer holdout ---
def nested_cv_xgboost_reduced(holdout_target, feature_list):
    print(f"--- Outer holdout: target {holdout_target} ---")
    train = df[df['target_id'] != holdout_target].reset_index(drop=True)
    test  = df[df['target_id'] == holdout_target].reset_index(drop=True)

    # compute class imbalance weight
    n_neg = (train['binder'] == 0).sum()
    n_pos = (train['binder'] == 1).sum()
    scale_pos_weight = n_neg / n_pos if n_pos > 0 else 1.0

    X_train = train[feature_list].values
    y_train = train['binder'].values
    groups  = train['target_id'].values

    pipe = Pipeline([
        ('imputer', SimpleImputer(strategy='median')),
        ('scaler',  StandardScaler()),
        ('clf',     XGBClassifier(
            objective='binary:logistic',
            eval_metric='aucpr',
            n_jobs=-1,
            use_label_encoder=False,
            tree_method='hist',
            n_estimators=150,
            scale_pos_weight=scale_pos_weight,
            verbosity=0
        ))
    ])

    # reduced grid: the "big four" + min_child_weight
    param_grid = {
        'clf__learning_rate':   [0.05, 0.1],
        'clf__max_depth':        [1, 3],
        'clf__subsample':        [0.8, 1.0],
        'clf__colsample_bytree': [0.5, 0.7],
        'clf__min_child_weight':[1, 5]
    }

    inner_cv = logo.split(X_train, y_train, groups)
    grid = GridSearchCV(
        estimator=pipe,
        param_grid=param_grid,
        scoring='average_precision',
        cv=inner_cv,
        verbose=1
    )

    grid.fit(X_train, y_train, groups=groups)
    print(" Best inner params:", grid.best_params_)
    print(f" Best inner AP: {grid.best_score_:.3f}")

    best_pipe = grid.best_estimator_

    X_test, y_test = test[feature_list].values, test['binder'].values
    y_scores = best_pipe.predict_proba(X_test)[:, 1]
    ap = average_precision_score(y_test, y_scores)
    metrics = compute_precision_metrics(y_test, y_scores)
    print(f" Holdout {holdout_target} → AP={ap:.3f}, metrics={metrics}\n")

    return holdout_target, ap, metrics

def run_all_holdouts_reduced():
    for tag in MODEL_TAGS:
        print(f"\n=== Tag: {tag} ===")
        # Build feature list
        if tag == 'all':
            # all model‐tag features, plus general & input
            
            features = [
            c for c in all_num_cols
            if c not in {"binder", "target_id"}
            ]
        else:
            # only this tag’s features
            other_tags = [t for t in MODEL_TAGS if t != tag]

            tag_feats = [
                c for c in all_num_cols
                if (
                    tag in c.lower()  # must contain this tag
                    and all(other not in c.lower() for other in other_tags)  # must contain none of the others
                    and c not in ['binder', 'target_id']
                )
            ]
            features = list(set(GENERAL_FEATURES + INPUT_FEATURES + tag_feats))

        results = {}
        for t in targets:
            hold, ap, metrics = nested_cv_xgboost_reduced(t, features)
            results[hold] = {'ap': ap, 'metrics': metrics}

        # summarize and save
        median_ap = np.median([v['ap'] for v in results.values()])
        print(f"Median AP for {tag}: {median_ap:.3f}")

        out_file = Path(OUTPUT_DIR) / f'{tag}_xgb_results.pkl'
        with open(out_file, 'wb') as f:
            pickle.dump(results, f)
        print(f"Saved results to {out_file}")

if __name__ == '__main__':
    run_all_holdouts_reduced()

