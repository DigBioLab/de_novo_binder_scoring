import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, precision_recall_curve
from sklearn.model_selection import LeaveOneGroupOut
import pickle
from pathlib import Path

# --- Config ---
DATA_PATH = '../data/prepared_training_dataset.csv'
RANKED_FEATURES_PATH = '../data/ranked_features_with_top50_interactions.csv'
OUTPUT_DIR = '../data/lr_greedy_selection_outputs'
C_VALUES = [0.01,0.1,1]
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# --- Load data ---
df = pd.read_csv(DATA_PATH)
df['binder'] = df['binder'].astype(int)
targets = df['target_id'].unique().tolist()
ranked_df = pd.read_csv(RANKED_FEATURES_PATH)


GENERAL_FEATURES = ['camsol_score', 'A_length', 'B_length', 'spatialPPI_poi']
MODEL_TAGS = ['af3', 'af2', 'boltz1', 'colab']


all_num_cols = df.select_dtypes(include=np.number).columns.tolist()
INPUT_FEATURES = [
    c for c in all_num_cols
    if ('input' in c.lower()) and not any(tag in c.lower() for tag in MODEL_TAGS)
]

interaction_feats = ranked_df[ranked_df['is_interaction']]['feature'].unique()

def add_interactions(df, interaction_list):
    new_cols = {}
    for feat in interaction_list:
        if feat in df.columns or '*' not in feat:
            continue
        f1, f2 = feat.split('*')
        if f1 in df.columns and f2 in df.columns:
            new_cols[feat] = df[f1] * df[f2]
    if new_cols:
        df = pd.concat([df, pd.DataFrame(new_cols)], axis=1)
    return df

df = add_interactions(df, interaction_feats)

def get_top_features(tag):
    df_tag = ranked_df[ranked_df['tag'] == tag]
    top_single = df_tag[df_tag['is_interaction'] == False].sort_values('median_AP', ascending=False).head(1)
    top_inter = df_tag[df_tag['is_interaction'] == True].sort_values('median_AP', ascending=False).head(1)
    return top_single['feature'].values[0], top_inter['feature'].values[0]

logo = LeaveOneGroupOut()

def compute_precision_metrics(y_true, y_scores):
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    f1_scores = 2 * precision * recall / (precision + recall + 1e-8)
    max_f1 = np.max(f1_scores)

    if np.any(recall >= 0.2):
        prec_at_0_2 = np.max(precision[recall >= 0.2])
    else:
        prec_at_0_2 = np.nan

    if np.any(recall >= 0.6):
        prec_at_0_6 = np.max(precision[recall >= 0.6])
    else:
        prec_at_0_6 = np.nan

    return {
        'prec_at_0.2': prec_at_0_2,
        'prec_at_0.6': prec_at_0_6,
        'max_f1': max_f1
    }

def nested_cv_baseline(holdout_target, feat, C_vals):
    inner = df[df['target_id'] != holdout_target].reset_index(drop=True)
    outer = df[df['target_id'] == holdout_target]
    best_C = None
    best_ap = -np.inf

    for C in C_vals:
        aps = []
        for tr, te in logo.split(inner, inner['binder'], inner['target_id']):
            X_tr, y_tr = inner.iloc[tr][[feat]], inner.iloc[tr]['binder']
            X_te, y_te = inner.iloc[te][[feat]], inner.iloc[te]['binder']
            if y_tr.nunique() < 2 or y_te.nunique() < 2:
                continue
            pipe = Pipeline([
                ('imp', SimpleImputer(strategy='median')),
                ('sc', StandardScaler()),
                ('clf', LogisticRegression(penalty='l1', solver='liblinear', C=C, class_weight='balanced'))
            ])
            pipe.fit(X_tr, y_tr)
            aps.append(average_precision_score(y_te, pipe.predict_proba(X_te)[:, 1]))
        median_ap = np.median(aps) if aps else -np.inf
        if median_ap > best_ap:
            best_ap = median_ap
            best_C = C

    pipe = Pipeline([
        ('imp', SimpleImputer(strategy='median')),
        ('sc', StandardScaler()),
        ('clf', LogisticRegression(penalty='l1', solver='liblinear', C=best_C, class_weight='balanced'))
    ])
    pipe.fit(inner[[feat]], inner['binder'])
    y_scores = pipe.predict_proba(outer[[feat]])[:, 1]
    ap = average_precision_score(outer['binder'], y_scores)
    metrics = compute_precision_metrics(outer['binder'], y_scores)
    return ap, metrics

from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut

def nested_cv_greedy(holdout_target, candidates, C_vals, max_feats=10, tol=0.005):
    """
    1) Greedy feature‐selection with fixed C=C_vals[0]
    2) Global C tuning via GridSearchCV on selected features
    3) Final evaluation on outer holdout
    """
    # split inner vs outer
    inner = df[df['target_id'] != holdout_target].reset_index(drop=True)
    outer = df[df['target_id'] == holdout_target].reset_index(drop=True)

    logo_inner = LeaveOneGroupOut()
    fixed_C = 0.01

    selected    = []
    track       = []   # [(feat, inner_median_ap)]
    outer_track = []   # [(feat, outer_ap_when_added)]

    best_score = -np.inf

    # --- Phase 1: greedy select up to max_feats using fixed_C ---
    while len(selected) < max_feats and candidates:
        feat_to_ap = {}
        for feat in candidates:
            aps = []
            for tr, te in logo_inner.split(inner, inner['binder'], inner['target_id']):
                X_tr = inner.iloc[tr][selected + [feat]] if selected else inner.iloc[tr][[feat]]
                X_te = inner.iloc[te][selected + [feat]] if selected else inner.iloc[te][[feat]]
                y_tr = inner.iloc[tr]['binder']
                y_te = inner.iloc[te]['binder']

                if y_tr.nunique() < 2 or y_te.nunique() < 2:
                    continue

                pipe = Pipeline([
                    ('imp', SimpleImputer(strategy='median')),
                    ('sc',  StandardScaler()),
                    ('clf', LogisticRegression(
                        penalty='l1',
                        solver='liblinear',
                        C=fixed_C,
                        class_weight='balanced'
                    ))
                ])
                pipe.fit(X_tr, y_tr)
                aps.append( average_precision_score(y_te, pipe.predict_proba(X_te)[:,1]) )

            feat_to_ap[feat] = np.median(aps) if aps else -np.inf

        # pick best feature this round
        best_feat, best_ap = max(feat_to_ap.items(), key=lambda x: x[1])
        if best_ap <= best_score + tol:
            break

        selected.append(best_feat)
        track.append((best_feat, best_ap))
        best_score = best_ap
        candidates.remove(best_feat)

        # record outer‐fold performance immediately after adding
        pipe_o = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('sc',  StandardScaler()),
            ('clf', LogisticRegression(
                penalty='l1',
                solver='liblinear',
                C=fixed_C,
                class_weight='balanced'
            ))
        ])
        pipe_o.fit(inner[selected], inner['binder'])
        outer_ap = average_precision_score(
            outer['binder'],
            pipe_o.predict_proba(outer[selected])[:,1]
        )
        outer_track.append((best_feat, outer_ap))

    # --- Phase 2: global C‐tuning on the selected features ---
    if selected:
        pipe_final = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('sc',  StandardScaler()),
            ('clf', LogisticRegression(
                penalty='l1',
                solver='liblinear',
                class_weight='balanced'
            ))
        ])
        param_grid = {'clf__C': C_vals}

        X_inner = inner[selected].values
        y_inner = inner['binder'].values
        groups  = inner['target_id'].values

        grid = GridSearchCV(
            estimator=pipe_final,
            param_grid=param_grid,
            scoring='average_precision',
            cv=logo_inner,
            n_jobs=-1,
            verbose=0
        )
        grid.fit(X_inner, y_inner, groups=groups)
        best_C = grid.best_params_['clf__C']

        # final fit & evaluate on outer
        pipe_final.set_params(clf__C=best_C)
        pipe_final.fit(X_inner, y_inner)
        y_scores = pipe_final.predict_proba(outer[selected])[:,1]

        ap      = average_precision_score(outer['binder'], y_scores)
        metrics = compute_precision_metrics(outer['binder'], y_scores)
        coef_dict = dict(zip(selected, pipe_final.named_steps['clf'].coef_[0]))
    else:
        best_C = None
        ap = np.nan
        metrics = {'prec_at_0.2': np.nan, 'prec_at_0.6': np.nan, 'max_f1': np.nan}
        coef_dict = {}

    return (
        holdout_target,
        best_C,
        ap,
        selected,
        track,
        outer_track,
        coef_dict,
        metrics
    )


def run_combined(holdout_target, baseline_feat, candidates):
    baseline_ap, baseline_metrics = nested_cv_baseline(holdout_target, baseline_feat, C_VALUES)
    hold, best_c, ap, selected, track, outer_track, coefs, metrics = nested_cv_greedy(holdout_target, candidates, C_VALUES)
    return hold, baseline_ap, baseline_metrics, ap,best_c, selected, track,outer_track, coefs, metrics

def run_all_models():
    for mode in ['single', 'interaction', 'all_single', 'all_interaction']:
        print(f"\n=== Running mode: {mode} ===")

        if mode.startswith('all'):
            # Only run once for "all" modes
            tag = 'all'
            global_top_single = ranked_df[~ranked_df['is_interaction']].sort_values('median_AP', ascending=False).iloc[0]['feature']
            global_top_inter = ranked_df[ranked_df['is_interaction']].sort_values('median_AP', ascending=False).iloc[0]['feature']

            if mode == 'all_single':
                baseline_feat = global_top_single
                rmsd_feats = [c for c in df.select_dtypes(include=np.number).columns if 'RMSD' in c and c not in ranked_df['feature'].values]
                candidates = list(set(ranked_df[~ranked_df['is_interaction']]['feature'].tolist() + rmsd_feats)) + GENERAL_FEATURES + INPUT_FEATURES
            else:  # 'all_interaction'
                baseline_feat = global_top_inter
                rmsd_feats = [c for c in df.select_dtypes(include=np.number).columns if 'RMSD' in c and c not in ranked_df['feature'].values]
                candidates = list(set(ranked_df[ranked_df['is_interaction']]['feature'].tolist() + rmsd_feats)) + GENERAL_FEATURES + INPUT_FEATURES

            results = {}
            with ProcessPoolExecutor() as executor:
                futures = {
                    executor.submit(run_combined, hold, baseline_feat, candidates.copy()): hold
                    for hold in targets
                }
                for fut in as_completed(futures):
                    hold, baseline_ap, baseline_metrics, ap,best_c, selected, track, outer_track, coefs, greedy_metrics = fut.result()
                    results[hold] = {
                        'baseline_ap': baseline_ap,
                        'baseline_metrics': baseline_metrics,
                        'greedy_ap': ap,
                        'best_c': best_c,
                        'greedy_selected': selected,
                        'greedy_track': track,
                        'greedy_outer_track': outer_track,
                        'greedy_coefficients': coefs,
                        'greedy_metrics': greedy_metrics
                    }
                    print(f"Target {hold}: baseline {baseline_ap:.3f}, greedy {ap:.3f}")

            if mode == 'all_single': 
                out_path = Path(OUTPUT_DIR) / f"nested_greedy_single_all.pkl"
            else:
                out_path = Path(OUTPUT_DIR) / f"nested_greedy_interaction_all.pkl"
            with open(out_path, "wb") as f:
                pickle.dump(results, f)
            print(f"Saved to {out_path}")

        else:
            for tag in MODEL_TAGS:
                print(f"--- Model tag: {tag} ---")
                top_single, top_inter = get_top_features(tag)

                if mode == 'single':
                    baseline_feat = top_single
                    candidates = ranked_df[(ranked_df['tag'] == tag) & (~ranked_df['is_interaction'])]['feature'].tolist() + GENERAL_FEATURES + INPUT_FEATURES
                    print(len(candidates))
                    print(candidates)
                else:  # 'interaction'
                    baseline_feat = top_inter
                    candidates = ranked_df[(ranked_df['tag'] == tag) & (ranked_df['is_interaction'])]['feature'].tolist() + GENERAL_FEATURES + INPUT_FEATURES

                results = {}
                with ProcessPoolExecutor() as executor:
                    futures = {
                        executor.submit(run_combined, hold, baseline_feat, candidates.copy()): hold
                        for hold in targets
                    }
                    for fut in as_completed(futures):
                        hold, baseline_ap, baseline_metrics, ap,best_c, selected, track, outer_track, coefs, greedy_metrics = fut.result()
                        results[hold] = {
                            'baseline_ap': baseline_ap,
                            'baseline_metrics': baseline_metrics,
                            'greedy_ap': ap,
                            'best_c': best_c,
                            'greedy_selected': selected,
                            'greedy_track': track,
                            'greedy_outer_track': outer_track,
                            'greedy_coefficients': coefs,
                            'greedy_metrics': greedy_metrics
                        }
                        print(f"Target {hold}: baseline {baseline_ap:.3f}, greedy {ap:.3f}")

                out_path = Path(OUTPUT_DIR) / f"nested_greedy_{mode}_{tag}.pkl"
                with open(out_path, "wb") as f:
                    pickle.dump(results, f)
                print(f"Saved to {out_path}")

if __name__ == '__main__':
    run_all_models()