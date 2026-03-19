"""
train.py — Train Alzheimer's classifier from AIBL data and save artifacts.

Usage:
    python train.py --data path/to/aibl_merged_clean.csv
"""
import argparse
import os
import json

import numpy as np
import pandas as pd
import joblib

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

DROP_COLS = [
    "RID", "SITEID", "SITEID_r", "VISCODE", "VISCODE_r",
    "PTDOB", "APTESTDT", "EXAMDATE", "EXAMDATE_r",
]

OUT_DIR = os.path.join(os.path.dirname(__file__), "artifacts")


def load_and_prepare(data_csv: str):
    df = pd.read_csv(data_csv)

    # Fix labels: -4 → 0 (healthy)
    df["DXAD"] = df["DXAD"].replace({-4: 0}).astype(int)
    df = df.dropna(subset=["DXAD"]).copy()

    # Leakage-safe feature selection: drop IDs, dates, and all DX* columns
    diag_cols = [c for c in df.columns if c.startswith("DX")]
    exclude = set(DROP_COLS + diag_cols + ["DXAD"])
    feature_cols = [c for c in df.columns if c not in exclude]

    X = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(numeric_only=True))
    y = df["DXAD"].astype(int)

    # Save feature column order — the API needs this to align inputs
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "feature_cols.json"), "w") as f:
        json.dump(feature_cols, f, indent=2)

    # Save per-column medians for imputation at inference time
    medians = X.median().to_dict()
    with open(os.path.join(OUT_DIR, "feature_medians.json"), "w") as f:
        json.dump(medians, f, indent=2)

    return X, y, feature_cols


def train(data_csv: str):
    print(f"Loading data from: {data_csv}")
    X, y, feature_cols = load_and_prepare(data_csv)
    print(f"  {X.shape[0]} samples, {X.shape[1]} features: {feature_cols}")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=RANDOM_SEED, stratify=y
    )

    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc  = scaler.transform(X_test)

    # --- Logistic Regression ---
    logreg = LogisticRegression(max_iter=1000, random_state=RANDOM_SEED)
    logreg.fit(X_train_sc, y_train)
    lr_acc = accuracy_score(y_test, logreg.predict(X_test_sc))
    lr_auc = roc_auc_score(y_test, logreg.predict_proba(X_test_sc)[:, 1])
    print(f"\nLogistic Regression  →  Accuracy: {lr_acc:.3f}  ROC-AUC: {lr_auc:.3f}")
    print(classification_report(y_test, logreg.predict(X_test_sc), digits=3))

    # --- MLP Neural Net ---
    mlp = MLPClassifier(
        hidden_layer_sizes=(64, 32),
        activation="relu",
        solver="adam",
        alpha=1e-4,
        batch_size=32,
        learning_rate_init=1e-3,
        max_iter=200,
        random_state=RANDOM_SEED,
        early_stopping=True,
        n_iter_no_change=10,
        validation_fraction=0.1,
    )
    mlp.fit(X_train_sc, y_train)
    mlp_acc = accuracy_score(y_test, mlp.predict(X_test_sc))
    mlp_auc = roc_auc_score(y_test, mlp.predict_proba(X_test_sc)[:, 1])
    print(f"\nMLP Neural Net       →  Accuracy: {mlp_acc:.3f}  ROC-AUC: {mlp_auc:.3f}")
    print(classification_report(y_test, mlp.predict(X_test_sc), digits=3))

    # --- Save artifacts ---
    os.makedirs(OUT_DIR, exist_ok=True)
    joblib.dump(scaler, os.path.join(OUT_DIR, "scaler.joblib"))
    joblib.dump(logreg, os.path.join(OUT_DIR, "logreg.joblib"))
    joblib.dump(mlp,    os.path.join(OUT_DIR, "mlp.joblib"))

    # Save evaluation metrics
    metrics = {
        "logreg": {"accuracy": round(lr_acc, 4), "roc_auc": round(lr_auc, 4)},
        "mlp":    {"accuracy": round(mlp_acc, 4), "roc_auc": round(mlp_auc, 4)},
    }
    with open(os.path.join(OUT_DIR, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nArtifacts saved to: {OUT_DIR}/")
    print("  scaler.joblib, logreg.joblib, mlp.joblib")
    print("  feature_cols.json, feature_medians.json, metrics.json")
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data",
        default="./Data/aibl_merged_clean.csv",
        help="Path to aibl_merged_clean.csv",
    )
    args = parser.parse_args()
    train(args.data)
