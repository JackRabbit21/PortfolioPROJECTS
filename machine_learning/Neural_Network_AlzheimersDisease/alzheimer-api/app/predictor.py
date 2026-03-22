"""
predictor.py — Loads trained artifacts once at startup and exposes a
               single `predict()` function used by the API endpoints.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import joblib
import numpy as np

from app.config import ARTIFACTS_DIR


class _ModelRegistry:
    """Singleton that lazily loads and caches all model artifacts."""

    def __init__(self, artifacts_dir: Path):
        self._dir = Path(artifacts_dir)
        self._loaded = False

    def _load(self):
        if self._loaded:
            return
        self.scaler  = joblib.load(self._dir / "scaler.joblib")
        self.logreg  = joblib.load(self._dir / "logreg.joblib")
        self.mlp     = joblib.load(self._dir / "mlp.joblib")

        with open(self._dir / "feature_cols.json") as f:
            self.feature_cols: list[str] = json.load(f)

        with open(self._dir / "feature_medians.json") as f:
            self.medians: dict[str, float] = json.load(f)

        with open(self._dir / "metrics.json") as f:
            self.metrics: dict = json.load(f)

        self._loaded = True

    # ── Public API ────────────────────────────────────────────────────────────

    def predict(
        self,
        features: dict,
        model_name: Literal["mlp", "logreg"] = "mlp",
    ) -> tuple[int, float]:
        """
        Given a dict of raw feature values (may be partial / None),
        returns (prediction_label, probability_of_AD).
        """
        self._load()

        # Build feature vector in training column order, imputing with medians
        row = []
        for col in self.feature_cols:
            val = features.get(col)
            if val is None or (isinstance(val, float) and np.isnan(val)):
                val = self.medians.get(col, 0.0)
            row.append(float(val))

        X = np.array(row, dtype=float).reshape(1, -1)
        X_sc = self.scaler.transform(X)

        clf = self.mlp if model_name == "mlp" else self.logreg
        pred  = int(clf.predict(X_sc)[0])
        proba = float(clf.predict_proba(X_sc)[0, 1])

        return pred, proba

    def get_metrics(self) -> dict:
        self._load()
        return self.metrics


# Module-level singleton — shared across all requests
registry = _ModelRegistry(ARTIFACTS_DIR)
