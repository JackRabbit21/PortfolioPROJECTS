"""
config.py — Central configuration loaded from environment variables.
All secrets / paths can be overridden via a .env file or Docker env.
"""
import os
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR      = Path(__file__).resolve().parent.parent
ARTIFACTS_DIR = Path(os.getenv("ARTIFACTS_DIR", BASE_DIR / "artifacts"))
DATABASE_URL  = os.getenv("DATABASE_URL", f"sqlite:///{BASE_DIR}/alzheimer_api.db")

# ── Auth ──────────────────────────────────────────────────────────────────────
# In production, override these with strong secrets via environment variables.
API_USERNAME = os.getenv("API_USERNAME", "admin")
API_PASSWORD = os.getenv("API_PASSWORD", "secret")

# ── App metadata ─────────────────────────────────────────────────────────────
APP_TITLE       = "Alzheimer's Disease Prediction API"
APP_DESCRIPTION = (
    "REST API serving an MLP and Logistic Regression classifier trained on the "
    "AIBL (Australian Imaging, Biomarker & Lifestyle) dataset. "
    "Endpoints require HTTP Basic Authentication."
)
APP_VERSION = "1.0.0"
