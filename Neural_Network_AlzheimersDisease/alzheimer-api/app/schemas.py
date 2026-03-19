"""
schemas.py — Pydantic request / response models.

Feature descriptions are taken from the AIBL data dictionary:
  PTGENDER    – Patient gender (1=Male, 2=Female)
  APGEN1/2    – APOE genotype alleles (2, 3, or 4)
  APRECEIVE   – APOE test received flag
  APAMBTEMP   – Ambient temperature at blood draw
  CDGLOBAL    – Clinical Dementia Rating global score
  LIMMTOTAL   – Logical Memory Immediate Recall total
  LDELTOTAL   – Logical Memory Delayed Recall total
  AXT117      – Thyroid stimulating hormone (TSH)
  BAT126      – Vitamin B12
  HMT3/7/13/40/100/102 – Haematology panel values
  RCT6/11/20/392       – Clinical chemistry panel values
"""
from __future__ import annotations

from datetime import datetime
from typing import Literal, Optional

from pydantic import BaseModel, Field


# ── Request ───────────────────────────────────────────────────────────────────

class PatientFeatures(BaseModel):
    """
    Input features for an Alzheimer's Disease prediction.
    All fields are optional — missing values are imputed with the
    training-set median for that feature.
    """
    PTGENDER:  Optional[float] = Field(None, description="Patient gender (1=Male, 2=Female)")
    APGEN1:    Optional[float] = Field(None, description="APOE allele 1 (2, 3, or 4)")
    APGEN2:    Optional[float] = Field(None, description="APOE allele 2 (2, 3, or 4)")
    APRECEIVE: Optional[float] = Field(None, description="APOE test received flag")
    APAMBTEMP: Optional[float] = Field(None, description="Ambient temperature at blood draw")
    CDGLOBAL:  Optional[float] = Field(None, description="CDR global score (0=normal, 3=severe)")
    LIMMTOTAL: Optional[float] = Field(None, description="Logical Memory Immediate Recall total")
    LDELTOTAL: Optional[float] = Field(None, description="Logical Memory Delayed Recall total")
    AXT117:    Optional[float] = Field(None, description="TSH (thyroid stimulating hormone)")
    BAT126:    Optional[float] = Field(None, description="Vitamin B12 level")
    HMT3:      Optional[float] = Field(None, description="Haematology panel — HMT3")
    HMT7:      Optional[float] = Field(None, description="Haematology panel — HMT7")
    HMT13:     Optional[float] = Field(None, description="Haematology panel — HMT13")
    HMT40:     Optional[float] = Field(None, description="Haematology panel — HMT40")
    HMT100:    Optional[float] = Field(None, description="Haematology panel — HMT100")
    HMT102:    Optional[float] = Field(None, description="Haematology panel — HMT102")
    RCT6:      Optional[float] = Field(None, description="Clinical chemistry — RCT6")
    RCT11:     Optional[float] = Field(None, description="Clinical chemistry — RCT11")
    RCT20:     Optional[float] = Field(None, description="Clinical chemistry — RCT20")
    RCT392:    Optional[float] = Field(None, description="Clinical chemistry — RCT392")

    model_config = {
        "json_schema_extra": {
            "example": {
                "PTGENDER": 1,
                "APGEN1": 3,
                "APGEN2": 3,
                "CDGLOBAL": 0.0,
                "LIMMTOTAL": 13,
                "LDELTOTAL": 12,
            }
        }
    }


class PredictRequest(BaseModel):
    """Top-level prediction request."""
    features: PatientFeatures
    model: Literal["mlp", "logreg"] = Field(
        "mlp",
        description="Which classifier to use: 'mlp' (neural net) or 'logreg' (logistic regression).",
    )

    model_config = {
        "json_schema_extra": {
            "example": {
                "features": {
                    "PTGENDER": 1,
                    "APGEN1": 3,
                    "APGEN2": 4,
                    "CDGLOBAL": 1.0,
                    "LIMMTOTAL": 8,
                    "LDELTOTAL": 5,
                },
                "model": "mlp",
            }
        }
    }


# ── Response ─────────────────────────────────────────────────────────────────

class PredictionResponse(BaseModel):
    """Prediction result returned by POST /predict."""
    id:           int
    model:        str
    prediction:   int   = Field(..., description="0 = No AD, 1 = AD")
    probability:  float = Field(..., description="P(Alzheimer's Disease = 1)")
    label:        str   = Field(..., description="Human-readable label")
    created_at:   datetime


class PredictionListItem(BaseModel):
    """Summary row used in GET /predictions."""
    id:          int
    model_name:  str
    prediction:  int
    probability: float
    created_at:  datetime

    model_config = {"from_attributes": True}


class ModelMetrics(BaseModel):
    """Evaluation metrics for a single model."""
    accuracy: float
    roc_auc:  float


class HealthResponse(BaseModel):
    """Response from GET /health."""
    status:  str
    version: str
    models:  dict[str, ModelMetrics]
