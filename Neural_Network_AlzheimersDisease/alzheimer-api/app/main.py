"""
main.py — FastAPI application entry point.

Endpoints
─────────
POST   /predict                  Submit patient features → get AD prediction
GET    /predictions              List all stored predictions (paginated)
GET    /predictions/{id}         Retrieve a single prediction by ID
DELETE /predictions/{id}         Delete a stored prediction record
GET    /health                   Service health + model metrics (no auth)
GET    /docs                     Auto-generated Swagger UI (no auth)
"""
from __future__ import annotations

import json
from datetime import datetime
from typing import List, Optional

from fastapi import Depends, FastAPI, HTTPException, Query, status
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session

from app.auth import verify_credentials
from app.config import APP_DESCRIPTION, APP_TITLE, APP_VERSION
from app.database import PredictionLog, create_tables, get_db
from app.predictor import registry
from app.schemas import (
    HealthResponse,
    ModelMetrics,
    PredictRequest,
    PredictionListItem,
    PredictionResponse,
)

# ── App setup ─────────────────────────────────────────────────────────────────

app = FastAPI(
    title=APP_TITLE,
    description=APP_DESCRIPTION,
    version=APP_VERSION,
    docs_url="/docs",
    redoc_url="/redoc",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def startup_event():
    create_tables()
    # Warm up the model registry so the first request isn't slow
    registry.get_metrics()


# ── Health ─────────────────────────────────────────────────────────────────────

@app.get(
    "/health",
    response_model=HealthResponse,
    summary="Service health check",
    tags=["System"],
)
def health():
    """
    Returns service status and model evaluation metrics.
    No authentication required.
    """
    raw_metrics = registry.get_metrics()
    return HealthResponse(
        status="ok",
        version=APP_VERSION,
        models={
            name: ModelMetrics(**vals)
            for name, vals in raw_metrics.items()
        },
    )


# ── Predictions ───────────────────────────────────────────────────────────────

@app.post(
    "/predict",
    response_model=PredictionResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Submit patient features and get an AD prediction",
    tags=["Predictions"],
)
def predict(
    body: PredictRequest,
    db: Session = Depends(get_db),
    username: str = Depends(verify_credentials),
):
    """
    Accepts a `features` dict (all fields optional — missing values are
    median-imputed) and an optional `model` selector (`"mlp"` or `"logreg"`).

    Returns:
    - **prediction** — `0` (No Alzheimer's) or `1` (Alzheimer's Disease)
    - **probability** — model confidence P(AD = 1)
    - **label** — human-readable diagnosis string
    - **id** — database record ID for future retrieval
    """
    features_dict = body.features.model_dump()
    pred, proba = registry.predict(features_dict, model_name=body.model)

    log = PredictionLog(
        model_name   = body.model,
        input_json   = json.dumps(features_dict),
        prediction   = pred,
        probability  = round(proba, 6),
        requested_by = username,
    )
    db.add(log)
    db.commit()
    db.refresh(log)

    return PredictionResponse(
        id          = log.id,
        model       = log.model_name,
        prediction  = pred,
        probability = round(proba, 4),
        label       = "Alzheimer's Disease" if pred == 1 else "No Alzheimer's Disease",
        created_at  = log.created_at,
    )


@app.get(
    "/predictions",
    response_model=List[PredictionListItem],
    summary="List stored predictions",
    tags=["Predictions"],
)
def list_predictions(
    skip:       int = Query(0,   ge=0,   description="Number of records to skip"),
    limit:      int = Query(20,  ge=1,  le=200, description="Max records to return"),
    model_name: Optional[str] = Query(None, description="Filter by model: 'mlp' or 'logreg'"),
    db: Session = Depends(get_db),
    _: str      = Depends(verify_credentials),
):
    """
    Returns a paginated list of all stored predictions, newest first.
    Optionally filter by model name.
    """
    q = db.query(PredictionLog).order_by(PredictionLog.created_at.desc())
    if model_name:
        q = q.filter(PredictionLog.model_name == model_name)
    records = q.offset(skip).limit(limit).all()
    return [
        PredictionListItem(
            id          = r.id,
            model_name  = r.model_name,
            prediction  = r.prediction,
            probability = round(r.probability, 4),
            created_at  = r.created_at,
        )
        for r in records
    ]


@app.get(
    "/predictions/{prediction_id}",
    response_model=PredictionResponse,
    summary="Retrieve a single prediction by ID",
    tags=["Predictions"],
)
def get_prediction(
    prediction_id: int,
    db: Session = Depends(get_db),
    _:  str     = Depends(verify_credentials),
):
    """Fetches a single prediction record including the original input features."""
    record = db.query(PredictionLog).filter(PredictionLog.id == prediction_id).first()
    if not record:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Prediction {prediction_id} not found.",
        )
    return PredictionResponse(
        id          = record.id,
        model       = record.model_name,
        prediction  = record.prediction,
        probability = round(record.probability, 4),
        label       = "Alzheimer's Disease" if record.prediction == 1 else "No Alzheimer's Disease",
        created_at  = record.created_at,
    )


@app.delete(
    "/predictions/{prediction_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    summary="Delete a prediction record",
    tags=["Predictions"],
)
def delete_prediction(
    prediction_id: int,
    db: Session = Depends(get_db),
    _:  str     = Depends(verify_credentials),
):
    """Removes a prediction log entry from the database."""
    record = db.query(PredictionLog).filter(PredictionLog.id == prediction_id).first()
    if not record:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Prediction {prediction_id} not found.",
        )
    db.delete(record)
    db.commit()
