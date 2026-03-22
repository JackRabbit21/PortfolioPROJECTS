# ── Stage 1: base image ───────────────────────────────────────────────────────
FROM python:3.11-slim AS base

WORKDIR /app

# Install dependencies first (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ── Stage 2: copy application code ────────────────────────────────────────────
COPY app/       ./app/
COPY artifacts/ ./artifacts/
COPY train.py   ./train.py

# ── Runtime configuration ──────────────────────────────────────────────────────
ENV ARTIFACTS_DIR=/app/artifacts
ENV DATABASE_URL=sqlite:////app/data/alzheimer_api.db
# Override these in docker-compose or at runtime:
ENV API_USERNAME=admin
ENV API_PASSWORD=secret

# Create a persistent volume mount point for the SQLite database
RUN mkdir -p /app/data

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
