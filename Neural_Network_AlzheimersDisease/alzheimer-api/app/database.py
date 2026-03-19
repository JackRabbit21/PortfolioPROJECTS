"""
database.py — SQLAlchemy setup and ORM model for prediction logs.
"""
from datetime import datetime

from sqlalchemy import (
    Column, DateTime, Float, Integer, String, Text, create_engine
)
from sqlalchemy.orm import DeclarativeBase, sessionmaker

from app.config import DATABASE_URL

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False},  # needed for SQLite
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


class Base(DeclarativeBase):
    pass


class PredictionLog(Base):
    """Persists every prediction request for auditability."""
    __tablename__ = "prediction_logs"

    id            = Column(Integer, primary_key=True, index=True)
    created_at    = Column(DateTime, default=datetime.utcnow, nullable=False)
    model_name    = Column(String(32), nullable=False)         # "mlp" | "logreg"
    input_json    = Column(Text, nullable=False)               # raw request features
    prediction    = Column(Integer, nullable=False)            # 0 or 1
    probability   = Column(Float,   nullable=False)            # P(AD=1)
    requested_by  = Column(String(128), nullable=True)         # username from auth


def create_tables():
    Base.metadata.create_all(bind=engine)


def get_db():
    """FastAPI dependency — yields a DB session, closes it on exit."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
