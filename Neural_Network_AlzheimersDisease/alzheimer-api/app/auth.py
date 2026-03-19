"""
auth.py — HTTP Basic Authentication dependency for FastAPI.
"""
import secrets

from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBasic, HTTPBasicCredentials

from app.config import API_PASSWORD, API_USERNAME

security = HTTPBasic()


def verify_credentials(credentials: HTTPBasicCredentials = Depends(security)) -> str:
    """
    Validates username and password using constant-time comparison to prevent
    timing attacks. Returns the verified username on success.
    """
    correct_username = secrets.compare_digest(
        credentials.username.encode("utf-8"),
        API_USERNAME.encode("utf-8"),
    )
    correct_password = secrets.compare_digest(
        credentials.password.encode("utf-8"),
        API_PASSWORD.encode("utf-8"),
    )
    if not (correct_username and correct_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid username or password.",
            headers={"WWW-Authenticate": "Basic"},
        )
    return credentials.username
