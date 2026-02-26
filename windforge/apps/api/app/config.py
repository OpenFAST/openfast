"""Application configuration via environment variables with Pydantic Settings."""

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """WindForge API configuration.

    All values can be overridden via environment variables or a .env file.
    """

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=True,
    )

    # --- Database ---
    DATABASE_URL: str = "postgresql+asyncpg://windforge:windforge@localhost:5432/windforge"

    # --- Authentication / JWT ---
    SECRET_KEY: str = "CHANGE-ME-in-production-use-openssl-rand-hex-32"
    JWT_ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24  # 24 hours

    # --- OpenFAST ---
    OPENFAST_LIB_PATH: str = "/usr/local/lib/libopenfastlib.so"
    OPENFAST_WORK_DIR: str = "/tmp/windforge/work"

    # --- CORS ---
    CORS_ORIGINS: list[str] = [
        "http://localhost:3000",
        "http://localhost:5173",
    ]


settings = Settings()
