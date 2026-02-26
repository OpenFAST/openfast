"""WindForge API — FastAPI application entry point."""

import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.config import settings
from app.database import engine
from app.routers import auth, blades, controllers, projects, towers, websocket
from app.routers.simulations import dlc_router, router as simulations_router

logger = logging.getLogger("windforge")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler — startup and shutdown logic."""

    # ---- startup ----
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    )
    logger.info("WindForge API starting up")

    # Ensure the OpenFAST work directory exists
    work_dir = Path(settings.OPENFAST_WORK_DIR)
    work_dir.mkdir(parents=True, exist_ok=True)
    logger.info("OpenFAST work directory: %s", work_dir)

    yield

    # ---- shutdown ----
    logger.info("WindForge API shutting down")
    await engine.dispose()


app = FastAPI(
    title="WindForge API",
    description="Wind turbine design SaaS powered by OpenFAST",
    version="0.1.0",
    lifespan=lifespan,
)

# ---------------------------------------------------------------------------
# Middleware
# ---------------------------------------------------------------------------
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------------
# Routers
# ---------------------------------------------------------------------------
app.include_router(auth.router, prefix="/api/v1")
app.include_router(projects.router, prefix="/api/v1")
app.include_router(towers.router, prefix="/api/v1")
app.include_router(blades.router, prefix="/api/v1")
app.include_router(controllers.router, prefix="/api/v1")
app.include_router(simulations_router, prefix="/api/v1")
app.include_router(dlc_router, prefix="/api/v1")
app.include_router(websocket.router)


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------
@app.get("/", tags=["health"])
async def root():
    """Root health-check endpoint."""
    return {
        "service": "windforge-api",
        "status": "healthy",
        "version": "0.1.0",
    }


@app.get("/health", tags=["health"])
async def health():
    """Detailed health check."""
    return {
        "service": "windforge-api",
        "status": "healthy",
        "version": "0.1.0",
        "database": settings.DATABASE_URL.split("@")[-1] if "@" in settings.DATABASE_URL else "configured",
        "openfast_lib": settings.OPENFAST_LIB_PATH,
    }
