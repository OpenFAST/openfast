"""Import all ORM models so Alembic and SQLAlchemy metadata discover them."""

from app.models.components import Airfoil, Blade, Controller, Tower, TurbineModel  # noqa: F401
from app.models.project import Project  # noqa: F401
from app.models.simulation import (  # noqa: F401
    DLCDefinition,
    ResultsDEL,
    ResultsExtreme,
    ResultsStatistics,
    Simulation,
    SimulationCase,
)
from app.models.user import Organization, User  # noqa: F401
