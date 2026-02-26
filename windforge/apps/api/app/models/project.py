"""Project ORM model — the top-level container for a turbine design."""

import enum
import uuid
from datetime import datetime

from sqlalchemy import DateTime, Enum, Float, ForeignKey, Integer, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, mapped_column, relationship

from app.database import Base


class ProjectStatus(str, enum.Enum):
    DRAFT = "draft"
    ACTIVE = "active"
    ARCHIVED = "archived"


class Project(Base):
    __tablename__ = "projects"

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True), primary_key=True, default=uuid.uuid4
    )
    org_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True), ForeignKey("organizations.id", ondelete="CASCADE"), nullable=False
    )
    created_by: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"), nullable=True
    )
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str | None] = mapped_column(Text, nullable=True)

    # IEC wind class parameters
    wind_class: Mapped[str | None] = mapped_column(String(10), nullable=True)  # I, II, III, S
    turbulence_class: Mapped[str | None] = mapped_column(String(5), nullable=True)  # A, B, C
    v_ref: Mapped[float | None] = mapped_column(Float, nullable=True)  # m/s
    i_ref: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Turbine-level parameters
    rated_power: Mapped[float | None] = mapped_column(Float, nullable=True)  # kW
    rotor_diameter: Mapped[float | None] = mapped_column(Float, nullable=True)  # m
    hub_height: Mapped[float | None] = mapped_column(Float, nullable=True)  # m
    num_blades: Mapped[int | None] = mapped_column(Integer, nullable=True, default=3)
    cut_in_speed: Mapped[float | None] = mapped_column(Float, nullable=True)  # m/s
    cut_out_speed: Mapped[float | None] = mapped_column(Float, nullable=True)  # m/s
    rated_speed: Mapped[float | None] = mapped_column(Float, nullable=True)  # m/s

    # Simulation defaults
    dt: Mapped[float | None] = mapped_column(Float, nullable=True, default=0.005)  # s
    t_max: Mapped[float | None] = mapped_column(Float, nullable=True, default=660.0)  # s

    status: Mapped[ProjectStatus] = mapped_column(
        Enum(ProjectStatus, name="project_status", create_constraint=True),
        default=ProjectStatus.DRAFT,
        nullable=False,
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False
    )

    # relationships
    organization: Mapped["Organization"] = relationship(  # noqa: F821
        "Organization", back_populates="projects"
    )
    creator: Mapped["User"] = relationship("User", lazy="selectin")  # noqa: F821
    towers: Mapped[list["Tower"]] = relationship(  # noqa: F821
        "Tower", back_populates="project", cascade="all, delete-orphan", lazy="selectin"
    )
    blades: Mapped[list["Blade"]] = relationship(  # noqa: F821
        "Blade", back_populates="project", cascade="all, delete-orphan", lazy="selectin"
    )
    controllers: Mapped[list["Controller"]] = relationship(  # noqa: F821
        "Controller", back_populates="project", cascade="all, delete-orphan", lazy="selectin"
    )
    turbine_models: Mapped[list["TurbineModel"]] = relationship(  # noqa: F821
        "TurbineModel", back_populates="project", cascade="all, delete-orphan", lazy="selectin"
    )
    simulations: Mapped[list["Simulation"]] = relationship(  # noqa: F821
        "Simulation", back_populates="project", cascade="all, delete-orphan", lazy="selectin"
    )

    def __repr__(self) -> str:
        return f"<Project id={self.id} name={self.name!r}>"
