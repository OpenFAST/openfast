"""Simulation, DLC, case, and results schemas."""

from datetime import datetime
from uuid import UUID

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# DLC case specification (nested inside DLCDefinition)
# ---------------------------------------------------------------------------
class DLCCaseSpec(BaseModel):
    """Specification for a single DLC within a definition."""

    dlc_number: str = Field(..., description="e.g. '1.1', '1.3', '6.1'")
    wind_speeds: list[float] = Field(..., description="List of wind speeds (m/s)")
    seeds: int = Field(default=6, ge=1, description="Number of random seeds")
    yaw_misalignments: list[float] = Field(
        default=[0.0], description="Yaw misalignment angles (deg)"
    )


class TurbSimParamsSchema(BaseModel):
    """TurbSim configuration parameters."""

    turbulence_model: str = Field(default="IECKAI", description="e.g. IECKAI, IECVKM, GP_LLJ")
    iec_standard: str = Field(default="1-ED3")
    iec_turbc: str = Field(default="B", description="IEC turbulence category")
    grid_height: float = Field(default=150.0, gt=0)
    grid_width: float = Field(default=150.0, gt=0)
    num_grid_z: int = Field(default=25, ge=3)
    num_grid_y: int = Field(default=25, ge=3)
    time_step: float = Field(default=0.05, gt=0)
    analysis_time: float = Field(default=660.0, gt=0)
    ref_height: float = Field(default=90.0, gt=0)


# ---------------------------------------------------------------------------
# DLC Definition
# ---------------------------------------------------------------------------
class DLCDefinitionCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    turbine_model_id: UUID
    dlc_cases: list[DLCCaseSpec] | None = None
    turbsim_params: TurbSimParamsSchema | None = None


class DLCDefinitionUpdate(BaseModel):
    name: str | None = None
    dlc_cases: list[DLCCaseSpec] | None = None
    turbsim_params: TurbSimParamsSchema | None = None


class DLCDefinitionResponse(BaseModel):
    id: UUID
    project_id: UUID
    turbine_model_id: UUID
    name: str
    dlc_cases: list[DLCCaseSpec] | None = None
    turbsim_params: TurbSimParamsSchema | None = None
    total_case_count: int
    status: str
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
class SimulationCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    dlc_definition_id: UUID
    turbine_model_id: UUID


class SimulationResponse(BaseModel):
    id: UUID
    project_id: UUID
    dlc_definition_id: UUID
    turbine_model_id: UUID
    name: str
    status: str
    total_cases: int
    completed_cases: int
    failed_cases: int
    agent_id: str | None = None
    started_at: datetime | None = None
    completed_at: datetime | None = None
    created_at: datetime

    model_config = {"from_attributes": True}


class SimulationWithProgress(SimulationResponse):
    """Extended response with computed progress fields."""

    progress_percent: float = 0.0
    estimated_remaining_seconds: float | None = None


# ---------------------------------------------------------------------------
# SimulationCase
# ---------------------------------------------------------------------------
class SimulationCaseResponse(BaseModel):
    id: UUID
    simulation_id: UUID
    dlc_number: str
    wind_speed: float
    seed_number: int
    yaw_misalignment: float
    wind_field_path: str | None = None
    input_files: dict | None = None
    status: str
    progress_percent: float
    current_time: float
    error_message: str | None = None
    started_at: datetime | None = None
    completed_at: datetime | None = None
    wall_time_seconds: float | None = None
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
class ChannelStatistic(BaseModel):
    """Statistics for a single output channel."""

    min: float
    max: float
    mean: float
    std: float
    abs_max: float | None = None
    integrated: float | None = None


class ResultsStatisticsResponse(BaseModel):
    id: UUID
    simulation_case_id: UUID
    dlc_number: str
    wind_speed: float
    channel_statistics: dict | None = None
    created_at: datetime

    model_config = {"from_attributes": True}


class ResultsDELResponse(BaseModel):
    id: UUID
    simulation_id: UUID
    del_results: dict | None = None
    m_exponent: float
    n_equivalent: float
    created_at: datetime

    model_config = {"from_attributes": True}


class ResultsExtremeResponse(BaseModel):
    id: UUID
    simulation_id: UUID
    extreme_loads: dict | None = None
    created_at: datetime

    model_config = {"from_attributes": True}
