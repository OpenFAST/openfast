"""Simulation endpoints: CRUD, execution control, results retrieval."""

from datetime import datetime, timezone
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.database import get_db
from app.dependencies import get_current_user
from app.models.project import Project
from app.models.simulation import (
    CaseStatus,
    DLCDefinition,
    ResultsDEL,
    ResultsExtreme,
    ResultsStatistics,
    Simulation,
    SimulationCase,
    SimulationStatus,
)
from app.models.user import User
from app.schemas.simulation import (
    DLCDefinitionCreate,
    DLCDefinitionResponse,
    ResultsDELResponse,
    ResultsExtremeResponse,
    ResultsStatisticsResponse,
    SimulationCaseResponse,
    SimulationCreate,
    SimulationResponse,
    SimulationWithProgress,
)

router = APIRouter(prefix="/projects/{project_id}/simulations", tags=["simulations"])


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
async def _verify_project(project_id: UUID, org_id: UUID, db: AsyncSession) -> Project:
    result = await db.execute(
        select(Project).where(Project.id == project_id, Project.org_id == org_id)
    )
    project = result.scalar_one_or_none()
    if project is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")
    return project


async def _get_simulation_or_404(
    sim_id: UUID, project_id: UUID, db: AsyncSession
) -> Simulation:
    result = await db.execute(
        select(Simulation).where(
            Simulation.id == sim_id, Simulation.project_id == project_id
        )
    )
    sim = result.scalar_one_or_none()
    if sim is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Simulation not found"
        )
    return sim


def _compute_progress(sim: Simulation) -> SimulationWithProgress:
    """Build an extended response with progress metrics."""
    total = max(sim.total_cases, 1)
    completed = sim.completed_cases + sim.failed_cases
    pct = (completed / total) * 100.0

    est_remaining: float | None = None
    if sim.started_at and completed > 0 and pct < 100.0:
        elapsed = (datetime.now(timezone.utc) - sim.started_at).total_seconds()
        per_case = elapsed / completed
        remaining_cases = total - completed
        est_remaining = per_case * remaining_cases

    return SimulationWithProgress(
        id=sim.id,
        project_id=sim.project_id,
        dlc_definition_id=sim.dlc_definition_id,
        turbine_model_id=sim.turbine_model_id,
        name=sim.name,
        status=sim.status.value,
        total_cases=sim.total_cases,
        completed_cases=sim.completed_cases,
        failed_cases=sim.failed_cases,
        agent_id=sim.agent_id,
        started_at=sim.started_at,
        completed_at=sim.completed_at,
        created_at=sim.created_at,
        progress_percent=round(pct, 2),
        estimated_remaining_seconds=round(est_remaining, 1) if est_remaining else None,
    )


def _expand_dlc_cases(dlc_def: DLCDefinition) -> list[dict]:
    """Expand the DLC case matrix into individual simulation cases."""
    cases: list[dict] = []
    for spec in (dlc_def.dlc_cases or []):
        dlc_number = spec.get("dlc_number", "1.1") if isinstance(spec, dict) else spec.dlc_number
        wind_speeds = spec.get("wind_speeds", []) if isinstance(spec, dict) else spec.wind_speeds
        seeds = spec.get("seeds", 6) if isinstance(spec, dict) else spec.seeds
        yaw_misalignments = (
            spec.get("yaw_misalignments", [0.0])
            if isinstance(spec, dict)
            else spec.yaw_misalignments
        )

        for ws in wind_speeds:
            for seed in range(1, seeds + 1):
                for yaw in yaw_misalignments:
                    cases.append(
                        {
                            "dlc_number": dlc_number,
                            "wind_speed": ws,
                            "seed_number": seed,
                            "yaw_misalignment": yaw,
                        }
                    )
    return cases


# ---------------------------------------------------------------------------
# DLC Definition endpoints (nested under simulations for convenience)
# ---------------------------------------------------------------------------
dlc_router = APIRouter(prefix="/projects/{project_id}/dlc-definitions", tags=["dlc"])


@dlc_router.post("/", response_model=DLCDefinitionResponse, status_code=status.HTTP_201_CREATED)
async def create_dlc_definition(
    project_id: UUID,
    body: DLCDefinitionCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)

    # Serialize nested schemas
    dlc_cases_data = None
    if body.dlc_cases:
        dlc_cases_data = [c.model_dump() for c in body.dlc_cases]

    turbsim_data = body.turbsim_params.model_dump() if body.turbsim_params else None

    # Count total cases
    total = 0
    for spec in (body.dlc_cases or []):
        total += len(spec.wind_speeds) * spec.seeds * len(spec.yaw_misalignments)

    dlc = DLCDefinition(
        project_id=project_id,
        turbine_model_id=body.turbine_model_id,
        name=body.name,
        dlc_cases=dlc_cases_data,
        turbsim_params=turbsim_data,
        total_case_count=total,
    )
    db.add(dlc)
    await db.flush()
    await db.refresh(dlc)
    return DLCDefinitionResponse.model_validate(dlc)


@dlc_router.get("/", response_model=list[DLCDefinitionResponse])
async def list_dlc_definitions(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    result = await db.execute(
        select(DLCDefinition)
        .where(DLCDefinition.project_id == project_id)
        .order_by(DLCDefinition.created_at.desc())
    )
    return [DLCDefinitionResponse.model_validate(d) for d in result.scalars().all()]


@dlc_router.get("/{dlc_id}", response_model=DLCDefinitionResponse)
async def get_dlc_definition(
    project_id: UUID,
    dlc_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    result = await db.execute(
        select(DLCDefinition).where(
            DLCDefinition.id == dlc_id, DLCDefinition.project_id == project_id
        )
    )
    dlc = result.scalar_one_or_none()
    if dlc is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="DLC definition not found")
    return DLCDefinitionResponse.model_validate(dlc)


# ---------------------------------------------------------------------------
# Simulation endpoints
# ---------------------------------------------------------------------------

@router.get("/", response_model=list[SimulationResponse])
async def list_simulations(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    result = await db.execute(
        select(Simulation)
        .where(Simulation.project_id == project_id)
        .order_by(Simulation.created_at.desc())
    )
    return [SimulationResponse.model_validate(s) for s in result.scalars().all()]


@router.post("/", response_model=SimulationWithProgress, status_code=status.HTTP_201_CREATED)
async def create_simulation(
    project_id: UUID,
    body: SimulationCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    """Create a simulation batch with individual cases expanded from the DLC definition."""
    await _verify_project(project_id, current_user.org_id, db)

    # Fetch the DLC definition
    result = await db.execute(
        select(DLCDefinition).where(
            DLCDefinition.id == body.dlc_definition_id,
            DLCDefinition.project_id == project_id,
        )
    )
    dlc_def = result.scalar_one_or_none()
    if dlc_def is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="DLC definition not found"
        )

    # Expand cases
    expanded = _expand_dlc_cases(dlc_def)

    sim = Simulation(
        project_id=project_id,
        dlc_definition_id=body.dlc_definition_id,
        turbine_model_id=body.turbine_model_id,
        name=body.name,
        total_cases=len(expanded),
    )
    db.add(sim)
    await db.flush()

    # Bulk-insert simulation cases
    for case_data in expanded:
        case = SimulationCase(simulation_id=sim.id, **case_data)
        db.add(case)

    await db.flush()
    await db.refresh(sim)
    return _compute_progress(sim)


@router.get("/{simulation_id}", response_model=SimulationWithProgress)
async def get_simulation(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    sim = await _get_simulation_or_404(simulation_id, project_id, db)
    return _compute_progress(sim)


@router.post("/{simulation_id}/start", response_model=SimulationWithProgress)
async def start_simulation(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    """Mark a simulation as ready to run (agent picks it up)."""
    await _verify_project(project_id, current_user.org_id, db)
    sim = await _get_simulation_or_404(simulation_id, project_id, db)

    if sim.status != SimulationStatus.PENDING:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"Cannot start simulation in '{sim.status.value}' state",
        )

    sim.status = SimulationStatus.GENERATING_WIND
    sim.started_at = datetime.now(timezone.utc)
    await db.flush()
    await db.refresh(sim)
    return _compute_progress(sim)


@router.post("/{simulation_id}/cancel", response_model=SimulationWithProgress)
async def cancel_simulation(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    """Cancel a running or pending simulation."""
    await _verify_project(project_id, current_user.org_id, db)
    sim = await _get_simulation_or_404(simulation_id, project_id, db)

    if sim.status in (SimulationStatus.COMPLETED, SimulationStatus.CANCELLED):
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"Cannot cancel simulation in '{sim.status.value}' state",
        )

    sim.status = SimulationStatus.CANCELLED
    sim.completed_at = datetime.now(timezone.utc)

    # Cancel any pending / running cases
    result = await db.execute(
        select(SimulationCase).where(
            SimulationCase.simulation_id == sim.id,
            SimulationCase.status.in_([CaseStatus.PENDING, CaseStatus.RUNNING, CaseStatus.GENERATING_WIND]),
        )
    )
    for case in result.scalars().all():
        case.status = CaseStatus.CANCELLED

    await db.flush()
    await db.refresh(sim)
    return _compute_progress(sim)


@router.get("/{simulation_id}/cases", response_model=list[SimulationCaseResponse])
async def list_cases(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    await _get_simulation_or_404(simulation_id, project_id, db)

    result = await db.execute(
        select(SimulationCase)
        .where(SimulationCase.simulation_id == simulation_id)
        .order_by(SimulationCase.wind_speed, SimulationCase.seed_number)
    )
    return [SimulationCaseResponse.model_validate(c) for c in result.scalars().all()]


@router.get(
    "/{simulation_id}/results/statistics",
    response_model=list[ResultsStatisticsResponse],
)
async def get_statistics(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    await _get_simulation_or_404(simulation_id, project_id, db)

    result = await db.execute(
        select(ResultsStatistics)
        .where(ResultsStatistics.simulation_id == simulation_id)
        .order_by(ResultsStatistics.wind_speed)
    )
    return [ResultsStatisticsResponse.model_validate(r) for r in result.scalars().all()]


@router.get("/{simulation_id}/results/del", response_model=list[ResultsDELResponse])
async def get_del_results(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    await _get_simulation_or_404(simulation_id, project_id, db)

    result = await db.execute(
        select(ResultsDEL).where(ResultsDEL.simulation_id == simulation_id)
    )
    return [ResultsDELResponse.model_validate(r) for r in result.scalars().all()]


@router.get(
    "/{simulation_id}/results/extreme",
    response_model=list[ResultsExtremeResponse],
)
async def get_extreme_results(
    project_id: UUID,
    simulation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    await _get_simulation_or_404(simulation_id, project_id, db)

    result = await db.execute(
        select(ResultsExtreme).where(ResultsExtreme.simulation_id == simulation_id)
    )
    return [ResultsExtremeResponse.model_validate(r) for r in result.scalars().all()]
