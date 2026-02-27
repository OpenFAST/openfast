"""Turbine Model CRUD endpoints."""

from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.database import get_db
from app.dependencies import get_current_user
from app.models.components import TurbineModel
from app.models.project import Project
from app.models.user import User
from app.schemas.components import TurbineModelCreate, TurbineModelResponse, TurbineModelUpdate

router = APIRouter(prefix="/projects/{project_id}/turbine-models", tags=["turbine-models"])


async def _verify_project(project_id: UUID, org_id: UUID, db: AsyncSession) -> Project:
    result = await db.execute(
        select(Project).where(Project.id == project_id, Project.org_id == org_id)
    )
    project = result.scalar_one_or_none()
    if project is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")
    return project


async def _get_model_or_404(model_id: UUID, project_id: UUID, db: AsyncSession) -> TurbineModel:
    result = await db.execute(
        select(TurbineModel).where(TurbineModel.id == model_id, TurbineModel.project_id == project_id)
    )
    model = result.scalar_one_or_none()
    if model is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Turbine model not found")
    return model


@router.get("", response_model=list[TurbineModelResponse])
async def list_turbine_models(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    result = await db.execute(
        select(TurbineModel)
        .where(TurbineModel.project_id == project_id)
        .order_by(TurbineModel.created_at.desc())
    )
    return [TurbineModelResponse.model_validate(m) for m in result.scalars().all()]


@router.post("", response_model=TurbineModelResponse, status_code=status.HTTP_201_CREATED)
async def create_turbine_model(
    project_id: UUID,
    body: TurbineModelCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    model = TurbineModel(project_id=project_id, **body.model_dump(exclude_unset=True))
    db.add(model)
    await db.flush()
    await db.refresh(model)
    return TurbineModelResponse.model_validate(model)


@router.get("/{model_id}", response_model=TurbineModelResponse)
async def get_turbine_model(
    project_id: UUID,
    model_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    model = await _get_model_or_404(model_id, project_id, db)
    return TurbineModelResponse.model_validate(model)


@router.put("/{model_id}", response_model=TurbineModelResponse)
async def update_turbine_model(
    project_id: UUID,
    model_id: UUID,
    body: TurbineModelUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    model = await _get_model_or_404(model_id, project_id, db)

    for field, value in body.model_dump(exclude_unset=True).items():
        setattr(model, field, value)

    model.version += 1
    await db.flush()
    await db.refresh(model)
    return TurbineModelResponse.model_validate(model)


@router.delete("/{model_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_turbine_model(
    project_id: UUID,
    model_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_db),
):
    await _verify_project(project_id, current_user.org_id, db)
    model = await _get_model_or_404(model_id, project_id, db)
    await db.delete(model)
    await db.flush()
