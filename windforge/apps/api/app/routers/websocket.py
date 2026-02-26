"""WebSocket endpoint for real-time simulation updates.

Protocol
--------
- Client connects to  ``/ws/{simulation_id}?token=<jwt>``
- After authentication the server pushes JSON messages whenever a simulation
  case makes progress, completes, or fails.
- The client may send JSON commands back:
  ``{"action": "subscribe", "channels": ["progress", "logs"]}``
  ``{"action": "cancel_case", "case_id": "<uuid>"}``
"""

from __future__ import annotations

import asyncio
import json
import logging
from datetime import datetime, timezone
from uuid import UUID

from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from jose import JWTError, jwt
from sqlalchemy import select

from app.config import settings
from app.database import async_session_factory
from app.models.simulation import CaseStatus, Simulation, SimulationCase

logger = logging.getLogger("windforge.ws")

router = APIRouter(tags=["websocket"])

# In-memory registry of connected WebSocket clients per simulation.
# In production this would be backed by Redis pub/sub or similar.
_connections: dict[UUID, set[WebSocket]] = {}

# Global event bus — agents push updates here, the WS loop forwards them.
# Keyed by simulation_id.  Each value is an asyncio.Queue of dicts.
_event_queues: dict[UUID, asyncio.Queue[dict]] = {}


def get_event_queue(simulation_id: UUID) -> asyncio.Queue[dict]:
    """Return (or create) the event queue for a given simulation."""
    if simulation_id not in _event_queues:
        _event_queues[simulation_id] = asyncio.Queue(maxsize=10_000)
    return _event_queues[simulation_id]


async def publish_event(simulation_id: UUID, event: dict) -> None:
    """Push an event into the simulation's queue so connected clients receive it."""
    q = get_event_queue(simulation_id)
    try:
        q.put_nowait(event)
    except asyncio.QueueFull:
        logger.warning("Event queue full for simulation %s — dropping event", simulation_id)


async def _authenticate(token: str) -> UUID | None:
    """Validate a JWT token and return the user_id, or None on failure."""
    try:
        payload = jwt.decode(token, settings.SECRET_KEY, algorithms=[settings.JWT_ALGORITHM])
        user_id_str: str | None = payload.get("sub")
        if user_id_str is None:
            return None
        return UUID(user_id_str)
    except (JWTError, ValueError):
        return None


async def _verify_simulation_access(
    simulation_id: UUID, user_id: UUID
) -> bool:
    """Check that the user's org owns the simulation."""
    async with async_session_factory() as db:
        result = await db.execute(
            select(Simulation).where(Simulation.id == simulation_id)
        )
        sim = result.scalar_one_or_none()
        if sim is None:
            return False

        # Minimal access check — the user must belong to the same org as the
        # project that owns the simulation.
        from app.models.project import Project
        from app.models.user import User

        proj_result = await db.execute(
            select(Project).where(Project.id == sim.project_id)
        )
        project = proj_result.scalar_one_or_none()
        if project is None:
            return False

        user_result = await db.execute(select(User).where(User.id == user_id))
        user = user_result.scalar_one_or_none()
        if user is None:
            return False

        return user.org_id == project.org_id


async def _forward_events(
    ws: WebSocket,
    simulation_id: UUID,
    subscribed_channels: set[str],
) -> None:
    """Background task that reads from the event queue and sends to the client."""
    q = get_event_queue(simulation_id)
    while True:
        event = await q.get()
        channel = event.get("channel", "progress")
        if subscribed_channels and channel not in subscribed_channels:
            continue
        try:
            await ws.send_json(event)
        except Exception:
            break


async def _handle_client_message(
    data: dict,
    subscribed_channels: set[str],
    simulation_id: UUID,
) -> dict | None:
    """Process a client command and return an optional response."""
    action = data.get("action")

    if action == "subscribe":
        channels = data.get("channels", [])
        if isinstance(channels, list):
            subscribed_channels.clear()
            subscribed_channels.update(channels)
        return {"type": "subscribed", "channels": list(subscribed_channels)}

    if action == "cancel_case":
        case_id_str = data.get("case_id")
        if not case_id_str:
            return {"type": "error", "message": "case_id required"}
        try:
            case_id = UUID(case_id_str)
        except ValueError:
            return {"type": "error", "message": "Invalid case_id"}

        async with async_session_factory() as db:
            result = await db.execute(
                select(SimulationCase).where(
                    SimulationCase.id == case_id,
                    SimulationCase.simulation_id == simulation_id,
                )
            )
            case = result.scalar_one_or_none()
            if case is None:
                return {"type": "error", "message": "Case not found"}
            if case.status in (CaseStatus.COMPLETED, CaseStatus.CANCELLED):
                return {
                    "type": "error",
                    "message": f"Cannot cancel case in '{case.status.value}' state",
                }
            case.status = CaseStatus.CANCELLED
            case.completed_at = datetime.now(timezone.utc)
            await db.commit()

        # Publish cancellation event
        await publish_event(
            simulation_id,
            {
                "channel": "progress",
                "type": "case_cancelled",
                "case_id": str(case_id),
                "timestamp": datetime.now(timezone.utc).isoformat(),
            },
        )
        return {"type": "case_cancelled", "case_id": str(case_id)}

    return {"type": "error", "message": f"Unknown action: {action}"}


@router.websocket("/ws/{simulation_id}")
async def websocket_endpoint(websocket: WebSocket, simulation_id: UUID, token: str = ""):
    """WebSocket endpoint for real-time simulation monitoring."""

    # --- authenticate ---
    if not token:
        await websocket.close(code=4001, reason="Missing authentication token")
        return

    user_id = await _authenticate(token)
    if user_id is None:
        await websocket.close(code=4001, reason="Invalid authentication token")
        return

    # --- verify access ---
    has_access = await _verify_simulation_access(simulation_id, user_id)
    if not has_access:
        await websocket.close(code=4003, reason="Access denied")
        return

    await websocket.accept()

    # Register connection
    if simulation_id not in _connections:
        _connections[simulation_id] = set()
    _connections[simulation_id].add(websocket)

    subscribed_channels: set[str] = {"progress", "logs"}

    # Start background task to forward events
    forward_task = asyncio.create_task(
        _forward_events(websocket, simulation_id, subscribed_channels)
    )

    try:
        # Send initial connection confirmation
        await websocket.send_json(
            {
                "type": "connected",
                "simulation_id": str(simulation_id),
                "channels": list(subscribed_channels),
            }
        )

        # Listen for client commands
        while True:
            raw = await websocket.receive_text()
            try:
                data = json.loads(raw)
            except json.JSONDecodeError:
                await websocket.send_json(
                    {"type": "error", "message": "Invalid JSON"}
                )
                continue

            response = await _handle_client_message(
                data, subscribed_channels, simulation_id
            )
            if response:
                await websocket.send_json(response)

    except WebSocketDisconnect:
        logger.info("Client disconnected from simulation %s", simulation_id)
    except Exception:
        logger.exception("WebSocket error for simulation %s", simulation_id)
    finally:
        forward_task.cancel()
        _connections[simulation_id].discard(websocket)
        if not _connections[simulation_id]:
            del _connections[simulation_id]
            # Clean up empty event queues
            _event_queues.pop(simulation_id, None)
