"""WebSocket connection manager with reconnection, heartbeat, and message dispatch."""

from __future__ import annotations

import asyncio
import json
import logging
import random
import time
from typing import Any

import websockets
import websockets.exceptions

from agent.config import AgentConfig

logger = logging.getLogger("windforge.agent.connection")


class ConnectionManager:
    """Manages the WebSocket lifecycle between the agent and the API server.

    Responsibilities
    ----------------
    * Connect to ``ws://{api_url}/ws/agent`` with JWT auth.
    * Reconnect with exponential back-off + jitter on transient failures.
    * Send periodic heartbeats so the server knows we are alive.
    * Dispatch incoming messages to the :class:`SimulationOrchestrator`.
    * Provide an ``async send(msg)`` method for the orchestrator to push
      progress / results back to the server.
    """

    def __init__(self, config: AgentConfig) -> None:
        self._config = config
        self._ws: websockets.WebSocketClientProtocol | None = None
        self._shutdown_event = asyncio.Event()
        self._orchestrator: Any = None  # set lazily to avoid circular import
        self._consecutive_failures = 0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def request_shutdown(self) -> None:
        """Signal the run-loop to shut down gracefully."""
        self._shutdown_event.set()

    async def run(self) -> None:
        """Main entry-point: connect, dispatch, reconnect until shutdown."""
        # Late import to avoid circular dependency at module level.
        from agent.executor import SimulationOrchestrator

        self._orchestrator = SimulationOrchestrator(
            config=self._config,
            send_fn=self.send,
        )

        while not self._shutdown_event.is_set():
            try:
                await self._connect_and_listen()
            except (
                OSError,
                websockets.exceptions.WebSocketException,
                ConnectionError,
            ) as exc:
                if self._shutdown_event.is_set():
                    break
                self._consecutive_failures += 1
                delay = self._backoff_delay()
                logger.warning(
                    "Connection lost (%s). Reconnecting in %.1fs "
                    "(attempt %d) ...",
                    exc,
                    delay,
                    self._consecutive_failures,
                )
                max_retries = self._config.reconnect_max_retries
                if max_retries >= 0 and self._consecutive_failures > max_retries:
                    logger.error(
                        "Exceeded max reconnection attempts (%d). Giving up.",
                        max_retries,
                    )
                    break
                await self._interruptible_sleep(delay)

        # Graceful teardown
        await self._orchestrator.cancel_all()
        if self._ws is not None:
            await self._ws.close()
            self._ws = None
        logger.info("Connection manager stopped.")

    async def send(self, message: dict[str, Any]) -> None:
        """Serialise *message* to JSON and send over the WebSocket.

        Silently drops messages when the socket is not connected -- the
        server will reconcile state on the next heartbeat cycle.
        """
        if self._ws is None:
            logger.debug("Cannot send -- no active WebSocket connection.")
            return
        try:
            payload = json.dumps(message)
            await self._ws.send(payload)
            logger.debug("Sent: %s", message.get("type", "?"))
        except websockets.exceptions.WebSocketException as exc:
            logger.warning("Send failed: %s", exc)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    async def _connect_and_listen(self) -> None:
        """Open the WebSocket and run the dispatch + heartbeat loops."""
        headers: dict[str, str] = {}
        if self._config.token:
            headers["Authorization"] = f"Bearer {self._config.token}"

        logger.info("Connecting to %s ...", self._config.ws_url)

        async with websockets.connect(
            self._config.ws_url,
            additional_headers=headers,
            ping_interval=None,  # we manage our own heartbeat
            close_timeout=10,
            max_size=50 * 1024 * 1024,  # 50 MiB for large input payloads
        ) as ws:
            self._ws = ws
            self._consecutive_failures = 0
            logger.info("Connected to API server.")

            # Announce readiness
            await self.send({
                "type": "agent_ready",
                "max_workers": self._config.max_workers,
                "timestamp": time.time(),
            })

            # Run heartbeat and message listener concurrently
            heartbeat_task = asyncio.create_task(
                self._heartbeat_loop(), name="heartbeat"
            )
            listener_task = asyncio.create_task(
                self._listen_loop(ws), name="listener"
            )

            done, pending = await asyncio.wait(
                {heartbeat_task, listener_task, self._shutdown_waiter()},
                return_when=asyncio.FIRST_COMPLETED,
            )

            # Cancel whichever is still running
            for task in pending:
                task.cancel()
                try:
                    await task
                except asyncio.CancelledError:
                    pass

            # Re-raise real exceptions (not shutdown / normal close)
            for task in done:
                if task.get_name() == "shutdown_waiter":
                    continue
                exc = task.exception()
                if exc is not None:
                    raise exc

        self._ws = None

    async def _shutdown_waiter(self) -> None:
        """Await the shutdown event so ``asyncio.wait`` can notice it."""
        await self._shutdown_event.wait()

    async def _heartbeat_loop(self) -> None:
        """Send heartbeat pings at the configured interval."""
        interval = self._config.heartbeat_interval
        while not self._shutdown_event.is_set():
            await self.send({
                "type": "heartbeat",
                "timestamp": time.time(),
                "active_cases": (
                    self._orchestrator.active_case_count
                    if self._orchestrator
                    else 0
                ),
            })
            await self._interruptible_sleep(interval)

    async def _listen_loop(
        self, ws: websockets.WebSocketClientProtocol
    ) -> None:
        """Receive and dispatch messages from the server."""
        async for raw in ws:
            if self._shutdown_event.is_set():
                break
            try:
                message = json.loads(raw)
            except json.JSONDecodeError:
                logger.warning("Received non-JSON message, ignoring.")
                continue

            await self._dispatch(message)

    async def _dispatch(self, message: dict[str, Any]) -> None:
        """Route an incoming message to the correct handler."""
        msg_type = message.get("type", "")
        logger.debug("Received message: %s", msg_type)

        match msg_type:
            case "job_assign":
                await self._handle_job_assign(message)
            case "cancel_case":
                await self._handle_cancel_case(message)
            case "cancel_simulation":
                await self._handle_cancel_simulation(message)
            case "ping":
                await self.send({"type": "pong", "timestamp": time.time()})
            case _:
                logger.warning("Unknown message type: %s", msg_type)

    # ------------------------------------------------------------------
    # Message handlers
    # ------------------------------------------------------------------

    async def _handle_job_assign(self, message: dict[str, Any]) -> None:
        """Process a job_assign message by forwarding to the orchestrator."""
        simulation_id: str = message.get("simulation_id", "")
        cases: list[dict[str, Any]] = message.get("cases", [])
        if not simulation_id or not cases:
            logger.error("Invalid job_assign message -- missing fields.")
            return

        logger.info(
            "Received job: simulation=%s  cases=%d",
            simulation_id,
            len(cases),
        )
        # Fire-and-forget -- the orchestrator manages its own lifecycle.
        asyncio.create_task(
            self._orchestrator.run_simulation(simulation_id, cases),
            name=f"sim-{simulation_id}",
        )

    async def _handle_cancel_case(self, message: dict[str, Any]) -> None:
        """Cancel a single running case."""
        case_id: str = message.get("case_id", "")
        if case_id:
            logger.info("Cancel requested for case %s", case_id)
            await self._orchestrator.cancel_case(case_id)

    async def _handle_cancel_simulation(self, message: dict[str, Any]) -> None:
        """Cancel all cases belonging to a simulation."""
        simulation_id: str = message.get("simulation_id", "")
        if simulation_id:
            logger.info("Cancel requested for simulation %s", simulation_id)
            await self._orchestrator.cancel_simulation(simulation_id)

    # ------------------------------------------------------------------
    # Back-off
    # ------------------------------------------------------------------

    def _backoff_delay(self) -> float:
        """Compute reconnection delay with exponential back-off and jitter."""
        base = self._config.reconnect_base_delay
        cap = self._config.reconnect_max_delay
        exp_delay = min(base * (2 ** (self._consecutive_failures - 1)), cap)
        # Full-jitter: uniform [0, exp_delay]
        return random.uniform(0, exp_delay)

    async def _interruptible_sleep(self, seconds: float) -> None:
        """Sleep that can be cut short by the shutdown event."""
        try:
            await asyncio.wait_for(
                self._shutdown_event.wait(), timeout=seconds
            )
        except asyncio.TimeoutError:
            pass
