"""Simulation orchestrator and per-case executor.

The :class:`SimulationOrchestrator` is the bridge between the WebSocket
connection layer and the subprocess-based OpenFAST runners.  It receives
simulation jobs (batches of cases), manages the two-phase execution
(TurbSim then OpenFAST), and streams progress back to the API server.
"""

from __future__ import annotations

import asyncio
import base64
import logging
import os
import shutil
import time
import uuid
from concurrent.futures import Future, ProcessPoolExecutor
from multiprocessing import Queue
from pathlib import Path
from typing import Any, Callable, Coroutine

from agent.config import AgentConfig
from agent.openfast_runner import (
    ProgressUpdate,
    RunResult,
    run_openfast_case,
    run_turbsim,
)
from agent.output_processor import process_output

logger = logging.getLogger("windforge.agent.executor")

# Type alias for the async send function provided by the connection manager.
SendFn = Callable[[dict[str, Any]], Coroutine[Any, Any, None]]


# ======================================================================
# Single-case executor
# ======================================================================

class CaseExecutor:
    """Prepares, runs, and post-processes a single OpenFAST case.

    Each case is described by a dict from the ``job_assign`` message.  The
    expected schema is::

        {
            "case_id": "uuid",
            "simulation_id": "uuid",
            "name": "case_01",
            "main_fst": "Test.fst",
            "input_files": {
                "Test.fst": "<base64 content>",
                "Test_ElastoDyn.dat": "<base64 content>",
                ...
            },
            "wind": {
                "type": "turbsim",          # or "steady", "uniform"
                "turbsim_input": "<base64>" # only if type == turbsim
            },
            "compute_fatigue": false
        }

    Parameters
    ----------
    case_spec : dict
        Case specification from the server.
    config : AgentConfig
        Agent configuration.
    send_fn : SendFn
        Coroutine used to push messages back over the WebSocket.
    """

    def __init__(
        self,
        case_spec: dict[str, Any],
        config: AgentConfig,
        send_fn: SendFn,
    ) -> None:
        self.case_id: str = case_spec["case_id"]
        self.simulation_id: str = case_spec.get("simulation_id", "")
        self.case_name: str = case_spec.get("name", self.case_id[:8])
        self.main_fst: str = case_spec.get("main_fst", "")
        self.input_files: dict[str, str] = case_spec.get("input_files", {})
        self.wind_spec: dict[str, Any] = case_spec.get("wind", {})
        self.compute_fatigue: bool = case_spec.get("compute_fatigue", False)

        self._config = config
        self._send = send_fn
        self._cancelled = False

        # Per-case working directory
        slug = f"{self.case_name}_{uuid.uuid4().hex[:8]}"
        self._work_dir = config.resolved_work_dir / self.simulation_id / slug

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def cancel(self) -> None:
        """Mark this case as cancelled (best-effort)."""
        self._cancelled = True

    async def execute(self, pool: ProcessPoolExecutor) -> dict[str, Any] | None:
        """Run the full lifecycle of a single case.

        Returns
        -------
        dict or None
            The processed output statistics, or ``None`` on cancellation /
            error (error messages are sent to the server directly).
        """
        t0 = time.monotonic()

        try:
            # 1. Write input files to disk
            self._write_input_files()

            # 2. Notify server that we have started
            await self._send({
                "type": "case_started",
                "case_id": self.case_id,
                "simulation_id": self.simulation_id,
                "timestamp": time.time(),
            })

            # 3. TurbSim phase (if needed)
            if self.wind_spec.get("type") == "turbsim":
                await self._run_turbsim(pool)

            if self._cancelled:
                return None

            # 4. Run OpenFAST
            result = await self._run_openfast(pool)

            if not result.success:
                await self._send({
                    "type": "case_error",
                    "case_id": self.case_id,
                    "simulation_id": self.simulation_id,
                    "error": result.error,
                    "timestamp": time.time(),
                })
                return None

            # 5. Process outputs
            stats: dict[str, Any] = {}
            if result.output_file:
                stats = process_output(
                    result.output_file,
                    compute_fatigue=self.compute_fatigue,
                )

            elapsed = time.monotonic() - t0

            await self._send({
                "type": "case_complete",
                "case_id": self.case_id,
                "simulation_id": self.simulation_id,
                "elapsed": elapsed,
                "total_steps": result.total_steps,
                "statistics": stats,
                "timestamp": time.time(),
            })

            return stats

        except asyncio.CancelledError:
            logger.info("Case %s cancelled.", self.case_id)
            return None

        except Exception as exc:
            logger.exception("Case %s failed with unexpected error.", self.case_id)
            await self._send({
                "type": "case_error",
                "case_id": self.case_id,
                "simulation_id": self.simulation_id,
                "error": str(exc),
                "timestamp": time.time(),
            })
            return None

        finally:
            if self._config.cleanup_temp and self._work_dir.exists():
                try:
                    shutil.rmtree(self._work_dir)
                except OSError as exc:
                    logger.warning(
                        "Failed to clean up %s: %s", self._work_dir, exc
                    )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _write_input_files(self) -> None:
        """Decode base64 input files and write them to the work directory."""
        self._work_dir.mkdir(parents=True, exist_ok=True)
        for rel_path, b64_content in self.input_files.items():
            dest = self._work_dir / rel_path
            dest.parent.mkdir(parents=True, exist_ok=True)
            content = base64.b64decode(b64_content)
            dest.write_bytes(content)

        logger.debug(
            "Wrote %d input files to %s", len(self.input_files), self._work_dir
        )

    async def _run_turbsim(self, pool: ProcessPoolExecutor) -> None:
        """Run TurbSim in a subprocess and wait for completion."""
        ts_input_b64 = self.wind_spec.get("turbsim_input", "")
        if not ts_input_b64:
            raise ValueError(
                "Wind type is 'turbsim' but no turbsim_input provided."
            )

        ts_input_file = self._work_dir / "TurbSim.inp"
        ts_input_file.write_bytes(base64.b64decode(ts_input_b64))

        loop = asyncio.get_running_loop()
        bts_path = await loop.run_in_executor(
            pool,
            run_turbsim,
            str(ts_input_file),
            self._config.turbsim_exe,
            None,  # no progress queue for TurbSim in this flow
        )
        logger.info("TurbSim wind file ready: %s", bts_path)

    async def _run_openfast(self, pool: ProcessPoolExecutor) -> RunResult:
        """Run OpenFAST in a subprocess, streaming progress updates."""
        loop = asyncio.get_running_loop()
        progress_queue: Queue[ProgressUpdate | RunResult] = Queue()

        # Submit the work to the process pool
        future: Future[RunResult] = pool.submit(
            run_openfast_case,
            lib_path=self._config.openfast_lib_path,
            input_dir=str(self._work_dir),
            main_fst_file=self.main_fst,
            progress_queue=progress_queue,
            progress_interval=0.5,
        )

        # Poll the progress queue while waiting for completion
        while not future.done():
            await asyncio.sleep(0.1)

            # Drain all queued progress updates
            while not progress_queue.empty():
                try:
                    update = progress_queue.get_nowait()
                except Exception:
                    break

                if isinstance(update, ProgressUpdate):
                    await self._send({
                        "type": "case_progress",
                        "case_id": self.case_id,
                        "simulation_id": self.simulation_id,
                        "percent": round(update.percent, 2),
                        "current_time": round(update.current_time, 4),
                        "total_time": round(update.total_time, 4),
                        "timestamp": time.time(),
                    })

            if self._cancelled:
                future.cancel()
                return RunResult(success=False, error="Cancelled by user")

        # Drain any remaining updates
        while not progress_queue.empty():
            try:
                progress_queue.get_nowait()
            except Exception:
                break

        return future.result()


# ======================================================================
# Simulation orchestrator
# ======================================================================

class SimulationOrchestrator:
    """Manages batch execution of simulation jobs.

    Responsibilities
    ----------------
    * Receive simulation jobs from the connection manager.
    * Phase 1: generate all TurbSim wind fields in parallel.
    * Phase 2: run all OpenFAST cases in parallel (up to worker limit).
    * Track overall progress and handle cancellation.
    * Report results back via the ``send_fn`` callback.
    """

    def __init__(self, config: AgentConfig, send_fn: SendFn) -> None:
        self._config = config
        self._send = send_fn

        # The process pool is shared across all simulations.  Each worker
        # is a separate OS process, which is necessary because the OpenFAST
        # Fortran library is not thread-safe.
        self._pool = ProcessPoolExecutor(
            max_workers=config.max_workers,
            mp_context=None,  # use default start method (fork/spawn)
        )

        # Track active simulations and cases
        self._active_cases: dict[str, CaseExecutor] = {}
        self._simulation_tasks: dict[str, asyncio.Task[None]] = {}

    @property
    def active_case_count(self) -> int:
        """Number of cases currently executing."""
        return len(self._active_cases)

    # ------------------------------------------------------------------
    # Public API called by the connection manager
    # ------------------------------------------------------------------

    async def run_simulation(
        self,
        simulation_id: str,
        cases: list[dict[str, Any]],
    ) -> None:
        """Execute a batch of cases for *simulation_id*.

        This method is called as a fire-and-forget ``asyncio.Task`` from
        the connection manager.
        """
        logger.info(
            "Starting simulation %s with %d cases (max %d workers).",
            simulation_id,
            len(cases),
            self._config.max_workers,
        )

        # Inject simulation_id into each case spec if not already present
        for case in cases:
            case.setdefault("simulation_id", simulation_id)

        # ---- Phase 1: identify TurbSim cases and generate wind files ----
        turbsim_cases = [
            c for c in cases
            if c.get("wind", {}).get("type") == "turbsim"
        ]
        if turbsim_cases:
            logger.info(
                "Phase 1: generating %d TurbSim wind fields.",
                len(turbsim_cases),
            )
            # TurbSim generation is handled inside CaseExecutor, but we
            # create the executors up front so we can parallelise.

        # ---- Phase 2: run all OpenFAST cases ----------------------------
        semaphore = asyncio.Semaphore(self._config.max_workers)

        async def _run_with_semaphore(case_spec: dict[str, Any]) -> dict[str, Any] | None:
            async with semaphore:
                executor = CaseExecutor(case_spec, self._config, self._send)
                self._active_cases[executor.case_id] = executor
                try:
                    return await executor.execute(self._pool)
                finally:
                    self._active_cases.pop(executor.case_id, None)

        tasks = [
            asyncio.create_task(
                _run_with_semaphore(case),
                name=f"case-{case.get('case_id', 'unknown')}",
            )
            for case in cases
        ]

        completed = 0
        failed = 0
        results: list[dict[str, Any] | None] = []

        for coro in asyncio.as_completed(tasks):
            try:
                result = await coro
                results.append(result)
                if result is not None:
                    completed += 1
                else:
                    failed += 1
            except Exception as exc:
                logger.exception("Case task raised: %s", exc)
                failed += 1
                results.append(None)

        # ---- Report simulation completion --------------------------------
        await self._send({
            "type": "simulation_complete",
            "simulation_id": simulation_id,
            "total_cases": len(cases),
            "completed": completed,
            "failed": failed,
            "timestamp": time.time(),
        })

        logger.info(
            "Simulation %s complete: %d/%d succeeded, %d failed.",
            simulation_id,
            completed,
            len(cases),
            failed,
        )

    async def cancel_case(self, case_id: str) -> None:
        """Cancel a single active case."""
        executor = self._active_cases.get(case_id)
        if executor:
            executor.cancel()
            logger.info("Cancellation signal sent to case %s.", case_id)
        else:
            logger.warning("Case %s not found among active cases.", case_id)

    async def cancel_simulation(self, simulation_id: str) -> None:
        """Cancel all cases belonging to *simulation_id*."""
        cancelled = 0
        for case_id, executor in list(self._active_cases.items()):
            if executor.simulation_id == simulation_id:
                executor.cancel()
                cancelled += 1
        logger.info(
            "Cancelled %d cases for simulation %s.", cancelled, simulation_id
        )

    async def cancel_all(self) -> None:
        """Cancel every active case (used during shutdown)."""
        for executor in list(self._active_cases.values()):
            executor.cancel()
        self._pool.shutdown(wait=False, cancel_futures=True)
        logger.info("All cases cancelled; process pool shut down.")
