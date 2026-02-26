"""Agent configuration dataclass and validation."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _default_workers() -> int:
    """Return number of CPUs, falling back to 1."""
    return os.cpu_count() or 1


@dataclass(frozen=True, slots=True)
class AgentConfig:
    """Immutable configuration for the WindForge agent.

    Attributes
    ----------
    api_url : str
        Base URL of the WindForge API server (HTTP scheme -- the WebSocket
        endpoint is derived automatically).
    openfast_lib_path : str
        Absolute path to the OpenFAST shared library
        (e.g. ``libopenfastlib.so``, ``libopenfastlib.dylib``, or
        ``openfastlib.dll``).
    turbsim_exe : str
        Name or path of the TurbSim executable.  If a bare name is given it
        must be on ``$PATH``.
    max_workers : int
        Maximum number of OpenFAST cases to run concurrently.  Each case
        runs in its own subprocess because the Fortran library is *not*
        thread-safe.
    work_dir : Path
        Root directory for temporary simulation files.  A sub-directory is
        created per job.
    cleanup_temp : bool
        Whether to remove temporary files after each case completes.
    heartbeat_interval : float
        Seconds between WebSocket heartbeat pings.
    reconnect_max_retries : int
        Maximum number of consecutive reconnection attempts before giving up.
        Set to ``-1`` for unlimited retries.
    reconnect_base_delay : float
        Base delay in seconds for the first reconnection attempt.  Subsequent
        attempts use exponential back-off with jitter.
    reconnect_max_delay : float
        Cap on the exponential back-off delay.
    token : str
        JWT authentication token sent in the WebSocket handshake.
    """

    api_url: str = "http://localhost:8000"
    openfast_lib_path: str = ""
    turbsim_exe: str = "turbsim"
    max_workers: int = field(default_factory=_default_workers)
    work_dir: Path = field(default_factory=lambda: Path("./windforge_work"))
    cleanup_temp: bool = True
    heartbeat_interval: float = 30.0
    reconnect_max_retries: int = -1  # unlimited
    reconnect_base_delay: float = 1.0
    reconnect_max_delay: float = 60.0
    token: str = ""

    def __post_init__(self) -> None:
        # Validate that critical paths are sane when they are provided.
        if self.openfast_lib_path:
            lib = Path(self.openfast_lib_path)
            if not lib.is_file():
                raise FileNotFoundError(
                    f"OpenFAST shared library not found: {lib}"
                )

    # ------------------------------------------------------------------
    # Derived helpers
    # ------------------------------------------------------------------

    @property
    def ws_url(self) -> str:
        """WebSocket URL derived from *api_url*."""
        base = self.api_url.rstrip("/")
        scheme = "wss" if base.startswith("https") else "ws"
        host = base.split("://", 1)[-1]
        return f"{scheme}://{host}/ws/agent"

    @property
    def resolved_work_dir(self) -> Path:
        """Return *work_dir* as an absolute, existing directory."""
        d = self.work_dir.resolve()
        d.mkdir(parents=True, exist_ok=True)
        return d
