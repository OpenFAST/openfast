"""Command-line interface for the WindForge agent."""

from __future__ import annotations

import asyncio
import logging
import os
import signal
import sys
from pathlib import Path

import click

from agent import __version__
from agent.config import AgentConfig


def _configure_logging(verbose: bool) -> None:
    """Set up structured logging to stderr."""
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)-8s] %(name)s: %(message)s"
    logging.basicConfig(level=level, format=fmt, stream=sys.stderr)
    # Quiet down noisy third-party loggers
    logging.getLogger("websockets").setLevel(logging.WARNING)


@click.group()
@click.version_option(version=__version__, prog_name="windforge-agent")
def cli() -> None:
    """WindForge Agent -- local OpenFAST simulation executor."""


@cli.command()
@click.option(
    "--api-url",
    default="http://localhost:8000",
    envvar="WINDFORGE_API_URL",
    show_default=True,
    help="WindForge API server URL.",
)
@click.option(
    "--openfast-lib",
    required=True,
    envvar="WINDFORGE_OPENFAST_LIB",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Path to the OpenFAST shared library (.so / .dylib / .dll).",
)
@click.option(
    "--workers",
    default=os.cpu_count() or 1,
    show_default=True,
    type=click.IntRange(min=1),
    help="Max parallel OpenFAST cases (each runs in a subprocess).",
)
@click.option(
    "--work-dir",
    default="./windforge_work",
    show_default=True,
    type=click.Path(file_okay=False, resolve_path=True),
    help="Working directory for temporary simulation files.",
)
@click.option(
    "--token",
    default="",
    envvar="WINDFORGE_TOKEN",
    help="JWT authentication token.",
)
@click.option(
    "--turbsim-exe",
    default="turbsim",
    show_default=True,
    envvar="WINDFORGE_TURBSIM_EXE",
    help="Name or path to the TurbSim executable.",
)
@click.option(
    "--no-cleanup",
    is_flag=True,
    default=False,
    help="Keep temporary simulation files after completion.",
)
@click.option("-v", "--verbose", is_flag=True, default=False, help="Enable debug logging.")
def connect(
    api_url: str,
    openfast_lib: str,
    workers: int,
    work_dir: str,
    token: str,
    turbsim_exe: str,
    no_cleanup: bool,
    verbose: bool,
) -> None:
    """Connect to the WindForge API and start processing simulation jobs."""
    _configure_logging(verbose)
    logger = logging.getLogger("windforge.agent.cli")

    config = AgentConfig(
        api_url=api_url,
        openfast_lib_path=openfast_lib,
        max_workers=workers,
        work_dir=Path(work_dir),
        token=token,
        turbsim_exe=turbsim_exe,
        cleanup_temp=not no_cleanup,
    )

    logger.info("WindForge Agent v%s starting", __version__)
    logger.info("API server : %s", config.api_url)
    logger.info("OpenFAST lib: %s", config.openfast_lib_path)
    logger.info("Workers     : %d", config.max_workers)
    logger.info("Work dir    : %s", config.resolved_work_dir)

    # Import here to keep CLI startup fast and allow validation above to fail early.
    from agent.connection import ConnectionManager

    manager = ConnectionManager(config)

    # ------------------------------------------------------------------
    # Graceful shutdown via SIGINT / SIGTERM
    # ------------------------------------------------------------------
    loop = asyncio.new_event_loop()

    def _request_shutdown(sig: signal.Signals) -> None:
        logger.info("Received %s -- shutting down gracefully", sig.name)
        loop.call_soon_threadsafe(manager.request_shutdown)

    for sig in (signal.SIGINT, signal.SIGTERM):
        try:
            loop.add_signal_handler(sig, _request_shutdown, sig)
        except NotImplementedError:
            # add_signal_handler is not available on Windows
            signal.signal(sig, lambda s, _f: _request_shutdown(signal.Signals(s)))

    try:
        loop.run_until_complete(manager.run())
    except KeyboardInterrupt:
        logger.info("Interrupted -- shutting down")
    finally:
        # Allow pending tasks to complete
        pending = asyncio.all_tasks(loop)
        if pending:
            loop.run_until_complete(asyncio.gather(*pending, return_exceptions=True))
        loop.run_until_complete(loop.shutdown_asyncgens())
        loop.close()
        logger.info("Agent stopped.")


def main() -> None:
    """Package entry-point."""
    cli()


if __name__ == "__main__":
    main()
