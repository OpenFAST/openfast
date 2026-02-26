"""Isolated OpenFAST and TurbSim runners designed to execute in subprocesses.

Every function in this module is intended to be called via
:class:`concurrent.futures.ProcessPoolExecutor` so that the non-thread-safe
OpenFAST Fortran library is loaded in a dedicated address space.

The public functions communicate progress back to the parent process through
a :class:`multiprocessing.Queue`.
"""

from __future__ import annotations

import ctypes
import logging
import math
import os
import subprocess
import sys
import time
from ctypes import (
    POINTER,
    byref,
    c_bool,
    c_char,
    c_double,
    c_float,
    c_int,
    create_string_buffer,
)
from dataclasses import dataclass, field
from multiprocessing import Queue
from pathlib import Path
from typing import Any

logger = logging.getLogger("windforge.agent.runner")

# OpenFAST shared-library constants (must match FAST_Library.f90)
INTF_STR_LEN = 1025
NUM_FIXED_INPUTS = 51


# ------------------------------------------------------------------
# Data containers for inter-process communication
# ------------------------------------------------------------------

@dataclass(slots=True)
class ProgressUpdate:
    """Lightweight progress payload sent through the queue."""

    current_time: float
    total_time: float
    percent: float
    step: int
    total_steps: int


@dataclass(slots=True)
class RunResult:
    """Final result payload sent through the queue on completion."""

    success: bool
    output_file: str = ""
    channel_names: list[str] = field(default_factory=list)
    error: str = ""
    elapsed: float = 0.0
    total_steps: int = 0


# ------------------------------------------------------------------
# OpenFAST runner (runs inside a subprocess)
# ------------------------------------------------------------------

def run_openfast_case(
    lib_path: str,
    input_dir: str,
    main_fst_file: str,
    progress_queue: Queue[ProgressUpdate | RunResult] | None = None,
    progress_interval: float = 0.25,
) -> RunResult:
    """Run a single OpenFAST case using the shared library via ctypes.

    This function is designed to be the *target* of a
    :class:`ProcessPoolExecutor` submission.  It loads the Fortran shared
    library into the current process, runs the simulation, and returns a
    :class:`RunResult`.

    Parameters
    ----------
    lib_path : str
        Absolute path to the OpenFAST shared library.
    input_dir : str
        Directory containing all input files for this case.
    main_fst_file : str
        Name of the primary ``.fst`` file (relative to *input_dir*).
    progress_queue : multiprocessing.Queue, optional
        Queue for streaming :class:`ProgressUpdate` messages back to the
        parent process.  ``None`` disables progress reporting.
    progress_interval : float
        Minimum seconds between progress updates to avoid flooding.

    Returns
    -------
    RunResult
        Outcome of the simulation.
    """
    t_start = time.monotonic()
    fst_path = str(Path(input_dir) / main_fst_file)

    # Ensure CWD is the input directory so that relative paths in .fst
    # files resolve correctly.
    original_cwd = os.getcwd()
    os.chdir(input_dir)

    try:
        result = _run_openfast_inner(
            lib_path, fst_path, progress_queue, progress_interval
        )
        result.elapsed = time.monotonic() - t_start
        return result
    except Exception as exc:
        elapsed = time.monotonic() - t_start
        return RunResult(success=False, error=str(exc), elapsed=elapsed)
    finally:
        os.chdir(original_cwd)


def _run_openfast_inner(
    lib_path: str,
    fst_path: str,
    progress_queue: Queue[ProgressUpdate | RunResult] | None,
    progress_interval: float,
) -> RunResult:
    """Core simulation loop -- separated for clarity."""
    lib = ctypes.CDLL(lib_path)

    # ---- Bind Fortran routines -----------------------------------------
    lib.FAST_AllocateTurbines.argtypes = [
        POINTER(c_int), POINTER(c_int), POINTER(c_char),
    ]
    lib.FAST_AllocateTurbines.restype = c_int

    lib.FAST_Sizes.argtypes = [
        POINTER(c_int), POINTER(c_char),
        POINTER(c_int), POINTER(c_int),
        POINTER(c_double), POINTER(c_double), POINTER(c_double),
        POINTER(c_int), POINTER(c_char), POINTER(c_char),
        POINTER(c_double), POINTER(c_double),
    ]
    lib.FAST_Sizes.restype = c_int

    lib.FAST_Start.argtypes = [
        POINTER(c_int), POINTER(c_int), POINTER(c_int),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int), POINTER(c_char),
    ]
    lib.FAST_Start.restype = c_int

    lib.FAST_Update.argtypes = [
        POINTER(c_int), POINTER(c_int), POINTER(c_int),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_bool), POINTER(c_int), POINTER(c_char),
    ]
    lib.FAST_Update.restype = c_int

    lib.FAST_End.argtypes = [POINTER(c_int), POINTER(c_bool)]
    lib.FAST_End.restype = c_int

    lib.FAST_DeallocateTurbines.argtypes = [POINTER(c_int), POINTER(c_char)]
    lib.FAST_DeallocateTurbines.restype = c_int

    # ---- Constants & buffers -------------------------------------------
    n_turbines = c_int(1)
    i_turb = c_int(0)
    abort_error_level = c_int(4)  # ErrID_Fatal
    num_outs = c_int(0)
    dt = c_double(0.0)
    dt_out = c_double(0.0)
    t_max = c_double(0.0)
    num_inputs = c_int(NUM_FIXED_INPUTS)
    inp_array = (c_double * NUM_FIXED_INPUTS)(0.0)
    end_early = c_bool(False)

    err_stat = c_int(0)
    err_msg = create_string_buffer(INTF_STR_LEN)
    input_file_buf = create_string_buffer(
        os.path.abspath(fst_path).encode("utf-8")
    )
    channel_names_buf = create_string_buffer(20 * 4000)

    def _fatal(status: c_int) -> bool:
        return status.value >= abort_error_level.value

    def _cleanup() -> None:
        try:
            lib.FAST_End(byref(i_turb), byref(c_bool(False)))
            lib.FAST_DeallocateTurbines(byref(err_stat), err_msg)
        except Exception:
            pass

    # ---- FAST_AllocateTurbines -----------------------------------------
    lib.FAST_AllocateTurbines(
        byref(n_turbines), byref(err_stat), err_msg,
    )
    if _fatal(err_stat):
        raise RuntimeError(
            f"FAST_AllocateTurbines failed ({err_stat.value}): "
            f"{err_msg.value.decode(errors='replace')}"
        )

    # ---- FAST_Sizes ----------------------------------------------------
    lib.FAST_Sizes(
        byref(i_turb), input_file_buf,
        byref(abort_error_level), byref(num_outs),
        byref(dt), byref(dt_out), byref(t_max),
        byref(err_stat), err_msg, channel_names_buf,
        None, None,  # optional TMax, InitInpAry
    )
    if _fatal(err_stat):
        _cleanup()
        raise RuntimeError(
            f"FAST_Sizes failed ({err_stat.value}): "
            f"{err_msg.value.decode(errors='replace')}"
        )

    raw_names = channel_names_buf.value.split()
    channel_names = [n.decode("utf-8") for n in raw_names] if raw_names else []

    total_time_steps = math.ceil(t_max.value / dt.value) + 1
    total_output_steps = math.ceil(t_max.value / dt_out.value) + 1
    output_frequency = round(dt_out.value / dt.value)

    # Allocate output buffer  (rows = output steps, cols = num_outs)
    import numpy as np

    output_values = np.zeros(
        (total_output_steps, num_outs.value), dtype=c_double, order="C"
    )

    # ---- FAST_Start ----------------------------------------------------
    lib.FAST_Start(
        byref(i_turb), byref(num_inputs), byref(num_outs),
        byref(inp_array),
        output_values[0].ctypes.data_as(POINTER(c_double)),
        byref(err_stat), err_msg,
    )
    if _fatal(err_stat):
        _cleanup()
        raise RuntimeError(
            f"FAST_Start failed ({err_stat.value}): "
            f"{err_msg.value.decode(errors='replace')}"
        )

    # ---- Time-stepping loop --------------------------------------------
    i_out = 1
    last_progress_time = time.monotonic()

    for i in range(1, total_time_steps):
        lib.FAST_Update(
            byref(i_turb), byref(num_inputs), byref(num_outs),
            byref(inp_array),
            output_values[i_out].ctypes.data_as(POINTER(c_double)),
            byref(end_early), byref(err_stat), err_msg,
        )

        if i % output_frequency == 0:
            i_out += 1

        if _fatal(err_stat):
            _cleanup()
            raise RuntimeError(
                f"FAST_Update failed at step {i} ({err_stat.value}): "
                f"{err_msg.value.decode(errors='replace')}"
            )

        if end_early.value:
            break

        # Stream progress
        if progress_queue is not None:
            now = time.monotonic()
            if now - last_progress_time >= progress_interval:
                current_time = i * dt.value
                pct = min(current_time / t_max.value * 100.0, 100.0)
                progress_queue.put(ProgressUpdate(
                    current_time=current_time,
                    total_time=t_max.value,
                    percent=pct,
                    step=i,
                    total_steps=total_time_steps,
                ))
                last_progress_time = now

    # ---- Cleanup -------------------------------------------------------
    _cleanup()

    # Determine output file path (.out or .outb)
    fst_stem = Path(fst_path).stem
    fst_dir = Path(fst_path).parent
    output_file = ""
    for suffix in (".outb", ".out"):
        candidate = fst_dir / (fst_stem + suffix)
        if candidate.is_file():
            output_file = str(candidate)
            break

    return RunResult(
        success=True,
        output_file=output_file,
        channel_names=channel_names,
        total_steps=total_time_steps,
    )


# ------------------------------------------------------------------
# TurbSim runner (subprocess, NOT ctypes)
# ------------------------------------------------------------------

def run_turbsim(
    input_file: str,
    turbsim_exe: str = "turbsim",
    progress_queue: Queue[dict[str, Any]] | None = None,
) -> str:
    """Run TurbSim to generate a wind field ``.bts`` file.

    Parameters
    ----------
    input_file : str
        Path to the TurbSim ``.inp`` file.
    turbsim_exe : str
        Executable name or full path.
    progress_queue : multiprocessing.Queue, optional
        Queue for progress updates (TurbSim progress is coarse-grained).

    Returns
    -------
    str
        Path to the generated ``.bts`` file.

    Raises
    ------
    RuntimeError
        If TurbSim exits with a non-zero return code.
    FileNotFoundError
        If the expected ``.bts`` file is not produced.
    """
    input_path = Path(input_file).resolve()
    working_dir = input_path.parent

    logger.info("Running TurbSim: %s", input_path.name)

    proc = subprocess.Popen(
        [turbsim_exe, str(input_path)],
        cwd=str(working_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    output_lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        stripped = line.strip()
        if stripped:
            output_lines.append(stripped)
            logger.debug("TurbSim: %s", stripped)
            # Attempt to parse progress from TurbSim output
            if progress_queue is not None and "%" in stripped:
                try:
                    pct_str = stripped.split("%")[0].split()[-1]
                    pct = float(pct_str)
                    progress_queue.put({
                        "type": "turbsim_progress",
                        "percent": pct,
                        "file": input_path.name,
                    })
                except (ValueError, IndexError):
                    pass

    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(
            f"TurbSim failed (exit {returncode}):\n"
            + "\n".join(output_lines[-20:])
        )

    # Find the generated .bts file
    bts_file = input_path.with_suffix(".bts")
    if not bts_file.is_file():
        # Some TurbSim configurations use a different output name; search
        # for any .bts file in the directory.
        bts_files = list(working_dir.glob("*.bts"))
        if not bts_files:
            raise FileNotFoundError(
                f"TurbSim did not produce a .bts file in {working_dir}"
            )
        bts_file = max(bts_files, key=lambda p: p.stat().st_mtime)

    logger.info("TurbSim complete: %s", bts_file.name)
    return str(bts_file)
