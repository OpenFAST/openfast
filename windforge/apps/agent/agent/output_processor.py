"""Post-simulation output processing.

Reads OpenFAST ``.out`` / ``.outb`` files (ported from the API's
``output_reader.py``), computes per-channel statistics, and optionally
computes Damage Equivalent Loads (DEL) for fatigue channels.
"""

from __future__ import annotations

import logging
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger("windforge.agent.output")


# ======================================================================
# Output file reader  (mirrors API's output_reader.py)
# ======================================================================

@dataclass(slots=True)
class OutputData:
    """Container for parsed OpenFAST output data."""

    channel_names: list[str]
    channel_units: list[str]
    data: np.ndarray  # shape (timesteps, channels)
    description: str
    filename: str

    @property
    def num_channels(self) -> int:
        return len(self.channel_names)

    @property
    def num_timesteps(self) -> int:
        return self.data.shape[0]

    @property
    def time(self) -> np.ndarray:
        return self.data[:, 0]

    @property
    def dt(self) -> float:
        if self.num_timesteps < 2:
            return 0.0
        return float(self.time[1] - self.time[0])


# Binary format IDs
_FMTID_WITH_TIME = 1
_FMTID_WITHOUT_TIME = 2
_FMTID_NO_COMPRESS = 3
_FMTID_CHAN_LEN_IN = 4


def _fread(fid, n: int, dtype: str) -> tuple:
    """Read *n* elements of *dtype* from a binary file."""
    fmt_map = {
        "uint8": ("B", 1),
        "int16": ("h", 2),
        "int32": ("i", 4),
        "float32": ("f", 4),
        "float64": ("d", 8),
    }
    fmt_char, nbytes = fmt_map[dtype]
    raw = fid.read(nbytes * n)
    actual_n = len(raw) // nbytes
    if actual_n == 0:
        return ()
    return struct.unpack(fmt_char * actual_n, raw[: actual_n * nbytes])


def load_output(filename: str | Path) -> OutputData:
    """Load an OpenFAST output file (ASCII ``.out`` or binary ``.outb``).

    Parameters
    ----------
    filename : str or Path
        Path to the output file.

    Returns
    -------
    OutputData
        Parsed output data.
    """
    filepath = Path(filename)
    if not filepath.is_file():
        raise FileNotFoundError(f"Output file not found: {filepath}")

    if filepath.suffix.lower() == ".outb":
        return _load_binary(filepath)

    try:
        with open(filepath) as f:
            f.readline()
        return _load_ascii(filepath)
    except UnicodeDecodeError:
        return _load_binary(filepath)


def _load_ascii(filepath: Path) -> OutputData:
    """Parse an ASCII ``.out`` file."""
    with open(filepath) as f:
        header_lines = [f.readline() for _ in range(8)]
        description = header_lines[4].strip()
        channel_names = header_lines[6].split()
        channel_units = [u.strip("()") for u in header_lines[7].split()]
        data_lines = f.readlines()

    if not data_lines:
        data = np.empty((0, len(channel_names)), dtype=np.float64)
    else:
        data = np.array(
            [line.split() for line in data_lines if line.strip()],
            dtype=np.float64,
        )

    return OutputData(
        channel_names=channel_names,
        channel_units=channel_units,
        data=data,
        description=description,
        filename=filepath.stem,
    )


def _load_binary(filepath: Path) -> OutputData:
    """Parse a binary ``.outb`` file supporting all four format IDs."""
    with open(filepath, "rb") as fid:
        file_id = _fread(fid, 1, "int16")[0]

        len_name = (
            _fread(fid, 1, "int16")[0]
            if file_id == _FMTID_CHAN_LEN_IN
            else 10
        )

        num_out_chans = _fread(fid, 1, "int32")[0]
        num_timesteps = _fread(fid, 1, "int32")[0]

        if file_id == _FMTID_WITH_TIME:
            time_scl = _fread(fid, 1, "float64")[0]
            time_off = _fread(fid, 1, "float64")[0]
        else:
            time_out1 = _fread(fid, 1, "float64")[0]
            time_incr = _fread(fid, 1, "float64")[0]

        if file_id != _FMTID_NO_COMPRESS:
            col_scl = np.array(_fread(fid, num_out_chans, "float32"))
            col_off = np.array(_fread(fid, num_out_chans, "float32"))

        len_desc = _fread(fid, 1, "int32")[0]
        desc_ascii = _fread(fid, len_desc, "uint8")
        description = "".join(map(chr, desc_ascii)).strip()

        channel_names: list[str] = []
        for _ in range(num_out_chans + 1):
            name_ascii = _fread(fid, len_name, "uint8")
            channel_names.append("".join(map(chr, name_ascii)).strip())

        channel_units: list[str] = []
        for _ in range(num_out_chans + 1):
            unit_ascii = _fread(fid, len_name, "uint8")
            channel_units.append(
                "".join(map(chr, unit_ascii)).strip().strip("()")
            )

        if file_id == _FMTID_WITH_TIME:
            packed_time = np.array(
                _fread(fid, num_timesteps, "int32"), dtype=np.float64
            )

        num_pts = num_timesteps * num_out_chans
        if file_id == _FMTID_NO_COMPRESS:
            packed_data = _fread(fid, num_pts, "float64")
        else:
            packed_data = _fread(fid, num_pts, "int16")

    raw = np.array(packed_data).reshape(num_timesteps, num_out_chans)

    if file_id == _FMTID_NO_COMPRESS:
        data_channels = raw
    else:
        data_channels = (raw - col_off) / col_scl

    if file_id == _FMTID_WITH_TIME:
        time_vec = (packed_time - time_off) / time_scl
    else:
        time_vec = time_out1 + time_incr * np.arange(num_timesteps)

    data = np.column_stack([time_vec.reshape(-1, 1), data_channels])

    return OutputData(
        channel_names=channel_names,
        channel_units=channel_units,
        data=data,
        description=description,
        filename=filepath.stem,
    )


# ======================================================================
# Statistics
# ======================================================================

@dataclass(slots=True)
class ChannelStats:
    """Per-channel summary statistics."""

    name: str
    unit: str
    min: float
    max: float
    mean: float
    std: float
    abs_max: float


def compute_statistics(output: OutputData) -> list[ChannelStats]:
    """Compute per-channel statistics for all output channels.

    Parameters
    ----------
    output : OutputData
        Parsed output from :func:`load_output`.

    Returns
    -------
    list[ChannelStats]
        One entry per channel (including Time).
    """
    stats: list[ChannelStats] = []
    for idx, (name, unit) in enumerate(
        zip(output.channel_names, output.channel_units)
    ):
        col = output.data[:, idx]
        stats.append(ChannelStats(
            name=name,
            unit=unit,
            min=float(np.nanmin(col)),
            max=float(np.nanmax(col)),
            mean=float(np.nanmean(col)),
            std=float(np.nanstd(col)),
            abs_max=float(np.nanmax(np.abs(col))),
        ))
    return stats


# ======================================================================
# DEL (Damage Equivalent Load) -- optional fatigue analysis
# ======================================================================

# Common Woehler exponents for wind-turbine materials.
_DEFAULT_WOEHLER: dict[str, float] = {
    "steel": 4.0,
    "composite": 10.0,
    "default": 4.0,
}

# Channel-name pattern to Woehler exponent mapping.  Blade root / tower
# base moments are typically composite / steel respectively.
_FATIGUE_CHANNEL_PATTERNS: dict[str, float] = {
    "RootMxc": 10.0,  # blade root flap
    "RootMyc": 10.0,  # blade root edge
    "RootMzc": 10.0,
    "TwrBsMxt": 4.0,   # tower base
    "TwrBsMyt": 4.0,
    "YawBrMxp": 4.0,
    "YawBrMyp": 4.0,
}


def compute_del(
    output: OutputData,
    channels: list[str] | None = None,
    woehler: float | None = None,
    n_eq: float = 1e7,
    t_lifetime: float = 600.0,
) -> dict[str, float]:
    """Compute Damage Equivalent Loads for selected fatigue channels.

    Uses the ``fatpack`` library for rainflow counting.  If ``fatpack`` is
    not installed the function returns an empty dict with a warning.

    Parameters
    ----------
    output : OutputData
        Parsed simulation output.
    channels : list[str], optional
        Channel names to analyse.  If ``None``, all channels matching
        known fatigue patterns are used.
    woehler : float, optional
        Woehler (S-N) exponent.  If ``None``, the exponent is inferred
        from the channel name using ``_FATIGUE_CHANNEL_PATTERNS``.
    n_eq : float
        Number of equivalent cycles (default 10^7).
    t_lifetime : float
        Duration of the time series in seconds.  Defaults to 600 s
        (10-minute IEC standard).

    Returns
    -------
    dict[str, float]
        Mapping ``channel_name -> DEL``.
    """
    try:
        import fatpack  # type: ignore[import-untyped]
    except ImportError:
        logger.warning(
            "fatpack not installed -- skipping DEL computation. "
            "Install with: pip install windforge-agent[fatigue]"
        )
        return {}

    if channels is None:
        channels = [
            name
            for name in output.channel_names
            if any(pat in name for pat in _FATIGUE_CHANNEL_PATTERNS)
        ]

    results: dict[str, float] = {}
    name_upper_map = {n.upper(): i for i, n in enumerate(output.channel_names)}

    for ch_name in channels:
        idx = name_upper_map.get(ch_name.upper())
        if idx is None:
            logger.warning("DEL channel '%s' not found -- skipping.", ch_name)
            continue

        signal = output.data[:, idx]

        # Determine Woehler exponent
        m = woehler
        if m is None:
            m = _DEFAULT_WOEHLER["default"]
            for pat, exp in _FATIGUE_CHANNEL_PATTERNS.items():
                if pat in ch_name:
                    m = exp
                    break

        try:
            # Rainflow counting
            reversals, _ = fatpack.find_reversals(signal)
            cycles = fatpack.find_rainflow_ranges(reversals)

            if len(cycles) == 0:
                results[ch_name] = 0.0
                continue

            # DEL = (sum(S_i^m) / n_eq)^(1/m)
            del_value = float(
                (np.sum(cycles ** m) / n_eq) ** (1.0 / m)
            )
            results[ch_name] = del_value
        except Exception as exc:
            logger.warning("DEL computation failed for '%s': %s", ch_name, exc)
            results[ch_name] = float("nan")

    return results


# ======================================================================
# High-level convenience
# ======================================================================

def process_output(
    output_file: str,
    compute_fatigue: bool = False,
) -> dict[str, Any]:
    """Load an output file, compute statistics, and optionally DEL.

    Returns a JSON-serialisable dict suitable for sending back to the
    API server.
    """
    output = load_output(output_file)
    stats = compute_statistics(output)

    result: dict[str, Any] = {
        "filename": output.filename,
        "num_timesteps": output.num_timesteps,
        "num_channels": output.num_channels,
        "dt": output.dt,
        "channels": [
            {
                "name": s.name,
                "unit": s.unit,
                "min": s.min,
                "max": s.max,
                "mean": s.mean,
                "std": s.std,
                "abs_max": s.abs_max,
            }
            for s in stats
        ],
    }

    if compute_fatigue:
        del_results = compute_del(output)
        result["del"] = {k: v for k, v in del_results.items()}

    return result
