"""
OpenFAST output file reader.

Ported from reg_tests/lib/fast_io.py with improvements for type safety,
error handling, and documentation. Handles both ASCII (.out) and binary
(.outb) output files from OpenFAST simulations.

Binary format IDs:
  1 = FileFmtID_WithTime       (compressed int16 with packed time)
  2 = FileFmtID_WithoutTime    (compressed int16, time from increment)
  3 = FileFmtID_NoCompressWithoutTime (uncompressed float64)
  4 = FileFmtID_ChanLen_In     (variable channel name length, compressed)
"""

from __future__ import annotations

import os
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import numpy as np


@dataclass
class OutputData:
    """Container for OpenFAST simulation output data.

    Attributes
    ----------
    channel_names : list[str]
        Names of all output channels (including Time).
    channel_units : list[str]
        Units for each channel (parentheses stripped).
    data : np.ndarray
        2D array of shape (num_timesteps, num_channels).
        First column is always the time channel.
    description : str
        Description string from the output file header.
    filename : str
        Base name of the source file (without extension).
    """
    channel_names: list[str]
    channel_units: list[str]
    data: np.ndarray
    description: str
    filename: str

    @property
    def num_channels(self) -> int:
        """Number of output channels including time."""
        return len(self.channel_names)

    @property
    def num_timesteps(self) -> int:
        """Number of time steps in the data."""
        return self.data.shape[0]

    @property
    def time(self) -> np.ndarray:
        """Time channel array."""
        return self.data[:, 0]

    @property
    def dt(self) -> float:
        """Time step between data points."""
        if self.num_timesteps < 2:
            return 0.0
        return float(self.time[1] - self.time[0])

    def channel(self, name: str) -> np.ndarray:
        """Get data for a specific channel by name.

        Parameters
        ----------
        name : str
            Channel name (case-insensitive).

        Returns
        -------
        np.ndarray
            1D array of channel data.

        Raises
        ------
        KeyError
            If the channel name is not found.
        """
        name_upper = name.upper()
        for idx, ch_name in enumerate(self.channel_names):
            if ch_name.upper() == name_upper:
                return self.data[:, idx]
        raise KeyError(
            f"Channel '{name}' not found. Available channels: "
            f"{', '.join(self.channel_names[:10])}..."
        )

    def channels_dict(self) -> dict[str, np.ndarray]:
        """Return all channels as a dictionary.

        Returns
        -------
        dict[str, np.ndarray]
            Mapping from channel name to data array.
        """
        return {
            name: self.data[:, idx]
            for idx, name in enumerate(self.channel_names)
        }


class OutputReader:
    """Reader for OpenFAST ASCII (.out) and binary (.outb) output files.

    This class provides methods to load simulation results from OpenFAST
    output files. Both ASCII and binary formats are supported with
    automatic detection.

    Usage
    -----
    >>> reader = OutputReader()
    >>> result = reader.load("simulation.outb")
    >>> print(result.channel_names)
    >>> blade_moment = result.channel("RootMxc1")
    """

    # Binary format identifiers used by OpenFAST
    _FMTID_WITH_TIME = 1
    _FMTID_WITHOUT_TIME = 2
    _FMTID_NO_COMPRESS = 3
    _FMTID_CHAN_LEN_IN = 4

    def load(self, filename: Union[str, Path]) -> OutputData:
        """Load an OpenFAST output file (ASCII or binary).

        Automatically detects the file format based on extension and
        content. Files ending in '.outb' are treated as binary; files
        ending in '.out' are treated as ASCII. If a '.out' file cannot
        be read as text, it falls back to binary parsing.

        Parameters
        ----------
        filename : str or Path
            Path to the OpenFAST output file.

        Returns
        -------
        OutputData
            Parsed output data with channel names, units, and values.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        ValueError
            If the file format cannot be determined or is corrupt.
        """
        filepath = Path(filename)
        if not filepath.is_file():
            raise FileNotFoundError(f"Output file not found: {filepath}")

        if filepath.suffix.lower() == ".outb":
            return self._load_binary(filepath)

        # Try ASCII first for .out files
        try:
            with open(filepath, "r") as f:
                f.readline()  # Test readability
            return self._load_ascii(filepath)
        except UnicodeDecodeError:
            return self._load_binary(filepath)

    def _load_ascii(self, filepath: Path) -> OutputData:
        """Load an ASCII .out file.

        The expected format is:
        - Lines 1-4: Header comment lines
        - Line 5: Description
        - Line 6: Empty or comment
        - Line 7: Channel names (space/tab delimited)
        - Line 8: Channel units (space/tab delimited, in parentheses)
        - Lines 9+: Numeric data

        Parameters
        ----------
        filepath : Path
            Path to the ASCII output file.

        Returns
        -------
        OutputData
            Parsed output data.
        """
        with open(filepath, "r") as f:
            header_lines = [f.readline() for _ in range(8)]

            description = header_lines[4].strip()
            channel_names = header_lines[6].split()
            # Strip parentheses from units
            channel_units = [
                unit.strip("()")
                for unit in header_lines[7].split()
            ]

            # Read all data lines into a numpy array
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

    def _load_binary(self, filepath: Path) -> OutputData:
        """Load a binary .outb file.

        Supports all four OpenFAST binary format IDs:
        - ID 1: Compressed (int16) with packed time channel
        - ID 2: Compressed (int16) with time from increment
        - ID 3: Uncompressed (float64) without time channel
        - ID 4: Variable channel name length, compressed (int16)

        Parameters
        ----------
        filepath : Path
            Path to the binary output file.

        Returns
        -------
        OutputData
            Parsed output data.

        Raises
        ------
        ValueError
            If the binary format is unrecognized or data is corrupt.
        """
        with open(filepath, "rb") as fid:
            # Read format ID
            file_id = self._fread(fid, 1, "int16")[0]

            # Determine channel name length
            if file_id == self._FMTID_CHAN_LEN_IN:
                len_name = self._fread(fid, 1, "int16")[0]
            else:
                len_name = 10  # Default OpenFAST channel name width

            # Read dimensions
            num_out_chans = self._fread(fid, 1, "int32")[0]
            num_timesteps = self._fread(fid, 1, "int32")[0]

            # Read time scaling/offset
            if file_id == self._FMTID_WITH_TIME:
                time_scl = self._fread(fid, 1, "float64")[0]
                time_off = self._fread(fid, 1, "float64")[0]
            else:
                time_out1 = self._fread(fid, 1, "float64")[0]
                time_incr = self._fread(fid, 1, "float64")[0]

            # Read channel scaling (not present for uncompressed format)
            if file_id != self._FMTID_NO_COMPRESS:
                col_scl = np.array(
                    self._fread(fid, num_out_chans, "float32")
                )
                col_off = np.array(
                    self._fread(fid, num_out_chans, "float32")
                )

            # Read description string
            len_desc = self._fread(fid, 1, "int32")[0]
            desc_ascii = self._fread(fid, len_desc, "uint8")
            description = "".join(map(chr, desc_ascii)).strip()

            # Read channel names
            channel_names: list[str] = []
            for _ in range(num_out_chans + 1):
                name_ascii = self._fread(fid, len_name, "uint8")
                channel_names.append("".join(map(chr, name_ascii)).strip())

            # Read channel units
            channel_units: list[str] = []
            for _ in range(num_out_chans + 1):
                unit_ascii = self._fread(fid, len_name, "uint8")
                raw_unit = "".join(map(chr, unit_ascii)).strip()
                channel_units.append(raw_unit.strip("()"))

            # Read time data (only for format ID 1)
            if file_id == self._FMTID_WITH_TIME:
                packed_time = np.array(
                    self._fread(fid, num_timesteps, "int32"),
                    dtype=np.float64,
                )
                if len(packed_time) < num_timesteps:
                    raise ValueError(
                        f"Incomplete time data: read {len(packed_time)} "
                        f"of {num_timesteps} values"
                    )

            # Read channel data
            num_pts = num_timesteps * num_out_chans
            if file_id == self._FMTID_NO_COMPRESS:
                packed_data = self._fread(fid, num_pts, "float64")
            else:
                packed_data = self._fread(fid, num_pts, "int16")

            if len(packed_data) < num_pts:
                raise ValueError(
                    f"Incomplete data: read {len(packed_data)} "
                    f"of {num_pts} values"
                )

        # Reshape and scale data
        raw = np.array(packed_data).reshape(num_timesteps, num_out_chans)

        if file_id == self._FMTID_NO_COMPRESS:
            data_channels = raw
        else:
            # Apply scaling: real_value = (packed - offset) / scale
            data_channels = (raw - col_off) / col_scl

        # Reconstruct time vector
        if file_id == self._FMTID_WITH_TIME:
            time = (packed_time - time_off) / time_scl
        else:
            time = time_out1 + time_incr * np.arange(num_timesteps)

        # Combine time and channel data
        data = np.column_stack([time.reshape(-1, 1), data_channels])

        return OutputData(
            channel_names=channel_names,
            channel_units=channel_units,
            data=data,
            description=description,
            filename=filepath.stem,
        )

    @staticmethod
    def _fread(fid, n: int, dtype: str) -> tuple:
        """Read n elements of the given type from a binary file.

        Parameters
        ----------
        fid : file-like
            Open binary file handle.
        n : int
            Number of elements to read.
        dtype : str
            Data type: 'uint8', 'int16', 'int32', 'float32', 'float64'.

        Returns
        -------
        tuple
            The unpacked values.
        """
        fmt_map = {
            "uint8": ("B", 1),
            "int16": ("h", 2),
            "int32": ("i", 4),
            "float32": ("f", 4),
            "float64": ("d", 8),
        }
        fmt_char, nbytes = fmt_map[dtype]
        raw = fid.read(nbytes * n)
        if len(raw) < nbytes * n:
            # Return what we could read
            actual_n = len(raw) // nbytes
            if actual_n == 0:
                return ()
            return struct.unpack(fmt_char * actual_n, raw[: actual_n * nbytes])
        return struct.unpack(fmt_char * n, raw)
