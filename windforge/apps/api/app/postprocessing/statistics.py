"""
Time series statistics calculator for OpenFAST simulation results.

Computes standard descriptive statistics (min, max, mean, std, abs_max)
for each output channel, and provides aggregation across multiple
simulation cases for reporting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np


@dataclass
class ChannelStatistics:
    """Statistics for a single output channel.

    Attributes
    ----------
    channel_name : str
        Name of the output channel.
    channel_unit : str
        Unit of the output channel.
    minimum : float
        Minimum value across all time steps.
    maximum : float
        Maximum value across all time steps.
    mean : float
        Arithmetic mean across all time steps.
    std : float
        Standard deviation across all time steps.
    abs_max : float
        Maximum of the absolute value (max(|min|, |max|)).
    time_of_min : float
        Time at which the minimum occurred (s).
    time_of_max : float
        Time at which the maximum occurred (s).
    """
    channel_name: str
    channel_unit: str = ""
    minimum: float = 0.0
    maximum: float = 0.0
    mean: float = 0.0
    std: float = 0.0
    abs_max: float = 0.0
    time_of_min: float = 0.0
    time_of_max: float = 0.0


@dataclass
class CaseStatistics:
    """Aggregated statistics for a single simulation case.

    Attributes
    ----------
    case_id : str
        Identifier for the simulation case.
    channels : dict[str, ChannelStatistics]
        Statistics per channel, keyed by channel name.
    """
    case_id: str
    channels: dict[str, ChannelStatistics] = field(default_factory=dict)


class StatisticsCalculator:
    """Calculates descriptive statistics for OpenFAST time series data.

    This calculator computes min, max, mean, std, and absolute max
    for each channel in a simulation output, optionally skipping
    an initial transient period.

    Usage
    -----
    >>> calc = StatisticsCalculator(t_start=30.0)
    >>> stats = calc.calculate(data, time, channel_names, channel_units)
    >>> print(stats["RootMxc1"].maximum)
    """

    def __init__(self, t_start: float = 0.0) -> None:
        """Initialize the statistics calculator.

        Parameters
        ----------
        t_start : float
            Time (s) at which to begin statistics computation,
            allowing for transient spinup to be excluded.
        """
        self.t_start = t_start

    def calculate(
        self,
        data: np.ndarray,
        time: np.ndarray,
        channel_names: list[str],
        channel_units: Optional[list[str]] = None,
        case_id: str = "",
    ) -> CaseStatistics:
        """Calculate statistics for all channels in a time series.

        Parameters
        ----------
        data : np.ndarray
            2D array of shape (num_timesteps, num_channels).
            Does NOT include time as the first column; time is
            provided separately.
        time : np.ndarray
            1D array of time values (s), same length as data rows.
        channel_names : list[str]
            Names for each data column.
        channel_units : list[str], optional
            Units for each channel. If None, empty strings are used.
        case_id : str
            Identifier for this case.

        Returns
        -------
        CaseStatistics
            Statistics for every channel.
        """
        if channel_units is None:
            channel_units = [""] * len(channel_names)

        # Apply t_start filter
        mask = time >= self.t_start
        if not np.any(mask):
            # If no data after t_start, use all data
            mask = np.ones(len(time), dtype=bool)

        filtered_time = time[mask]
        filtered_data = data[mask]

        result = CaseStatistics(case_id=case_id)

        for idx, (name, unit) in enumerate(zip(channel_names, channel_units)):
            if idx >= filtered_data.shape[1]:
                break

            col = filtered_data[:, idx]

            if len(col) == 0:
                result.channels[name] = ChannelStatistics(
                    channel_name=name, channel_unit=unit
                )
                continue

            min_val = float(np.nanmin(col))
            max_val = float(np.nanmax(col))
            mean_val = float(np.nanmean(col))
            std_val = float(np.nanstd(col))
            abs_max_val = max(abs(min_val), abs(max_val))

            min_idx = int(np.nanargmin(col))
            max_idx = int(np.nanargmax(col))
            time_of_min = float(filtered_time[min_idx])
            time_of_max = float(filtered_time[max_idx])

            result.channels[name] = ChannelStatistics(
                channel_name=name,
                channel_unit=unit,
                minimum=min_val,
                maximum=max_val,
                mean=mean_val,
                std=std_val,
                abs_max=abs_max_val,
                time_of_min=time_of_min,
                time_of_max=time_of_max,
            )

        return result

    def calculate_from_output(
        self,
        data: np.ndarray,
        channel_names: list[str],
        channel_units: Optional[list[str]] = None,
        case_id: str = "",
    ) -> CaseStatistics:
        """Calculate statistics from a full OpenFAST output array.

        This convenience method handles the common case where the first
        column of the data array is the time channel.

        Parameters
        ----------
        data : np.ndarray
            2D array where column 0 is time and remaining columns
            are output channels.
        channel_names : list[str]
            Channel names including "Time" as the first entry.
        channel_units : list[str], optional
            Channel units including the time unit.
        case_id : str
            Identifier for this case.

        Returns
        -------
        CaseStatistics
            Statistics for every non-time channel.
        """
        time = data[:, 0]
        ch_data = data[:, 1:]
        ch_names = channel_names[1:]  # Skip "Time"
        ch_units = channel_units[1:] if channel_units else None

        return self.calculate(ch_data, time, ch_names, ch_units, case_id)

    @staticmethod
    def aggregate_across_cases(
        case_stats_list: list[CaseStatistics],
    ) -> dict[str, dict[str, float]]:
        """Aggregate statistics across multiple simulation cases.

        For each channel, computes the overall min, max, mean, and
        max standard deviation across all cases.

        Parameters
        ----------
        case_stats_list : list[CaseStatistics]
            Statistics from multiple cases.

        Returns
        -------
        dict[str, dict[str, float]]
            Aggregated statistics per channel:
            {channel_name: {
                "min": overall minimum,
                "max": overall maximum,
                "mean": average of means,
                "std_max": maximum std deviation,
                "abs_max": overall absolute maximum,
                "n_cases": number of cases with this channel,
            }}
        """
        if not case_stats_list:
            return {}

        # Collect all channel names
        all_channels: set[str] = set()
        for cs in case_stats_list:
            all_channels.update(cs.channels.keys())

        result: dict[str, dict[str, float]] = {}

        for ch_name in sorted(all_channels):
            mins: list[float] = []
            maxs: list[float] = []
            means: list[float] = []
            stds: list[float] = []
            abs_maxs: list[float] = []

            for cs in case_stats_list:
                if ch_name in cs.channels:
                    ch = cs.channels[ch_name]
                    mins.append(ch.minimum)
                    maxs.append(ch.maximum)
                    means.append(ch.mean)
                    stds.append(ch.std)
                    abs_maxs.append(ch.abs_max)

            if not mins:
                continue

            result[ch_name] = {
                "min": min(mins),
                "max": max(maxs),
                "mean": float(np.mean(means)),
                "std_max": max(stds),
                "abs_max": max(abs_maxs),
                "n_cases": len(mins),
            }

        return result
