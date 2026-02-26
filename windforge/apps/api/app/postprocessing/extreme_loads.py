"""
Extreme load extractor for IEC ultimate load analysis.

Extracts characteristic extreme loads from simulation results across
all ultimate DLC cases, applies partial safety factors from the
governing DLC, and produces design extreme loads per IEC 61400-1.

The extreme load extraction process:
  1. Find the global maximum and minimum for each channel across all
     ultimate DLC cases.
  2. Record the source DLC, time of occurrence, and wind conditions.
  3. Apply the appropriate partial safety factor (gamma_f * gamma_n)
     to produce design loads.
  4. Identify the governing DLC for each channel.

References
----------
IEC 61400-1:2019, Section 7.6 (Ultimate strength analysis)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class ExtremeValue:
    """A single extreme load value with metadata.

    Attributes
    ----------
    characteristic : float
        Characteristic (unfactored) extreme value.
    design : float
        Design (factored) extreme value = characteristic * safety_factor.
    safety_factor : float
        Applied partial safety factor (gamma_f * gamma_n).
    source_dlc : str
        DLC number that produced this extreme.
    source_case_id : str
        Specific case ID that produced this extreme.
    time : float
        Time of occurrence within the simulation (s).
    wind_speed : float
        Hub-height wind speed for the governing case (m/s).
    seed : int
        Random seed of the governing case.
    """
    characteristic: float
    design: float
    safety_factor: float
    source_dlc: str
    source_case_id: str
    time: float
    wind_speed: float = 0.0
    seed: int = 0


@dataclass
class ChannelExtremes:
    """Extreme values (max and min) for a single channel.

    Attributes
    ----------
    channel_name : str
        Name of the load channel.
    channel_unit : str
        Unit of the channel.
    max_extreme : ExtremeValue
        Maximum extreme value with metadata.
    min_extreme : ExtremeValue
        Minimum extreme value with metadata.
    """
    channel_name: str
    channel_unit: str
    max_extreme: ExtremeValue
    min_extreme: ExtremeValue


@dataclass
class SimulationResult:
    """Container for a single simulation case result.

    Attributes
    ----------
    case_id : str
        Unique case identifier.
    dlc_number : str
        DLC number (e.g., "1.1").
    wind_speed : float
        Hub-height mean wind speed (m/s).
    seed : int
        Random seed number.
    safety_factor : float
        Partial safety factor for this DLC.
    time : np.ndarray
        Time channel array (s).
    data : dict[str, np.ndarray]
        Channel data keyed by channel name.
    channel_units : dict[str, str]
        Channel units keyed by channel name.
    """
    case_id: str
    dlc_number: str
    wind_speed: float
    seed: int
    safety_factor: float
    time: np.ndarray
    data: dict[str, np.ndarray]
    channel_units: dict[str, str] = field(default_factory=dict)


class ExtremeLoadExtractor:
    """Extracts characteristic and design extreme loads from simulations.

    Processes results from multiple ultimate DLC cases to find the
    governing extremes for each load channel and applies partial
    safety factors.

    Usage
    -----
    >>> extractor = ExtremeLoadExtractor(t_start=30.0)
    >>> results = [SimulationResult(...), SimulationResult(...)]
    >>> channels = ["RootMxc1", "TwrBsMxt", "YawBrMxp"]
    >>> extremes = extractor.extract_extremes(results, channels)
    >>> print(extremes["RootMxc1"].max_extreme.design)
    """

    def __init__(
        self,
        t_start: float = 0.0,
        consequence_factor: float = 1.0,
    ) -> None:
        """Initialize the extreme load extractor.

        Parameters
        ----------
        t_start : float
            Time (s) after which to begin searching for extremes
            (skip initial transient).
        consequence_factor : float
            Consequence of failure factor gamma_n (default 1.0).
        """
        self.t_start = t_start
        self.consequence_factor = consequence_factor

    def extract_extremes(
        self,
        simulation_results: list[SimulationResult],
        channels: Optional[list[str]] = None,
    ) -> dict[str, ChannelExtremes]:
        """Find characteristic and design extreme loads across all cases.

        For each specified channel (or all channels if none specified),
        finds the global maximum and minimum values across all
        simulation results, records the governing DLC and time, and
        applies the appropriate partial safety factor.

        Parameters
        ----------
        simulation_results : list[SimulationResult]
            Results from all ultimate DLC simulation cases.
        channels : list[str], optional
            Specific channels to extract extremes for. If None, all
            channels present in any result are processed.

        Returns
        -------
        dict[str, ChannelExtremes]
            Extreme values per channel, keyed by channel name.
        """
        if not simulation_results:
            return {}

        # Determine channels to process
        if channels is None:
            channels_set: set[str] = set()
            for result in simulation_results:
                channels_set.update(result.data.keys())
            channels = sorted(channels_set)

        extremes: dict[str, ChannelExtremes] = {}

        for ch_name in channels:
            max_val = -np.inf
            min_val = np.inf
            max_info: Optional[dict] = None
            min_info: Optional[dict] = None

            for result in simulation_results:
                if ch_name not in result.data:
                    continue

                signal = result.data[ch_name]
                time = result.time

                # Apply t_start filter
                mask = time >= self.t_start
                if not np.any(mask):
                    mask = np.ones(len(time), dtype=bool)

                filtered_signal = signal[mask]
                filtered_time = time[mask]

                if len(filtered_signal) == 0:
                    continue

                local_max_idx = int(np.nanargmax(filtered_signal))
                local_min_idx = int(np.nanargmin(filtered_signal))
                local_max = float(filtered_signal[local_max_idx])
                local_min = float(filtered_signal[local_min_idx])

                if local_max > max_val:
                    max_val = local_max
                    max_info = {
                        "characteristic": local_max,
                        "safety_factor": result.safety_factor * self.consequence_factor,
                        "source_dlc": result.dlc_number,
                        "source_case_id": result.case_id,
                        "time": float(filtered_time[local_max_idx]),
                        "wind_speed": result.wind_speed,
                        "seed": result.seed,
                    }

                if local_min < min_val:
                    min_val = local_min
                    min_info = {
                        "characteristic": local_min,
                        "safety_factor": result.safety_factor * self.consequence_factor,
                        "source_dlc": result.dlc_number,
                        "source_case_id": result.case_id,
                        "time": float(filtered_time[local_min_idx]),
                        "wind_speed": result.wind_speed,
                        "seed": result.seed,
                    }

            if max_info is None or min_info is None:
                continue

            # Apply safety factors to get design values
            max_design = max_info["characteristic"] * max_info["safety_factor"]
            min_design = min_info["characteristic"] * min_info["safety_factor"]

            ch_unit = ""
            for result in simulation_results:
                if ch_name in result.channel_units:
                    ch_unit = result.channel_units[ch_name]
                    break

            extremes[ch_name] = ChannelExtremes(
                channel_name=ch_name,
                channel_unit=ch_unit,
                max_extreme=ExtremeValue(
                    characteristic=max_info["characteristic"],
                    design=max_design,
                    safety_factor=max_info["safety_factor"],
                    source_dlc=max_info["source_dlc"],
                    source_case_id=max_info["source_case_id"],
                    time=max_info["time"],
                    wind_speed=max_info["wind_speed"],
                    seed=max_info["seed"],
                ),
                min_extreme=ExtremeValue(
                    characteristic=min_info["characteristic"],
                    design=min_design,
                    safety_factor=min_info["safety_factor"],
                    source_dlc=min_info["source_dlc"],
                    source_case_id=min_info["source_case_id"],
                    time=min_info["time"],
                    wind_speed=min_info["wind_speed"],
                    seed=min_info["seed"],
                ),
            )

        return extremes

    @staticmethod
    def to_summary_table(
        extremes: dict[str, ChannelExtremes],
    ) -> list[dict[str, object]]:
        """Convert extreme results to a flat summary table format.

        Parameters
        ----------
        extremes : dict[str, ChannelExtremes]
            Extreme values per channel.

        Returns
        -------
        list[dict[str, object]]
            List of dictionaries suitable for tabular display
            or DataFrame construction.
        """
        rows: list[dict[str, object]] = []

        for ch_name, ch_ext in sorted(extremes.items()):
            rows.append({
                "Channel": ch_name,
                "Unit": ch_ext.channel_unit,
                "Max_Characteristic": ch_ext.max_extreme.characteristic,
                "Max_Design": ch_ext.max_extreme.design,
                "Max_DLC": ch_ext.max_extreme.source_dlc,
                "Max_Case": ch_ext.max_extreme.source_case_id,
                "Max_Time_s": ch_ext.max_extreme.time,
                "Max_Vhub_mps": ch_ext.max_extreme.wind_speed,
                "Min_Characteristic": ch_ext.min_extreme.characteristic,
                "Min_Design": ch_ext.min_extreme.design,
                "Min_DLC": ch_ext.min_extreme.source_dlc,
                "Min_Case": ch_ext.min_extreme.source_case_id,
                "Min_Time_s": ch_ext.min_extreme.time,
                "Min_Vhub_mps": ch_ext.min_extreme.wind_speed,
                "Safety_Factor_Max": ch_ext.max_extreme.safety_factor,
                "Safety_Factor_Min": ch_ext.min_extreme.safety_factor,
            })

        return rows
