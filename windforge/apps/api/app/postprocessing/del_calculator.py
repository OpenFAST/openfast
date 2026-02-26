"""
Damage Equivalent Load (DEL) calculator with built-in rainflow counting.

Implements the 4-point rainflow counting algorithm from scratch (no
external dependencies beyond numpy) and calculates DELs per IEC
fatigue methodology.

The 4-point rainflow algorithm follows ASTM E1049-85 and is the
standard approach for extracting fatigue cycles from irregular
load time histories.

DEL calculation:
  DEL = ( sum(S_i^m * n_i) / N_eq )^(1/m)

where:
  S_i = stress/load range for cycle i
  n_i = count for cycle i (0.5 for half-cycles, 1.0 for full)
  m   = Woehler exponent (material-dependent)
  N_eq = number of equivalent cycles at the DEL level

References
----------
ASTM E1049-85: Standard Practices for Cycle Counting in Fatigue Analysis
IEC 61400-1:2019, Section 7.8 (Fatigue analysis)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class RainflowCycle:
    """A single rainflow cycle extracted from a load time history.

    Attributes
    ----------
    range : float
        Load range (peak-to-valley or valley-to-peak).
    mean : float
        Mean load level of the cycle.
    count : float
        Cycle count (1.0 for full cycle, 0.5 for half cycle).
    """
    range: float
    mean: float
    count: float


class DELCalculator:
    """Damage Equivalent Load calculator with integrated rainflow counting.

    Implements the 4-point rainflow counting method and DEL computation
    per IEC 61400-1 fatigue methodology.

    Usage
    -----
    >>> calc = DELCalculator()
    >>> cycles = calc.rainflow_count(signal)
    >>> del_val = calc.calculate_del(signal, dt=0.05, m_exponent=10.0)
    """

    def rainflow_count(
        self, signal: np.ndarray
    ) -> list[tuple[float, float]]:
        """Perform 4-point rainflow cycle counting on a load signal.

        Implements the ASTM E1049-85 4-point rainflow counting
        algorithm. The signal is first reduced to its peaks and
        valleys (turning points), then cycles are extracted using
        the 4-point method.

        Parameters
        ----------
        signal : np.ndarray
            1D array of load values (time series).

        Returns
        -------
        list[tuple[float, float]]
            List of (range, count) tuples. Range is the load
            amplitude of each cycle, count is 0.5 for half-cycles
            and 1.0 for full cycles.
        """
        if len(signal) < 3:
            return []

        # Step 1: Extract turning points (peaks and valleys)
        turning_points = self._extract_turning_points(signal)

        if len(turning_points) < 3:
            return []

        # Step 2: Apply 4-point rainflow counting
        cycles = self._four_point_rainflow(turning_points)

        # Step 3: Convert to (range, count) tuples
        return [(c.range, c.count) for c in cycles]

    def rainflow_count_detailed(
        self, signal: np.ndarray
    ) -> list[RainflowCycle]:
        """Perform rainflow counting with full cycle information.

        Like rainflow_count but returns RainflowCycle objects with
        range, mean, and count information.

        Parameters
        ----------
        signal : np.ndarray
            1D array of load values.

        Returns
        -------
        list[RainflowCycle]
            Detailed cycle information.
        """
        if len(signal) < 3:
            return []

        turning_points = self._extract_turning_points(signal)
        if len(turning_points) < 3:
            return []

        return self._four_point_rainflow(turning_points)

    def calculate_del(
        self,
        signal: np.ndarray,
        dt: float,
        m_exponent: float = 10.0,
        n_equivalent: float = 1.0e7,
        t_start: float = 0.0,
    ) -> float:
        """Calculate the Damage Equivalent Load for a single channel.

        DEL = ( sum(S_i^m * n_i) / N_eq )^(1/m)

        Parameters
        ----------
        signal : np.ndarray
            1D array of load time series values.
        dt : float
            Time step between samples (s).
        m_exponent : float
            Woehler curve exponent (S-N slope). Typical values:
            - m=3-5 for steel/welded joints
            - m=8-12 for fiberglass composites
            - m=10 default for blade root moments
        n_equivalent : float
            Number of equivalent cycles for the DEL reference.
            Default 1e7 corresponds to ~20 year lifetime at 1P.
        t_start : float
            Time (s) at which to begin analysis (skip transient).

        Returns
        -------
        float
            Damage Equivalent Load value. Returns 0.0 if no cycles
            are found.
        """
        # Skip transient period
        if t_start > 0.0 and dt > 0.0:
            n_skip = int(t_start / dt)
            signal = signal[n_skip:]

        if len(signal) < 3:
            return 0.0

        cycles = self.rainflow_count(signal)

        if not cycles:
            return 0.0

        # Sum S_i^m * n_i
        damage_sum = 0.0
        for load_range, count in cycles:
            damage_sum += (abs(load_range) ** m_exponent) * count

        if damage_sum <= 0.0:
            return 0.0

        # DEL = (damage_sum / N_eq)^(1/m)
        del_value = (damage_sum / n_equivalent) ** (1.0 / m_exponent)

        return float(del_value)

    @staticmethod
    def combine_del_across_cases(
        del_per_case: list[float],
        case_weights: list[float],
        m_exponent: float = 10.0,
    ) -> float:
        """Combine DELs from multiple cases using weighted summation.

        The combined DEL accounts for the probability of each wind
        speed bin (from the Weibull distribution) and sums the
        damage contributions:

          DEL_combined = ( sum(w_i * DEL_i^m) )^(1/m)

        Parameters
        ----------
        del_per_case : list[float]
            DEL value from each simulation case.
        case_weights : list[float]
            Probability weight for each case (must sum to ~1.0
            for proper normalization).
        m_exponent : float
            Woehler exponent (same as used for individual DELs).

        Returns
        -------
        float
            Combined damage equivalent load.
        """
        if not del_per_case or not case_weights:
            return 0.0

        if len(del_per_case) != len(case_weights):
            raise ValueError(
                f"Number of DEL values ({len(del_per_case)}) must match "
                f"number of weights ({len(case_weights)})"
            )

        weighted_damage = 0.0
        for del_val, weight in zip(del_per_case, case_weights):
            weighted_damage += weight * (abs(del_val) ** m_exponent)

        if weighted_damage <= 0.0:
            return 0.0

        return float(weighted_damage ** (1.0 / m_exponent))

    @staticmethod
    def _extract_turning_points(signal: np.ndarray) -> list[float]:
        """Extract peaks and valleys (turning points) from a signal.

        Removes intermediate points that lie on monotonic segments,
        retaining only the local extrema and the first/last points.

        Parameters
        ----------
        signal : np.ndarray
            1D array of load values.

        Returns
        -------
        list[float]
            Turning point values in order.
        """
        if len(signal) <= 2:
            return list(signal)

        turning_points: list[float] = [float(signal[0])]

        for i in range(1, len(signal) - 1):
            prev_val = signal[i - 1]
            curr_val = signal[i]
            next_val = signal[i + 1]

            # Check if current point is a local extremum
            if (curr_val >= prev_val and curr_val >= next_val) or \
               (curr_val <= prev_val and curr_val <= next_val):
                # Only add if different from previous turning point
                if curr_val != turning_points[-1]:
                    turning_points.append(float(curr_val))

        # Always include the last point
        last_val = float(signal[-1])
        if last_val != turning_points[-1]:
            turning_points.append(last_val)

        return turning_points

    @staticmethod
    def _four_point_rainflow(
        turning_points: list[float],
    ) -> list[RainflowCycle]:
        """Apply the 4-point rainflow counting algorithm.

        The 4-point method examines four consecutive turning points
        at a time. If the inner range is less than or equal to the
        outer range, a full cycle is extracted from the inner pair.

        After processing, any remaining turning points are counted
        as half-cycles using the 3-point method.

        Parameters
        ----------
        turning_points : list[float]
            Ordered list of turning point values.

        Returns
        -------
        list[RainflowCycle]
            Extracted cycles.
        """
        cycles: list[RainflowCycle] = []
        points = list(turning_points)  # working copy

        # Phase 1: 4-point method (extract full cycles)
        i = 0
        while len(points) >= 4:
            found = False
            i = 0
            while i <= len(points) - 4:
                s1 = points[i]
                s2 = points[i + 1]
                s3 = points[i + 2]
                s4 = points[i + 3]

                range_outer = abs(s1 - s4)
                range_inner = abs(s2 - s3)
                range_12 = abs(s1 - s2)
                range_34 = abs(s3 - s4)

                # Check: inner range <= both adjacent ranges
                if range_inner <= range_12 and range_inner <= range_34:
                    # Extract full cycle from inner pair
                    cycle_range = range_inner
                    cycle_mean = (s2 + s3) / 2.0
                    cycles.append(RainflowCycle(
                        range=cycle_range,
                        mean=cycle_mean,
                        count=1.0,
                    ))
                    # Remove inner two points
                    del points[i + 1: i + 3]
                    found = True
                    break
                else:
                    i += 1

            if not found:
                break

        # Phase 2: Count remaining turning points as half-cycles
        # Using the simple range-pair method for residuals
        for j in range(len(points) - 1):
            cycle_range = abs(points[j] - points[j + 1])
            cycle_mean = (points[j] + points[j + 1]) / 2.0
            if cycle_range > 0:
                cycles.append(RainflowCycle(
                    range=cycle_range,
                    mean=cycle_mean,
                    count=0.5,
                ))

        return cycles

    @staticmethod
    def compute_damage(
        cycles: list[tuple[float, float]],
        m_exponent: float,
        n_reference: float = 1.0e7,
        s_reference: float = 1.0,
    ) -> float:
        """Compute Miner's rule cumulative damage from cycle counts.

        D = sum( n_i / N_i )

        where N_i is the number of cycles to failure at range S_i,
        determined from the S-N curve: N = N_ref * (S_ref / S_i)^m

        Parameters
        ----------
        cycles : list[tuple[float, float]]
            List of (range, count) tuples from rainflow counting.
        m_exponent : float
            Woehler exponent.
        n_reference : float
            Reference number of cycles on the S-N curve.
        s_reference : float
            Reference stress/load range on the S-N curve.

        Returns
        -------
        float
            Cumulative damage ratio D. Failure when D >= 1.0.
        """
        damage = 0.0
        for load_range, count in cycles:
            if load_range <= 0:
                continue
            # N_i = N_ref * (S_ref / S_i)^m
            n_to_fail = n_reference * (s_reference / load_range) ** m_exponent
            if n_to_fail > 0:
                damage += count / n_to_fail
        return damage
