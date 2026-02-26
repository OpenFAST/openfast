"""
2D load envelope calculator for correlated load channel analysis.

Computes the convex hull of two correlated load channels across
multiple simulation cases, providing the load envelope used for
structural design. This is essential for cross-sectional design
where the combined loading from multiple directions must be
considered simultaneously.

Typical use cases:
  - Blade root Mx vs My envelope
  - Tower base Fx vs My envelope
  - Yaw bearing Mx vs My envelope

References
----------
IEC 61400-1:2019, Section 7.6 (Combined loading)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class LoadEnvelopeResult:
    """Result of a 2D load envelope calculation.

    Attributes
    ----------
    channel_x : str
        Name of the x-axis load channel.
    channel_y : str
        Name of the y-axis load channel.
    unit_x : str
        Unit of the x-axis channel.
    unit_y : str
        Unit of the y-axis channel.
    points : list[list[float]]
        All data points [[x, y], ...] from all cases.
    hull_points : list[list[float]]
        Convex hull vertices [[x, y], ...] in counter-clockwise order.
    per_dlc : dict[str, list[list[float]]]
        Data points grouped by DLC number.
    per_case : dict[str, list[list[float]]]
        Data points grouped by case ID.
    hull_area : float
        Area enclosed by the convex hull.
    """
    channel_x: str
    channel_y: str
    unit_x: str = ""
    unit_y: str = ""
    points: list[list[float]] = field(default_factory=list)
    hull_points: list[list[float]] = field(default_factory=list)
    per_dlc: dict[str, list[list[float]]] = field(default_factory=dict)
    per_case: dict[str, list[list[float]]] = field(default_factory=dict)
    hull_area: float = 0.0


class LoadEnvelopeCalculator:
    """Calculates 2D load envelopes from correlated load channels.

    The envelope is defined by the convex hull of all simultaneous
    (channel_x, channel_y) data points across all specified simulation
    cases. The convex hull is computed using the Graham scan algorithm,
    implemented without external computational geometry dependencies.

    Usage
    -----
    >>> calc = LoadEnvelopeCalculator(t_start=30.0)
    >>> result = calc.calculate_envelope(
    ...     channel_x_data={"case1": x_data_1, "case2": x_data_2},
    ...     channel_y_data={"case1": y_data_1, "case2": y_data_2},
    ...     case_dlc_map={"case1": "1.1", "case2": "1.3"},
    ...     time_data={"case1": time_1, "case2": time_2},
    ...     channel_x_name="RootMxc1",
    ...     channel_y_name="RootMyc1",
    ... )
    >>> print(f"Hull has {len(result.hull_points)} vertices")
    """

    def __init__(
        self,
        t_start: float = 0.0,
        downsample_factor: int = 1,
    ) -> None:
        """Initialize the load envelope calculator.

        Parameters
        ----------
        t_start : float
            Time (s) after which to begin analysis (skip transient).
        downsample_factor : int
            Factor by which to downsample data points for efficiency.
            1 = no downsampling (use all points).
        """
        self.t_start = t_start
        self.downsample_factor = max(1, downsample_factor)

    def calculate_envelope(
        self,
        channel_x_data: dict[str, np.ndarray],
        channel_y_data: dict[str, np.ndarray],
        case_dlc_map: Optional[dict[str, str]] = None,
        time_data: Optional[dict[str, np.ndarray]] = None,
        channel_x_name: str = "Channel_X",
        channel_y_name: str = "Channel_Y",
        unit_x: str = "",
        unit_y: str = "",
    ) -> LoadEnvelopeResult:
        """Calculate the 2D load envelope from correlated channels.

        Parameters
        ----------
        channel_x_data : dict[str, np.ndarray]
            X-channel data per case: {case_id: 1D array}.
        channel_y_data : dict[str, np.ndarray]
            Y-channel data per case: {case_id: 1D array}.
        case_dlc_map : dict[str, str], optional
            Mapping from case_id to DLC number for grouping.
        time_data : dict[str, np.ndarray], optional
            Time arrays per case for t_start filtering.
        channel_x_name : str
            Name of the x-axis channel.
        channel_y_name : str
            Name of the y-axis channel.
        unit_x : str
            Unit of the x channel.
        unit_y : str
            Unit of the y channel.

        Returns
        -------
        LoadEnvelopeResult
            Complete envelope results including hull vertices.
        """
        all_points: list[list[float]] = []
        per_dlc: dict[str, list[list[float]]] = {}
        per_case: dict[str, list[list[float]]] = {}

        for case_id in channel_x_data:
            if case_id not in channel_y_data:
                continue

            x_data = channel_x_data[case_id]
            y_data = channel_y_data[case_id]

            # Filter by t_start
            if time_data and case_id in time_data:
                time = time_data[case_id]
                mask = time >= self.t_start
                if not np.any(mask):
                    mask = np.ones(len(time), dtype=bool)
                x_data = x_data[mask]
                y_data = y_data[mask]

            # Ensure same length
            n = min(len(x_data), len(y_data))
            x_data = x_data[:n]
            y_data = y_data[:n]

            # Downsample
            if self.downsample_factor > 1:
                x_data = x_data[::self.downsample_factor]
                y_data = y_data[::self.downsample_factor]

            # Collect points
            case_points: list[list[float]] = []
            for x, y in zip(x_data, y_data):
                pt = [float(x), float(y)]
                all_points.append(pt)
                case_points.append(pt)

            per_case[case_id] = case_points

            # Group by DLC
            if case_dlc_map and case_id in case_dlc_map:
                dlc = case_dlc_map[case_id]
                if dlc not in per_dlc:
                    per_dlc[dlc] = []
                per_dlc[dlc].extend(case_points)

        # Compute convex hull
        if len(all_points) < 3:
            hull_points = list(all_points)
            hull_area = 0.0
        else:
            hull_points = self._convex_hull(all_points)
            hull_area = self._polygon_area(hull_points)

        return LoadEnvelopeResult(
            channel_x=channel_x_name,
            channel_y=channel_y_name,
            unit_x=unit_x,
            unit_y=unit_y,
            points=all_points,
            hull_points=hull_points,
            per_dlc=per_dlc,
            per_case=per_case,
            hull_area=hull_area,
        )

    @staticmethod
    def _convex_hull(points: list[list[float]]) -> list[list[float]]:
        """Compute the 2D convex hull using Graham scan.

        Parameters
        ----------
        points : list[list[float]]
            List of [x, y] points.

        Returns
        -------
        list[list[float]]
            Convex hull vertices in counter-clockwise order.
        """
        if len(points) < 3:
            return list(points)

        # Find the bottom-most point (lowest y, then leftmost x)
        pts = sorted(points, key=lambda p: (p[1], p[0]))
        pivot = pts[0]

        def polar_angle(p: list[float]) -> float:
            """Polar angle from pivot to point p."""
            dx = p[0] - pivot[0]
            dy = p[1] - pivot[1]
            import math
            return math.atan2(dy, dx)

        def distance_sq(p: list[float]) -> float:
            """Squared distance from pivot."""
            return (p[0] - pivot[0]) ** 2 + (p[1] - pivot[1]) ** 2

        # Sort points by polar angle from pivot
        remaining = pts[1:]
        remaining.sort(key=lambda p: (polar_angle(p), distance_sq(p)))

        # Remove duplicates at same angle (keep farthest)
        if not remaining:
            return [pivot]

        # Graham scan
        stack: list[list[float]] = [pivot]

        for pt in remaining:
            while len(stack) > 1:
                # Cross product of (stack[-2] -> stack[-1]) x (stack[-2] -> pt)
                o = stack[-2]
                a = stack[-1]
                cross = (a[0] - o[0]) * (pt[1] - o[1]) - \
                        (a[1] - o[1]) * (pt[0] - o[0])
                if cross <= 0:
                    stack.pop()
                else:
                    break
            stack.append(pt)

        return stack

    @staticmethod
    def _polygon_area(vertices: list[list[float]]) -> float:
        """Calculate the area of a polygon using the shoelace formula.

        Parameters
        ----------
        vertices : list[list[float]]
            Polygon vertices in order.

        Returns
        -------
        float
            Polygon area.
        """
        n = len(vertices)
        if n < 3:
            return 0.0

        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += vertices[i][0] * vertices[j][1]
            area -= vertices[j][0] * vertices[i][1]

        return abs(area) / 2.0
