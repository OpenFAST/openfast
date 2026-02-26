"""
HTML report generator for IEC load analysis results.

Generates comprehensive HTML reports with embedded Plotly charts and
styled tables. Reports cover the complete IEC load analysis workflow
including turbine description, DLC matrix, statistics, DEL results,
extreme loads, and load envelopes.

The generated HTML is self-contained (no external dependencies at
render time) with inline CSS and Plotly.js loaded from CDN.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

from .statistics import CaseStatistics
from .extreme_loads import ChannelExtremes
from .load_envelope import LoadEnvelopeResult


@dataclass
class TurbineDescription:
    """Turbine specification data for report display."""
    name: str = ""
    rated_power_kw: float = 0.0
    rotor_diameter: float = 0.0
    hub_height: float = 0.0
    num_blades: int = 3
    cut_in_speed: float = 0.0
    cut_out_speed: float = 0.0
    rated_speed: float = 0.0
    iec_class: str = ""


@dataclass
class DLCMatrixEntry:
    """Summary of a DLC in the analysis."""
    dlc_number: str
    description: str
    wind_model: str
    analysis_type: str
    num_cases: int
    wind_speed_range: str
    safety_factor: float


@dataclass
class DELResult:
    """DEL result for a single channel."""
    channel_name: str
    channel_unit: str
    del_value: float
    m_exponent: float
    n_equivalent: float


@dataclass
class ReportData:
    """Complete data package for report generation.

    Attributes
    ----------
    project_name : str
        Name of the project.
    turbine : TurbineDescription
        Turbine specification.
    dlc_matrix : list[DLCMatrixEntry]
        DLC matrix summary.
    statistics_summary : dict[str, dict[str, float]]
        Aggregated statistics per channel.
    del_results : list[DELResult]
        DEL values per channel.
    extreme_results : dict[str, ChannelExtremes]
        Extreme load results.
    load_envelopes : list[LoadEnvelopeResult]
        Load envelope results.
    total_cases : int
        Total number of simulation cases.
    total_simulation_hours : float
        Total simulated time in hours.
    """
    project_name: str = ""
    turbine: TurbineDescription = field(default_factory=TurbineDescription)
    dlc_matrix: list[DLCMatrixEntry] = field(default_factory=list)
    statistics_summary: dict[str, dict[str, float]] = field(default_factory=dict)
    del_results: list[DELResult] = field(default_factory=list)
    extreme_results: dict[str, ChannelExtremes] = field(default_factory=dict)
    load_envelopes: list[LoadEnvelopeResult] = field(default_factory=list)
    total_cases: int = 0
    total_simulation_hours: float = 0.0


class ReportGenerator:
    """Generates comprehensive HTML reports for IEC load analysis.

    The report includes:
    - Executive Summary
    - Turbine Description
    - DLC Matrix overview
    - Statistics Summary with bar charts
    - DEL Results with material exponent comparison
    - Extreme Loads table with governing DLCs
    - Load Envelope plots
    - Conclusions

    Usage
    -----
    >>> gen = ReportGenerator()
    >>> report_data = ReportData(
    ...     project_name="MyProject",
    ...     turbine=TurbineDescription(name="NREL 5MW"),
    ...     ...
    ... )
    >>> html = gen.generate_html_report(report_data)
    >>> with open("report.html", "w") as f:
    ...     f.write(html)
    """

    def generate_html_report(self, data: ReportData) -> str:
        """Generate a complete HTML report.

        Parameters
        ----------
        data : ReportData
            All data needed for the report.

        Returns
        -------
        str
            Complete HTML document string.
        """
        timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

        sections = [
            self._section_executive_summary(data, timestamp),
            self._section_turbine_description(data.turbine),
            self._section_dlc_matrix(data.dlc_matrix),
            self._section_statistics(data.statistics_summary),
            self._section_del_results(data.del_results),
            self._section_extreme_loads(data.extreme_results),
            self._section_load_envelopes(data.load_envelopes),
            self._section_conclusions(data),
        ]

        body_content = "\n".join(sections)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{_esc(data.project_name)} - IEC Load Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    {self._css()}
</head>
<body>
    <div class="container">
        <header>
            <h1>{_esc(data.project_name)}</h1>
            <h2>IEC 61400-1 Ed.4 Load Analysis Report</h2>
            <p class="timestamp">Generated: {timestamp} by WindForge</p>
        </header>
        <nav>
            <a href="#executive-summary">Executive Summary</a>
            <a href="#turbine-description">Turbine</a>
            <a href="#dlc-matrix">DLC Matrix</a>
            <a href="#statistics">Statistics</a>
            <a href="#del-results">DEL Results</a>
            <a href="#extreme-loads">Extreme Loads</a>
            <a href="#load-envelopes">Load Envelopes</a>
            <a href="#conclusions">Conclusions</a>
        </nav>
        {body_content}
        <footer>
            <p>Report generated by WindForge | IEC 61400-1:2019 (Ed.4) compliant analysis</p>
        </footer>
    </div>
</body>
</html>"""
        return html

    def _css(self) -> str:
        """Return the inline CSS stylesheet."""
        return """<style>
    * { margin: 0; padding: 0; box-sizing: border-box; }
    body {
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        line-height: 1.6;
        color: #333;
        background: #f5f5f5;
    }
    .container { max-width: 1200px; margin: 0 auto; padding: 20px; }
    header {
        background: linear-gradient(135deg, #1a365d, #2b6cb0);
        color: white;
        padding: 40px;
        border-radius: 8px 8px 0 0;
        margin-bottom: 0;
    }
    header h1 { font-size: 2em; margin-bottom: 8px; }
    header h2 { font-size: 1.2em; font-weight: 300; opacity: 0.9; }
    .timestamp { margin-top: 12px; font-size: 0.9em; opacity: 0.7; }
    nav {
        background: #2b6cb0;
        padding: 12px 40px;
        display: flex;
        gap: 20px;
        flex-wrap: wrap;
        border-radius: 0 0 8px 8px;
        margin-bottom: 30px;
    }
    nav a {
        color: white;
        text-decoration: none;
        font-size: 0.9em;
        opacity: 0.8;
        transition: opacity 0.2s;
    }
    nav a:hover { opacity: 1; }
    section {
        background: white;
        border-radius: 8px;
        padding: 30px;
        margin-bottom: 20px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    section h2 {
        font-size: 1.5em;
        color: #1a365d;
        border-bottom: 2px solid #e2e8f0;
        padding-bottom: 10px;
        margin-bottom: 20px;
    }
    section h3 {
        font-size: 1.1em;
        color: #2d3748;
        margin: 20px 0 10px;
    }
    table {
        width: 100%;
        border-collapse: collapse;
        margin: 15px 0;
        font-size: 0.9em;
    }
    th {
        background: #edf2f7;
        color: #1a365d;
        font-weight: 600;
        text-align: left;
        padding: 10px 12px;
        border-bottom: 2px solid #cbd5e0;
    }
    td {
        padding: 8px 12px;
        border-bottom: 1px solid #e2e8f0;
    }
    tr:hover td { background: #f7fafc; }
    .stat-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 15px;
        margin: 20px 0;
    }
    .stat-card {
        background: #f7fafc;
        border-radius: 6px;
        padding: 15px;
        border-left: 4px solid #3182ce;
    }
    .stat-card .label { font-size: 0.85em; color: #718096; }
    .stat-card .value { font-size: 1.4em; font-weight: 600; color: #1a365d; }
    .plot-container { margin: 20px 0; }
    .highlight-max { background: #fff5f5; color: #c53030; font-weight: 600; }
    .highlight-min { background: #f0fff4; color: #276749; font-weight: 600; }
    footer {
        text-align: center;
        padding: 20px;
        color: #718096;
        font-size: 0.85em;
    }
    @media print {
        nav { display: none; }
        section { break-inside: avoid; }
    }
</style>"""

    def _section_executive_summary(
        self, data: ReportData, timestamp: str
    ) -> str:
        """Generate the Executive Summary section."""
        return f"""
<section id="executive-summary">
    <h2>Executive Summary</h2>
    <div class="stat-grid">
        <div class="stat-card">
            <div class="label">Project</div>
            <div class="value">{_esc(data.project_name)}</div>
        </div>
        <div class="stat-card">
            <div class="label">Turbine</div>
            <div class="value">{_esc(data.turbine.name)}</div>
        </div>
        <div class="stat-card">
            <div class="label">IEC Class</div>
            <div class="value">{_esc(data.turbine.iec_class)}</div>
        </div>
        <div class="stat-card">
            <div class="label">Total Cases</div>
            <div class="value">{data.total_cases}</div>
        </div>
        <div class="stat-card">
            <div class="label">Simulation Time</div>
            <div class="value">{data.total_simulation_hours:.1f} hrs</div>
        </div>
        <div class="stat-card">
            <div class="label">DLCs Analyzed</div>
            <div class="value">{len(data.dlc_matrix)}</div>
        </div>
    </div>
    <p>This report presents the results of the IEC 61400-1 Ed.4 design load
    analysis performed using OpenFAST aeroelastic simulations. The analysis
    covers both ultimate and fatigue load assessments across the design
    load case matrix.</p>
</section>"""

    def _section_turbine_description(self, turbine: TurbineDescription) -> str:
        """Generate the Turbine Description section."""
        return f"""
<section id="turbine-description">
    <h2>Turbine Description</h2>
    <table>
        <tr><td><strong>Turbine Model</strong></td><td>{_esc(turbine.name)}</td></tr>
        <tr><td><strong>Rated Power</strong></td><td>{turbine.rated_power_kw:.0f} kW</td></tr>
        <tr><td><strong>Rotor Diameter</strong></td><td>{turbine.rotor_diameter:.1f} m</td></tr>
        <tr><td><strong>Hub Height</strong></td><td>{turbine.hub_height:.1f} m</td></tr>
        <tr><td><strong>Number of Blades</strong></td><td>{turbine.num_blades}</td></tr>
        <tr><td><strong>Cut-in Wind Speed</strong></td><td>{turbine.cut_in_speed:.1f} m/s</td></tr>
        <tr><td><strong>Rated Wind Speed</strong></td><td>{turbine.rated_speed:.1f} m/s</td></tr>
        <tr><td><strong>Cut-out Wind Speed</strong></td><td>{turbine.cut_out_speed:.1f} m/s</td></tr>
        <tr><td><strong>IEC Classification</strong></td><td>{_esc(turbine.iec_class)}</td></tr>
    </table>
</section>"""

    def _section_dlc_matrix(self, dlc_matrix: list[DLCMatrixEntry]) -> str:
        """Generate the DLC Matrix section."""
        if not dlc_matrix:
            return """
<section id="dlc-matrix">
    <h2>Design Load Case Matrix</h2>
    <p>No DLC data available.</p>
</section>"""

        rows = ""
        for dlc in dlc_matrix:
            rows += f"""
        <tr>
            <td><strong>{_esc(dlc.dlc_number)}</strong></td>
            <td>{_esc(dlc.description)}</td>
            <td>{_esc(dlc.wind_model)}</td>
            <td>{_esc(dlc.analysis_type)}</td>
            <td>{dlc.num_cases}</td>
            <td>{_esc(dlc.wind_speed_range)}</td>
            <td>{dlc.safety_factor:.2f}</td>
        </tr>"""

        return f"""
<section id="dlc-matrix">
    <h2>Design Load Case Matrix</h2>
    <p>Summary of analyzed design load cases per IEC 61400-1 Ed.4 Table 4.</p>
    <table>
        <thead>
            <tr>
                <th>DLC</th>
                <th>Description</th>
                <th>Wind Model</th>
                <th>Type</th>
                <th>Cases</th>
                <th>Wind Range</th>
                <th>gamma_f</th>
            </tr>
        </thead>
        <tbody>{rows}
        </tbody>
    </table>
</section>"""

    def _section_statistics(
        self, stats: dict[str, dict[str, float]]
    ) -> str:
        """Generate the Statistics Summary section with bar chart."""
        if not stats:
            return """
<section id="statistics">
    <h2>Statistics Summary</h2>
    <p>No statistics data available.</p>
</section>"""

        # Build table rows
        rows = ""
        channels = []
        abs_maxes = []
        for ch_name, ch_stats in sorted(stats.items()):
            channels.append(ch_name)
            abs_maxes.append(ch_stats.get("abs_max", 0.0))
            rows += f"""
        <tr>
            <td><strong>{_esc(ch_name)}</strong></td>
            <td>{ch_stats.get('min', 0.0):.4E}</td>
            <td>{ch_stats.get('max', 0.0):.4E}</td>
            <td>{ch_stats.get('mean', 0.0):.4E}</td>
            <td>{ch_stats.get('std_max', 0.0):.4E}</td>
            <td>{ch_stats.get('abs_max', 0.0):.4E}</td>
            <td>{int(ch_stats.get('n_cases', 0))}</td>
        </tr>"""

        # Generate Plotly bar chart data
        chart_data = json.dumps({
            "x": channels[:20],  # Limit to 20 channels for readability
            "y": abs_maxes[:20],
            "type": "bar",
            "marker": {"color": "#3182ce"},
        })

        chart_id = "stats-chart"

        return f"""
<section id="statistics">
    <h2>Statistics Summary</h2>
    <p>Aggregated statistics across all simulation cases. Values represent
    the envelope across all seeds and wind speeds.</p>
    <table>
        <thead>
            <tr>
                <th>Channel</th>
                <th>Min</th>
                <th>Max</th>
                <th>Mean</th>
                <th>Max Std</th>
                <th>Abs Max</th>
                <th>Cases</th>
            </tr>
        </thead>
        <tbody>{rows}
        </tbody>
    </table>
    <div class="plot-container">
        <h3>Absolute Maximum Values</h3>
        <div id="{chart_id}" style="height: 400px;"></div>
        <script>
            Plotly.newPlot('{chart_id}', [{chart_data}], {{
                title: 'Absolute Maximum per Channel',
                xaxis: {{tickangle: -45}},
                yaxis: {{title: 'Absolute Maximum'}},
                margin: {{b: 120}},
            }});
        </script>
    </div>
</section>"""

    def _section_del_results(self, del_results: list[DELResult]) -> str:
        """Generate the DEL Results section."""
        if not del_results:
            return """
<section id="del-results">
    <h2>Damage Equivalent Loads</h2>
    <p>No DEL results available.</p>
</section>"""

        rows = ""
        channels = []
        del_values = []
        for dr in del_results:
            channels.append(dr.channel_name)
            del_values.append(dr.del_value)
            rows += f"""
        <tr>
            <td><strong>{_esc(dr.channel_name)}</strong></td>
            <td>{_esc(dr.channel_unit)}</td>
            <td>{dr.del_value:.4E}</td>
            <td>{dr.m_exponent:.1f}</td>
            <td>{dr.n_equivalent:.2E}</td>
        </tr>"""

        chart_data = json.dumps({
            "x": channels[:20],
            "y": del_values[:20],
            "type": "bar",
            "marker": {"color": "#e53e3e"},
        })
        chart_id = "del-chart"

        return f"""
<section id="del-results">
    <h2>Damage Equivalent Loads</h2>
    <p>Lifetime DEL values computed using rainflow cycle counting and
    Woehler S-N curve methodology per IEC 61400-1 Ed.4.</p>
    <table>
        <thead>
            <tr>
                <th>Channel</th>
                <th>Unit</th>
                <th>DEL</th>
                <th>m Exponent</th>
                <th>N_eq</th>
            </tr>
        </thead>
        <tbody>{rows}
        </tbody>
    </table>
    <div class="plot-container">
        <h3>DEL Values by Channel</h3>
        <div id="{chart_id}" style="height: 400px;"></div>
        <script>
            Plotly.newPlot('{chart_id}', [{chart_data}], {{
                title: 'Damage Equivalent Loads',
                xaxis: {{tickangle: -45}},
                yaxis: {{title: 'DEL'}},
                margin: {{b: 120}},
            }});
        </script>
    </div>
</section>"""

    def _section_extreme_loads(
        self, extremes: dict[str, ChannelExtremes]
    ) -> str:
        """Generate the Extreme Loads section."""
        if not extremes:
            return """
<section id="extreme-loads">
    <h2>Extreme Loads</h2>
    <p>No extreme load results available.</p>
</section>"""

        rows = ""
        for ch_name, ch_ext in sorted(extremes.items()):
            me = ch_ext.max_extreme
            ie = ch_ext.min_extreme
            rows += f"""
        <tr>
            <td><strong>{_esc(ch_name)}</strong></td>
            <td>{_esc(ch_ext.channel_unit)}</td>
            <td class="highlight-max">{me.characteristic:.4E}</td>
            <td class="highlight-max">{me.design:.4E}</td>
            <td>{me.source_dlc}</td>
            <td>{me.wind_speed:.1f}</td>
            <td>{me.safety_factor:.2f}</td>
            <td class="highlight-min">{ie.characteristic:.4E}</td>
            <td class="highlight-min">{ie.design:.4E}</td>
            <td>{ie.source_dlc}</td>
        </tr>"""

        return f"""
<section id="extreme-loads">
    <h2>Extreme Loads</h2>
    <p>Characteristic and design extreme loads across all ultimate DLCs.
    Design values include partial safety factors per IEC 61400-1 Ed.4 Table 3.</p>
    <table>
        <thead>
            <tr>
                <th rowspan="2">Channel</th>
                <th rowspan="2">Unit</th>
                <th colspan="5">Maximum</th>
                <th colspan="3">Minimum</th>
            </tr>
            <tr>
                <th>Char.</th>
                <th>Design</th>
                <th>DLC</th>
                <th>V_hub</th>
                <th>gamma</th>
                <th>Char.</th>
                <th>Design</th>
                <th>DLC</th>
            </tr>
        </thead>
        <tbody>{rows}
        </tbody>
    </table>
</section>"""

    def _section_load_envelopes(
        self, envelopes: list[LoadEnvelopeResult]
    ) -> str:
        """Generate the Load Envelopes section with scatter plots."""
        if not envelopes:
            return """
<section id="load-envelopes">
    <h2>Load Envelopes</h2>
    <p>No load envelope results available.</p>
</section>"""

        plots_html = ""
        for idx, env in enumerate(envelopes):
            chart_id = f"envelope-chart-{idx}"

            # Prepare data traces
            traces = []

            # DLC-grouped scatter points
            for dlc_num, pts in sorted(env.per_dlc.items()):
                if not pts:
                    continue
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                traces.append({
                    "x": xs[:500],  # Limit points for performance
                    "y": ys[:500],
                    "mode": "markers",
                    "type": "scatter",
                    "name": f"DLC {dlc_num}",
                    "marker": {"size": 3, "opacity": 0.5},
                })

            # Convex hull outline
            if env.hull_points:
                hull_x = [p[0] for p in env.hull_points]
                hull_y = [p[1] for p in env.hull_points]
                # Close the hull
                hull_x.append(hull_x[0])
                hull_y.append(hull_y[0])
                traces.append({
                    "x": hull_x,
                    "y": hull_y,
                    "mode": "lines",
                    "type": "scatter",
                    "name": "Convex Hull",
                    "line": {"color": "red", "width": 2},
                })

            traces_json = json.dumps(traces)

            plots_html += f"""
    <div class="plot-container">
        <h3>{_esc(env.channel_x)} vs {_esc(env.channel_y)}</h3>
        <p>Hull area: {env.hull_area:.4E} {_esc(env.unit_x)}*{_esc(env.unit_y)}</p>
        <div id="{chart_id}" style="height: 500px;"></div>
        <script>
            Plotly.newPlot('{chart_id}', {traces_json}, {{
                title: '{_esc(env.channel_x)} vs {_esc(env.channel_y)} Load Envelope',
                xaxis: {{title: '{_esc(env.channel_x)} ({_esc(env.unit_x)})'}},
                yaxis: {{title: '{_esc(env.channel_y)} ({_esc(env.unit_y)})'}},
                showlegend: true,
            }});
        </script>
    </div>"""

        return f"""
<section id="load-envelopes">
    <h2>Load Envelopes</h2>
    <p>2D load envelopes showing correlated loading between pairs of
    channels. The convex hull defines the design envelope.</p>
    {plots_html}
</section>"""

    def _section_conclusions(self, data: ReportData) -> str:
        """Generate the Conclusions section."""
        n_extreme = len(data.extreme_results)
        n_del = len(data.del_results)
        n_envelope = len(data.load_envelopes)

        return f"""
<section id="conclusions">
    <h2>Conclusions</h2>
    <p>The IEC 61400-1 Ed.4 load analysis has been completed with the
    following scope:</p>
    <ul style="margin: 15px 0; padding-left: 25px;">
        <li>{len(data.dlc_matrix)} design load cases analyzed</li>
        <li>{data.total_cases} total simulation cases executed</li>
        <li>{data.total_simulation_hours:.1f} hours of simulation time</li>
        <li>{n_extreme} channels assessed for extreme loads</li>
        <li>{n_del} channels assessed for fatigue (DEL)</li>
        <li>{n_envelope} load envelope pairs computed</li>
    </ul>
    <p>All results should be reviewed against the turbine design limits
    and material strength values. The partial safety factors applied
    follow IEC 61400-1:2019 (Ed.4) Tables 3 and 5.</p>
    <p><strong>Note:</strong> This report is generated automatically by
    WindForge and should be reviewed by a qualified structural engineer
    before use in certification or design decisions.</p>
</section>"""


def _esc(text: str) -> str:
    """Escape HTML special characters.

    Parameters
    ----------
    text : str
        Raw text to escape.

    Returns
    -------
    str
        HTML-safe text.
    """
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#x27;")
    )
