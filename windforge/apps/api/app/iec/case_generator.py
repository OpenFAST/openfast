"""
IEC DLC case expansion engine.

Takes a user-selected set of DLCs with configuration and expands them
into a complete list of simulation case specifications. Each case spec
fully defines one OpenFAST simulation run with all wind/turbulence
parameters resolved.

The expansion handles:
  - Wind speed binning (Vin to Vout in user-defined steps)
  - Turbulent seed generation (multiple seeds per wind speed bin)
  - TurbSim parameter computation with correct turbulence model
  - Yaw misalignment application
  - Special wind speeds (Vr, Vout, Ve50, Ve1)
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Optional

from .dlc_matrix import (
    DLCTemplate,
    IEC_DLC_TABLE,
    WindModel,
    AnalysisType,
    get_dlc_by_number,
)
from .wind_conditions import (
    get_iref,
    get_vref,
    ntm_sigma,
    etm_sigma,
    ewm_ve50,
    ewm_ve1,
    wind_speed_bin_probability,
)
from .safety_factors import get_load_factor


@dataclass
class DLCSelection:
    """User selection and configuration for a single DLC.

    Attributes
    ----------
    dlc_number : str
        DLC identifier (e.g., "1.1", "6.2").
    enabled : bool
        Whether this DLC is enabled for analysis.
    num_seeds : Optional[int]
        Override for number of seeds (None = use DLC default).
    wind_speed_step : float
        Wind speed bin width (m/s) for operating range DLCs.
    yaw_misalignment_override : Optional[float]
        Override for yaw misalignment angle (degrees).
    custom_wind_speeds : Optional[list[float]]
        If provided, use these exact wind speeds instead of binning.
    """
    dlc_number: str
    enabled: bool = True
    num_seeds: Optional[int] = None
    wind_speed_step: float = 2.0
    yaw_misalignment_override: Optional[float] = None
    custom_wind_speeds: Optional[list[float]] = None


@dataclass
class ProjectConfig:
    """Project-level configuration for case generation.

    Attributes
    ----------
    name : str
        Project name used for case ID generation.
    turbine_class : int
        IEC wind turbine class (1, 2, or 3).
    turbulence_class : str
        IEC turbulence category ("A", "B", or "C").
    v_in : float
        Cut-in wind speed (m/s).
    v_out : float
        Cut-out wind speed (m/s).
    v_rated : float
        Rated wind speed (m/s).
    hub_height : float
        Hub height above ground (m).
    rotor_diameter : float
        Rotor diameter (m).
    simulation_time : float
        Simulation time for each case (s).
    dt : float
        Simulation time step (s).
    initial_rotor_speed : float
        Initial rotor speed (rpm).
    base_seed : int
        Base random seed for reproducibility.
    """
    name: str = "WindForge"
    turbine_class: int = 1
    turbulence_class: str = "B"
    v_in: float = 3.0
    v_out: float = 25.0
    v_rated: float = 11.4
    hub_height: float = 90.0
    rotor_diameter: float = 126.0
    simulation_time: float = 630.0
    dt: float = 0.005
    initial_rotor_speed: float = 12.1
    base_seed: int = 12345


@dataclass
class SimulationCaseSpec:
    """Fully specified simulation case ready for file generation.

    Contains all parameters needed to generate OpenFAST input files
    for a single simulation run.

    Attributes
    ----------
    case_id : str
        Unique identifier for this case.
    dlc_number : str
        Source DLC number.
    dlc_description : str
        Human-readable DLC description.
    wind_speed : float
        Hub-height mean wind speed (m/s).
    seed_number : int
        Turbulent seed number.
    random_seed : int
        Actual random seed value for TurbSim.
    yaw_misalignment : float
        Yaw misalignment angle (degrees).
    simulation_time : float
        Total simulation time (s).
    dt : float
        Simulation time step (s).
    initial_rotor_speed : float
        Initial rotor speed (rpm).
    wind_type : int
        InflowWind wind type (1=steady, 3=TurbSim).
    turbulence_model : str
        TurbSim turbulence model name.
    iec_wind_type : str
        IEC wind type string for TurbSim.
    turbulence_class : str
        IEC turbulence class ("A", "B", "C").
    sigma_1 : float
        Target longitudinal turbulence std dev (m/s).
    analysis_type : str
        "ULTIMATE" or "FATIGUE".
    safety_factor : float
        Partial safety factor for loads (gamma_f).
    wind_speed_probability : float
        Weibull probability weight for this wind speed bin.
    """
    case_id: str
    dlc_number: str
    dlc_description: str
    wind_speed: float
    seed_number: int
    random_seed: int
    yaw_misalignment: float
    simulation_time: float
    dt: float
    initial_rotor_speed: float
    wind_type: int
    turbulence_model: str
    iec_wind_type: str
    turbulence_class: str
    sigma_1: float
    analysis_type: str
    safety_factor: float
    wind_speed_probability: float = 1.0


class CaseGenerator:
    """Expands DLC configurations into fully specified simulation cases.

    This class takes a list of user-selected DLCs and project
    configuration and produces a complete list of simulation case
    specifications ready for file generation and execution.

    Usage
    -----
    >>> generator = CaseGenerator()
    >>> dlc_config = [
    ...     DLCSelection(dlc_number="1.1", num_seeds=6),
    ...     DLCSelection(dlc_number="1.2", num_seeds=6),
    ... ]
    >>> project = ProjectConfig(turbine_class=1, turbulence_class="B")
    >>> cases = generator.expand_dlc_matrix(dlc_config, project)
    >>> print(f"Generated {len(cases)} simulation cases")
    """

    def expand_dlc_matrix(
        self,
        dlc_config: list[DLCSelection],
        project: ProjectConfig,
    ) -> list[SimulationCaseSpec]:
        """Expand the DLC matrix into individual simulation cases.

        For each enabled DLC:
        1. Determine wind speed bins based on operating range or
           special wind speeds.
        2. Generate random seeds for each wind speed bin.
        3. Compute turbulence parameters (sigma, turbulence model).
        4. Apply yaw misalignment.
        5. Calculate Weibull probability weights for fatigue DLCs.

        Parameters
        ----------
        dlc_config : list[DLCSelection]
            User-selected DLCs with configuration overrides.
        project : ProjectConfig
            Project-level settings.

        Returns
        -------
        list[SimulationCaseSpec]
            Complete list of simulation case specifications.
        """
        cases: list[SimulationCaseSpec] = []
        v_ref = get_vref(project.turbine_class)
        i_ref = get_iref(project.turbulence_class)
        v_ave = 0.2 * v_ref

        for dlc_sel in dlc_config:
            if not dlc_sel.enabled:
                continue

            dlc = get_dlc_by_number(dlc_sel.dlc_number)
            if dlc is None:
                continue

            # Determine wind speeds for this DLC
            wind_speeds = self._resolve_wind_speeds(
                dlc, dlc_sel, project, v_ref
            )

            # Number of seeds
            num_seeds = dlc_sel.num_seeds or dlc.default_num_seeds

            # Yaw misalignment
            yaw = dlc_sel.yaw_misalignment_override
            if yaw is None:
                yaw = dlc.yaw_misalignment

            # Wind type (turbulent vs deterministic)
            wind_type = self._resolve_wind_type(dlc)

            # Turbulence model and IEC wind type for TurbSim
            turb_model = self._resolve_turbulence_model(dlc)
            iec_wt = self._resolve_iec_wind_type(dlc, project.turbine_class)

            for v_hub in wind_speeds:
                # Calculate turbulence sigma
                sigma_1 = self._compute_sigma(
                    dlc, v_hub, i_ref, v_ref
                )

                # Wind speed bin probability (for fatigue weighting)
                if dlc.analysis_type == AnalysisType.FATIGUE:
                    bin_prob = wind_speed_bin_probability(
                        max(v_hub - dlc_sel.wind_speed_step / 2.0, 0.0),
                        v_hub + dlc_sel.wind_speed_step / 2.0,
                        v_ave,
                    )
                else:
                    bin_prob = 1.0

                for seed_idx in range(1, num_seeds + 1):
                    # Generate deterministic random seed from base
                    rand_seed = self._generate_seed(
                        project.base_seed, dlc.number, v_hub, seed_idx
                    )

                    case_id = self._make_case_id(
                        project.name, dlc.number, v_hub, seed_idx
                    )

                    cases.append(SimulationCaseSpec(
                        case_id=case_id,
                        dlc_number=dlc.number,
                        dlc_description=dlc.description,
                        wind_speed=v_hub,
                        seed_number=seed_idx,
                        random_seed=rand_seed,
                        yaw_misalignment=yaw,
                        simulation_time=project.simulation_time,
                        dt=project.dt,
                        initial_rotor_speed=project.initial_rotor_speed,
                        wind_type=wind_type,
                        turbulence_model=turb_model,
                        iec_wind_type=iec_wt,
                        turbulence_class=project.turbulence_class,
                        sigma_1=sigma_1,
                        analysis_type=dlc.analysis_type.name,
                        safety_factor=dlc.default_safety_factor,
                        wind_speed_probability=bin_prob,
                    ))

        return cases

    def _resolve_wind_speeds(
        self,
        dlc: DLCTemplate,
        dlc_sel: DLCSelection,
        project: ProjectConfig,
        v_ref: float,
    ) -> list[float]:
        """Determine the wind speed values for a DLC.

        Parameters
        ----------
        dlc : DLCTemplate
            DLC template with wind speed range specification.
        dlc_sel : DLCSelection
            User selection with possible custom wind speeds.
        project : ProjectConfig
            Project config with Vin, Vout, Vrated.
        v_ref : float
            Reference wind speed.

        Returns
        -------
        list[float]
            List of wind speed values (m/s) for this DLC.
        """
        if dlc_sel.custom_wind_speeds:
            return sorted(dlc_sel.custom_wind_speeds)

        v_in = project.v_in
        v_out = project.v_out
        v_rated = project.v_rated
        step = dlc_sel.wind_speed_step

        if dlc.wind_speed_range == "Vin_to_Vout":
            # Full operating range in steps
            speeds = []
            v = v_in
            while v <= v_out + 0.01:
                speeds.append(round(v, 1))
                v += step
            if speeds[-1] < v_out:
                speeds.append(v_out)
            return speeds

        elif dlc.wind_speed_range == "Vr-2_Vr_Vr+2":
            return [v_rated - 2.0, v_rated, v_rated + 2.0]

        elif dlc.wind_speed_range == "Vr-2_Vr_Vr+2_Vout":
            return [v_rated - 2.0, v_rated, v_rated + 2.0, v_out]

        elif dlc.wind_speed_range == "Vin_Vr-2_Vr_Vr+2":
            return [v_in, v_rated - 2.0, v_rated, v_rated + 2.0]

        elif dlc.wind_speed_range == "Vin_Vr_Vout":
            return [v_in, v_rated, v_out]

        elif dlc.wind_speed_range == "Ve50":
            return [ewm_ve50(v_ref)]

        elif dlc.wind_speed_range == "Ve1":
            return [ewm_ve1(v_ref)]

        elif dlc.wind_speed_range == "Vmaint":
            # Maintenance wind speed typically 15 m/s
            return [15.0]

        else:
            # Default: rated wind speed
            return [v_rated]

    def _resolve_wind_type(self, dlc: DLCTemplate) -> int:
        """Determine the InflowWind wind type for a DLC.

        Parameters
        ----------
        dlc : DLCTemplate
            DLC template.

        Returns
        -------
        int
            Wind type: 1=steady, 3=TurbSim full-field.
        """
        # Deterministic wind models use steady/transient approaches
        if dlc.wind_model in (
            WindModel.EOG, WindModel.EDC, WindModel.ECD,
            WindModel.EWS, WindModel.NWP,
        ):
            return 1  # Steady wind (deterministic events handled separately)
        else:
            return 3  # TurbSim full-field binary

    def _resolve_turbulence_model(self, dlc: DLCTemplate) -> str:
        """Get the TurbSim turbulence model string for a DLC.

        Parameters
        ----------
        dlc : DLCTemplate
            DLC template.

        Returns
        -------
        str
            TurbSim turbulence model name.
        """
        if dlc.wind_model in (WindModel.NTM, WindModel.ETM):
            return "IECKAI"
        elif dlc.wind_model == WindModel.EWM:
            return "IECKAI"
        else:
            return "IECKAI"  # Default for deterministic also (not used)

    def _resolve_iec_wind_type(
        self, dlc: DLCTemplate, turbine_class: int
    ) -> str:
        """Get the IEC wind type string for TurbSim.

        Parameters
        ----------
        dlc : DLCTemplate
            DLC template.
        turbine_class : int
            IEC turbine class (1, 2, 3).

        Returns
        -------
        str
            IEC wind type string (e.g., "NTM", "1ETM", "1EWM50").
        """
        tc = str(turbine_class)

        if dlc.wind_model == WindModel.NTM:
            return "NTM"
        elif dlc.wind_model == WindModel.ETM:
            return f"{tc}ETM"
        elif dlc.wind_model == WindModel.EWM:
            if dlc.special_wind_speed and "Ve50" in dlc.special_wind_speed:
                return f"{tc}EWM50"
            else:
                return f"{tc}EWM1"
        else:
            return "NTM"

    def _compute_sigma(
        self,
        dlc: DLCTemplate,
        v_hub: float,
        i_ref: float,
        v_ref: float,
    ) -> float:
        """Compute the target turbulence standard deviation.

        Parameters
        ----------
        dlc : DLCTemplate
            DLC template specifying the wind model.
        v_hub : float
            Hub-height mean wind speed (m/s).
        i_ref : float
            Reference turbulence intensity.
        v_ref : float
            Reference wind speed (m/s).

        Returns
        -------
        float
            Target sigma_1 (m/s).
        """
        if dlc.wind_model == WindModel.NTM:
            return ntm_sigma(v_hub, i_ref)
        elif dlc.wind_model == WindModel.ETM:
            return etm_sigma(v_hub, i_ref, v_ref)
        elif dlc.wind_model == WindModel.EWM:
            # EWM uses NTM sigma for the turbulence in the wind field
            return ntm_sigma(v_hub, i_ref)
        else:
            # Deterministic models: sigma not directly used
            return ntm_sigma(v_hub, i_ref)

    def _generate_seed(
        self,
        base_seed: int,
        dlc_number: str,
        wind_speed: float,
        seed_idx: int,
    ) -> int:
        """Generate a deterministic random seed.

        Uses a hash-based approach to produce reproducible but
        well-distributed seeds across DLCs and wind speeds.

        Parameters
        ----------
        base_seed : int
            Project base seed.
        dlc_number : str
            DLC number string.
        wind_speed : float
            Wind speed value.
        seed_idx : int
            Seed index (1-based).

        Returns
        -------
        int
            Random seed value for TurbSim.
        """
        # Create a deterministic seed from the combination
        seed_str = f"{base_seed}_{dlc_number}_{wind_speed:.1f}_{seed_idx}"
        return abs(hash(seed_str)) % 2147483647 + 1

    def _make_case_id(
        self,
        project_name: str,
        dlc_number: str,
        wind_speed: float,
        seed_idx: int,
    ) -> str:
        """Generate a unique case identifier string.

        Parameters
        ----------
        project_name : str
            Project name.
        dlc_number : str
            DLC number.
        wind_speed : float
            Wind speed (m/s).
        seed_idx : int
            Seed index.

        Returns
        -------
        str
            Case ID string (e.g., "DLC1.1_V11.4_S001").
        """
        dlc_str = dlc_number.replace(".", "p")
        return f"DLC{dlc_str}_V{wind_speed:05.1f}_S{seed_idx:03d}"
