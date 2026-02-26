"""
Master OpenFAST input file generator.

Orchestrates all sub-generators to produce the complete set of input
files needed for a single OpenFAST simulation case. This is the main
entry point for the file generation pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from .fst_generator import FSTConfig, FSTGenerator
from .elastodyn_generator import (
    ElastoDynConfig,
    ElastoDynBladeConfig,
    ElastoDynTowerConfig,
    ElastoDynGenerator,
)
from .aerodyn_generator import (
    AeroDynConfig,
    AeroDynBladeConfig,
    AeroDynGenerator,
)
from .servodyn_generator import (
    ServoDynConfig,
    DISCONConfig,
    ServoDynGenerator,
)
from .inflowwind_generator import (
    InflowWindConfig,
    InflowWindGenerator,
    WindType,
)
from .turbsim_generator import TurbSimConfig, TurbSimGenerator


# ---------------------------------------------------------------------------
# Domain model interfaces (simplified stand-ins for full ORM models)
# ---------------------------------------------------------------------------

@dataclass
class TurbineModel:
    """Simplified turbine model data container.

    In a full application this would correspond to the ORM model
    from the database, with all turbine specification fields.
    """
    name: str = "NREL_5MW"
    num_blades: int = 3
    tip_radius: float = 63.0
    hub_radius: float = 1.5
    hub_height: float = 90.0
    tower_height: float = 87.6
    tower_base_height: float = 10.0
    rated_power_kw: float = 5000.0
    rated_wind_speed: float = 11.4
    rated_rotor_speed: float = 12.1
    cut_in_wind_speed: float = 3.0
    cut_out_wind_speed: float = 25.0
    gearbox_ratio: float = 97.0
    generator_efficiency: float = 94.4
    rotor_overhang: float = -5.0191
    shaft_tilt: float = -5.0
    precone: float = -2.5

    # Sub-configs (populated by the application layer)
    elastodyn_config: Optional[ElastoDynConfig] = None
    blade_config: Optional[ElastoDynBladeConfig] = None
    tower_config: Optional[ElastoDynTowerConfig] = None
    aerodyn_config: Optional[AeroDynConfig] = None
    aerodyn_blade_config: Optional[AeroDynBladeConfig] = None
    servodyn_config: Optional[ServoDynConfig] = None
    discon_config: Optional[DISCONConfig] = None


@dataclass
class SimulationCase:
    """A single simulation case specification.

    Combines wind conditions, controller settings, and simulation
    parameters for one OpenFAST run.
    """
    case_id: str = "DLC1.1_V11.4_S001"
    dlc_number: str = "1.1"
    wind_speed: float = 11.4
    seed_number: int = 1
    yaw_misalignment: float = 0.0
    simulation_time: float = 630.0
    dt: float = 0.005
    initial_rotor_speed: float = 12.1
    initial_blade_pitch: float = 0.0

    # Wind conditions
    wind_type: int = 3          # 1=steady, 3=TurbSim FF
    turbsim_config: Optional[TurbSimConfig] = None

    # Inflow
    inflowwind_config: Optional[InflowWindConfig] = None


@dataclass
class Project:
    """Project-level configuration."""
    name: str = "WindForge_Project"
    iec_class: str = "IB"
    turbine_class: int = 1     # IEC wind turbine class (I, II, III)
    turbulence_class: str = "B"  # A, B, or C


class OpenFASTFileGenerator:
    """Master generator that orchestrates all OpenFAST input file generation.

    This class coordinates the generation of all input files needed for a
    single OpenFAST simulation case, including the primary .fst file and
    all referenced module input files.

    Usage
    -----
    >>> generator = OpenFASTFileGenerator()
    >>> turbine = TurbineModel(name="NREL_5MW")
    >>> case = SimulationCase(wind_speed=11.4, seed_number=1)
    >>> project = Project(name="MyProject")
    >>> files = generator.generate_all(turbine, case, project)
    >>> for filename, content in files.items():
    ...     with open(filename, "w") as f:
    ...         f.write(content)
    """

    def __init__(self) -> None:
        """Initialize all sub-generators."""
        self._fst_gen = FSTGenerator()
        self._ed_gen = ElastoDynGenerator()
        self._ad_gen = AeroDynGenerator()
        self._srvd_gen = ServoDynGenerator()
        self._ifw_gen = InflowWindGenerator()
        self._ts_gen = TurbSimGenerator()

    def generate_all(
        self,
        turbine_model: TurbineModel,
        sim_case: SimulationCase,
        project: Project,
    ) -> dict[str, str]:
        """Generate all input files for a single simulation case.

        Produces a dictionary of {filename: content} for every input
        file needed by OpenFAST. The caller is responsible for writing
        these files to disk.

        Parameters
        ----------
        turbine_model : TurbineModel
            Complete turbine specification including structural and
            aerodynamic data.
        sim_case : SimulationCase
            Case-specific parameters (wind speed, seed, DLC, etc.).
        project : Project
            Project-level settings (IEC class, name).

        Returns
        -------
        dict[str, str]
            Mapping from filename to file content string for all
            generated input files.
        """
        files: dict[str, str] = {}
        root_name = sim_case.case_id

        # --- File naming convention ---
        ed_main_file = f"{root_name}_ElastoDyn.dat"
        ed_blade_file = f"{root_name}_ElastoDyn_Blade.dat"
        ed_tower_file = f"{root_name}_ElastoDyn_Tower.dat"
        ad_main_file = f"{root_name}_AeroDyn.dat"
        ad_blade_file = f"{root_name}_AeroDyn_Blade.dat"
        srvd_file = f"{root_name}_ServoDyn.dat"
        ifw_file = f"{root_name}_InflowWind.dat"
        discon_file = "DISCON.IN"
        ts_file = f"{root_name}_TurbSim.inp"

        # ----------------------------------------------------------------
        # 1. ElastoDyn files
        # ----------------------------------------------------------------
        ed_config = turbine_model.elastodyn_config or self._build_default_ed_config(
            turbine_model, sim_case
        )
        ed_config.blade_file = ed_blade_file
        ed_config.tower_file = ed_tower_file

        files[ed_main_file] = self._ed_gen.generate_main_file(ed_config)

        blade_config = turbine_model.blade_config or ElastoDynBladeConfig()
        files[ed_blade_file] = self._ed_gen.generate_blade_file(blade_config)

        tower_config = turbine_model.tower_config or ElastoDynTowerConfig()
        files[ed_tower_file] = self._ed_gen.generate_tower_file(tower_config)

        # ----------------------------------------------------------------
        # 2. AeroDyn files
        # ----------------------------------------------------------------
        ad_config = turbine_model.aerodyn_config or AeroDynConfig(
            num_bl=turbine_model.num_blades,
            ad_bl_file=ad_blade_file,
        )
        ad_config.ad_bl_file = ad_blade_file

        files[ad_main_file] = self._ad_gen.generate_main_file(ad_config)

        ad_blade_config = turbine_model.aerodyn_blade_config or AeroDynBladeConfig()
        files[ad_blade_file] = self._ad_gen.generate_blade_file(ad_blade_config)

        # ----------------------------------------------------------------
        # 3. ServoDyn files
        # ----------------------------------------------------------------
        srvd_config = turbine_model.servodyn_config or ServoDynConfig(
            dll_in_file=discon_file,
        )
        srvd_config.dll_in_file = discon_file
        files[srvd_file] = self._srvd_gen.generate_servodyn_file(srvd_config)

        discon_cfg = turbine_model.discon_config or DISCONConfig(
            we_blade_radius=turbine_model.tip_radius,
            we_gear_ratio=turbine_model.gearbox_ratio,
            vs_rated_gen_pwr=turbine_model.rated_power_kw * 1000.0,
            we_rated_pwr=turbine_model.rated_power_kw * 1000.0,
            we_rated_v=turbine_model.rated_wind_speed,
        )
        files[discon_file] = self._srvd_gen.generate_discon_in(discon_cfg)

        # ----------------------------------------------------------------
        # 4. InflowWind file
        # ----------------------------------------------------------------
        ifw_config = sim_case.inflowwind_config or self._build_default_ifw_config(
            sim_case, turbine_model
        )
        files[ifw_file] = self._ifw_gen.generate_inflowwind_file(ifw_config)

        # ----------------------------------------------------------------
        # 5. TurbSim input (only for turbulent cases)
        # ----------------------------------------------------------------
        if sim_case.wind_type == 3:
            ts_config = sim_case.turbsim_config or self._build_default_ts_config(
                sim_case, turbine_model, project
            )
            files[ts_file] = self._ts_gen.generate_turbsim_input(ts_config)

        # ----------------------------------------------------------------
        # 6. Primary .fst file (must reference all other files)
        # ----------------------------------------------------------------
        fst_config = FSTConfig(
            t_max=sim_case.simulation_time,
            dt=sim_case.dt,
            comp_elast=1,
            comp_inflow=1 if sim_case.wind_type > 0 else 0,
            comp_aero=2,
            comp_servo=1,
            comp_hydro=0,
            comp_sub=0,
            comp_mooring=0,
            comp_ice=0,
            air_dens=1.225,
            gravity=9.80665,
            ed_file=ed_main_file,
            inflow_file=ifw_file,
            aero_file=ad_main_file,
            servo_file=srvd_file,
        )
        files[f"{root_name}.fst"] = self._fst_gen.generate(
            fst_config, root_name=root_name
        )

        return files

    def _build_default_ed_config(
        self,
        turbine: TurbineModel,
        case: SimulationCase,
    ) -> ElastoDynConfig:
        """Build a default ElastoDyn configuration from turbine parameters.

        Parameters
        ----------
        turbine : TurbineModel
            Turbine specification.
        case : SimulationCase
            Simulation case parameters.

        Returns
        -------
        ElastoDynConfig
            Populated ElastoDyn configuration.
        """
        return ElastoDynConfig(
            num_bl=turbine.num_blades,
            tip_rad=turbine.tip_radius,
            hub_rad=turbine.hub_radius,
            tower_ht=turbine.tower_height,
            tower_bsht=turbine.tower_base_height,
            overhang=turbine.rotor_overhang,
            shft_tilt=turbine.shaft_tilt,
            pre_cone_1=turbine.precone,
            pre_cone_2=turbine.precone,
            pre_cone_3=turbine.precone,
            gb_ratio=turbine.gearbox_ratio,
            rot_speed=case.initial_rotor_speed,
            bl_pitch1=case.initial_blade_pitch,
            bl_pitch2=case.initial_blade_pitch,
            bl_pitch3=case.initial_blade_pitch,
        )

    def _build_default_ifw_config(
        self,
        case: SimulationCase,
        turbine: TurbineModel,
    ) -> InflowWindConfig:
        """Build a default InflowWind configuration.

        Parameters
        ----------
        case : SimulationCase
            Simulation case with wind type and speed.
        turbine : TurbineModel
            Turbine for reference height.

        Returns
        -------
        InflowWindConfig
            Populated InflowWind configuration.
        """
        if case.wind_type == 1:
            # Steady wind
            return InflowWindConfig(
                wind_type=WindType.STEADY,
                h_wind_speed=case.wind_speed,
                ref_ht=turbine.hub_height,
                propagation_dir=case.yaw_misalignment,
            )
        else:
            # TurbSim full-field
            ts_filename = f"{case.case_id}_TurbSim.bts"
            return InflowWindConfig(
                wind_type=WindType.TURBSIM_FF,
                turbsim_filename=ts_filename,
                ref_ht=turbine.hub_height,
                h_wind_speed=case.wind_speed,
                propagation_dir=case.yaw_misalignment,
            )

    def _build_default_ts_config(
        self,
        case: SimulationCase,
        turbine: TurbineModel,
        project: Project,
    ) -> TurbSimConfig:
        """Build a default TurbSim configuration.

        Parameters
        ----------
        case : SimulationCase
            Case with wind speed and seed.
        turbine : TurbineModel
            Turbine for hub height and rotor dimensions.
        project : Project
            Project for IEC class settings.

        Returns
        -------
        TurbSimConfig
            Populated TurbSim configuration.
        """
        from .turbsim_generator import TurbulenceModel

        # Grid dimensions should cover the rotor with some margin
        rotor_diameter = 2.0 * turbine.tip_radius
        grid_size = rotor_diameter * 1.2  # 20% margin

        return TurbSimConfig(
            rand_seed1=case.seed_number,
            hub_ht=turbine.hub_height,
            grid_height=grid_size,
            grid_width=grid_size,
            u_ref=case.wind_speed,
            ref_ht=turbine.hub_height,
            analysis_time=case.simulation_time + 30.0,  # Extra for spinup
            iec_turbc=project.turbulence_class,
            iec_wind_type="NTM",
            turb_model=TurbulenceModel.IECKAI,
        )
