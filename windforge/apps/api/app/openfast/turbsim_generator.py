"""
TurbSim input file generator.

Generates input files for the TurbSim stochastic turbulence simulator
used to create full-field wind files for OpenFAST simulations.

Supports IEC Kaimal (IECKAI) and IEC von Karman (IECVKM) spectral models,
as well as IEC wind types (NTM, ETM, EWM).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
import random


class TurbulenceModel(str, Enum):
    """IEC turbulence spectral model."""
    IECKAI = "IECKAI"   # IEC Kaimal
    IECVKM = "IECVKM"   # IEC von Karman


class IECWindType(str, Enum):
    """IEC wind type for TurbSim."""
    NTM = "NTM"     # Normal Turbulence Model
    xETM = "xETM"   # Extreme Turbulence Model (x = turbine class)
    xEWM1 = "xEWM1" # Extreme Wind Model - 1 year recurrence
    xEWM50 = "xEWM50" # Extreme Wind Model - 50 year recurrence


class IECTurbulenceClass(str, Enum):
    """IEC turbulence characteristic (category)."""
    A = "A"   # High turbulence (Iref = 0.16)
    B = "B"   # Medium turbulence (Iref = 0.14)
    C = "C"   # Low turbulence (Iref = 0.12)


@dataclass
class TurbSimConfig:
    """Configuration for the TurbSim input file."""

    # --- Runtime Options ---
    echo: bool = False
    rand_seed1: int = -1          # First random seed (-2147483648 to 2147483647). -1 = auto-generate
    rand_seed2: str = "RanLux"    # Second random seed: "RanLux" or "RNSNLW"
    wrap_files: bool = False

    # --- Turbine/Model Specifications ---
    num_grid_z: int = 31          # Vertical grid points
    num_grid_y: int = 31          # Horizontal grid points
    time_step: float = 0.05       # Time step (s)
    analysis_time: float = 660.0  # Length of analysis time series (s) (use 630s + 30s spinup)
    usable_time: str = "ALL"      # Usable time (s) or "ALL"
    hub_ht: float = 90.0          # Hub height (m)
    grid_height: float = 150.0    # Grid height (m)
    grid_width: float = 150.0     # Grid width (m)
    v_flow_ang: float = 0.0       # Vertical mean flow angle (deg)
    h_flow_ang: float = 0.0       # Horizontal mean flow angle (deg)

    # --- Meteorological Boundary Conditions ---
    turb_model: TurbulenceModel = TurbulenceModel.IECKAI
    usable_turb_model: str = "default"  # Usually same as TurbModel
    iec_standard: str = "1-ED3"    # IEC standard edition (1-ED2, 1-ED3, 1-ED4)
    iec_turbc: str = "B"           # IEC turbulence characteristic (A, B, C, or TI in %%)
    iec_wind_type: str = "NTM"     # IEC wind type (NTM, xETM, xEWM1, xEWM50)
    etm_c: str = "default"        # IEC ETM c parameter (m/s) or "default"
    wind_profile_type: str = "IEC" # Wind profile type (LOG, PL, IEC, JET, default)
    profile_file: str = "unused"   # Wind profile data file
    ref_ht: float = 90.0          # Reference height for wind speed (m)
    u_ref: float = 11.4           # Reference mean wind speed at RefHt (m/s)
    z_jet: str = "default"        # Jet height (m) or "default"
    pl_exp: str = "default"       # Power law exponent or "default"
    z0: str = "default"           # Surface roughness length (m) or "default"

    # --- Non-IEC Meteorological Parameters ---
    latitude: str = "default"
    rich_no: float = 0.05
    usr_vkm: str = "default"
    usr_vkm_l: str = "default"
    ustar_file: str = "default"
    ustar_profile: str = "default"

    # --- Spatial Coherence Parameters ---
    scy_lux: str = "default"      # u-component coherence scale
    incdc1: str = "default"
    incdc2: str = "default"
    coh_exp: str = "default"      # Coherence exponent

    # --- Coherent Turbulence Scaling Parameters ---
    ct_event_path: str = "unused"
    ct_event_file: str = "random"
    randomize: bool = True
    num_events: int = 1000
    ct_scaling: str = "default"

    # --- Output ---
    wr_bts: bool = True           # Write .bts binary output?
    wr_formatted: bool = False    # Write formatted output?
    wr_ats: bool = False          # Write .ats (AeroDyn) output?
    wr_bladed: bool = False       # Write Bladed-style output?
    wr_hawc: bool = False         # Write HAWC format output?
    wr_vtk: bool = False          # Write VTK output?
    wr_twrs: bool = False         # Write tower data?

    @property
    def actual_seed1(self) -> int:
        """Return the random seed, generating one if set to -1."""
        if self.rand_seed1 == -1:
            return random.randint(1, 2147483647)
        return self.rand_seed1


class TurbSimGenerator:
    """Generates TurbSim input files for wind field generation."""

    def generate_turbsim_input(self, config: TurbSimConfig) -> str:
        """Generate the TurbSim input file.

        Parameters
        ----------
        config : TurbSimConfig
            Complete TurbSim configuration.

        Returns
        -------
        str
            Complete TurbSim input file content.
        """
        lines: list[str] = []
        _a = lines.append
        _f = _flag

        _a("---------TurbSim v2.00.* Input File------------------------")
        _a("Generated by WindForge for IEC wind turbine load analysis")
        _a("")

        # --- Runtime Options ---
        _a("---------Runtime Options-----------------------------------")
        _a(f"{_f(config.echo):<14s}   Echo            - Echo input data to <RootName>.ech (flag)")
        _a(f"{config.actual_seed1:<14d}   RandSeed1       - First random seed  (-2147483648 to 2147483647)")
        _a(f'{"\"" + config.rand_seed2 + "\"":<14s}   RandSeed2       - Second random seed {{\"RanLux\" or \"RNSNLW\"}} (-)')
        _a(f"{_f(config.wrap_files):<14s}   WrBHHTP         - Output hub-height turbulence parameters in binary form?  (Generates RootName.bin) (flag)")
        _a(f"{_f(False):<14s}   WrFHHTP         - Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat) (flag)")
        _a(f"{_f(False):<14s}   WrADHH          - Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh) (flag)")
        _a(f"{_f(False):<14s}   WrADFF          - Output full-field time-series data in TurbSim/AeroDyn form? (Generates RootName.bts) (flag)")
        _a(f"{_f(config.wr_bts):<14s}   WrBLFF          - Output full-field time-series data in BLADED/AeroDyn form? (Generates RootName.wnd) (flag)")
        _a(f"{_f(config.wr_ats):<14s}   WrADTWR         - Output tower time-series data? (Generates RootName.twr) (flag)")
        _a(f"{_f(config.wr_formatted):<14s}   WrFMTFF         - Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, .v, .w) (flag)")
        _a(f"{_f(False):<14s}   WrACT           - Output coherent turbulence time step file? (Generates RootName.cts) (flag)")
        _a(f"{_f(config.wr_vtk):<14s}   Clockwise       - Clockwise rotation looking downwind? (used only for full-field binary files - notass with AeroDyn) (flag)")
        _a(f"{'0':<14s}   ScaleIEC        - Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale 2=use grid scale] (-)")
        _a("")

        # --- Turbine/Model Specifications ---
        _a("---------Turbine/Model Specifications-----------------------")
        _a(f"{config.num_grid_z:<14d}   NumGrid_Z       - Vertical grid-point matrix dimension (-)")
        _a(f"{config.num_grid_y:<14d}   NumGrid_Y       - Horizontal grid-point matrix dimension (-)")
        _a(f"{config.time_step:<14.4f}   TimeStep        - Time step [seconds]")
        _a(f"{config.analysis_time:<14.1f}   AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )")
        _a(f'{"\"" + config.usable_time + "\"":<14s}   UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is \"ALL\")')
        _a(f"{config.hub_ht:<14.1f}   HubHt           - Hub height [m] (should be > 0.5*GridHeight)")
        _a(f"{config.grid_height:<14.1f}   GridHeight      - Grid height [m]")
        _a(f"{config.grid_width:<14.1f}   GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+owner shaft) )")
        _a(f"{config.v_flow_ang:<14.1f}   VFlowAng        - Vertical mean flow (uptilt) angle [degrees]")
        _a(f"{config.h_flow_ang:<14.1f}   HFlowAng        - Horizontal mean flow (skew) angle [degrees]")
        _a("")

        # --- Meteorological Boundary Conditions ---
        _a("---------Meteorological Boundary Conditions-------------------")
        _a(f'{"\"" + config.turb_model.value + "\"":<14s}   TurbModel       - Turbulence model (\"IECKAI\",\"IECVKM\",\"GP_LLJ\",\"NWTCUP\",\"SMOOTH\",\"WF_UPW\",\"WF_07D\",\"WF_14D\",\"TIDAL\",\"API\",\"USRINP\",\"TIMESR\", or \"NONE\") (-)')
        _a(f'{"\"" + config.usable_turb_model + "\"":<14s}   UserHH_or_default - Name of the file that contains inputs for user-defined HH or default settings (-)')
        _a(f'{"\"" + config.iec_standard + "\"":<14s}   IECstandard     - Number of IEC 61400-x]standard edition number (x=1,2, or 3) with optional 61400-1 edition number (e.g. \"1-ED2\") (-)')
        _a(f'{"\"" + config.iec_turbc + "\"":<14s}   IECturbc        - IEC turbulence characteristic (\"A\", \"B\", \"C\" or the turbulence intensity in percent) (\"KHTEST\" option with NWTCUP model, not used for other models)')
        _a(f'{"\"" + config.iec_wind_type + "\"":<14s}   IEC_WindType    - IEC turbulence type (\"NTM\"=normal, \"xETM\"=extreme turbulence, \"xEWM1\"=extreme 1-yr wind, \"xEWM50\"=extreme 50-yr wind, where x=wind turbine class 1, 2, or 3) (-)')
        _a(f'{"\"" + config.etm_c + "\"":<14s}   ETMc            - IEC ETM \"c\" parameter [m/s] (or \"default\")')
        _a(f'{"\"" + config.wind_profile_type + "\"":<14s}   WindProfileType - Velocity profile type (\"LOG\";\"PL\"=power law;\"JET\";\"IEC\"=PL on rotor disk, LOG elsewhere; or \"default\")')
        _a(f'{"\"" + config.profile_file + "\"":<14s}   ProfileFile     - Name of the file that contains input profiles for WindProfileType=\"USR\" and target turbulence profiles for the UserHH or TIMESR TurbModels (-)')
        _a(f"{config.ref_ht:<14.1f}   RefHt           - Height of the reference speed (URef) [m]")
        _a(f"{config.u_ref:<14.2f}   URef            - Mean (total) wind speed at the reference height [m/s] (or \"default\" for JET type)")
        _a(f'{"\"" + config.z_jet + "\"":<14s}   ZJetMax         - Jet height [m] (used only for JET type, valid 70-490 m)')
        _a(f'{"\"" + config.pl_exp + "\"":<14s}   PLExp           - Power law exponent [-] (or \"default\")')
        _a(f'{"\"" + config.z0 + "\"":<14s}   Z0              - Surface roughness length [m] (or \"default\")')
        _a("")

        # --- Non-IEC Meteorological Parameters ---
        _a("---------Non-IEC Meteorological Boundary Conditions------------")
        _a(f'{"\"" + config.latitude + "\"":<14s}   Latitude        - Site latitude [degrees] (or \"default\")')
        _a(f"{config.rich_no:<14.2f}   RICH_NO         - Gradient Richardson number (-)")
        _a(f'{"\"" + config.usr_vkm + "\"":<14s}   UStar           - Friction or shear velocity [m/s] (or \"default\")')
        _a(f'{"\"" + config.usr_vkm_l + "\"":<14s}   ZI              - Mixing layer depth [m] (or \"default\")')
        _a(f'{"\"" + config.ustar_file + "\"":<14s}   PC_UW           - Hub mean uw Reynolds stress [(m/s)^2] (or \"default\")')
        _a(f'{"\"" + config.ustar_profile + "\"":<14s}   PC_UV           - Hub mean uv Reynolds stress [(m/s)^2] (or \"default\")')
        _a(f'{"\"default\"":<14s}   PC_VW           - Hub mean vw Reynolds stress [(m/s)^2] (or \"default\")')
        _a("")

        # --- Spatial Coherence Parameters ---
        _a("---------Spatial Coherence Parameters----------------------------")
        _a(f'{"\"" + config.scy_lux + "\"":<14s}   SCMod1          - u-component  coherence model (\"GENERAL\", \"IEC\", \"API\", \"NONE\", or \"default\") (-)')
        _a(f'{"\"" + config.incdc1 + "\"":<14s}   SCMod2          - v-component  coherence model (\"GENERAL\", \"IEC\", \"NONE\", or \"default\") (-)')
        _a(f'{"\"" + config.incdc2 + "\"":<14s}   SCMod3          - w-component  coherence model (\"GENERAL\", \"IEC\", \"NONE\", or \"default\") (-)')
        _a(f'{"\"" + config.coh_exp + "\"":<14s}   InCDec1         - u-component  coherence parameters for the GENERAL model [-, m^-1] (e.g. \"10.0  0.3e-3\" in quotes) (or \"default\")')
        _a(f'{"\"default\"":<14s}   InCDec2         - v-component  coherence parameters for the GENERAL model [-, m^-1] (e.g. \"10.0  0.3e-3\" in quotes) (or \"default\")')
        _a(f'{"\"default\"":<14s}   InCDec3         - w-component  coherence parameters for the GENERAL model [-, m^-1] (e.g. \"10.0  0.3e-3\" in quotes) (or \"default\")')
        _a(f'{"\"default\"":<14s}   CohExp          - Coherence exponent for general model [-] (or \"default\")')
        _a("")

        # --- Coherent Turbulence ---
        _a("---------Coherent Turbulence Scaling Parameters-------------------")
        _a(f'{"\"" + config.ct_event_path + "\"":<40s}   CTEventPath     - Name of the path where event data files are located')
        _a(f'{"\"" + config.ct_event_file + "\"":<14s}   CTEventFile     - Type of event files (\"LES\", \"DNS\", or \"RANDOM\") (-)')
        _a(f"{_f(config.randomize):<14s}   Randomize       - Randomize the disturbance scale and location? [true/false]")
        _a(f"{config.num_events:<14d}   NumEvents       - Number of discrete events to model (-)")
        _a(f'{"\"" + config.ct_scaling + "\"":<14s}   CTScaling       - Scaling factor for the turbulence [-] (or \"default\")')
        _a("")

        # --- Output ---
        # Note: The actual output format flags were written above in Runtime Options
        # per TurbSim's expected format.

        return "\n".join(lines) + "\n"


def _flag(value: bool) -> str:
    """Convert a boolean to OpenFAST flag format."""
    return "True" if value else "False"
