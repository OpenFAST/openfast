"""
AeroDyn v15 input file generator for OpenFAST.

Generates two files:
  1. AeroDyn primary input file (.dat)
  2. AeroDyn blade definition input file (.dat)

Format verified against:
  docs/source/user/aerodyn/examples/ad_primary_example.dat
  docs/source/user/aerodyn/examples/ad_blade_example.dat
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Data classes for blade aerodynamic stations
# ---------------------------------------------------------------------------

@dataclass
class AeroBladeStation:
    """A single AeroDyn blade station entry."""
    bl_spn: float       # Blade span position (m)
    bl_crv_ac: float    # Out-of-plane offset of aerodynamic center (m)
    bl_swp_ac: float    # In-plane offset of aerodynamic center (m)
    bl_crv_ang: float   # Angle of curvature (deg)
    bl_twist: float     # Blade twist (deg)
    bl_chord: float     # Blade chord length (m)
    bl_af_id: int       # Airfoil ID number (1-based index into AFNames)


@dataclass
class TowerAeroStation:
    """A single AeroDyn tower node entry."""
    twr_elev: float     # Tower elevation (m)
    twr_diam: float     # Tower diameter (m)
    twr_cd: float       # Tower drag coefficient (-)
    twr_ti: float       # Tower turbulence intensity used with TwrShadow=2 (-)
    twr_cb: float       # Tower buoyancy coefficient used with Buoyancy=True (-)


# ---------------------------------------------------------------------------
# Configuration data classes
# ---------------------------------------------------------------------------

@dataclass
class AeroDynBladeConfig:
    """Configuration for an AeroDyn blade definition file."""
    num_bl_nds: int = 19
    stations: list[AeroBladeStation] = field(default_factory=list)


@dataclass
class AeroDynConfig:
    """Configuration for the AeroDyn primary input file."""

    # --- General Options ---
    echo: bool = False
    dt_aero: str = "default"
    wake_mod: int = 1          # 0=none, 1=BEMT, 2=DBEMT, 3=OLAF
    af_aero_mod: int = 1       # 1=steady, 2=Beddoes-Leishman
    twr_potent: int = 0        # 0=none, 1=baseline, 2=Bak correction
    twr_shadow: int = 0        # 0=none, 1=Powles, 2=Eames
    twr_aero: bool = False
    frozen_wake: bool = False
    cavit_check: bool = False
    buoyancy: bool = False
    comp_aa: bool = False
    aa_input_file: str = "unused"

    # --- Environmental Conditions ---
    air_dens: str = "default"
    kin_visc: str = "default"
    spd_sound: str = "default"
    patm: str = "default"
    pvap: str = "default"

    # --- BEMT Options ---
    skew_mod: int = 1
    skew_mod_factor: str = "default"
    tip_loss: bool = True
    hub_loss: bool = True
    tan_ind: bool = True
    ai_drag: bool = True
    ti_drag: bool = True
    ind_toler: float = 1.0e-05
    max_iter: int = 100

    # --- DBEMT Options ---
    dbemt_mod: int = 2
    tau1_const: float = 4.0

    # --- OLAF Options ---
    olaf_input_file: str = "unused"

    # --- Unsteady Aero Options ---
    ua_mod: int = 4
    f_lookup: bool = True

    # --- Airfoil Information ---
    af_tab_mod: int = 1
    in_col_alfa: int = 1
    in_col_cl: int = 2
    in_col_cd: int = 3
    in_col_cm: int = 4
    in_col_cpmin: int = 0
    num_af_files: int = 1
    af_names: list[str] = field(default_factory=lambda: ["Airfoils/cylinder.dat"])

    # --- Blade Properties ---
    use_bl_cm: bool = True
    num_bl: int = 3
    ad_bl_file: str = "AeroDyn_blade.dat"

    # --- Hub Properties (Buoyancy) ---
    vol_hub: float = 0.0
    hub_cen_bx: float = 0.0

    # --- Nacelle Properties (Buoyancy) ---
    vol_nac: float = 0.0
    nac_cen_b: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # --- Tower Influence ---
    num_twr_nds: int = 0
    tower_stations: list[TowerAeroStation] = field(default_factory=list)

    # --- Outputs ---
    sum_print: bool = False
    n_bl_outs: int = 0
    bl_out_nd: list[int] = field(default_factory=list)
    n_tw_outs: int = 0
    tw_out_nd: list[int] = field(default_factory=list)
    out_list: list[str] = field(default_factory=lambda: [
        '"RtAeroCp"', '"RtAeroCt"', '"RtAeroCq"',
        '"RtTSR"', '"RtSpeed"', '"RtAeroFxh"', '"RtAeroFyh"', '"RtAeroFzh"',
        '"RtAeroMxh"', '"RtAeroMyh"', '"RtAeroMzh"',
        '"RtAeroPwr"',
    ])


class AeroDynGenerator:
    """Generates AeroDyn primary and blade definition input files."""

    def generate_main_file(self, config: AeroDynConfig) -> str:
        """Generate the AeroDyn primary input file.

        Parameters
        ----------
        config : AeroDynConfig
            Complete AeroDyn configuration parameters.

        Returns
        -------
        str
            Complete AeroDyn primary input file content.
        """
        lines: list[str] = []
        _a = lines.append
        _f = _flag

        _a("------- AERODYN v15 for OpenFAST INPUT FILE -----------------------------------------------")
        _a("Generated by WindForge - AeroDyn primary input file")

        # --- General Options ---
        _a("======  General Options  ============================================================================")
        _a(f"{_f(config.echo):<14s}   Echo               - Echo the input to \"<rootname>.AD.ech\"?  (flag)")
        _a(f'{"\"" + config.dt_aero + "\"":<14s}   DTAero             - Time interval for aerodynamic calculations {{or \"default\"}} (s)')
        _a(f"{config.wake_mod:<14d}   WakeMod            - Type of wake/induction model (switch) {{0=none, 1=BEMT, 2=DBEMT, 3=OLAF}} [WakeMod cannot be 2 or 3 when linearizing]")
        _a(f"{config.af_aero_mod:<14d}   AFAeroMod          - Type of blade airfoil aerodynamics model (switch) {{1=steady model, 2=Beddoes-Leishman unsteady model}} [AFAeroMod must be 1 when linearizing]")
        _a(f"{config.twr_potent:<14d}   TwrPotent          - Type of tower influence on wind based on potential flow around the tower (switch) {{0=none, 1=baseline potential flow, 2=potential flow with Bak correction}}")
        _a(f"{config.twr_shadow:<14d}   TwrShadow          - Calculate tower influence on wind based on downstream tower shadow (switch) {{0=none, 1=Powles model, 2=Eames model}}")
        _a(f"{_f(config.twr_aero):<14s}   TwrAero            - Calculate tower aerodynamic loads? (flag)")
        _a(f"{_f(config.frozen_wake):<14s}   FrozenWake         - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]")
        _a(f"{_f(config.cavit_check):<14s}   CavitCheck         - Perform cavitation check? (flag) [AFAeroMod must be 1 when CavitCheck=true]")
        _a(f"{_f(config.buoyancy):<14s}   Buoyancy           - Include buoyancy effects? (flag)")
        _a(f"{_f(config.comp_aa):<14s}   CompAA             - Flag to compute AeroAcoustics calculation [used only when WakeMod = 1 or 2]")
        _a(f'{"\"" + config.aa_input_file + "\"":<40s}   AA_InputFile       - AeroAcoustics input file [used only when CompAA=true]')

        # --- Environmental Conditions ---
        _a("======  Environmental Conditions  ===================================================================")
        _a(f'{"\"" + config.air_dens + "\"":<14s}   AirDens            - Air density (kg/m^3)')
        _a(f'{"\"" + config.kin_visc + "\"":<14s}   KinVisc            - Kinematic viscosity of working fluid (m^2/s)')
        _a(f'{"\"" + config.spd_sound + "\"":<14s}   SpdSound           - Speed of sound in working fluid (m/s)')
        _a(f'{"\"" + config.patm + "\"":<14s}   Patm               - Atmospheric pressure (Pa) [used only when CavitCheck=True]')
        _a(f'{"\"" + config.pvap + "\"":<14s}   Pvap               - Vapour pressure of working fluid (Pa) [used only when CavitCheck=True]')

        # --- BEMT Options ---
        _a("======  Blade-Element/Momentum Theory Options  ====================================================== [unused when WakeMod=0 or 3]")
        _a(f"{config.skew_mod:<14d}   SkewMod            - Type of skewed-wake correction model (switch) {{1=uncoupled, 2=Pitt/Peters, 3=coupled}} [unused when WakeMod=0 or 3]")
        _a(f'{"\"" + config.skew_mod_factor + "\"":<14s}   SkewModFactor      - Constant used in Pitt/Peters skewed wake model {{or \"default\" is 15/32*pi}} (-) [used only when SkewMod=2; unused when WakeMod=0 or 3]')
        _a(f"{_f(config.tip_loss):<14s}   TipLoss            - Use the Prandtl tip-loss model? (flag) [unused when WakeMod=0 or 3]")
        _a(f"{_f(config.hub_loss):<14s}   HubLoss            - Use the Prandtl hub-loss model? (flag) [unused when WakeMod=0 or 3]")
        _a(f"{_f(config.tan_ind):<14s}   TanInd             - Include tangential induction in BEMT calculations? (flag) [unused when WakeMod=0 or 3]")
        _a(f"{_f(config.ai_drag):<14s}   AIDrag             - Include the drag term in the axial-induction calculation? (flag) [unused when WakeMod=0 or 3]")
        _a(f"{_f(config.ti_drag):<14s}   TIDrag             - Include the drag term in the tangential-induction calculation? (flag) [unused when WakeMod=0,3 or TanInd=FALSE]")
        _a(f"{config.ind_toler:<14.1E}   IndToler           - Convergence tolerance for BEMT nonlinear solve residual equation {{or \"default\"}} (-) [unused when WakeMod=0 or 3]")
        _a(f"{config.max_iter:<14d}   MaxIter            - Maximum number of iteration steps (-) [unused when WakeMod=0]")

        # --- DBEMT Options ---
        _a("======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]")
        _a(f"{config.dbemt_mod:<14d}   DBEMT_Mod          - Type of dynamic BEMT (DBEMT) model {{1=constant tau1, 2=time-dependent tau1, 3=constant tau1 with continuous formulation}} (-) [used only when WakeMod=2]")
        _a(f"{config.tau1_const:<14.1f}   tau1_const         - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1 or 3]")

        # --- OLAF ---
        _a("======  OLAF -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when WakeMod=3]")
        _a(f'{"\"" + config.olaf_input_file + "\"":<40s}   OLAFInputFileName - Input file for OLAF [used only when WakeMod=3]')

        # --- Unsteady Aero Options ---
        _a("======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]")
        _a(f"{config.ua_mod:<14d}   UAMod              - Unsteady Aero Model Switch (switch) {{2=B-L Gonzalez, 3=B-L Minnema/Pierce, 4=B-L HGM 4-states, 5=B-L 5 states, 6=Oye, 7=Boeing-Vertol}} [used only when AFAeroMod=2]")
        _a(f"{_f(config.f_lookup):<14s}   FLookup            - Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]")

        # --- Airfoil Information ---
        _a("======  Airfoil Information =========================================================================")
        _a(f"{config.af_tab_mod:<14d}   AFTabMod           - Interpolation method for multiple airfoil tables {{1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp}} (-)")
        _a(f"{config.in_col_alfa:<14d}   InCol_Alfa         - The column in the airfoil tables that contains the angle of attack (-)")
        _a(f"{config.in_col_cl:<14d}   InCol_Cl           - The column in the airfoil tables that contains the lift coefficient (-)")
        _a(f"{config.in_col_cd:<14d}   InCol_Cd           - The column in the airfoil tables that contains the drag coefficient (-)")
        _a(f"{config.in_col_cm:<14d}   InCol_Cm           - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)")
        _a(f"{config.in_col_cpmin:<14d}   InCol_Cpmin        - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)")
        _a(f"{config.num_af_files:<14d}   NumAFfiles         - Number of airfoil files used (-)")
        for af_name in config.af_names:
            _a(f'"{ af_name }"')

        # --- Rotor/Blade Properties ---
        _a("======  Rotor/Blade Properties  =====================================================================")
        _a(f"{_f(config.use_bl_cm):<14s}   UseBlCm            - Include aerodynamic pitching moment in calculations?  (flag)")
        for i in range(config.num_bl):
            _a(f'{"\"" + config.ad_bl_file + "\"":<40s}   ADBlFile({i+1})      - Name of file containing distributed aerodynamic properties for Blade #{i+1} (-)')

        # --- Hub Properties ---
        _a("======  Hub Properties ============================================================================== [used only when Buoyancy=True]")
        _a(f"{config.vol_hub:<14.1f}   VolHub             - Hub volume (m^3)")
        _a(f"{config.hub_cen_bx:<14.1f}   HubCenBx           - Hub center of buoyancy x direction offset (m)")

        # --- Nacelle Properties ---
        _a("======  Nacelle Properties ========================================================================== [used only when Buoyancy=True]")
        _a(f"{config.vol_nac:<14.1f}   VolNac             - Nacelle volume (m^3)")
        nac_b_str = ", ".join(f"{v:.1f}" for v in config.nac_cen_b)
        _a(f"{nac_b_str:<40s}   NacCenB            - Position of nacelle center of buoyancy from yaw bearing in nacelle coordinates (m)")

        # --- Tower Influence ---
        _a("======  Tower Influence and Aerodynamics ============================================================ [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]")
        _a(f"{config.num_twr_nds:<14d}   NumTwrNds         - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]")
        if config.num_twr_nds > 0:
            _a("TwrElev        TwrDiam        TwrCd          TwrTI          TwrCb")
            _a("(m)            (m)            (-)            (-)            (-)")
            for ts in config.tower_stations:
                _a(f"{ts.twr_elev:14.7E}  {ts.twr_diam:14.7E}  {ts.twr_cd:14.7E}  {ts.twr_ti:14.7E}  {ts.twr_cb:14.1f}")

        # --- Outputs ---
        _a("======  Outputs  ====================================================================================")
        _a(f"{_f(config.sum_print):<14s}   SumPrint            - Generate a summary file listing input options and interpolated properties to \"<rootname>.AD.sum\"?  (flag)")
        _a(f"{config.n_bl_outs:<14d}   NBlOuts             - Number of blade node outputs [0 - 9] (-)")
        if config.bl_out_nd:
            bl_out_str = ", ".join(str(n) for n in config.bl_out_nd)
            _a(f"{bl_out_str:<40s}   BlOutNd             - Blade nodes whose values will be output  (-)")
        else:
            _a(f"{'1':<14s}   BlOutNd             - Blade nodes whose values will be output  (-)")
        _a(f"{config.n_tw_outs:<14d}   NTwOuts             - Number of tower node outputs [0 - 9]  (-)")
        if config.tw_out_nd:
            tw_out_str = ", ".join(str(n) for n in config.tw_out_nd)
            _a(f"{tw_out_str:<14s}   TwOutNd             - Tower nodes whose values will be output  (-)")
        else:
            _a(f"{'1':<14s}   TwOutNd             - Tower nodes whose values will be output  (-)")

        _a("                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)")
        for out in config.out_list:
            _a(out)
        _a('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)')
        _a("---------------------------------------------------------------------------------------")

        return "\n".join(lines) + "\n"

    def generate_blade_file(self, config: AeroDynBladeConfig) -> str:
        """Generate the AeroDyn blade definition input file.

        Parameters
        ----------
        config : AeroDynBladeConfig
            Blade aerodynamic station data.

        Returns
        -------
        str
            Complete AeroDyn blade definition file content.
        """
        lines: list[str] = []
        _a = lines.append

        _a("------- AERODYN v15.04.* BLADE DEFINITION INPUT FILE -------------------------------------")
        _a("Generated by WindForge - AeroDyn blade definition file")

        _a("======  Blade Properties =================================================================")
        _a(f"   {config.num_bl_nds:<15d}   NumBlNds           - Number of blade nodes used in the analysis (-)")
        _a("  BlSpn     BlCrvAC    BlSwpAC    BlCrvAng    BlTwist    BlChord    BlAFID")
        _a("  (m)       (m)        (m)        (deg)       (deg)      (m)        (-)")
        for s in config.stations:
            _a(
                f" {s.bl_spn:9.4f}"
                f"  {s.bl_crv_ac:9.4f}"
                f"  {s.bl_swp_ac:9.4f}"
                f"  {s.bl_crv_ang:10.4f}"
                f"  {s.bl_twist:10.4f}"
                f"  {s.bl_chord:10.4f}"
                f"  {s.bl_af_id:5d}"
            )

        return "\n".join(lines) + "\n"


def _flag(value: bool) -> str:
    """Convert a boolean to OpenFAST flag format."""
    return "True" if value else "False"
