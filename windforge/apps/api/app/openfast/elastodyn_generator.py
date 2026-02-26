"""
ElastoDyn input file generator for OpenFAST.

Generates three files:
  1. ElastoDyn primary input (.dat)
  2. ElastoDyn blade input (.dat)
  3. ElastoDyn tower input (.dat)

The primary file covers simulation control, DOFs, initial conditions,
turbine configuration, mass & inertia, drivetrain, and output lists.
Blade and tower files contain structural property tables and mode shapes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Data classes for blade / tower structural station tables
# ---------------------------------------------------------------------------

@dataclass
class BladeStation:
    """A single ElastoDyn blade station entry."""
    bl_fract: float       # Blade fractional radius (0..1)
    pitch_ax: float       # Pitch axis location (fraction of chord)
    strc_twist: float     # Structural twist (deg)
    b_mass_den: float     # Blade mass density (kg/m)
    flp_stff: float       # Flapwise stiffness EI_flap (N-m^2)
    edg_stff: float       # Edgewise stiffness EI_edge (N-m^2)


@dataclass
class TowerStation:
    """A single ElastoDyn tower station entry."""
    ht_fract: float       # Fractional height along tower (0..1)
    t_mass_den: float     # Tower mass density (kg/m)
    tw_fa_stif: float     # Tower fore-aft stiffness EI (N-m^2)
    tw_ss_stif: float     # Tower side-side stiffness EI (N-m^2)


# ---------------------------------------------------------------------------
# Configuration data classes
# ---------------------------------------------------------------------------

@dataclass
class ElastoDynBladeConfig:
    """Blade structural configuration for ElastoDyn blade file."""
    n_bl_inp_st: int = 19
    bld_flex_l: float = 61.5    # Blade flexible length (m)
    bld_flp_dmp: float = 0.477  # Blade flapwise damping (%)
    bld_edg_dmp: float = 0.477  # Blade edgewise damping (%)

    # Blade station table (BlFract, PitchAx, StrcTwst, BMassDen, FlpStff, EdgStff)
    stations: list[BladeStation] = field(default_factory=list)

    # Blade mode shape polynomial coefficients (6 terms: x^2..x^6)
    flp_mode_1: list[float] = field(
        default_factory=lambda: [0.0622, 1.7254, -3.2452, 4.7131, -2.2555]
    )
    flp_mode_2: list[float] = field(
        default_factory=lambda: [-0.5809, 1.2067, -15.5349, 29.7347, -13.8255]
    )
    edg_mode_1: list[float] = field(
        default_factory=lambda: [0.3627, 2.5337, -3.5772, 2.3760, -0.6952]
    )


@dataclass
class ElastoDynTowerConfig:
    """Tower structural configuration for ElastoDyn tower file."""
    n_tw_inp_st: int = 11
    twr_fa_dmp1: float = 1.0    # Tower 1st FA mode structural damping (%)
    twr_ss_dmp1: float = 1.0    # Tower 1st SS mode structural damping (%)
    twr_fa_dmp2: float = 1.0    # Tower 2nd FA mode structural damping (%)
    twr_ss_dmp2: float = 1.0    # Tower 2nd SS mode structural damping (%)

    # Tower station table (HtFract, TMassDen, TwFAStif, TwSSStif)
    stations: list[TowerStation] = field(default_factory=list)

    # Tower mode shape polynomial coefficients (5 terms: x^2..x^6)
    fa_mode_1: list[float] = field(
        default_factory=lambda: [0.7004, 2.1963, -5.6202, 6.2275, -2.5040]
    )
    fa_mode_2: list[float] = field(
        default_factory=lambda: [-26.0840, 73.8440, -78.5640, 34.1800, -2.3760]
    )
    ss_mode_1: list[float] = field(
        default_factory=lambda: [0.6360, 2.2124, -5.5836, 6.2433, -2.5081]
    )
    ss_mode_2: list[float] = field(
        default_factory=lambda: [-26.5340, 75.0040, -79.3660, 34.2560, -2.3600]
    )


@dataclass
class ElastoDynConfig:
    """Full ElastoDyn primary file configuration."""

    # --- Simulation Control ---
    echo: bool = False
    method: int = 3       # Integration method: 1=RK4, 2=AB4, 3=ABM4
    dt: str = "default"   # Integration time step (s) or "default"

    # --- Environmental Condition ---
    gravity: float = 9.80665

    # --- Degrees of Freedom ---
    flap_dof1: bool = True
    flap_dof2: bool = True
    edge_dof: bool = True
    drtr_dof: bool = True
    gen_dof: bool = True
    yaw_dof: bool = False
    tw_fa_dof1: bool = True
    tw_fa_dof2: bool = True
    tw_ss_dof1: bool = True
    tw_ss_dof2: bool = True
    ptfm_sg_dof: bool = False
    ptfm_sw_dof: bool = False
    ptfm_hv_dof: bool = False
    ptfm_r_dof: bool = False
    ptfm_p_dof: bool = False
    ptfm_y_dof: bool = False

    # --- Initial Conditions ---
    oopo_a: float = 0.0      # Out of plane blade 1 deflection (m)
    ipo_a: float = 0.0       # In plane blade 1 deflection (m)
    bl_pitch1: float = 0.0   # Blade 1 initial pitch (deg)
    bl_pitch2: float = 0.0   # Blade 2 initial pitch (deg)
    bl_pitch3: float = 0.0   # Blade 3 initial pitch (deg)
    teeter: float = 0.0      # Initial teeter angle (deg)
    azimuth: float = 0.0     # Initial azimuth (deg)
    rot_speed: float = 12.1  # Initial rotor speed (rpm)
    nac_yaw: float = 0.0     # Initial nacelle yaw (deg)
    ttdspfa: float = 0.0     # Initial fore-aft tower-top displacement (m)
    ttdspss: float = 0.0     # Initial side-side tower-top displacement (m)
    ptfm_surge: float = 0.0
    ptfm_sway: float = 0.0
    ptfm_heave: float = 0.0
    ptfm_roll: float = 0.0
    ptfm_pitch: float = 0.0
    ptfm_yaw: float = 0.0

    # --- Turbine Configuration ---
    num_bl: int = 3
    tip_rad: float = 63.0    # Tip radius (m)
    hub_rad: float = 1.5     # Hub radius (m)
    pre_cone_1: float = -2.5  # Blade 1 cone angle (deg)
    pre_cone_2: float = -2.5  # Blade 2 cone angle (deg)
    pre_cone_3: float = -2.5  # Blade 3 cone angle (deg)
    hub_cm: float = 0.0       # Hub CM location (m)
    undslng: float = 0.0      # Undersling (m)
    delta_3: float = 0.0      # Delta-3 angle for teetering (deg)
    azim_b1_up: float = 0.0   # Azimuth when blade 1 is Up (deg)
    overhang: float = -5.0191 # Rotor overhang (m) (negative downwind)
    shft_gag_l: float = 0.0   # Distance from hub to shaft strain gages (m)
    shft_tilt: float = -5.0   # Rotor shaft tilt angle (deg)
    nacelle_cm_xn: float = -1.0  # Nacelle CM fore-aft location (m)
    nacelle_cm_yn: float = 0.0   # Nacelle CM lateral location (m)
    nacelle_cm_zn: float = 1.0   # Nacelle CM vertical location (m)
    ncd_cmpn: float = 0.0     # Downwind distance from yaw bearing to nacelle CM (m)
    ncl_mpn: float = 0.0      # Lateral distance from yaw bearing to nacelle CM (m)
    ncu_mpn: float = 0.0      # Vertical distance from yaw bearing to nacelle CM (m)
    twr2shft: float = 1.96256  # Vertical distance from tower-top to rotor shaft (m)
    tower_ht: float = 87.6     # Tower height from ground/MSL (m)
    tower_bsht: float = 10.0   # Tower base height above ground/MSL (m)
    ptfm_cm_xt: float = 0.0
    ptfm_cm_yt: float = 0.0
    ptfm_cm_zt: float = -10.0
    ptfm_ref_zt: float = 0.0

    # --- Mass and Inertia ---
    tip_mass_1: float = 0.0   # Tip brake mass 1 (kg)
    tip_mass_2: float = 0.0
    tip_mass_3: float = 0.0
    hub_mass: float = 56780.0  # Hub mass (kg)
    hub_iner: float = 115926.0 # Hub inertia about rotor axis (kg-m^2)
    gen_iner: float = 534.116  # Generator inertia about HSS (kg-m^2)
    nac_mass: float = 240000.0 # Nacelle mass (kg)
    nac_yiner: float = 2607890.0  # Nacelle yaw inertia (kg-m^2)
    yaw_br_mass: float = 0.0  # Yaw bearing mass (kg)
    ptfm_mass: float = 0.0    # Platform mass (kg)
    ptfm_r_iner: float = 0.0  # Platform roll inertia (kg-m^2)
    ptfm_p_iner: float = 0.0  # Platform pitch inertia (kg-m^2)
    ptfm_y_iner: float = 0.0  # Platform yaw inertia (kg-m^2)

    # --- Drivetrain ---
    gb_ratio: float = 97.0    # Gearbox ratio
    gb_eff: float = 100.0     # Gearbox efficiency (%) -- NOT USED when using ServoDyn
    dt_tor_spr: float = 867637000.0  # Drivetrain torsional spring (N-m/rad)
    dt_tor_dmp: float = 6215000.0    # Drivetrain torsional damper (N-m/(rad/s))

    # --- Furling ---
    furling: bool = False

    # --- Tower File ---
    tower_file: str = "NRELOffs662hrBl_ElastoDyn_Tower.dat"

    # --- Output ---
    sum_print: bool = False
    out_file: int = 1
    tab_delim: bool = True
    out_fmt: str = "ES10.3E2"
    dec_fact: int = 1
    n_tw_gages: int = 0
    tw_gage_nds: list[int] = field(default_factory=list)
    n_bl_gages: int = 0
    bl_gage_nds: list[int] = field(default_factory=list)
    out_list: list[str] = field(default_factory=lambda: [
        '"BldPitch1"', '"BldPitch2"', '"BldPitch3"',
        '"Azimuth"', '"RotSpeed"', '"GenSpeed"',
        '"NacYaw"', '"OoPDefl1"', '"IPDefl1"', '"TTDspFA"', '"TTDspSS"',
        '"RootMxc1"', '"RootMyc1"', '"RootMzc1"',
        '"RootMxc2"', '"RootMyc2"', '"RootMzc2"',
        '"RootMxc3"', '"RootMyc3"', '"RootMzc3"',
        '"RotTorq"', '"LSSGagMya"', '"LSSGagMza"',
        '"YawBrFxp"', '"YawBrFyp"', '"YawBrFzp"',
        '"YawBrMxp"', '"YawBrMyp"', '"YawBrMzp"',
        '"TwrBsFxt"', '"TwrBsFyt"', '"TwrBsFzt"',
        '"TwrBsMxt"', '"TwrBsMyt"', '"TwrBsMzt"',
    ])

    # --- Blade files (referenced from primary ED file) ---
    blade_file: str = "NRELOffs662hrBl_ElastoDyn_Blade.dat"


class ElastoDynGenerator:
    """Generates ElastoDyn primary, blade, and tower input files."""

    def generate_main_file(self, config: ElastoDynConfig) -> str:
        """Generate the ElastoDyn primary input file.

        Parameters
        ----------
        config : ElastoDynConfig
            Full ElastoDyn configuration.

        Returns
        -------
        str
            Complete ElastoDyn primary input file content.
        """
        lines: list[str] = []
        _a = lines.append
        _f = _flag

        _a("------- ELASTODYN for OpenFAST INPUT FILE -------------------------------------------")
        _a("Generated by WindForge - ElastoDyn primary input file")

        # --- Simulation Control ---
        _a("---------------------- SIMULATION CONTROL --------------------------------------")
        _a(f"{_f(config.echo):<14s}   Echo            - Echo input data to <RootName>.ech (flag)")
        _a(f'{"\"" + config.method.__str__() + "\"":<14s}   Method          - Integration method: {{1: RK4, 2: AB4, or 3: ABM4}} (-)')
        _a(f'{"\"" + config.dt + "\"":<14s}   DT              - Integration time step (s)')

        # --- Environmental Condition ---
        _a("---------------------- ENVIRONMENTAL CONDITION ---------------------------------")
        _a(f"{config.gravity:<14.5f}   Gravity         - Gravitational acceleration (m/s^2)")

        # --- Degrees of Freedom ---
        _a("---------------------- DEGREES OF FREEDOM --------------------------------------")
        _a(f"{_f(config.flap_dof1):<14s}   FlapDOF1        - First flapwise blade mode DOF (flag)")
        _a(f"{_f(config.flap_dof2):<14s}   FlapDOF2        - Second flapwise blade mode DOF (flag)")
        _a(f"{_f(config.edge_dof):<14s}   EdgeDOF         - First edgewise blade mode DOF (flag)")
        _a(f"{_f(config.drtr_dof):<14s}   DrTrDOF         - Drivetrain rotational-flexibility DOF (flag)")
        _a(f"{_f(config.gen_dof):<14s}   GenDOF          - Generator DOF (flag)")
        _a(f"{_f(config.yaw_dof):<14s}   YawDOF          - Yaw DOF (flag)")
        _a(f"{_f(config.tw_fa_dof1):<14s}   TwFADOF1        - First fore-aft tower bending-mode DOF (flag)")
        _a(f"{_f(config.tw_fa_dof2):<14s}   TwFADOF2        - Second fore-aft tower bending-mode DOF (flag)")
        _a(f"{_f(config.tw_ss_dof1):<14s}   TwSSDOF1        - First side-to-side tower bending-mode DOF (flag)")
        _a(f"{_f(config.tw_ss_dof2):<14s}   TwSSDOF2        - Second side-to-side tower bending-mode DOF (flag)")
        _a(f"{_f(config.ptfm_sg_dof):<14s}   PtfmSgDOF       - Platform horizontal surge translation DOF (flag)")
        _a(f"{_f(config.ptfm_sw_dof):<14s}   PtfmSwDOF       - Platform horizontal sway translation DOF (flag)")
        _a(f"{_f(config.ptfm_hv_dof):<14s}   PtfmHvDOF       - Platform vertical heave translation DOF (flag)")
        _a(f"{_f(config.ptfm_r_dof):<14s}   PtfmRDOF        - Platform roll tilt rotation DOF (flag)")
        _a(f"{_f(config.ptfm_p_dof):<14s}   PtfmPDOF        - Platform pitch tilt rotation DOF (flag)")
        _a(f"{_f(config.ptfm_y_dof):<14s}   PtfmYDOF        - Platform yaw rotation DOF (flag)")

        # --- Initial Conditions ---
        _a("---------------------- INITIAL CONDITIONS --------------------------------------")
        _a(f"{config.oopo_a:<14.3f}   OoPDefl         - Initial out-of-plane blade-tip displacement (meters)")
        _a(f"{config.ipo_a:<14.3f}   IPDefl          - Initial in-plane blade-tip deflection (meters)")
        _a(f"{config.bl_pitch1:<14.3f}   BlPitch(1)      - Blade 1 initial pitch (degrees)")
        _a(f"{config.bl_pitch2:<14.3f}   BlPitch(2)      - Blade 2 initial pitch (degrees)")
        _a(f"{config.bl_pitch3:<14.3f}   BlPitch(3)      - Blade 3 initial pitch (degrees)")
        _a(f"{config.teeter:<14.3f}   Teeter          - Initial or fixed teeter angle (degrees) [unused for 3 blades]")
        _a(f"{config.azimuth:<14.3f}   Azimuth         - Initial azimuth angle for blade 1 (degrees)")
        _a(f"{config.rot_speed:<14.4f}   RotSpeed        - Initial or fixed rotor speed (rpm)")
        _a(f"{config.nac_yaw:<14.3f}   NacYaw          - Initial or fixed nacelle-yaw angle (degrees)")
        _a(f"{config.ttdspfa:<14.3f}   TTDspFA         - Initial fore-aft tower-top displacement (meters)")
        _a(f"{config.ttdspss:<14.3f}   TTDspSS         - Initial side-to-side tower-top displacement (meters)")
        _a(f"{config.ptfm_surge:<14.3f}   PtfmSurge       - Initial or fixed horizontal surge translational displacement of platform (meters)")
        _a(f"{config.ptfm_sway:<14.3f}   PtfmSway        - Initial or fixed horizontal sway translational displacement of platform (meters)")
        _a(f"{config.ptfm_heave:<14.3f}   PtfmHeave       - Initial or fixed vertical heave translational displacement of platform (meters)")
        _a(f"{config.ptfm_roll:<14.3f}   PtfmRoll        - Initial or fixed roll tilt rotational displacement of platform (degrees)")
        _a(f"{config.ptfm_pitch:<14.3f}   PtfmPitch       - Initial or fixed pitch tilt rotational displacement of platform (degrees)")
        _a(f"{config.ptfm_yaw:<14.3f}   PtfmYaw         - Initial or fixed yaw rotational displacement of platform (degrees)")

        # --- Turbine Configuration ---
        _a("---------------------- TURBINE CONFIGURATION -----------------------------------")
        _a(f"{config.num_bl:<14d}   NumBl           - Number of blades (-)")
        _a(f"{config.tip_rad:<14.4f}   TipRad          - The distance from the rotor apex to the blade tip (meters)")
        _a(f"{config.hub_rad:<14.4f}   HubRad          - The distance from the rotor apex to the blade root (meters)")
        _a(f"{config.pre_cone_1:<14.4f}   PreCone(1)      - Blade 1 cone angle (degrees)")
        _a(f"{config.pre_cone_2:<14.4f}   PreCone(2)      - Blade 2 cone angle (degrees)")
        _a(f"{config.pre_cone_3:<14.4f}   PreCone(3)      - Blade 3 cone angle (degrees)")
        _a(f"{config.hub_cm:<14.4f}   HubCM           - Distance from rotor apex to hub mass [positive downwind] (meters)")
        _a(f"{config.undslng:<14.4f}   UndSling        - Undersling length [distance from yaw axis to teeter pin or additional rotor mass] (meters)")
        _a(f"{config.delta_3:<14.4f}   Delta3          - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]")
        _a(f"{config.azim_b1_up:<14.4f}   AzimB1Up        - Azimuth value to use for I/O when blade 1 points up (degrees)")
        _a(f"{config.overhang:<14.4f}   OverHang        - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)")
        _a(f"{config.shft_gag_l:<14.4f}   ShftGagL        - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gage/gage [positive downwind] (meters)")
        _a(f"{config.shft_tilt:<14.4f}   ShftTilt        - Rotor shaft tilt angle (degrees)")
        _a(f"{config.nacelle_cm_xn:<14.4f}   NacCMxn         - Downwind distance from the tower-top to the nacelle CM (meters)")
        _a(f"{config.nacelle_cm_yn:<14.4f}   NacCMyn         - Lateral  distance from the tower-top to the nacelle CM (meters)")
        _a(f"{config.nacelle_cm_zn:<14.4f}   NacCMzn         - Vertical distance from the tower-top to the nacelle CM (meters)")
        _a(f"{config.ncd_cmpn:<14.4f}   NcIMUxn         - Downwind distance from the tower-top to the nacelle IMU (meters)")
        _a(f"{config.ncl_mpn:<14.4f}   NcIMUyn         - Lateral  distance from the tower-top to the nacelle IMU (meters)")
        _a(f"{config.ncu_mpn:<14.4f}   NcIMUzn         - Vertical distance from the tower-top to the nacelle IMU (meters)")
        _a(f"{config.twr2shft:<14.5f}   Twr2Shft        - Vertical distance from the tower-top to the rotor shaft (meters)")
        _a(f"{config.tower_ht:<14.4f}   TowerHt         - Height of tower above ground level [onshore] or MSL [offshore] (meters)")
        _a(f"{config.tower_bsht:<14.4f}   TowerBsHt       - Height of tower base above ground level [onshore] or MSL [offshore] (meters)")
        _a(f"{config.ptfm_cm_xt:<14.4f}   PtfmCMxt        - Downwind distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)")
        _a(f"{config.ptfm_cm_yt:<14.4f}   PtfmCMyt        - Lateral distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)")
        _a(f"{config.ptfm_cm_zt:<14.4f}   PtfmCMzt        - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)")
        _a(f"{config.ptfm_ref_zt:<14.4f}   PtfmRefzt       - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform reference point (meters)")

        # --- Mass and Inertia ---
        _a("---------------------- MASS AND INERTIA ----------------------------------------")
        _a(f"{config.tip_mass_1:<14.1f}   TipMass(1)      - Tip-brake mass, blade 1 (kg)")
        _a(f"{config.tip_mass_2:<14.1f}   TipMass(2)      - Tip-brake mass, blade 2 (kg)")
        _a(f"{config.tip_mass_3:<14.1f}   TipMass(3)      - Tip-brake mass, blade 3 (kg)")
        _a(f"{config.hub_mass:<14.1f}   HubMass         - Hub mass (kg)")
        _a(f"{config.hub_iner:<14.1f}   HubIner         - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)")
        _a(f"{config.gen_iner:<14.3f}   GenIner         - Generator inertia about HSS (kg m^2)")
        _a(f"{config.nac_mass:<14.1f}   NacMass         - Nacelle mass (kg)")
        _a(f"{config.nac_yiner:<14.1f}   NacYIner        - Nacelle inertia about yaw axis (kg m^2)")
        _a(f"{config.yaw_br_mass:<14.1f}   YawBrMass       - Yaw bearing mass (kg)")
        _a(f"{config.ptfm_mass:<14.1f}   PtfmMass        - Platform mass (kg)")
        _a(f"{config.ptfm_r_iner:<14.1f}   PtfmRIner       - Platform inertia for roll tilt rotation about the platform CM (kg m^2)")
        _a(f"{config.ptfm_p_iner:<14.1f}   PtfmPIner       - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)")
        _a(f"{config.ptfm_y_iner:<14.1f}   PtfmYIner       - Platform inertia for yaw rotation about the platform CM (kg m^2)")

        # --- Blade ---
        _a("---------------------- BLADE ---------------------------------------------------")
        _a(f"{config.num_bl:<14d}   BldNodes        - Number of blade nodes (per blade) used for analysis (-)")
        _a(f'{"\"" + config.blade_file + "\"":<40s}   BldFile(1)      - Name of file containing properties for blade 1 (quoted string)')
        _a(f'{"\"" + config.blade_file + "\"":<40s}   BldFile(2)      - Name of file containing properties for blade 2 (quoted string)')
        _a(f'{"\"" + config.blade_file + "\"":<40s}   BldFile(3)      - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]')

        # --- Rotor-Teeter ---
        _a("---------------------- ROTOR-TEETER --------------------------------------------")
        _a(f"{'0':<14s}   TeetMod         - Rotor-teeter spring/damper model {{0: none, 1: standard, 2: user-defined from routine UserTeet}} (switch) [unused for 3 blades]")
        _a(f"{'0.0':<14s}   TeetDmpP        - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetDmp         - Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetCDmp        - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetSStP        - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetHStP        - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetSSSp        - Rotor-teeter soft-stop linear-loss constant (N-m/rad) [used only for 2 blades and when TeetMod=1]")
        _a(f"{'0.0':<14s}   TeetHSSp        - Rotor-teeter hard-stop linear-loss constant (N-m/rad) [used only for 2 blades and when TeetMod=1]")

        # --- Drivetrain ---
        _a("---------------------- DRIVETRAIN ----------------------------------------------")
        _a(f"{config.gb_eff:<14.4f}   GBoxEff         - Gearbox efficiency (%%)")
        _a(f"{config.gb_ratio:<14.4f}   GBRatio         - Gearbox ratio (-)")
        _a(f"{config.dt_tor_spr:<14.1f}   DTTorSpr        - Drivetrain torsional spring (N-m/rad)")
        _a(f"{config.dt_tor_dmp:<14.1f}   DTTorDmp        - Drivetrain torsional damper (N-m/(rad/s))")

        # --- Furling ---
        _a("---------------------- FURLING -------------------------------------------------")
        _a(f"{_f(config.furling):<14s}   Furling         - Read in additional model properties for furling turbine (flag) [must currently be FALSE]")

        # --- Tower ---
        _a("---------------------- TOWER ---------------------------------------------------")
        _a(f"{20:<14d}   TwrNodes        - Number of tower nodes used for analysis (-)")
        _a(f'{"\"" + config.tower_file + "\"":<40s}   TwrFile         - Name of file containing tower properties (quoted string)')

        # --- Output ---
        _a("---------------------- OUTPUT --------------------------------------------------")
        _a(f"{_f(config.sum_print):<14s}   SumPrint        - Print summary data to <RootName>.sum (flag)")
        _a(f"{config.out_file:<14d}   OutFile         - Switch to determine where output will be placed: {{1: in module output file only; 2: in glue code output file only; 3: both}} (currently unused)")
        _a(f"{_f(config.tab_delim):<14s}   TabDelim        - Use tab delimiters in text tabular output file? (flag) (currently unused)")
        _a(f'{"\"" + config.out_fmt + "\"":<14s}   OutFmt          - Format used for text tabular output, excluding the time channel (quoted string) (currently unused)')
        _a(f"{config.dec_fact:<14d}   DecFact         - Decimation factor for tabular output {{1: output every time step}} (-) (currently unused)")
        _a(f"{config.n_tw_gages:<14d}   NTwGages        - Number of tower nodes that have strain gages for output [0 to 9] (-)")
        if config.tw_gage_nds:
            _a(f"{', '.join(str(n) for n in config.tw_gage_nds):<40s}   TwGagNd         - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]")
        else:
            _a(f"{'0':<40s}   TwGagNd         - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]")
        _a(f"{config.n_bl_gages:<14d}   NBlGages        - Number of blade nodes that have strain gages for output [0 to 9] (-)")
        if config.bl_gage_nds:
            _a(f"{', '.join(str(n) for n in config.bl_gage_nds):<40s}   BlGagNd         - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]")
        else:
            _a(f"{'0':<40s}   BlGagNd         - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]")

        # --- OutList ---
        _a("                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)")
        for out in config.out_list:
            _a(out)
        _a('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)')
        _a("---------------------------------------------------------------------------------------")

        return "\n".join(lines) + "\n"

    def generate_blade_file(self, config: ElastoDynBladeConfig) -> str:
        """Generate the ElastoDyn blade input file.

        Parameters
        ----------
        config : ElastoDynBladeConfig
            Blade structural properties and mode shapes.

        Returns
        -------
        str
            Complete ElastoDyn blade input file content.
        """
        lines: list[str] = []
        _a = lines.append

        _a("------- ELASTODYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------")
        _a("Generated by WindForge - ElastoDyn blade input file")

        # --- Blade Parameters ---
        _a("---------------------- BLADE PARAMETERS ----------------------------------------")
        _a(f"{config.n_bl_inp_st:<14d}   NBlInpSt        - Number of blade input stations (-)")
        _a(f"{config.bld_flex_l:<14.4f}   BldFlexL        - Blade flexible length (meters)")
        _a(f"{config.bld_flp_dmp:<14.4f}   BldFlpDmp       - Blade flapwise structural damping ratio (%%)")
        _a(f"{config.bld_edg_dmp:<14.4f}   BldEdgDmp       - Blade edgewise structural damping ratio (%%)")

        # --- Blade Adjustment Factors ---
        _a("---------------------- BLADE ADJUSTMENT FACTORS --------------------------------")
        _a(f"{'1.0':<14s}   FlStTunr(1)     - Blade flapwise modal stiffness tuner, 1st mode (-)")
        _a(f"{'1.0':<14s}   FlStTunr(2)     - Blade flapwise modal stiffness tuner, 2nd mode (-)")
        _a(f"{'1.0':<14s}   AdjBlMs         - Factor to adjust blade mass density (-)")
        _a(f"{'1.0':<14s}   AdjFlSt         - Factor to adjust blade flap stiffness (-)")
        _a(f"{'1.0':<14s}   AdjEdSt         - Factor to adjust blade edge stiffness (-)")

        # --- Distributed Blade Properties ---
        _a("---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------")
        _a("    BlFract      PitchAx      StrcTwst       BMassDen        FlpStff        EdgStff")
        _a("      (-)          (-)          (deg)         (kg/m)         (Nm^2)         (Nm^2)")
        for s in config.stations:
            _a(f"  {s.bl_fract:10.6f}  {s.pitch_ax:10.6f}  {s.strc_twist:12.4f}  {s.b_mass_den:13.4f}  {s.flp_stff:14.4E}  {s.edg_stff:14.4E}")

        # --- Blade Mode Shapes ---
        _a("---------------------- BLADE MODE SHAPES ---------------------------------------")
        _a(f"{config.flp_mode_1[0]:<14.4f}   BldFl1Sh(2) - Flap mode 1, coeff of x^2")
        _a(f"{config.flp_mode_1[1]:<14.4f}   BldFl1Sh(3) -            , coeff of x^3")
        _a(f"{config.flp_mode_1[2]:<14.4f}   BldFl1Sh(4) -            , coeff of x^4")
        _a(f"{config.flp_mode_1[3]:<14.4f}   BldFl1Sh(5) -            , coeff of x^5")
        _a(f"{config.flp_mode_1[4]:<14.4f}   BldFl1Sh(6) -            , coeff of x^6")

        _a(f"{config.flp_mode_2[0]:<14.4f}   BldFl2Sh(2) - Flap mode 2, coeff of x^2")
        _a(f"{config.flp_mode_2[1]:<14.4f}   BldFl2Sh(3) -            , coeff of x^3")
        _a(f"{config.flp_mode_2[2]:<14.4f}   BldFl2Sh(4) -            , coeff of x^4")
        _a(f"{config.flp_mode_2[3]:<14.4f}   BldFl2Sh(5) -            , coeff of x^5")
        _a(f"{config.flp_mode_2[4]:<14.4f}   BldFl2Sh(6) -            , coeff of x^6")

        _a(f"{config.edg_mode_1[0]:<14.4f}   BldEdgSh(2) - Edge mode 1, coeff of x^2")
        _a(f"{config.edg_mode_1[1]:<14.4f}   BldEdgSh(3) -            , coeff of x^3")
        _a(f"{config.edg_mode_1[2]:<14.4f}   BldEdgSh(4) -            , coeff of x^4")
        _a(f"{config.edg_mode_1[3]:<14.4f}   BldEdgSh(5) -            , coeff of x^5")
        _a(f"{config.edg_mode_1[4]:<14.4f}   BldEdgSh(6) -            , coeff of x^6")

        return "\n".join(lines) + "\n"

    def generate_tower_file(self, config: ElastoDynTowerConfig) -> str:
        """Generate the ElastoDyn tower input file.

        Parameters
        ----------
        config : ElastoDynTowerConfig
            Tower structural properties and mode shapes.

        Returns
        -------
        str
            Complete ElastoDyn tower input file content.
        """
        lines: list[str] = []
        _a = lines.append

        _a("------- ELASTODYN V1.00.* TOWER INPUT FILE -------------------------------------")
        _a("Generated by WindForge - ElastoDyn tower input file")

        # --- Tower Parameters ---
        _a("---------------------- TOWER PARAMETERS ----------------------------------------")
        _a(f"{config.n_tw_inp_st:<14d}   NTwInpSt        - Number of input stations to specify tower geometry")
        _a(f"{config.twr_fa_dmp1:<14.4f}   TwrFADmp(1)     - Tower 1st fore-aft mode structural damping ratio (%%)")
        _a(f"{config.twr_ss_dmp1:<14.4f}   TwrSSDmp(1)     - Tower 1st side-to-side mode structural damping ratio (%%)")
        _a(f"{config.twr_fa_dmp2:<14.4f}   TwrFADmp(2)     - Tower 2nd fore-aft mode structural damping ratio (%%)")
        _a(f"{config.twr_ss_dmp2:<14.4f}   TwrSSDmp(2)     - Tower 2nd side-to-side mode structural damping ratio (%%)")

        # --- Tower Adjustment Factors ---
        _a("---------------------- TOWER ADJUSTMNET FACTORS --------------------------------")
        _a(f"{'1.0':<14s}   FAStTunr(1)     - Tower fore-aft modal stiffness tuner, 1st mode (-)")
        _a(f"{'1.0':<14s}   FAStTunr(2)     - Tower fore-aft modal stiffness tuner, 2nd mode (-)")
        _a(f"{'1.0':<14s}   SSStTunr(1)     - Tower side-to-side stiffness tuner, 1st mode (-)")
        _a(f"{'1.0':<14s}   SSStTunr(2)     - Tower side-to-side stiffness tuner, 2nd mode (-)")
        _a(f"{'1.0':<14s}   AdjTwMa         - Factor to adjust tower mass density (-)")
        _a(f"{'1.0':<14s}   AdjFASt         - Factor to adjust tower fore-aft stiffness (-)")
        _a(f"{'1.0':<14s}   AdjSSSt         - Factor to adjust tower side-to-side stiffness (-)")

        # --- Distributed Tower Properties ---
        _a("---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------")
        _a("    HtFract       TMassDen         TwFAStif         TwSSStif")
        _a("      (-)          (kg/m)           (Nm^2)           (Nm^2)")
        for s in config.stations:
            _a(f"  {s.ht_fract:10.6f}  {s.t_mass_den:13.4f}  {s.tw_fa_stif:16.4E}  {s.tw_ss_stif:16.4E}")

        # --- Tower Mode Shapes ---
        _a("---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------")
        _a(f"{config.fa_mode_1[0]:<14.4f}   TwFAM1Sh(2) - Mode 1, coeff of x^2")
        _a(f"{config.fa_mode_1[1]:<14.4f}   TwFAM1Sh(3) -         coeff of x^3")
        _a(f"{config.fa_mode_1[2]:<14.4f}   TwFAM1Sh(4) -         coeff of x^4")
        _a(f"{config.fa_mode_1[3]:<14.4f}   TwFAM1Sh(5) -         coeff of x^5")
        _a(f"{config.fa_mode_1[4]:<14.4f}   TwFAM1Sh(6) -         coeff of x^6")

        _a(f"{config.fa_mode_2[0]:<14.4f}   TwFAM2Sh(2) - Mode 2, coeff of x^2")
        _a(f"{config.fa_mode_2[1]:<14.4f}   TwFAM2Sh(3) -         coeff of x^3")
        _a(f"{config.fa_mode_2[2]:<14.4f}   TwFAM2Sh(4) -         coeff of x^4")
        _a(f"{config.fa_mode_2[3]:<14.4f}   TwFAM2Sh(5) -         coeff of x^5")
        _a(f"{config.fa_mode_2[4]:<14.4f}   TwFAM2Sh(6) -         coeff of x^6")

        _a("---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------")
        _a(f"{config.ss_mode_1[0]:<14.4f}   TwSSM1Sh(2) - Mode 1, coeff of x^2")
        _a(f"{config.ss_mode_1[1]:<14.4f}   TwSSM1Sh(3) -         coeff of x^3")
        _a(f"{config.ss_mode_1[2]:<14.4f}   TwSSM1Sh(4) -         coeff of x^4")
        _a(f"{config.ss_mode_1[3]:<14.4f}   TwSSM1Sh(5) -         coeff of x^5")
        _a(f"{config.ss_mode_1[4]:<14.4f}   TwSSM1Sh(6) -         coeff of x^6")

        _a(f"{config.ss_mode_2[0]:<14.4f}   TwSSM2Sh(2) - Mode 2, coeff of x^2")
        _a(f"{config.ss_mode_2[1]:<14.4f}   TwSSM2Sh(3) -         coeff of x^3")
        _a(f"{config.ss_mode_2[2]:<14.4f}   TwSSM2Sh(4) -         coeff of x^4")
        _a(f"{config.ss_mode_2[3]:<14.4f}   TwSSM2Sh(5) -         coeff of x^5")
        _a(f"{config.ss_mode_2[4]:<14.4f}   TwSSM2Sh(6) -         coeff of x^6")

        return "\n".join(lines) + "\n"


def _flag(value: bool) -> str:
    """Convert a boolean to OpenFAST flag format."""
    return "True" if value else "False"
