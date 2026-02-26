"""
ServoDyn input file generator for OpenFAST.

Generates two files:
  1. ServoDyn primary input file (.dat)
  2. DISCON.IN controller parameter file (for ROSCO or Bladed-style DLLs)

Format verified against docs/source/user/servodyn/input.rst.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ServoDynConfig:
    """Configuration for the ServoDyn primary input file."""

    # --- Simulation Control ---
    echo: bool = False
    dt: str = "default"

    # --- Pitch Control ---
    pc_mode: int = 5           # 0=none, 3=user routine, 4=Simulink, 5=Bladed DLL
    tpc_on: float = 0.0
    tpit_man_s_1: float = 9999.9
    tpit_man_s_2: float = 9999.9
    tpit_man_s_3: float = 9999.9
    pit_man_rat_1: float = 2.0
    pit_man_rat_2: float = 2.0
    pit_man_rat_3: float = 2.0
    bl_pitch_f_1: float = 0.0
    bl_pitch_f_2: float = 0.0
    bl_pitch_f_3: float = 0.0

    # --- Generator and Torque Control ---
    vs_contrl: int = 5         # 0=none, 1=simple VS, 3=user routine, 4=Simulink, 5=Bladed DLL
    gen_model: int = 1
    gen_eff: float = 94.4
    gen_ti_str: bool = True
    gen_ti_stp: bool = True
    spd_gen_on: float = 9999.9
    tim_gen_on: float = 0.0
    tim_gen_of: float = 9999.9

    # --- Simple Variable-Speed Torque Control ---
    vs_rt_gn_sp: float = 9999.9
    vs_rt_tq: float = 9999.9
    vs_rgn2_k: float = 0.0
    vs_sl_pc: float = 10.0

    # --- Simple Induction Generator ---
    sig_sl_pc: float = 1.0
    sig_sy_sp: float = 1500.0
    sig_rt_tq: float = 9999.9
    sig_port: float = 2.0

    # --- Thevenin-Equivalent Induction Generator ---
    tec_freq: float = 60.0
    tec_npol: int = 2
    tec_sres: float = 0.0
    tec_rres: float = 0.0
    tec_vll: float = 0.0
    tec_slr: float = 0.0
    tec_rlr: float = 0.0
    tec_mr: float = 0.0

    # --- HSS Brake ---
    hss_br_mode: int = 0      # 0=none, 1=simple, 3=user routine, 5=Bladed DLL
    thss_br_dp: float = 9999.9
    hss_br_dt: float = 9999.9
    hss_br_tqf: float = 9999.9

    # --- Nacelle-Yaw Control ---
    yc_mode: int = 0           # 0=none, 3=user routine, 4=Simulink, 5=Bladed DLL
    tyc_on: float = 0.0
    yaw_neut: float = 0.0
    yaw_spr: float = 0.0
    yaw_damp: float = 0.0
    tyaw_man_s: float = 9999.9
    yaw_man_rat: float = 0.25
    nac_yaw_f: float = 0.0

    # --- Aerodynamic Flow Control ---
    afc_mode: int = 0
    afc_mean: float = 0.0
    afc_amp: float = 0.0
    afc_phase: float = 0.0

    # --- Cable Control ---
    cc_mode: int = 0

    # --- Structural Control ---
    num_b_stc: int = 0
    b_stc_files: list[str] = field(default_factory=list)
    num_n_stc: int = 0
    n_stc_files: list[str] = field(default_factory=list)
    num_t_stc: int = 0
    t_stc_files: list[str] = field(default_factory=list)
    num_s_stc: int = 0
    s_stc_files: list[str] = field(default_factory=list)

    # --- Bladed Interface ---
    dll_file_name: str = "libdiscon.so"
    dll_in_file: str = "DISCON.IN"
    dll_proc_name: str = "DISCON"
    dll_dt: str = "default"
    dll_ramp: bool = False
    bp_cutoff: float = 9999.9
    nac_yaw_north: float = 0.0
    ptch_cntrl: int = 1       # 0=collective, 1=individual
    ptch_set_pnt: float = 0.0
    ptch_min: float = 0.0
    ptch_max: float = 90.0
    ptch_rate_min: float = -8.0
    ptch_rate_max: float = 8.0
    gain_om: float = 0.0
    gen_spd_min_om: float = 0.0
    gen_spd_max_om: float = 9999.9
    gen_spd_dem: float = 9999.9
    gen_trq_dem: float = 9999.9
    gen_pwr_dem: float = 9999.9

    # --- Bladed Interface Torque-Speed LUT ---
    dll_num_trq: int = 0

    # --- Output ---
    sum_print: bool = False
    out_file: int = 1
    tab_delim: bool = True
    out_fmt: str = "ES10.3E2"
    t_start: float = 0.0
    out_list: list[str] = field(default_factory=lambda: [
        '"GenPwr"', '"GenTq"', '"BldPitch1"',
    ])


@dataclass
class DISCONConfig:
    """Configuration for a ROSCO-style DISCON.IN controller parameter file."""

    # --- Logging ---
    log_file: int = 1          # 0=no log, 1=write to file, 2=write to screen
    log_level: int = 1         # 0=minimal, 1=standard, 2=verbose

    # --- Controller Flags ---
    f_lp_corner_freq: float = 1.5    # Corner frequency for LP filters (rad/s)
    f_fl_hp_corner_freq: float = 0.0  # Corner frequency for HP filter on floating feedback (rad/s)
    f_notify_filter: int = 0

    # --- Filter Parameters ---
    f_lp_type: int = 1         # 1=first-order, 2=second-order
    f_lp_damping: float = 0.7

    # --- Pitch Control ---
    pc_control_mode: int = 1   # 0=power, 1=speed reference
    pc_gain_schedule_mode: int = 1
    pc_min_pitch: float = 0.0
    pc_max_pitch: float = 1.5708  # 90 deg in rad
    pc_min_rat: float = -0.1396   # Min pitch rate (rad/s)
    pc_max_rat: float = 0.1396    # Max pitch rate (rad/s)
    pc_ref_speed: float = 1.2671  # Rated rotor speed (rad/s)
    pc_fine_pitch: float = 0.0
    pc_switch: float = 0.0

    # --- Pitch Controller Gains ---
    pc_gs_n: int = 5
    pc_gs_angles: list[float] = field(default_factory=lambda: [0.0588, 0.0967, 0.1283, 0.1573, 0.1856])
    pc_gs_kp: list[float] = field(default_factory=lambda: [-0.0159, -0.0183, -0.0203, -0.0222, -0.024])
    pc_gs_ki: list[float] = field(default_factory=lambda: [-0.0068, -0.0078, -0.0087, -0.0095, -0.0103])

    # --- Torque Control ---
    vs_control_mode: int = 2   # 0=power ref, 1=simple VS, 2=tip-speed-ratio tracking
    vs_gen_eff: float = 0.944
    vs_arith_bi_num: float = 0.0
    vs_rated_gen_pwr: float = 5000000.0
    vs_max_rat: float = 15000.0
    vs_max_tq: float = 47402.91
    vs_min_tq: float = 0.0
    vs_min_om_spd: float = 0.7854  # Min generator speed (rad/s)
    vs_rated_om_spd: float = 1.2671 # Rated generator speed (rad/s)
    vs_rgn2_k: float = 2.3323      # Region 2 gain (Nm/(rad/s)^2)
    vs_ref_spd: float = 1.2671
    vs_tsr: float = 7.55           # Optimal tip-speed ratio

    # --- Setpoint Smoother ---
    ss_vsgain: float = 1.0
    ss_pcgain: float = 0.001

    # --- Wind Speed Estimator ---
    we_blade_radius: float = 63.0
    we_cp_n: int = 1
    we_cp: list[float] = field(default_factory=lambda: [0.0])
    we_gamma: float = 0.0
    we_gear_ratio: float = 97.0
    we_rlspd: float = 1.2671
    we_rated_pwr: float = 5000000.0
    we_min_pitch: float = 0.0
    we_rated_v: float = 11.4
    we_rated_omega: float = 1.2671

    # --- Floating Feedback ---
    fl_mode: int = 0

    # --- Yaw Control ---
    y_control_mode: int = 0
    y_me_max_rate: float = 0.0

    # --- Tower Fore-aft Damping ---
    fa_k_constant: float = 0.0

    # --- Minimum Pitch Schedule ---
    ps_bld_pitch_min: float = 0.0

    # --- Shutdown ---
    sd_max_pitch: float = 1.5708
    sd_time_limit: float = 9999.9


class ServoDynGenerator:
    """Generates ServoDyn primary and DISCON.IN controller input files."""

    def generate_servodyn_file(self, config: ServoDynConfig) -> str:
        """Generate the ServoDyn primary input file.

        Parameters
        ----------
        config : ServoDynConfig
            Complete ServoDyn configuration.

        Returns
        -------
        str
            Complete ServoDyn primary input file content.
        """
        lines: list[str] = []
        _a = lines.append
        _f = _flag

        _a("------- SERVODYN v1.05.* INPUT FILE ------------------------------------------------")
        _a("Generated by WindForge - ServoDyn primary input file")

        # --- Simulation Control ---
        _a("---------------------- SIMULATION CONTROL -----------------------------------------")
        _a(f"{_f(config.echo):<14s}   Echo         - Echo input data to <RootName>.ech (flag)")
        _a(f'{"\"" + config.dt + "\"":<14s}   DT           - Communication interval for controllers (s) (or \"default\")')

        # --- Pitch Control ---
        _a("---------------------- PITCH CONTROL -----------------------------------------------")
        _a(f"{config.pc_mode:<14d}   PCMode       - Pitch control mode {{0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")
        _a(f"{config.tpc_on:<14.1f}   TPCOn        - Time to enable active pitch control [unused when PCMode = 0] (s)")
        _a(f"{config.tpit_man_s_1:<14.1f}   TPitManS(1)  - Time to start override pitch maneuver for blade 1 and target pitch angle for pitch controller setpoints (s)")
        _a(f"{config.tpit_man_s_2:<14.1f}   TPitManS(2)  - Time to start override pitch maneuver for blade 2 and target pitch angle for pitch controller setpoints (s)")
        _a(f"{config.tpit_man_s_3:<14.1f}   TPitManS(3)  - Time to start override pitch maneuver for blade 3 and target pitch angle for pitch controller setpoints (s) [unused for 2 blades]")
        _a(f"{config.pit_man_rat_1:<14.4f}   PitManRat(1) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 1 (deg/s)")
        _a(f"{config.pit_man_rat_2:<14.4f}   PitManRat(2) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 2 (deg/s)")
        _a(f"{config.pit_man_rat_3:<14.4f}   PitManRat(3) - Pitch rate at which override pitch maneuver heads toward final pitch angle for blade 3 (deg/s) [unused for 2 blades]")
        _a(f"{config.bl_pitch_f_1:<14.4f}   BlPitchF(1)  - Blade 1 final pitch for override pitch maneuvers (degrees)")
        _a(f"{config.bl_pitch_f_2:<14.4f}   BlPitchF(2)  - Blade 2 final pitch for override pitch maneuvers (degrees)")
        _a(f"{config.bl_pitch_f_3:<14.4f}   BlPitchF(3)  - Blade 3 final pitch for override pitch maneuvers (degrees) [unused for 2 blades]")

        # --- Generator and Torque Control ---
        _a("---------------------- GENERATOR AND TORQUE CONTROL --------------------------------")
        _a(f"{config.vs_contrl:<14d}   VSContrl     - Variable-speed control mode {{0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")
        _a(f"{config.gen_model:<14d}   GenModel     - Generator model {{1: simple, 2: Thevenin, 3: user-defined from routine UserGen}} (switch) [used only when VSContrl=0]")
        _a(f"{config.gen_eff:<14.2f}   GenEff       - Generator efficiency [ignored by the Thevenin and user-defined generator models] (%%)")
        _a(f"{_f(config.gen_ti_str):<14s}   GenTiStr     - Method to start the generator {{T: timed using TimGenOn, F: generator speed using SpdGenOn}} (flag)")
        _a(f"{_f(config.gen_ti_stp):<14s}   GenTiStp     - Method to stop the generator {{T: timed using TimGenOf, F: when generator power = 0}} (flag)")
        _a(f"{config.spd_gen_on:<14.1f}   SpdGenOn     - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]")
        _a(f"{config.tim_gen_on:<14.1f}   TimGenOn     - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]")
        _a(f"{config.tim_gen_of:<14.1f}   TimGenOf     - Time to turn off the generator (s) [used only when GenTiStp=True]")

        # --- Simple Variable-Speed Torque Control ---
        _a("---------------------- SIMPLE VARIABLE-SPEED TORQUE CONTROL ------------------------")
        _a(f"{config.vs_rt_gn_sp:<14.1f}   VS_RtGnSp    - Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]")
        _a(f"{config.vs_rt_tq:<14.1f}   VS_RtTq      - Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]")
        _a(f"{config.vs_rgn2_k:<14.6f}   VS_Rgn2K     - Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]")
        _a(f"{config.vs_sl_pc:<14.2f}   VS_SlPc      - Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%%)[used only when VSContrl=1]")

        # --- Simple Induction Generator ---
        _a("---------------------- SIMPLE INDUCTION GENERATOR ----------------------------------")
        _a(f"{config.sig_sl_pc:<14.2f}   SIG_SlPc     - Rated generator slip percentage (%%)[used only when VSContrl=0 and GenModel=1]")
        _a(f"{config.sig_sy_sp:<14.1f}   SIG_SySp     - Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]")
        _a(f"{config.sig_rt_tq:<14.1f}   SIG_RtTq     - Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]")
        _a(f"{config.sig_port:<14.2f}   SIG_PORt     - Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]")

        # --- Thevenin-Equivalent Induction Generator ---
        _a("---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR ---------------------")
        _a(f"{config.tec_freq:<14.1f}   TEC_Freq     - Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_npol:<14d}   TEC_NPol     - Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_sres:<14.4f}   TEC_SRes     - Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_rres:<14.4f}   TEC_RRes     - Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_vll:<14.4f}   TEC_VLL      - Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_slr:<14.4f}   TEC_SLR      - Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_rlr:<14.4f}   TEC_RLR      - Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]")
        _a(f"{config.tec_mr:<14.4f}   TEC_MR       - Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]")

        # --- HSS Brake ---
        _a("---------------------- HIGH-SPEED SHAFT BRAKE --------------------------------------")
        _a(f"{config.hss_br_mode:<14d}   HSSBrMode    - HSS brake model {{0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")
        _a(f"{config.thss_br_dp:<14.1f}   THSSBrDp     - Time to initiate deployment of the HSS brake (s)")
        _a(f"{config.hss_br_dt:<14.1f}   HSSBrDT      - Time for HSS-brake to reach full deployment once initiated (used only when HSSBrMode=1) (s)")
        _a(f"{config.hss_br_tqf:<14.1f}   HSSBrTqF     - Fully deployed HSS-brake torque (N-m)")

        # --- Nacelle-Yaw Control ---
        _a("---------------------- NACELLE-YAW CONTROL -----------------------------------------")
        _a(f"{config.yc_mode:<14d}   YCMode       - Yaw control mode {{0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")
        _a(f"{config.tyc_on:<14.1f}   TYCOn        - Time to enable active yaw control [unused when YCMode=0] (s)")
        _a(f"{config.yaw_neut:<14.4f}   YawNeut      - Neutral yaw position--yaw spring force is zero at this yaw (degrees)")
        _a(f"{config.yaw_spr:<14.1f}   YawSpr       - Nacelle-yaw spring constant (N-m/rad)")
        _a(f"{config.yaw_damp:<14.1f}   YawDamp      - Nacelle-yaw damping constant (N-m/(rad/s))")
        _a(f"{config.tyaw_man_s:<14.1f}   TYawManS     - Time to start override yaw maneuver and target yaw angle (s)")
        _a(f"{config.yaw_man_rat:<14.4f}   YawManRat    - Yaw maneuver rate (in absolute value) (deg/s)")
        _a(f"{config.nac_yaw_f:<14.4f}   NacYawF      - Final yaw angle for override yaw maneuvers (degrees)")

        # --- Aerodynamic Flow Control ---
        _a("---------------------- AERODYNAMIC FLOW CONTROL ------------------------------------")
        _a(f"{config.afc_mode:<14d}   AfCmode      - Airfoil control mode {{0: none, 1: sine wave cycle, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")
        _a(f"{config.afc_mean:<14.4f}   AfC_Mean     - Mean level for cosine cycling or steady value (-) [used only with AfCmode==1]")
        _a(f"{config.afc_amp:<14.4f}   AfC_Amp      - Amplitude for cosine cycling of flap signal (-) [used only with AfCmode==1]")
        _a(f"{config.afc_phase:<14.4f}   AfC_Phase    - Phase relative to the blade azimuth (0 is vertical) for cosine cycling of flap signal (degrees) [used only with AfCmode==1]")

        # --- Cable Control ---
        _a("---------------------- CABLE CONTROL -----------------------------------------------")
        _a(f"{config.cc_mode:<14d}   CCmode       - Cable control mode {{0: none, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}} (switch)")

        # --- Structural Control ---
        _a("---------------------- STRUCTURAL CONTROL ------------------------------------------")
        _a(f"{config.num_b_stc:<14d}   NumBStC      - Number of blade structural controllers (integer)")
        if config.b_stc_files:
            _a(f'{"  ".join(f\'"{f}\'' for f in config.b_stc_files):<40s}   BStCfiles    - Name of the files for blade structural controllers (quoted strings) [unused when NumBStC==0]')
        else:
            _a(f'{"\"unused\"":<40s}   BStCfiles    - Name of the files for blade structural controllers (quoted strings) [unused when NumBStC==0]')
        _a(f"{config.num_n_stc:<14d}   NumNStC      - Number of nacelle structural controllers (integer)")
        if config.n_stc_files:
            _a(f'{"  ".join(f\'"{f}\'' for f in config.n_stc_files):<40s}   NStCfiles    - Name of the files for nacelle structural controllers (quoted strings) [unused when NumNStC==0]')
        else:
            _a(f'{"\"unused\"":<40s}   NStCfiles    - Name of the files for nacelle structural controllers (quoted strings) [unused when NumNStC==0]')
        _a(f"{config.num_t_stc:<14d}   NumTStC      - Number of tower structural controllers (integer)")
        if config.t_stc_files:
            _a(f'{"  ".join(f\'"{f}\'' for f in config.t_stc_files):<40s}   TStCfiles    - Name of the files for tower structural controllers (quoted strings) [unused when NumTStC==0]')
        else:
            _a(f'{"\"unused\"":<40s}   TStCfiles    - Name of the files for tower structural controllers (quoted strings) [unused when NumTStC==0]')
        _a(f"{config.num_s_stc:<14d}   NumSStC      - Number of substructure structural controllers (integer)")
        if config.s_stc_files:
            _a(f'{"  ".join(f\'"{f}\'' for f in config.s_stc_files):<40s}   SStCfiles    - Name of the files for substructure structural controllers (quoted strings) [unused when NumSStC==0]')
        else:
            _a(f'{"\"unused\"":<40s}   SStCfiles    - Name of the files for substructure structural controllers (quoted strings) [unused when NumSStC==0]')

        # --- Bladed Interface ---
        _a("---------------------- BLADED INTERFACE ---------------------------------------- [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f'{"\"" + config.dll_file_name + "\"":<40s}   DLL_FileName - Name/location of the dynamic library {{.dll [Windows] or .so [Linux]}} in the Bladed-DLL format (-) [used with PCMode=5, VSContrl=5, or YCMode=5]')
        _a(f'{"\"" + config.dll_in_file + "\"":<40s}   DLL_InFile   - Name of input file sent to the DLL (-) [used with PCMode=5, VSContrl=5, or YCMode=5]')
        _a(f'{"\"" + config.dll_proc_name + "\"":<40s}   DLL_ProcName - Name of procedure in DLL to be called (-) [case sensitive; used with PCMode=5, VSContrl=5, or YCMode=5]')
        _a(f'{"\"" + config.dll_dt + "\"":<14s}   DLL_DT       - Communication interval for dynamic library (s) (or \"default\") [used with PCMode=5, VSContrl=5, or YCMode=5]')
        _a(f"{_f(config.dll_ramp):<14s}   DLL_Ramp     - Whether a linear ramp should be used between DLL_DT time steps [introduces time shift when true] (flag) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.bp_cutoff:<14.1f}   BPCutoff     - Cutoff frequency for low-pass filter on blade pitch from DLL (Hz) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.nac_yaw_north:<14.1f}   NacYaw_North - Reference yaw angle of the nacelle when the upwind end points due North (deg) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_cntrl:<14d}   Ptch_Cntrl   - Record 28: Use individual pitch control {{0: collective pitch; 1: individual pitch control}} (switch) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_set_pnt:<14.4f}   Ptch_SetPnt  - Record  5: Below-rated pitch angle set-point (deg) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_min:<14.4f}   Ptch_Min     - Record  6: Minimum pitch angle (deg) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_max:<14.4f}   Ptch_Max     - Record  7: Maximum pitch angle (deg) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_rate_min:<14.4f}   PtchRate_Min - Record  8: Minimum pitch rate (most negative value allowed) (deg/s) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.ptch_rate_max:<14.4f}   PtchRate_Max - Record  9: Maximum pitch rate  (deg/s) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gain_om:<14.4f}   Gain_OM      - Record 16: Optimal mode gain (N-m/(rad/s)^2) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gen_spd_min_om:<14.4f}   GenSpd_MinOM - Record 17: Minimum generator speed (rpm) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gen_spd_max_om:<14.4f}   GenSpd_MaxOM - Record 18: Optimal mode maximum speed (rpm) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gen_spd_dem:<14.4f}   GenSpd_Dem   - Record 19: Demanded generator speed above rated (rpm) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gen_trq_dem:<14.4f}   GenTrq_Dem   - Record 22: Demanded generator torque above rated (N-m) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a(f"{config.gen_pwr_dem:<14.1f}   GenPwr_Dem   - Record 13: Demanded power (W) [used with PCMode=5, VSContrl=5, or YCMode=5]")

        # --- Bladed Interface Torque-Speed LUT ---
        _a("---------------------- BLADED INTERFACE TORQUE-SPEED LOOK-UP TABLE -----------------")
        _a(f"{config.dll_num_trq:<14d}   DLL_NumTrq   - Record 26: No. of points in torque-speed look-up table {{0 = none and use the optimal mode parameters; nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0}} (-) [used with PCMode=5, VSContrl=5, or YCMode=5]")
        _a("GenSpd_TLU   GenTrq_TLU")
        _a("(rpm)        (N-m)")

        # --- Output ---
        _a("---------------------- OUTPUT --------------------------------------------------")
        _a(f"{_f(config.sum_print):<14s}   SumPrint     - Print summary data to <RootName>.sum (flag) (currently unused)")
        _a(f"{config.out_file:<14d}   OutFile      - Switch to determine where output will be placed: {{1: in module output file only; 2: in glue code output file only; 3: both}} (currently unused)")
        _a(f"{_f(config.tab_delim):<14s}   TabDelim     - Use tab delimiters in text tabular output file? (flag) (currently unused)")
        _a(f'{"\"" + config.out_fmt + "\"":<14s}   OutFmt       - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)')
        _a(f"{config.t_start:<14.1f}   TStart       - Time to begin tabular output (s) (currently unused)")
        _a("                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)")
        for out in config.out_list:
            _a(out)
        _a('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)')
        _a("---------------------------------------------------------------------------------------")

        return "\n".join(lines) + "\n"

    def generate_discon_in(self, config: DISCONConfig) -> str:
        """Generate a ROSCO-style DISCON.IN controller parameter file.

        Parameters
        ----------
        config : DISCONConfig
            ROSCO controller configuration parameters.

        Returns
        -------
        str
            Complete DISCON.IN file content.
        """
        lines: list[str] = []
        _a = lines.append

        _a("! Controller parameter input file for the WindForge ROSCO controller")
        _a("! Generated by WindForge")
        _a("")

        _a("!------- DEBUG ------------------------------------------------------------")
        _a(f"{config.log_file:<20d}! LoggingLevel        - {{0: write no debug files, 1: write standard output .dbg-file, 2: write standard output .dbg2-file with more}}")
        _a(f"{config.log_level:<20d}! LoggingLevel2       - {{0: minimal, 1: standard, 2: verbose}}")
        _a("")

        _a("!------- CONTROLLER FLAGS -------------------------------------------------")
        _a(f"{config.pc_control_mode:<20d}! PC_ControlMode      - Blade pitch control mode {{0: No pitch, fix to fine pitch, 1: active PI blade pitch control}}")
        _a(f"{config.vs_control_mode:<20d}! VS_ControlMode      - Generator torque control mode in above rated conditions {{0: constant power, 1: constant torque, 2: TSR tracking with constant power, 3: TSR tracking with constant torque}}")
        _a(f"{config.y_control_mode:<20d}! Y_ControlMode       - Yaw control mode {{0: no yaw control, 1: yaw rate control, 2: yaw-by-IPC}}")
        _a(f"{config.fl_mode:<20d}! FL_ControlMode      - Floating specific control mode {{0: no floating feedback, 1: nacelle velocity feedback, 2: nacelle & platform velocity feedback}}")
        _a("")

        _a("!------- FILTERS ----------------------------------------------------------")
        _a(f"{config.f_lp_corner_freq:<20.5f}! F_LPFCornerFreq     - Corner frequency (-3dB point) in the low-pass filters (rad/s)")
        _a(f"{config.f_fl_hp_corner_freq:<20.5f}! F_FlpCornerFreq     - Corner frequency (-3dB point) in the HP filter on floating feedback (rad/s)")
        _a(f"{config.f_lp_type:<20d}! F_LPFType           - {{1: first-order low-pass filter, 2: second-order low-pass filter}}")
        _a(f"{config.f_lp_damping:<20.5f}! F_LPFDamping        - Damping coefficient [used only when F_LPFType=2, ignored otherwise]")
        _a(f"{config.f_notify_filter:<20d}! F_NotchFilterF      - Notch filter mode {{0: disable, 1: enable}}")
        _a("")

        _a("!------- BLADE PITCH CONTROL ----------------------------------------------")
        _a(f"{config.pc_gs_n:<20d}! PC_GS_n             - Amount of gain-scheduling table entries (-)")
        pc_gs_angles_str = "  ".join(f"{v:.6f}" for v in config.pc_gs_angles)
        _a(f"{pc_gs_angles_str}       ! PC_GS_angles        - Gain-scheduling table: pitch angles (rad)")
        pc_gs_kp_str = "  ".join(f"{v:.6f}" for v in config.pc_gs_kp)
        _a(f"{pc_gs_kp_str}       ! PC_GS_KP            - Gain-scheduling table: pitch controller kp gains (s)")
        pc_gs_ki_str = "  ".join(f"{v:.6f}" for v in config.pc_gs_ki)
        _a(f"{pc_gs_ki_str}       ! PC_GS_KI            - Gain-scheduling table: pitch controller ki gains (-)")
        _a(f"{config.pc_max_pitch:<20.5f}! PC_MaxPit            - Maximum physical pitch limit (rad)")
        _a(f"{config.pc_min_pitch:<20.5f}! PC_MinPit            - Minimum physical pitch limit (rad)")
        _a(f"{config.pc_max_rat:<20.5f}! PC_MaxRat            - Maximum pitch rate (in absolute value) (rad/s)")
        _a(f"{config.pc_min_rat:<20.5f}! PC_MinRat            - Minimum pitch rate (in absolute value, most negative) (rad/s)")
        _a(f"{config.pc_ref_speed:<20.5f}! PC_RefSpd            - Desired (reference) HSS speed for pitch controller (rad/s)")
        _a(f"{config.pc_fine_pitch:<20.5f}! PC_FinePit           - Record 5: Below-rated pitch angle set-point (rad)")
        _a(f"{config.pc_switch:<20.5f}! PC_Switch            - Angle above lowest minimum pitch angle for switch (rad)")
        _a("")

        _a("!------- TORQUE CONTROL ---------------------------------------------------")
        _a(f"{config.vs_gen_eff:<20.5f}! VS_GenEff            - Generator efficiency (mechanical to electrical) (-)")
        _a(f"{config.vs_max_rat:<20.5f}! VS_MaxRat            - Maximum torque rate (in absolute value) (N-m/s)")
        _a(f"{config.vs_max_tq:<20.5f}! VS_MaxTq             - Maximum generator torque in Region 3 (HSS side) (N-m)")
        _a(f"{config.vs_min_tq:<20.5f}! VS_MinTq             - Minimum generator torque (N-m)")
        _a(f"{config.vs_min_om_spd:<20.5f}! VS_MinOMSpd          - Minimum generator speed (rad/s)")
        _a(f"{config.vs_rated_om_spd:<20.5f}! VS_RtOMSpd           - Rated generator speed (rad/s)")
        _a(f"{config.vs_rgn2_k:<20.5f}! VS_Rgn2K             - Generator torque constant in Region 2 (HSS side) (N-m/(rad/s)^2)")
        _a(f"{config.vs_ref_spd:<20.5f}! VS_RefSpd            - Rated generator speed (rad/s)")
        _a(f"{config.vs_rated_gen_pwr:<20.1f}! VS_RtPwr             - Wind turbine rated power (W)")
        _a(f"{config.vs_tsr:<20.5f}! VS_TSRopt            - Power-maximizing controller desired TSR (-)")
        _a("")

        _a("!------- SETPOINT SMOOTHER ------------------------------------------------")
        _a(f"{config.ss_vsgain:<20.5f}! SS_VSGain            - Variable speed torque controller setpoint smoother gain bias percentage [%%]")
        _a(f"{config.ss_pcgain:<20.5f}! SS_PCGain            - Collective pitch controller setpoint smoother gain bias percentage [%%]")
        _a("")

        _a("!------- WIND SPEED ESTIMATOR ---------------------------------------------")
        _a(f"{config.we_blade_radius:<20.5f}! WE_BladeRadius       - Blade length (distance from hub center to blade tip) (m)")
        _a(f"{config.we_cp_n:<20d}! WE_CP_n              - Amount of parameters in the Cp array (-)")
        we_cp_str = "  ".join(f"{v:.5f}" for v in config.we_cp)
        _a(f"{we_cp_str:<20s}  ! WE_CP                - Parameters that define the parameterized CP(lambda) function (-)")
        _a(f"{config.we_gamma:<20.5f}! WE_Gamma             - Adaption gain of the wind speed estimator algorithm (m/rad)")
        _a(f"{config.we_gear_ratio:<20.5f}! WE_GearboxRatio      - Gearbox ratio (-, used only for determination of generator speed and target determine generator speed)")
        _a(f"{config.we_rlspd:<20.5f}! WE_Rlspd             - Rated generator speed (rad/s)")
        _a(f"{config.we_rated_pwr:<20.1f}! WE_RatedPwr          - Rated power (W)")
        _a(f"{config.we_min_pitch:<20.5f}! WE_MinPitch          - Minimum pitch angle (rad)")
        _a(f"{config.we_rated_v:<20.5f}! WE_RatedV            - Rated wind speed (m/s)")
        _a(f"{config.we_rated_omega:<20.5f}! WE_RatedOmega        - Rated generator speed (rad/s)")
        _a("")

        _a("!------- MINIMUM PITCH SATURATION -----------------------------------------")
        _a(f"{config.ps_bld_pitch_min:<20.5f}! PS_BldPitchMin       - Minimum blade pitch angle (rad)")
        _a("")

        _a("!------- SHUTDOWN ---------------------------------------------------------")
        _a(f"{config.sd_max_pitch:<20.5f}! SD_MaxPit            - Maximum blade pitch angle to initiate shutdown (rad)")
        _a(f"{config.sd_time_limit:<20.1f}! SD_TimeLimit         - Maximum time to complete shutdown (s)")
        _a("")

        return "\n".join(lines) + "\n"


def _flag(value: bool) -> str:
    """Convert a boolean to OpenFAST flag format."""
    return "True" if value else "False"
