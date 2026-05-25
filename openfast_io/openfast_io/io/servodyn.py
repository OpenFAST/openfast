"""
ServoDynIO  – read / write ServoDyn, StC, DISCON, and spd_trq files.

Produces the dict structure expected by ``fst_vt``:

    {
        'ServoDyn': { ... },
        'BStC': [ {stc_dict}, ... ],
        'NStC': [ {stc_dict}, ... ],
        'TStC': [ {stc_dict}, ... ],
        'SStC': [ {stc_dict}, ... ],
        'DISCON_in': { ... } | None,
        'spd_trq':   { ... } | None,
    }
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Callable

from .base import ModuleIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
)


# ---------------------------------------------------------------------------
# tiny helpers copied from FAST_writer (avoids coupling to the writer class)
# ---------------------------------------------------------------------------

def _float_default_out(val: Any) -> str:
    if isinstance(val, float):
        return f'{val:<22}'
    return f'{val!s:<22}'


def _int_default_out(val: Any) -> str:
    if isinstance(val, int):
        return f'{val:<22d}'
    return f'{val!s:<22}'


def _get_outlist(outlist_dict: dict, keys: list) -> list:
    """Return the list-of-lists for *keys* in an outlist dict."""
    out = []
    for key in keys:
        if key in outlist_dict:
            out.append(outlist_dict[key])
    return out


# ---------------------------------------------------------------------------
# Optional ROSCO imports
# ---------------------------------------------------------------------------
try:
    from rosco.toolbox.utilities import read_DISCON, load_from_txt
    from rosco.toolbox import turbine as ROSCO_turbine
    from rosco.toolbox import utilities as ROSCO_utilities
    _ROSCO = True
except ImportError:
    _ROSCO = False


class ServoDynIO(ModuleIO):
    """Read / write ServoDyn v1.05+ input files and sub-files."""

    # ------------------------------------------------------------------
    # read
    # ------------------------------------------------------------------
    def read(
        self,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        read_outlist_fn: Optional[Callable] = None,
        path2dll: Optional[str] = None,
        servo_file_rel: str = '',
    ) -> dict:
        """Read a ServoDyn input file and all referenced sub-files.

        Parameters
        ----------
        file_path : str
            Path to the ServoDyn ``.dat`` file.
        base_dir : str
            Base directory for resolving relative paths.
        outlist : dict, optional
            Pre-allocated outlist dict; a 'ServoDyn' key will be populated.
        read_outlist_fn : callable, optional
            ``read_outlist(f, key)`` helper that reads the OutList section.
        path2dll : str, optional
            Override path for DLL_FileName.
        servo_file_rel : str
            Relative path from ``base_dir`` to the ServoDyn file (used to
            resolve StC file paths the same way the legacy reader does).

        Returns
        -------
        dict with keys ``ServoDyn``, ``BStC``, ``NStC``, ``TStC``, ``SStC``,
        ``DISCON_in`` (or absent), ``spd_trq`` (or absent).
        """
        sd = {}
        sd_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        with open(sd_file) as f:
            f.readline()  # header 1
            f.readline()  # header 2

            # -- Simulation Control --
            f.readline()
            sd['Echo'] = bool_read(f.readline().split()[0])
            sd['DT'] = float_read(f.readline().split()[0])

            # -- Pitch Control --
            f.readline()
            sd['PCMode']       = int(f.readline().split()[0])
            sd['TPCOn']        = float_read(f.readline().split()[0])
            sd['PitNeut(1)']   = float_read(f.readline().split()[0])
            sd['PitNeut(2)']   = float_read(f.readline().split()[0])
            sd['PitNeut(3)']   = float_read(f.readline().split()[0])
            sd['PitSpr(1)']    = float_read(f.readline().split()[0])
            sd['PitSpr(2)']    = float_read(f.readline().split()[0])
            sd['PitSpr(3)']    = float_read(f.readline().split()[0])
            sd['PitDamp(1)']   = float_read(f.readline().split()[0])
            sd['PitDamp(2)']   = float_read(f.readline().split()[0])
            sd['PitDamp(3)']   = float_read(f.readline().split()[0])
            sd['TPitManS(1)']  = float_read(f.readline().split()[0])
            sd['TPitManS(2)']  = float_read(f.readline().split()[0])
            sd['TPitManS(3)']  = float_read(f.readline().split()[0])
            sd['PitManRat(1)'] = float_read(f.readline().split()[0])
            sd['PitManRat(2)'] = float_read(f.readline().split()[0])
            sd['PitManRat(3)'] = float_read(f.readline().split()[0])
            sd['BlPitchF(1)']  = float_read(f.readline().split()[0])
            sd['BlPitchF(2)']  = float_read(f.readline().split()[0])
            sd['BlPitchF(3)']  = float_read(f.readline().split()[0])

            # -- Generator and Torque Control --
            f.readline()
            sd['VSContrl'] = int(f.readline().split()[0])
            sd['GenModel'] = int(f.readline().split()[0])
            sd['GenEff']   = float_read(f.readline().split()[0])
            sd['GenTiStr'] = bool_read(f.readline().split()[0])
            sd['GenTiStp'] = bool_read(f.readline().split()[0])
            sd['SpdGenOn'] = float_read(f.readline().split()[0])
            sd['TimGenOn'] = float_read(f.readline().split()[0])
            sd['TimGenOf'] = float_read(f.readline().split()[0])

            # -- Simple Variable-Speed Torque Control --
            f.readline()
            sd['VS_RtGnSp'] = float_read(f.readline().split()[0])
            sd['VS_RtTq']   = float_read(f.readline().split()[0])
            sd['VS_Rgn2K']  = float_read(f.readline().split()[0])
            sd['VS_SlPc']   = float_read(f.readline().split()[0])

            # -- Simple Induction Generator --
            f.readline()
            sd['SIG_SlPc'] = float_read(f.readline().split()[0])
            sd['SIG_SySp'] = float_read(f.readline().split()[0])
            sd['SIG_RtTq'] = float_read(f.readline().split()[0])
            sd['SIG_PORt'] = float_read(f.readline().split()[0])

            # -- Thevenin-Equivalent Induction Generator --
            f.readline()
            sd['TEC_Freq'] = float_read(f.readline().split()[0])
            sd['TEC_NPol'] = int(f.readline().split()[0])
            sd['TEC_SRes'] = float_read(f.readline().split()[0])
            sd['TEC_RRes'] = float_read(f.readline().split()[0])
            sd['TEC_VLL']  = float_read(f.readline().split()[0])
            sd['TEC_SLR']  = float_read(f.readline().split()[0])
            sd['TEC_RLR']  = float_read(f.readline().split()[0])
            sd['TEC_MR']   = float_read(f.readline().split()[0])

            # -- High-Speed Shaft Brake --
            f.readline()
            sd['HSSBrMode'] = int(f.readline().split()[0])
            sd['THSSBrDp']  = float_read(f.readline().split()[0])
            sd['HSSBrDT']   = float_read(f.readline().split()[0])
            sd['HSSBrTqF']  = float_read(f.readline().split()[0])

            # -- Nacelle-Yaw Control --
            f.readline()
            sd['YCMode']    = int(f.readline().split()[0])
            sd['TYCOn']     = float_read(f.readline().split()[0])
            sd['YawNeut']   = float_read(f.readline().split()[0])
            sd['YawSpr']    = float_read(f.readline().split()[0])
            sd['YawDamp']   = float_read(f.readline().split()[0])
            sd['TYawManS']  = float_read(f.readline().split()[0])
            sd['YawManRat'] = float_read(f.readline().split()[0])
            sd['NacYawF']   = float_read(f.readline().split()[0])

            # -- Aero Flow Control --
            f.readline()
            sd['AfCmode']   = int(f.readline().split()[0])
            sd['AfC_Mean']  = float_read(f.readline().split()[0])
            sd['AfC_Amp']   = float_read(f.readline().split()[0])
            sd['AfC_Phase'] = float_read(f.readline().split()[0])

            # -- Structural Control --
            f.readline()
            sd['NumBStC']   = int(f.readline().split()[0])
            sd['BStCfiles'] = read_array(f, sd['NumBStC'], array_type=str)
            sd['NumNStC']   = int(f.readline().split()[0])
            sd['NStCfiles'] = read_array(f, sd['NumNStC'], array_type=str)
            sd['NumTStC']   = int(f.readline().split()[0])
            sd['TStCfiles'] = read_array(f, sd['NumTStC'], array_type=str)
            sd['NumSStC']   = int(f.readline().split()[0])
            sd['SStCfiles'] = read_array(f, sd['NumSStC'], array_type=str)

            # -- Cable Control --
            f.readline()
            sd['CCmode'] = int(f.readline().split()[0])

            # -- Bladed Interface --
            f.readline()
            if not path2dll:
                sd['DLL_FileName'] = os.path.abspath(
                    os.path.normpath(os.path.join(os.path.split(sd_file)[0],
                                                  quoted_read(f.readline().split()[0]))))
            else:
                f.readline()
                sd['DLL_FileName'] = path2dll

            sd['DLL_InFile']   = os.path.abspath(
                os.path.normpath(os.path.join(os.path.split(sd_file)[0],
                                              quoted_read(f.readline().split()[0]))))
            sd['DLL_ProcName'] = quoted_read(f.readline().split()[0])
            dll_dt_line = f.readline().split()[0]
            try:
                sd['DLL_DT'] = float_read(dll_dt_line)
            except Exception:
                sd['DLL_DT'] = dll_dt_line[1:-1]
            sd['DLL_Ramp']     = bool_read(f.readline().split()[0])
            sd['BPCutoff']     = float_read(f.readline().split()[0])
            sd['NacYaw_North'] = float_read(f.readline().split()[0])
            sd['Ptch_Cntrl']   = int(f.readline().split()[0])
            sd['Ptch_SetPnt']  = float_read(f.readline().split()[0])
            sd['Ptch_Min']     = float_read(f.readline().split()[0])
            sd['Ptch_Max']     = float_read(f.readline().split()[0])
            sd['PtchRate_Min'] = float_read(f.readline().split()[0])
            sd['PtchRate_Max'] = float_read(f.readline().split()[0])
            sd['Gain_OM']      = float_read(f.readline().split()[0])
            sd['GenSpd_MinOM'] = float_read(f.readline().split()[0])
            sd['GenSpd_MaxOM'] = float_read(f.readline().split()[0])
            sd['GenSpd_Dem']   = float_read(f.readline().split()[0])
            sd['GenTrq_Dem']   = float_read(f.readline().split()[0])
            sd['GenPwr_Dem']   = float_read(f.readline().split()[0])

            f.readline()  # section header

            sd['DLL_NumTrq'] = int(f.readline().split()[0])
            f.readline()  # column names
            f.readline()  # units
            sd['GenSpd_TLU'] = [None] * sd['DLL_NumTrq']
            sd['GenTrq_TLU'] = [None] * sd['DLL_NumTrq']
            for i in range(sd['DLL_NumTrq']):
                data = f.readline().split()
                sd['GenSpd_TLU'][i] = float_read(data[0])
                sd['GenTrq_TLU'][i] = float_read(data[1])

            # -- Output --
            f.readline()
            sd['SumPrint'] = bool_read(f.readline().split()[0])
            sd['OutFile']  = int(f.readline().split()[0])
            sd['TabDelim'] = bool_read(f.readline().split()[0])
            sd['OutFmt']   = quoted_read(f.readline().split()[0])
            sd['TStart']   = float_read(f.readline().split()[0])

            # -- Outlist --
            f.readline()
            if read_outlist_fn is not None:
                read_outlist_fn(f, 'ServoDyn')

        # ------------------------------------------------------------------
        # Build result dict
        # ------------------------------------------------------------------
        result: Dict[str, Any] = {'ServoDyn': sd}

        # Read StC files  (paths are relative to the ServoDyn file directory)
        svd_dir = os.path.dirname(servo_file_rel) if servo_file_rel else ''

        result['BStC'] = []
        for fn in sd['BStCfiles']:
            result['BStC'].append(self._read_stc(fn, base_dir, svd_dir))
        result['NStC'] = []
        for fn in sd['NStCfiles']:
            result['NStC'].append(self._read_stc(fn, base_dir, svd_dir))
        result['TStC'] = []
        for fn in sd['TStCfiles']:
            result['TStC'].append(self._read_stc(fn, base_dir, svd_dir))
        result['SStC'] = []
        for fn in sd['SStCfiles']:
            result['SStC'].append(self._read_stc(fn, base_dir, svd_dir))

        # Read DISCON_in (ROSCO)
        if _ROSCO:
            discon = self._read_discon(sd, base_dir)
            if discon is not None:
                result['DISCON_in'] = discon

        # Read spd_trq
        if sd['VSContrl'] == 3:
            result['spd_trq'] = self._read_spd_trq('spd_trq.dat', base_dir)

        return result

    # ------------------------------------------------------------------
    # write
    # ------------------------------------------------------------------
    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        run_dir: str = '',
        naming_out: str = '',
    ) -> None:
        """Write a ServoDyn input file and referenced sub-files.

        Parameters
        ----------
        data
            Dict with keys ``ServoDyn``, ``BStC``, ``NStC``, etc.
        file_path
            Target path for the ServoDyn ``.dat`` file.
        base_dir
            Base directory for output (defaults to dirname of *file_path*).
        outlist
            Outlist dict – 'ServoDyn' key used.
        run_dir
            Run directory used for writing sub-files (DISCON, StC, spd_trq).
            Defaults to ``os.path.dirname(file_path)``.
        naming_out
            Naming prefix for generated sub-files.
        """
        sd = data['ServoDyn']
        out_dir = run_dir or os.path.dirname(file_path) or '.'

        with open(file_path, 'w') as f:
            f.write('------- SERVODYN INPUT FILE --------------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sd['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['DT'], 'DT', '- Communication interval for controllers (s) (or "default")\n'))
            f.write('---------------------- PITCH CONTROL -------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['PCMode'], 'PCMode', '- Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TPCOn'], 'TPCOn', '- Time to enable active pitch control (s) [unused when PCMode=0]\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'PitNeut({idx})'], f'PitNeut({idx})', f'- Blade {idx} neutral pitch position--pitch spring moment is zero at this pitch (degrees)\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'PitSpr({idx})'], f'PitSpr({idx})', f'- Blade {idx} pitch spring constant (N-m/rad)\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'PitDamp({idx})'], f'PitDamp({idx})', f'- Blade {idx} pitch damping constant (N-m/(rad/s))\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'TPitManS({idx})'], f'TPitManS({idx})', f'- Time to start override pitch maneuver for blade {idx} and end standard pitch control (s)\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'PitManRat({idx})'], f'PitManRat({idx})', f'- Pitch rate at which override pitch maneuver heads toward final pitch angle for blade {idx} (deg/s)\n'))
            for idx in range(1, 4):
                f.write('{:<22} {:<11} {:}'.format(sd[f'BlPitchF({idx})'], f'BlPitchF({idx})', f'- Blade {idx} final pitch for pitch maneuvers (degrees)\n'))

            f.write('---------------------- GENERATOR AND TORQUE CONTROL ----------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['VSContrl'], 'VSContrl', '- Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenModel'], 'GenModel', '- Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenEff'], 'GenEff', '- Generator efficiency [ignored by the Thevenin and user-defined generator models] (%)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['GenTiStr'], 'GenTiStr', '- Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['GenTiStp'], 'GenTiStp', '- Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['SpdGenOn'], 'SpdGenOn', '- Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TimGenOn'], 'TimGenOn', '- Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TimGenOf'], 'TimGenOf', '- Time to turn off the generator (s) [used only when GenTiStp=True]\n'))

            f.write('---------------------- SIMPLE VARIABLE-SPEED TORQUE CONTROL --------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['VS_RtGnSp'], 'VS_RtGnSp', '- Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['VS_RtTq'], 'VS_RtTq', '- Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['VS_Rgn2K'], 'VS_Rgn2K', '- Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['VS_SlPc'], 'VS_SlPc', '- Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]\n'))

            f.write('---------------------- SIMPLE INDUCTION GENERATOR ------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['SIG_SlPc'], 'SIG_SlPc', '- Rated generator slip percentage (%) [used only when VSContrl=0 and GenModel=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['SIG_SySp'], 'SIG_SySp', '- Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['SIG_RtTq'], 'SIG_RtTq', '- Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['SIG_PORt'], 'SIG_PORt', '- Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]\n'))

            f.write('---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_Freq'], 'TEC_Freq', '- Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_NPol'], 'TEC_NPol', '- Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_SRes'], 'TEC_SRes', '- Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_RRes'], 'TEC_RRes', '- Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_VLL'], 'TEC_VLL', '- Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_SLR'], 'TEC_SLR', '- Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_RLR'], 'TEC_RLR', '- Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TEC_MR'], 'TEC_MR', '- Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n'))

            f.write('---------------------- HIGH-SPEED SHAFT BRAKE ----------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['HSSBrMode'], 'HSSBrMode', '- HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['THSSBrDp'], 'THSSBrDp', '- Time to initiate deployment of the HSS brake (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['HSSBrDT'], 'HSSBrDT', '- Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['HSSBrTqF'], 'HSSBrTqF', '- Fully deployed HSS-brake torque (N-m)\n'))

            f.write('---------------------- NACELLE-YAW CONTROL -------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['YCMode'], 'YCMode', '- Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TYCOn'], 'TYCOn', '- Time to enable active yaw control (s) [unused when YCMode=0]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['YawNeut'], 'YawNeut', '- Neutral yaw position--yaw spring force is zero at this yaw (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['YawSpr'], 'YawSpr', '- Nacelle-yaw spring constant (N-m/rad)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['YawDamp'], 'YawDamp', '- Nacelle-yaw damping constant (N-m/(rad/s))\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TYawManS'], 'TYawManS', '- Time to start override yaw maneuver and end standard yaw control (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['YawManRat'], 'YawManRat', '- Yaw maneuver rate (in absolute value) (deg/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['NacYawF'], 'NacYawF', '- Final yaw angle for override yaw maneuvers (degrees)\n'))

            f.write('---------------------- Aerodynamic Flow Control -------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['AfCmode'], 'AfCmode', '- Airfoil control mode {0: none, 1: cosine wave cycle, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['AfC_Mean'], 'AfC_Mean', '- Mean level for cosine cycling or steady value (-) [used only with AfCmode==1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['AfC_Amp'], 'AfC_Amp', '- Amplitude for for cosine cycling of flap signal (AfC = AfC_Amp*cos(Azimuth+phase)+AfC_mean) (-) [used only with AfCmode==1]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['AfC_Phase'], 'AfC_phase', '- Phase relative to the blade azimuth (0 is vertical) for for cosine cycling of flap signal (deg) [used only with AfCmode==1]\n'))

            f.write('---------------------- STRUCTURAL CONTROL ---------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['NumBStC'], 'NumBStC', '- Number of blade structural controllers (integer)\n'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(sd['BStCfiles']) + '"', 'BStCfiles', '- Name of the files for blade structural controllers (quoted strings) [unused when NumBStC==0]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['NumNStC'], 'NumNStC', '- Number of nacelle structural controllers (integer)\n'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(sd['NStCfiles']) + '"', 'NStCfiles', '- Name of the files for nacelle structural controllers (quoted strings) [unused when NumNStC==0]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['NumTStC'], 'NumTStC', '- Number of tower structural controllers (integer)\n'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(sd['TStCfiles']) + '"', 'TStCfiles', '- Name of the files for tower structural controllers (quoted strings) [unused when NumTStC==0]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['NumSStC'], 'NumSStC', '- Number of substructure structural controllers (integer)\n'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + '" "'.join(sd['SStCfiles']) + '"', 'SStCfiles', '- Name of the files for substructure structural controllers (quoted strings) [unused when NumSStC==0]\n'))

            f.write('---------------------- CABLE CONTROL ---------------------------------------- \n')
            f.write('{:<22} {:<11} {:}'.format(sd['CCmode'], 'CCmode', '- Cable control mode {0: none, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)\n'))

            f.write('---------------------- BLADED INTERFACE ---------------------------------------- [used only with Bladed Interface]\n')
            f.write('{:<22} {:<11} {:}'.format('"'+sd['DLL_FileName']+'"', 'DLL_FileName', '- Name/location of the dynamic library {.dll [Windows] or .so [Linux]} in the Bladed-DLL format (-) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format('"'+sd['DLL_InFile']+'"', 'DLL_InFile', '- Name of input file sent to the DLL (-) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format('"'+sd['DLL_ProcName']+'"', 'DLL_ProcName', '- Name of procedure in DLL to be called (-) [case sensitive; used only with DLL Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['DLL_DT'], 'DLL_DT', '- Communication interval for dynamic library (s) (or "default") [used only with Bladed Interface]\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['DLL_Ramp'], 'DLL_Ramp', '- Whether a linear ramp should be used between DLL_DT time steps [introduces time shift when true] (flag) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['BPCutoff'], 'BPCutoff', '- Cutoff frequency for low-pass filter on blade pitch from DLL (Hz) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['NacYaw_North'], 'NacYaw_North', '- Reference yaw angle of the nacelle when the upwind end points due North (deg) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['Ptch_Cntrl'], 'Ptch_Cntrl', '- Record 28: Use individual pitch control {0: collective pitch; 1: individual pitch control} (switch) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['Ptch_SetPnt'], 'Ptch_SetPnt', '- Record  5: Below-rated pitch angle set-point (deg) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['Ptch_Min'], 'Ptch_Min', '- Record  6: Minimum pitch angle (deg) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['Ptch_Max'], 'Ptch_Max', '- Record  7: Maximum pitch angle (deg) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['PtchRate_Min'], 'PtchRate_Min', '- Record  8: Minimum pitch rate (most negative value allowed) (deg/s) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['PtchRate_Max'], 'PtchRate_Max', '- Record  9: Maximum pitch rate  (deg/s) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['Gain_OM'], 'Gain_OM', '- Record 16: Optimal mode gain (Nm/(rad/s)^2) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenSpd_MinOM'], 'GenSpd_MinOM', '- Record 17: Minimum generator speed (rpm) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenSpd_MaxOM'], 'GenSpd_MaxOM', '- Record 18: Optimal mode maximum speed (rpm) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenSpd_Dem'], 'GenSpd_Dem', '- Record 19: Demanded generator speed above rated (rpm) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenTrq_Dem'], 'GenTrq_Dem', '- Record 22: Demanded generator torque above rated (Nm) [used only with Bladed Interface]\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['GenPwr_Dem'], 'GenPwr_Dem', '- Record 13: Demanded power (W) [used only with Bladed Interface]\n'))

            f.write('---------------------- BLADED INTERFACE TORQUE-SPEED LOOK-UP TABLE -------------\n')
            f.write('{:<22} {:<11} {:}'.format(sd['DLL_NumTrq'], 'DLL_NumTrq', '- Record 26: No. of points in torque-speed look-up table {0 = none and use the optimal mode parameters; nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0} (-) [used only with Bladed Interface]\n'))
            f.write('{:<22}\t{:<22}\n'.format('GenSpd_TLU', 'GenTrq_TLU'))
            f.write('{:<22}\t{:<22}\n'.format('(rpm)', '(Nm)'))
            for i in range(sd['DLL_NumTrq']):
                f.write('{:<22}\t{:<22}\n'.format(sd['GenSpd_TLU'][i], sd['GenTrq_TLU'][i]))

            f.write('---------------------- OUTPUT --------------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sd['SumPrint'], 'SumPrint', '- Print summary data to <RootName>.sum (flag) (currently unused)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['OutFile'], 'OutFile', '- Switch to determine where output will be placed: {1: in module output file only; 2: in glue code output file only; 3: both} (currently unused)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) (currently unused)\n'))
            f.write('{:<22} {:<11} {:}'.format('"'+sd['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string) (currently unused)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['TStart'], 'TStart', '- Time to begin tabular output (s) (currently unused)\n'))
            f.write('              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

            if outlist is not None:
                ol = _get_outlist(outlist, ['ServoDyn'])
                for channel_list in ol:
                    for ch in channel_list:
                        f.write('"' + ch + '"\n')

            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
            f.write('---------------------------------------------------------------------------------------\n')

        # Write StC files
        for i, stc in enumerate(data.get('BStC', [])):
            self._write_stc(stc, sd['BStCfiles'][i], file_path, out_dir)
        for i, stc in enumerate(data.get('NStC', [])):
            self._write_stc(stc, sd['NStCfiles'][i], file_path, out_dir)
        for i, stc in enumerate(data.get('TStC', [])):
            self._write_stc(stc, sd['TStCfiles'][i], file_path, out_dir)
        for i, stc in enumerate(data.get('SStC', [])):
            self._write_stc(stc, sd['SStCfiles'][i], file_path, out_dir)

        # Write spd_trq
        if sd['VSContrl'] == 3 and 'spd_trq' in data:
            self._write_spd_trq(data['spd_trq'], out_dir)

    # ==================================================================
    # Private helpers
    # ==================================================================

    # ---- StC reader --------------------------------------------------
    @staticmethod
    def _read_stc(filename: str, base_dir: str, svd_dir: str) -> dict:
        """Read a single Structural-Controller input file."""
        stc: Dict[str, Any] = {}

        with open(os.path.join(base_dir, svd_dir, filename)) as f:
            f.readline()  # header 1
            f.readline()  # header 2
            f.readline()  # sim control header
            stc['Echo'] = bool_read(f.readline().split()[0])

            # DOF
            f.readline()
            stc['StC_DOF_MODE'] = int_read(f.readline().split()[0])
            stc['StC_X_DOF']    = bool_read(f.readline().split()[0])
            stc['StC_Y_DOF']    = bool_read(f.readline().split()[0])
            stc['StC_Z_DOF']    = bool_read(f.readline().split()[0])

            # Location
            f.readline()
            stc['StC_P_X'] = float_read(f.readline().split()[0])
            stc['StC_P_Y'] = float_read(f.readline().split()[0])
            stc['StC_P_Z'] = float_read(f.readline().split()[0])

            # Initial conditions
            f.readline()
            stc['StC_X_DSP']   = float_read(f.readline().split()[0])
            stc['StC_Y_DSP']   = float_read(f.readline().split()[0])
            stc['StC_Z_DSP']   = float_read(f.readline().split()[0])
            stc['StC_Z_PreLd'] = f.readline().split()[0]

            # Configuration
            f.readline()
            stc['StC_X_PSP'] = float_read(f.readline().split()[0])
            stc['StC_X_NSP'] = float_read(f.readline().split()[0])
            stc['StC_Y_PSP'] = float_read(f.readline().split()[0])
            stc['StC_Y_NSP'] = float_read(f.readline().split()[0])
            stc['StC_Z_PSP'] = float_read(f.readline().split()[0])
            stc['StC_Z_NSP'] = float_read(f.readline().split()[0])

            # Mass, stiffness, damping
            f.readline()
            stc['StC_X_M']    = float_read(f.readline().split()[0])
            stc['StC_Y_M']    = float_read(f.readline().split()[0])
            stc['StC_Z_M']    = float_read(f.readline().split()[0])
            stc['StC_Omni_M'] = float_read(f.readline().split()[0])
            stc['StC_X_K']    = float_read(f.readline().split()[0])
            stc['StC_Y_K']    = float_read(f.readline().split()[0])
            stc['StC_Z_K']    = float_read(f.readline().split()[0])
            stc['StC_X_C']    = float_read(f.readline().split()[0])
            stc['StC_Y_C']    = float_read(f.readline().split()[0])
            stc['StC_Z_C']    = float_read(f.readline().split()[0])
            stc['StC_X_KS']   = float_read(f.readline().split()[0])
            stc['StC_Y_KS']   = float_read(f.readline().split()[0])
            stc['StC_Z_KS']   = float_read(f.readline().split()[0])
            stc['StC_X_CS']   = float_read(f.readline().split()[0])
            stc['StC_Y_CS']   = float_read(f.readline().split()[0])
            stc['StC_Z_CS']   = float_read(f.readline().split()[0])

            # User-defined spring forces
            f.readline()
            stc['Use_F_TBL'] = bool_read(f.readline().split()[0])
            stc['NKInpSt']   = int_read(f.readline().split()[0])

            table: Dict[str, list] = {}
            table['X']   = [None] * stc['NKInpSt']
            table['F_X'] = [None] * stc['NKInpSt']
            table['Y']   = [None] * stc['NKInpSt']
            table['F_Y'] = [None] * stc['NKInpSt']
            table['Z']   = [None] * stc['NKInpSt']
            table['F_Z'] = [None] * stc['NKInpSt']

            f.readline()  # section header
            f.readline()  # col names
            f.readline()  # units
            for i in range(stc['NKInpSt']):
                ln = f.readline().split()
                table['X'][i]   = float(ln[0])
                table['F_X'][i] = float(ln[1])
                table['Y'][i]   = float(ln[2])
                table['F_Y'][i] = float(ln[3])
                table['Z'][i]   = float(ln[4])
                table['F_Z'][i] = float(ln[5])
            stc['SpringForceTable'] = table

            # Control
            f.readline()
            stc['StC_CMODE']    = int_read(f.readline().split()[0])
            stc['StC_CChan']    = int_read(f.readline().split()[0])
            stc['StC_SA_MODE']  = int_read(f.readline().split()[0])
            stc['StC_X_C_LOW']  = float_read(f.readline().split()[0])
            stc['StC_X_C_HIGH'] = float_read(f.readline().split()[0])
            stc['StC_Y_C_HIGH'] = float_read(f.readline().split()[0])
            stc['StC_Y_C_LOW']  = float_read(f.readline().split()[0])
            stc['StC_Z_C_HIGH'] = float_read(f.readline().split()[0])
            stc['StC_Z_C_LOW']  = float_read(f.readline().split()[0])
            stc['StC_X_C_BRAKE'] = float_read(f.readline().split()[0])
            stc['StC_Y_C_BRAKE'] = float_read(f.readline().split()[0])
            stc['StC_Z_C_BRAKE'] = float_read(f.readline().split()[0])

            # TLCD
            f.readline()
            stc['L_X']             = float_read(f.readline().split()[0])
            stc['B_X']             = float_read(f.readline().split()[0])
            stc['area_X']          = float_read(f.readline().split()[0])
            stc['area_ratio_X']    = float_read(f.readline().split()[0])
            stc['headLossCoeff_X'] = float_read(f.readline().split()[0])
            stc['rho_X']           = float_read(f.readline().split()[0])
            stc['L_Y']             = float_read(f.readline().split()[0])
            stc['B_Y']             = float_read(f.readline().split()[0])
            stc['area_Y']          = float_read(f.readline().split()[0])
            stc['area_ratio_Y']    = float_read(f.readline().split()[0])
            stc['headLossCoeff_Y'] = float_read(f.readline().split()[0])
            stc['rho_Y']           = float_read(f.readline().split()[0])

            # Prescribed time series
            f.readline()
            stc['PrescribedForcesCoord'] = int_read(f.readline().split()[0])
            stc['PrescribedForcesFile']  = os.path.join(base_dir, quoted_read(f.readline().split()[0]))
            f.readline()  # trailing line

        return stc

    # ---- StC writer --------------------------------------------------
    @staticmethod
    def _write_stc(stc: dict, filename: str, sd_file_path: str, out_dir: str) -> None:
        """Write a single StC input file."""
        sd_dir = os.path.dirname(sd_file_path)
        stc_file = os.path.join(sd_dir, filename)

        with open(stc_file, 'w') as f:
            f.write('------- STRUCTURAL CONTROL (StC) INPUT FILE ----------------------------\n')
            f.write('Generated with OpenFAST_IO\n')

            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(stc['Echo'], 'Echo', '- Echo input data to "<rootname>.SD.ech" (flag)\n'))

            f.write('---------------------- StC DEGREES OF FREEDOM ----------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_DOF_MODE'], 'StC_DOF_MODE', '- DOF mode (switch) {0: No StC or TLCD DOF; 1: StC_X_DOF, StC_Y_DOF, and/or StC_Z_DOF (three independent StC DOFs); 2: StC_XY_DOF (Omni-Directional StC); 3: StC_XYZ_DOF (Omni-Directional StC); 5: TLCD; 6: Prescribed force/moment time series; 7: Force determined by external DLL}\n'))
            f.write('{!s:<22} {:<11} {:}'.format(stc['StC_X_DOF'], 'StC_X_DOF', '- DOF on or off for StC X (flag) [Used only when StC_DOF_MODE=1]\n'))
            f.write('{!s:<22} {:<11} {:}'.format(stc['StC_Y_DOF'], 'StC_Y_DOF', '- DOF on or off for StC Y (flag) [Used only when StC_DOF_MODE=1]\n'))
            f.write('{!s:<22} {:<11} {:}'.format(stc['StC_Z_DOF'], 'StC_Z_DOF', '- DOF on or off for StC Z (flag) [Used only when StC_DOF_MODE=1]\n'))

            f.write('---------------------- StC LOCATION ------------------------------------------- [relative to the reference origin of component attached to]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_P_X'], 'StC_P_X', '- At rest X position of StC (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_P_Y'], 'StC_P_Y', '- At rest Y position of StC (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_P_Z'], 'StC_P_Z', '- At rest Z position of StC (m)\n'))

            f.write('---------------------- StC INITIAL CONDITIONS --------------------------------- [used only when StC_DOF_MODE=1, 2, or 3]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_DSP'], 'StC_X_DSP', '- StC X initial displacement (m) [relative to at rest position]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_DSP'], 'StC_Y_DSP', '- StC Y initial displacement (m) [relative to at rest position]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_DSP'], 'StC_Z_DSP', '- StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{!s:<22} {:<11} {:}'.format(stc['StC_Z_PreLd'], 'StC_Z_PreLd', '- StC Z pre-load (N) {"gravity" to offset for gravity load; "none" or 0 to turn off} [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))

            f.write('---------------------- StC CONFIGURATION -------------------------------------- [used only when StC_DOF_MODE=1, 2, or 3]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_PSP'], 'StC_X_PSP', '- Positive stop position (maximum X mass displacement) (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_NSP'], 'StC_X_NSP', '- Negative stop position (minimum X mass displacement) (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_PSP'], 'StC_Y_PSP', '- Positive stop position (maximum Y mass displacement) (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_NSP'], 'StC_Y_NSP', '- Negative stop position (minimum Y mass displacement) (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_PSP'], 'StC_Z_PSP', '- Positive stop position (maximum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_NSP'], 'StC_Z_NSP', '- Negative stop position (minimum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))

            f.write('---------------------- StC MASS, STIFFNESS, & DAMPING ------------------------- [used only when StC_DOF_MODE=1, 2, or 3]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_M'], 'StC_X_M', '- StC X mass (kg) [used only when StC_DOF_MODE=1 and StC_X_DOF=TRUE]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_M'], 'StC_Y_M', '- StC Y mass (kg) [used only when StC_DOF_MODE=1 and StC_Y_DOF=TRUE]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_M'], 'StC_Z_M', '- StC Z mass (kg) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Omni_M'], 'StC_Omni_M', '- StC omni mass (kg) [used only when StC_DOF_MODE=2 or 3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_K'], 'StC_X_K', '- StC X stiffness (N/m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_K'], 'StC_Y_K', '- StC Y stiffness (N/m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_K'], 'StC_Z_K', '- StC Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_C'], 'StC_X_C', '- StC X damping (N/(m/s))\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_C'], 'StC_Y_C', '- StC Y damping (N/(m/s))\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_C'], 'StC_Z_C', '- StC Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_KS'], 'StC_X_KS', '- Stop spring X stiffness (N/m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_KS'], 'StC_Y_KS', '- Stop spring Y stiffness (N/m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_KS'], 'StC_Z_KS', '- Stop spring Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_CS'], 'StC_X_CS', '- Stop spring X damping (N/(m/s))\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_CS'], 'StC_Y_CS', '- Stop spring Y damping (N/(m/s))\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_CS'], 'StC_Z_CS', '- Stop spring Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))

            f.write('---------------------- StC USER-DEFINED SPRING FORCES ------------------------- [used only when StC_DOF_MODE=1, 2, or 3]\n')
            f.write('{!s:<22} {:<11} {:}'.format(stc['Use_F_TBL'], 'Use_F_TBL', '- Use spring force from user-defined table (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['NKInpSt'], 'NKInpSt', '- Number of spring force input stations\n'))

            f.write('---------------------- StC SPRING FORCES TABLE -------------------------------- [used only when StC_DOF_MODE=1, 2, or 3]\n')
            f.write('X                F_X               Y              F_Y              Z              F_Z\n')
            f.write('(m)               (N)              (m)             (N)             (m)             (N)\n')
            table = stc['SpringForceTable']
            for x, f_x, y, f_y, z, f_z in zip(table['X'], table['F_X'], table['Y'], table['F_Y'], table['Z'], table['F_Z']):
                row = [x, f_x, y, f_y, z, f_z]
                f.write(' '.join(['{: 2.8e}'.format(val) for val in row]) + '\n')

            f.write('---------------------- StructUserProp CONTROL -------------------------------------------- [used only when StC_DOF_MODE=1, 2, 3, or 7]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['StC_CMODE'], 'StC_CMODE', '- Control mode (switch) {0:none; 1: Semi-Active Control Mode; 3: Active Control Mode through user subroutine; 4: Active Control Mode through Simulink (not available); 5: Active Control Mode through Bladed interface}\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_CChan'], 'StC_CChan', '- Control channel group (1:10) for stiffness and damping (StC_[XYZ]_K, StC_[XYZ]_C, and StC_[XYZ]_Brake) (specify additional channels for blade instances of StC active control -- one channel per blade) [used only when StC_DOF_MODE=1, 2, 3, or 7, and StC_CMODE=4 or 5]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_SA_MODE'], 'StC_SA_MODE', '- Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse velocity-based ground hook control; 3: displacement-based ground hook control 4: Phase difference Algorithm with Friction Force 5: Phase difference Algorithm with Damping Force} (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_C_HIGH'], 'StC_X_C_HIGH', '- StC X high damping for ground hook control\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_C_LOW'], 'StC_X_C_LOW', '- StC X low damping for ground hook control\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_C_HIGH'], 'StC_Y_C_HIGH', '- StC Y high damping for ground hook control\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_C_LOW'], 'StC_Y_C_LOW', '- StC Y low damping for ground hook control\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_C_HIGH'], 'StC_Z_C_HIGH', '- StC Z high damping for ground hook control [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_C_LOW'], 'StC_Z_C_LOW', '- StC Z low damping for ground hook control  [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_X_C_BRAKE'], 'StC_X_C_BRAKE', '- StC X high damping for braking the StC (Don\'t use it now. should be zero)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Y_C_BRAKE'], 'StC_Y_C_BRAKE', '- StC Y high damping for braking the StC (Don\'t use it now. should be zero)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['StC_Z_C_BRAKE'], 'StC_Z_C_BRAKE', '- StC Z high damping for braking the StC (Don\'t use it now. should be zero) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE or when StC_DOF_MODE=3]\n'))

            f.write('---------------------- TLCD --------------------------------------------------- [used only when StC_DOF_MODE=5]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['L_X'], 'L_X', '- X TLCD total length (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['B_X'], 'B_X', '- X TLCD horizontal length (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['area_X'], 'area_X', '- X TLCD cross-sectional area of vertical column (m^2)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['area_ratio_X'], 'area_ratio_X', '- X TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['headLossCoeff_X'], 'headLossCoeff_X', '- X TLCD head loss coeff (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['rho_X'], 'rho_X', '- X TLCD liquid density (kg/m^3)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['L_Y'], 'L_Y', '- Y TLCD total length (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['B_Y'], 'B_Y', '- Y TLCD horizontal length (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['area_Y'], 'area_Y', '- Y TLCD cross-sectional area of vertical column (m^2)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['area_ratio_Y'], 'area_ratio_Y', '- Y TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['headLossCoeff_Y'], 'headLossCoeff_Y', '- Y TLCD head loss coeff (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(stc['rho_Y'], 'rho_Y', '- Y TLCD liquid density (kg/m^3)\n'))

            f.write('---------------------- PRESCRIBED TIME SERIES --------------------------------- [used only when StC_DOF_MODE=6]\n')
            f.write('{:<22} {:<11} {:}'.format(stc['PrescribedForcesCoord'], 'PrescribedForcesCoord', '- Prescribed forces are in global or local coordinates (switch) {1: global; 2: local}\n'))
            f.write('{!s:<22} {:<11} {:}'.format(stc['PrescribedForcesFile'], 'PrescribedForcesFile', '- Time series force and moment (7 columns of time, FX, FY, FZ, MX, MY, MZ)\n'))
            f.write('-------------------------------------------------------------------------------\n')

    # ---- DISCON reader -----------------------------------------------
    @staticmethod
    def _read_discon(sd: dict, base_dir: str) -> Optional[dict]:
        """Read ROSCO DISCON input file if it exists."""
        if not _ROSCO:
            return None

        discon_in_file = os.path.normpath(
            os.path.join(base_dir, sd['DLL_InFile']))

        if not os.path.exists(discon_in_file):
            return None

        discon = read_DISCON(discon_in_file)

        # Additional filename parsing
        discon_dir = os.path.dirname(discon_in_file)
        discon['PerfFileName'] = os.path.abspath(
            os.path.join(discon_dir, discon['PerfFileName']))

        # Try to read rotor performance data
        try:
            pitch_vector, tsr_vector, Cp_table, Ct_table, Cq_table = \
                load_from_txt(discon['PerfFileName'])
            RotorPerformance = ROSCO_turbine.RotorPerformance
            discon['Cp'] = RotorPerformance(Cp_table, pitch_vector, tsr_vector)
            discon['Ct'] = RotorPerformance(Ct_table, pitch_vector, tsr_vector)
            discon['Cq'] = RotorPerformance(Cq_table, pitch_vector, tsr_vector)
            discon['Cp_pitch_initial_rad'] = pitch_vector
            discon['Cp_TSR_initial']       = tsr_vector
            discon['Cp_table'] = Cp_table
            discon['Ct_table'] = Ct_table
            discon['Cq_table'] = Cq_table
        except Exception:
            print('WARNING: Cp table not loaded!')

        discon['v_rated'] = 1.0
        return discon

    # ---- spd_trq reader ----------------------------------------------
    @staticmethod
    def _read_spd_trq(filename: str, base_dir: str) -> dict:
        """Read speed-torque look-up table."""
        spd_trq: Dict[str, Any] = {}
        filepath = os.path.normpath(os.path.join(base_dir, filename))
        with open(filepath) as f:
            spd_trq['header'] = f.readline()
            data = f.readlines()
            spd_trq['RPM']    = [float(line.split()[0]) for line in data]
            spd_trq['Torque'] = [float(line.split()[1]) for line in data]
        return spd_trq

    # ---- spd_trq writer ----------------------------------------------
    @staticmethod
    def _write_spd_trq(spd_trq: dict, out_dir: str) -> None:
        """Write speed-torque look-up table."""
        filepath = os.path.join(out_dir, 'spd_trq.dat')
        with open(filepath, 'w') as f:
            f.write('{:}'.format(spd_trq['header'], '\n'))
            for i in range(len(spd_trq['RPM'])):
                f.write('{:<22f} {:<22f} {:}'.format(
                    spd_trq['RPM'][i], spd_trq['Torque'][i], '\n'))
