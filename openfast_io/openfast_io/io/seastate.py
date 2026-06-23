"""
SeaStateIO -- read / write SeaState input files.

Produces ``{'SeaState': ss}``
"""
from __future__ import annotations

import os
from typing import Any, Dict, Optional, Callable

from .base import ModuleIO
from ..outlist import emit_outlist
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
)


def _get_outlist(outlist_dict: dict, keys: list) -> list:
    out = []
    for key in keys:
        if key in outlist_dict:
            out.append(outlist_dict[key])
    return out


class SeaStateIO(ModuleIO):
    """Read / write SeaState input files."""

    # ------------------------------------------------------------------
    def read(
        self,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        read_outlist_fn: Optional[Callable] = None,
        **kwargs,
    ) -> dict:
        ss: Dict[str, Any] = {}
        ss_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        f = open(ss_file)
        f.readline()
        f.readline()

        ss['Echo'] = bool_read(f.readline().split()[0])

        # ENVIRONMENTAL CONDITIONS
        f.readline()
        ss['WtrDens'] = float_read(f.readline().split()[0])
        ss['WtrDpth'] = float_read(f.readline().split()[0])
        ss['MSL2SWL'] = float_read(f.readline().split()[0])

        # SPATIAL DISCRETIZATION
        f.readline()
        ss['X_HalfWidth'] = float_read(f.readline().split()[0])
        ss['Y_HalfWidth'] = float_read(f.readline().split()[0])
        ss['Z_Depth']     = float_read(f.readline().split()[0])
        ss['NX']          = int_read(f.readline().split()[0])
        ss['NY']          = int_read(f.readline().split()[0])
        ss['NZ']          = int_read(f.readline().split()[0])

        # WAVES
        f.readline()
        ss['WaveMod']       = int_read(f.readline().split()[0])
        ss['WaveStMod']     = int_read(f.readline().split()[0])
        ss['WvCrntMod']     = int_read(f.readline().split()[0])
        ss['WaveTMax']      = float_read(f.readline().split()[0])
        ss['WaveDT']        = float_read(f.readline().split()[0])
        ss['WaveHs']        = float_read(f.readline().split()[0])
        ss['WaveTp']        = float_read(f.readline().split()[0])
        ss['WavePkShp']     = float_read(f.readline().split()[0])
        ss['WvLowCOff']     = float_read(f.readline().split()[0])
        ss['WvHiCOff']      = float_read(f.readline().split()[0])
        ss['WaveDir']       = float_read(f.readline().split()[0])
        ss['WaveDirMod']    = int_read(f.readline().split()[0])
        ss['WaveDirSpread'] = float_read(f.readline().split()[0])
        ss['WaveNDir']      = int_read(f.readline().split()[0])
        ss['WaveDirRange']  = float_read(f.readline().split()[0])
        ss['WaveSeed1']     = int_read(f.readline().split()[0])
        ss['WaveSeed2']     = int_read(f.readline().split()[0])
        ss['WaveNDAmp']     = bool_read(f.readline().split()[0])
        ss['WvKinFile']     = quoted_read(f.readline().split()[0])

        # 2ND-ORDER WAVES
        f.readline()
        ss['WvDiffQTF']  = bool_read(f.readline().split()[0])
        ss['WvSumQTF']   = bool_read(f.readline().split()[0])
        ss['WvLowCOffD'] = float_read(f.readline().split()[0])
        ss['WvHiCOffD']  = float_read(f.readline().split()[0])
        ss['WvLowCOffS'] = float_read(f.readline().split()[0])
        ss['WvHiCOffS']  = float_read(f.readline().split()[0])

        # CONSTRAINED WAVE
        f.readline()
        ss['ConstWaveMod'] = int_read(f.readline().split()[0])
        ss['CrestHmax']    = float_read(f.readline().split()[0])
        ss['CrestTime']    = float_read(f.readline().split()[0])
        ss['CrestXi']      = float_read(f.readline().split()[0])
        ss['CrestYi']      = float_read(f.readline().split()[0])

        # CURRENT
        f.readline()
        ss['CurrMod']   = int_read(f.readline().split()[0])
        ss['CurrSSV0']  = float_read(f.readline().split()[0])
        ss['CurrSSDir'] = float_read(f.readline().split()[0])
        ss['CurrNSRef'] = float_read(f.readline().split()[0])
        ss['CurrNSV0']  = float_read(f.readline().split()[0])
        ss['CurrNSDir'] = float_read(f.readline().split()[0])
        ss['CurrDIV']   = float_read(f.readline().split()[0])
        ss['CurrDIDir'] = float_read(f.readline().split()[0])

        # MacCamy-Fuchs
        f.readline()
        ss['MCFD'] = float_read(f.readline().split()[0])

        # OUTPUT
        f.readline()
        ss['SeaStSum']   = bool_read(f.readline().split()[0])
        ss['OutSwtch']   = int_read(f.readline().split()[0])
        ss['OutFmt']     = quoted_read(f.readline().split()[0])
        ss['OutSFmt']    = quoted_read(f.readline().split()[0])
        ss['NWaveElev']  = int_read(f.readline().split()[0])
        ss['WaveElevxi'] = [float_read(idx.strip()) for idx in f.readline().split('WaveElevxi')[0].replace(',', ' ').split()]
        ss['WaveElevyi'] = [float_read(idx.strip()) for idx in f.readline().split('WaveElevyi')[0].replace(',', ' ').split()]
        ss['NWaveKin']   = int_read(f.readline().split()[0])
        if ss['NWaveKin']:
            ss['WaveKinxi'] = [float_read(idx.strip()) for idx in f.readline().split('WaveKinxi')[0].replace(',', ' ').split()]
            ss['WaveKinyi'] = [float_read(idx.strip()) for idx in f.readline().split('WaveKinyi')[0].replace(',', ' ').split()]
            ss['WaveKinzi'] = [float_read(idx.strip()) for idx in f.readline().split('WaveKinzi')[0].replace(',', ' ').split()]
        else:
            [f.readline() for _ in range(3)]
            ss['WaveKinxi'] = [0]
            ss['WaveKinyi'] = [0]
            ss['WaveKinzi'] = [0]

        # Outlist  (legacy uses read_outlist_freeForm)
        f.readline()
        if read_outlist_fn is not None:
            read_outlist_fn(f, 'SeaState')

        f.close()
        return {'SeaState': ss}

    # ------------------------------------------------------------------
    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        **kwargs,
    ) -> None:
        ss = data['SeaState']

        with open(file_path, 'w') as f:
            f.write('------- SeaState Input File --------------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('{!s:<22} {:<11} {:}'.format(ss['Echo'], 'Echo', '- Echo the input file data (flag)\n'))

            f.write('---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(ss['WtrDens'], 'WtrDens', '- Water density (kg/m^3)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WtrDpth'], 'WtrDpth', '- Water depth (m) relative to MSL\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['MSL2SWL'], 'MSL2SWL', '- Offset between SWL and MSL (m)\n'))

            f.write('---------------------- SPATIAL DISCRETIZATION ---------------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(ss['X_HalfWidth'], 'X_HalfWidth', '- Half-width of X domain (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['Y_HalfWidth'], 'Y_HalfWidth', '- Half-width of Y domain (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['Z_Depth'], 'Z_Depth', '- Depth of Z domain (m)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['NX'], 'NX', '- Number of nodes in half X-direction (-)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['NY'], 'NY', '- Number of nodes in half Y-direction (-)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['NZ'], 'NZ', '- Number of nodes in Z direction (-)\n'))

            f.write('---------------------- WAVES ---------------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(ss['WaveMod'], 'WaveMod', '- Incident wave kinematics model (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['WaveStMod'], 'WaveStMod', '- Wave stretching model (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['WvCrntMod'], 'WvCrntMod', '- Wave-current model (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveTMax'], 'WaveTMax', '- Analysis time for wave calculations (sec)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveDT'], 'WaveDT', '- Time step for wave calculations (sec)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveHs'], 'WaveHs', '- Significant wave height (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveTp'], 'WaveTp', '- Peak-spectral period (sec)\n'))
            # WavePkShp special handling
            wpks = ss['WavePkShp']
            if isinstance(wpks, float) and wpks == 0.0:
                wpks = 'Default'
            f.write('{:<22} {:<11} {:}'.format(wpks, 'WavePkShp', '- Peak-shape parameter (-) or DEFAULT\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvLowCOff'], 'WvLowCOff', '- Low cut-off frequency (rad/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvHiCOff'], 'WvHiCOff', '- High cut-off frequency (rad/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveDir'], 'WaveDir', '- Wave propagation heading (deg)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['WaveDirMod'], 'WaveDirMod', '- Directional spreading function (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveDirSpread'], 'WaveDirSpread', '- Wave direction spreading coeff (-)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['WaveNDir'], 'WaveNDir', '- Number of wave directions (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WaveDirRange'], 'WaveDirRange', '- Range of wave directions (deg)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['WaveSeed1'], 'WaveSeed(1)', '- First random seed (-)\n'))
            try:
                seed2 = int(ss['WaveSeed2'])
                f.write('{:<22d} {:<11} {:}'.format(ss['WaveSeed2'], 'WaveSeed(2)', '- Second random seed (-)\n'))
            except (ValueError, TypeError):
                f.write('{!s:<22} {:<11} {:}'.format(ss['WaveSeed2'], 'WaveSeed(2)', '- Second random seed (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ss['WaveNDAmp'], 'WaveNDAmp', '- Normally distributed amplitudes (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + ss['WvKinFile'] + '"', 'WvKinFile', '- External wave data file root name\n'))

            f.write('---------------------- 2ND-ORDER WAVES -----------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ss['WvDiffQTF'], 'WvDiffQTF', '- Difference-frequency 2nd-order (flag)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ss['WvSumQTF'], 'WvSumQTF', '- Summation-frequency 2nd-order (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvLowCOffD'], 'WvLowCOffD', '- Low freq cutoff for diff (rad/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvHiCOffD'], 'WvHiCOffD', '- High freq cutoff for diff (rad/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvLowCOffS'], 'WvLowCOffS', '- Low freq cutoff for sum (rad/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['WvHiCOffS'], 'WvHiCOffS', '- High freq cutoff for sum (rad/s)\n'))

            f.write('---------------------- CONSTRAINED WAVES ----------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(ss['ConstWaveMod'], 'ConstWaveMod', '- Constrained wave model (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CrestHmax'], 'CrestHmax', '- Crest height (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CrestTime'], 'CrestTime', '- Time of crest (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CrestXi'], 'CrestXi', '- X-position of crest (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CrestYi'], 'CrestYi', '- Y-position of crest (m)\n'))

            f.write('---------------------- CURRENT -------------------------------------------------\n')
            f.write('{:<22d} {:<11} {:}'.format(ss['CurrMod'], 'CurrMod', '- Current profile model (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrSSV0'], 'CurrSSV0', '- Sub-surface current velocity (m/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrSSDir'], 'CurrSSDir', '- Sub-surface current heading (deg)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrNSRef'], 'CurrNSRef', '- Near-surface reference depth (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrNSV0'], 'CurrNSV0', '- Near-surface current velocity (m/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrNSDir'], 'CurrNSDir', '- Near-surface current heading (deg)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrDIV'], 'CurrDIV', '- Depth-independent current velocity (m/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(ss['CurrDIDir'], 'CurrDIDir', '- Depth-independent current heading (deg)\n'))

            f.write('---------------------- MacCamy-Fuchs Diffraction Model -------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(ss['MCFD'], 'MCFD', '- MacCamy-Fuchs member radius\n'))

            f.write('---------------------- OUTPUT --------------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ss['SeaStSum'], 'SeaStSum', '- Output a summary file [flag]\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['OutSwtch'], 'OutSwtch', '- Output channels (switch)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ss['OutFmt'], 'OutFmt', '- Output format\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ss['OutSFmt'], 'OutSFmt', '- Output format for headers\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['NWaveElev'], 'NWaveElev', '- Number of wave elevation output points (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{float(v):f}' for v in ss['WaveElevxi']]), 'WaveElevxi', '- xi-coordinates for wave elevation output (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(", ".join([f'{float(v):f}' for v in ss['WaveElevyi']]), 'WaveElevyi', '- yi-coordinates for wave elevation output (m)\n'))
            f.write('{:<22d} {:<11} {:}'.format(ss['NWaveKin'], 'NWaveKin', '- Number of wave kinematics output points (-)\n'))

            if ss['NWaveKin'] > 0:
                f.write('{:<22} {:<11} {:}'.format(", ".join([f'{v:f}' for v in ss['WaveKinxi']]), 'WaveKinxi', '- xi-coordinates for wave kinematics (m)\n'))
                f.write('{:<22} {:<11} {:}'.format(", ".join([f'{v:f}' for v in ss['WaveKinyi']]), 'WaveKinyi', '- yi-coordinates for wave kinematics (m)\n'))
                f.write('{:<22} {:<11} {:}'.format(", ".join([f'{v:f}' for v in ss['WaveKinzi']]), 'WaveKinzi', '- zi-coordinates for wave kinematics (m)\n'))
            else:
                f.write('{:<11} {:}'.format('WaveKinxi', '- xi-coordinates for wave kinematics (m)\n'))
                f.write('{:<11} {:}'.format('WaveKinyi', '- yi-coordinates for wave kinematics (m)\n'))
                f.write('{:<11} {:}'.format('WaveKinzi', '- zi-coordinates for wave kinematics (m)\n'))

            f.write('---------------------- OUTPUT CHANNELS -----------------------------------------\n')
            if outlist is not None:
                emit_outlist(f, outlist, 'SeaState')
            f.write('END of output channels and end of file.\n')
