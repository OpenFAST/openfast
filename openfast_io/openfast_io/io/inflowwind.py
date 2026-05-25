"""InflowWind module IO — reads and writes InflowWind input files.

Extracted from FAST_reader.py and FAST_writer.py.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from .base import ModuleIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
)


def _get_outlist(outlist_dict, channel_list):
    def loop_dict(vartree, outlist_i):
        for var in vartree.keys():
            if isinstance(vartree[var], dict):
                loop_dict(vartree[var], outlist_i)
            else:
                if vartree[var]:
                    outlist_i.append(var)
        return outlist_i

    if not channel_list:
        channel_list = outlist_dict.keys()
    outlist = []
    for var in channel_list:
        var = var.replace(' ', '')
        outlist_i = loop_dict(outlist_dict[var], [])
        if outlist_i:
            outlist.append(sorted(outlist_i))
    return outlist


class InflowWindIO(ModuleIO):
    """Reads and writes InflowWind v3.01+ input files.

    read() returns::

        {'InflowWind': { ... all params ... }}

    write() accepts data with key 'InflowWind'.
    """

    def read(self, file_path: Path, base_dir: Path, *,
             outlist: dict | None = None,
             read_outlist_fn=None) -> dict:
        ifw = {}
        file_path = str(file_path)
        inflow_dir = os.path.dirname(file_path)

        f = open(file_path)
        f.readline(); f.readline(); f.readline()

        # Header
        ifw['Echo']           = bool_read(f.readline().split()[0])
        ifw['WindType']       = int(f.readline().split()[0])
        ifw['PropagationDir'] = float_read(f.readline().split()[0])
        ifw['VFlowAng']       = float_read(f.readline().split()[0])
        ifw['VelInterpCubic'] = bool_read(f.readline().split()[0])
        ifw['NWindVel']       = int(f.readline().split()[0])
        ifw['WindVxiList']    = [idx.strip() for idx in f.readline().split('WindVxiList')[0].split(',')]
        ifw['WindVyiList']    = [idx.strip() for idx in f.readline().split('WindVyiList')[0].split(',')]
        ifw['WindVziList']    = [idx.strip() for idx in f.readline().split('WindVziList')[0].split(',')]

        # Steady Wind
        f.readline()
        ifw['HWindSpeed'] = float_read(f.readline().split()[0])
        ifw['RefHt']      = float_read(f.readline().split()[0])
        ifw['PLExp']      = float_read(f.readline().split()[0])

        # Uniform Wind
        f.readline()
        ifw['FileName_Uni'] = os.path.join(inflow_dir, quoted_read(f.readline().split()[0]))
        ifw['RefHt_Uni']    = float_read(f.readline().split()[0])
        ifw['RefLength']    = float_read(f.readline().split()[0])

        # TurbSim FF
        f.readline()
        ifw['FileName_BTS'] = os.path.join(inflow_dir, quoted_read(f.readline().split()[0]))

        # Bladed FF
        f.readline()
        ifw['FileNameRoot'] = os.path.join(inflow_dir, quoted_read(f.readline().split()[0]))
        ifw['TowerFile']    = bool_read(f.readline().split()[0])

        # HAWC
        f.readline()
        ifw['FileName_u']   = os.path.normpath(os.path.join(inflow_dir, quoted_read(f.readline().split()[0])))
        ifw['FileName_v']   = os.path.normpath(os.path.join(inflow_dir, quoted_read(f.readline().split()[0])))
        ifw['FileName_w']   = os.path.normpath(os.path.join(inflow_dir, quoted_read(f.readline().split()[0])))
        ifw['nx']           = int(f.readline().split()[0])
        ifw['ny']           = int(f.readline().split()[0])
        ifw['nz']           = int(f.readline().split()[0])
        ifw['dx']           = float_read(f.readline().split()[0])
        ifw['dy']           = float_read(f.readline().split()[0])
        ifw['dz']           = float_read(f.readline().split()[0])
        ifw['RefHt_Hawc']   = float_read(f.readline().split()[0])

        # HAWC scaling
        f.readline()
        ifw['ScaleMethod']  = int(f.readline().split()[0])
        ifw['SFx']          = float_read(f.readline().split()[0])
        ifw['SFy']          = float_read(f.readline().split()[0])
        ifw['SFz']          = float_read(f.readline().split()[0])
        ifw['SigmaFx']      = float_read(f.readline().split()[0])
        ifw['SigmaFy']      = float_read(f.readline().split()[0])
        ifw['SigmaFz']      = float_read(f.readline().split()[0])

        # HAWC mean profile
        f.readline()
        ifw['URef']         = float_read(f.readline().split()[0])
        ifw['WindProfile']  = int(f.readline().split()[0])
        ifw['PLExp_Hawc']   = float_read(f.readline().split()[0])
        ifw['Z0']           = float_read(f.readline().split()[0])
        ifw['XOffset']      = float_read(f.readline().split()[0])

        # LIDAR
        f.readline()
        ifw['SensorType']          = int(f.readline().split()[0])
        ifw['NumPulseGate']        = int(f.readline().split()[0])
        ifw['PulseSpacing']        = float_read(f.readline().split()[0])
        ifw['NumBeam']             = int(f.readline().split()[0])
        ifw['FocalDistanceX']      = [idx.strip() for idx in f.readline().split('FocalDistanceX')[0].split(',')]
        ifw['FocalDistanceY']      = [idx.strip() for idx in f.readline().split('FocalDistanceY')[0].split(',')]
        ifw['FocalDistanceZ']      = [idx.strip() for idx in f.readline().split('FocalDistanceZ')[0].split(',')]
        ifw['RotorApexOffsetPos']  = [idx.strip() for idx in f.readline().split('RotorApexOffsetPos')[0].split(',')]
        ifw['URefLid']             = float_read(f.readline().split()[0])
        ifw['MeasurementInterval'] = float_read(f.readline().split()[0])
        ifw['LidRadialVel']        = bool_read(f.readline().split()[0])
        ifw['ConsiderHubMotion']   = int(f.readline().split()[0])

        # Output
        f.readline()
        ifw['SumPrint'] = bool_read(f.readline().split()[0])

        # OutList
        f.readline()
        if read_outlist_fn is not None and outlist is not None:
            read_outlist_fn(f, 'InflowWind')
        else:
            line = f.readline()
            while line and 'END' not in line.split('!')[0].upper()[:3]:
                line = f.readline()

        f.close()
        return {'InflowWind': ifw}

    def write(self, data: dict, file_path: Path, base_dir: Path, *,
              naming_out: str = 'openfast', outlist: dict | None = None) -> None:
        ifw = data['InflowWind']
        f = open(str(file_path), 'w')

        f.write('------- InflowWind INPUT FILE -------------------------------------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------------------------------------------------------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(ifw['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['WindType'], 'WindType', '- switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['PropagationDir'], 'PropagationDir', '- Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees) (not used for native Bladed format WindType=7)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['VFlowAng'], 'VFlowAng', '- Upflow angle (degrees) (not used for native Bladed format WindType=7)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ifw['VelInterpCubic'], 'VelInterpCubic', '- Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['NWindVel'], 'NWindVel', '- Number of points to output the wind velocity    (0 to 9)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['WindVxiList'], dtype=str)), 'WindVxiList', '- List of coordinates in the inertial X direction (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['WindVyiList'], dtype=str)), 'WindVyiList', '- List of coordinates in the inertial Y direction (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['WindVziList'], dtype=str)), 'WindVziList', '- List of coordinates in the inertial Z direction (m)\n'))
        f.write('================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================\n')
        f.write('{:<22} {:<11} {:}'.format(ifw['HWindSpeed'], 'HWindSpeed', '- Horizontal wind speed                            (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['RefHt'], 'RefHt', '- Reference height for horizontal wind speed      (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['PLExp'], 'PLExp', '- Power law exponent                              (-)\n'))
        f.write('================== Parameters for Uniform wind file   [used only for WindType = 2] ============================\n')
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileName_Uni'] + '"', 'FileName_Uni', '- Filename of time series data for uniform wind field.      (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['RefHt_Uni'], 'RefHt_Uni', '- Reference height for horizontal wind speed                (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['RefLength'], 'RefLength', '- Reference length for linear horizontal and vertical sheer (-)\n'))
        f.write('================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============\n')
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileName_BTS'] + '"', 'FileName_BTS', '- Name of the Full field wind file to use (.bts)\n'))
        f.write('================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4 or WindType = 7] =========\n')
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileNameRoot'] + '"', 'FileNameRoot', '- WindType=4: Rootname of the full-field wind file to use (.wnd, .sum); WindType=7: name of the intermediate file with wind scaling values\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ifw['TowerFile'], 'TowerFile', '- Have tower file (.twr) (flag) ignored when WindType = 7\n'))
        f.write('================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================\n')
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileName_u'] + '"', 'FileName_u', '- name of the file containing the u-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileName_v'] + '"', 'FileName_v', '- name of the file containing the v-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ifw['FileName_w'] + '"', 'FileName_w', '- name of the file containing the w-component fluctuating wind (.bin)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['nx'], 'nx', '- number of grids in the x direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['ny'], 'ny', '- number of grids in the y direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['nz'], 'nz', '- number of grids in the z direction (in the 3 files above) (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['dx'], 'dx', '- distance (in meters) between points in the x direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['dy'], 'dy', '- distance (in meters) between points in the y direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['dz'], 'dz', '- distance (in meters) between points in the z direction    (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['RefHt_Hawc'], 'RefHt_Hawc', '- reference height; the height (in meters) of the vertical center of the grid (m)\n'))
        f.write('-------------   Scaling parameters for turbulence   ---------------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ifw['ScaleMethod'], 'ScaleMethod', '- Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SFx'], 'SFx', '- Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SFy'], 'SFy', '- Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SFz'], 'SFz', '- Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SigmaFx'], 'SigmaFx', '- Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SigmaFy'], 'SigmaFy', '- Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['SigmaFz'], 'SigmaFz', '- Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]\n'))
        f.write('-------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ifw['URef'], 'URef', '- Mean u-component wind speed at the reference height (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['WindProfile'], 'WindProfile', '- Wind profile type (0=constant;1=logarithmic,2=power law)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['PLExp_Hawc'], 'PLExp_Hawc', '- Power law exponent (-) (used for PL wind profile type only)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['Z0'], 'Z0', '- Surface roughness length (m) (used for LG wind profile type only)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['XOffset'], 'XOffset', '- Initial offset in +x direction (shift of wind box) (-)\n'))
        f.write('-------------   LIDAR Parameters   --------------------------------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ifw['SensorType'], 'SensorType', '- Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['NumPulseGate'], 'NumPulseGate', '- Number of lidar measurement gates (used when SensorType = 3)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['PulseSpacing'], 'PulseSpacing', '- Distance between range gates (m) (used when SensorType = 3)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['NumBeam'], 'NumBeam', '- Number of lidar measurement beams (0-5)(used when SensorType = 1)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['FocalDistanceX'], dtype=str)), 'FocalDistanceX', '- Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['FocalDistanceY'], dtype=str)), 'FocalDistanceY', '- Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['FocalDistanceZ'], dtype=str)), 'FocalDistanceZ', '- Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ifw['RotorApexOffsetPos'], dtype=str)), 'RotorApexOffsetPos', '- Offset of the lidar from hub height (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['URefLid'], 'URefLid', '- Reference average wind speed for the lidar[m/s]\n'))
        f.write('{:<22} {:<11} {:}'.format(ifw['MeasurementInterval'], 'MeasurementInterval', '- Time between each measurement [s]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ifw['LidRadialVel'], 'LidRadialVel', "- TRUE => return radial component, FALSE => return 'x' direction estimate\n"))
        f.write('{:<22} {:<11} {:}'.format(ifw['ConsiderHubMotion'], 'ConsiderHubMotion', "- Flag whether to consider the hub motion's impact on Lidar measurements\n"))
        f.write('====================== OUTPUT ==================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(ifw['SumPrint'], 'SumPrint', '- Print summary data to <RootName>.IfW.sum (flag)\n'))
        f.write('OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        if outlist is not None:
            ol = _get_outlist(outlist, ['InflowWind'])
            for channel_list in ol:
                for ch in channel_list:
                    f.write('"' + ch + '"\n')

        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
        f.write('---------------------------------------------------------------------------------------\n')

        f.flush()
        os.fsync(f)
        f.close()
