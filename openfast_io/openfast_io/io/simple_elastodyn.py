"""
SimpleElastoDynIO – read / write Simplified ElastoDyn (SED) input files.

Produces ``{'SimpleElastoDyn': sed}``
"""
from __future__ import annotations

import os
from typing import Any, Callable, Dict, Optional

from .base import ModuleIO
from ..parsing import bool_read, float_read, int_read


class SimpleElastoDynIO(ModuleIO):
    """Read / write SimpleElastoDyn input files."""

    # ------------------------------------------------------------------
    # READ
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
        sed: Dict[str, Any] = {}
        f = open(file_path)

        f.readline()
        f.readline()
        f.readline()
        sed['Echo'] = bool_read(f.readline().split()[0])
        sed['IntMethod'] = int_read(f.readline().split()[0])
        sed['DT'] = float_read(f.readline().split()[0])

        # Degrees of Freedom
        f.readline()
        sed['GenDOF'] = bool_read(f.readline().split()[0])
        sed['YawDOF'] = bool_read(f.readline().split()[0])

        # Initial Conditions
        f.readline()
        sed['Azimuth'] = float_read(f.readline().split()[0])
        sed['BlPitch'] = float_read(f.readline().split()[0])
        sed['RotSpeed'] = float_read(f.readline().split()[0])
        sed['NacYaw'] = float_read(f.readline().split()[0])
        sed['PtfmPitch'] = float_read(f.readline().split()[0])

        # Turbine Configuration
        f.readline()
        sed['NumBl'] = int_read(f.readline().split()[0])
        sed['TipRad'] = float_read(f.readline().split()[0])
        sed['HubRad'] = float_read(f.readline().split()[0])
        sed['PreCone'] = float_read(f.readline().split()[0])
        sed['OverHang'] = float_read(f.readline().split()[0])
        sed['ShftTilt'] = float_read(f.readline().split()[0])
        sed['Twr2Shft'] = float_read(f.readline().split()[0])
        sed['TowerHt'] = float_read(f.readline().split()[0])

        # Mass and Inertia
        f.readline()
        sed['RotIner'] = float_read(f.readline().split()[0])
        sed['GenIner'] = float_read(f.readline().split()[0])

        # Drivetrain
        f.readline()
        sed['GBoxRatio'] = float_read(f.readline().split()[0])

        # Output
        f.readline()
        f.readline()

        # Read output list — route into the shared registry (mirrors baseline)
        if read_outlist_fn is not None and outlist is not None:
            read_outlist_fn(f, 'SimpleElastoDyn')
        else:
            sed['_outlist'] = {}
            data = f.readline()
            while data.split().__len__() == 0:
                data = f.readline()
            while data.split()[0] != 'END':
                if data.find('"') >= 0:
                    channels = data.split('"')
                    channel_list = channels[1].split(',')
                else:
                    row_string = data.split(',')
                    if len(row_string) == 1:
                        channel_list = row_string[0].split('\n')[0]
                    else:
                        channel_list = row_string
                if isinstance(channel_list, list):
                    for ch in channel_list:
                        ch = ch.strip()
                        if ch:
                            sed['_outlist'][ch] = True
                else:
                    ch = channel_list.strip()
                    if ch:
                        sed['_outlist'][ch] = True
                data = f.readline()

        f.close()
        return {'SimpleElastoDyn': sed}

    # ------------------------------------------------------------------
    # WRITE
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
        sed = data['SimpleElastoDyn']
        with open(file_path, 'w') as f:
            f.write('------- SIMPLIFIED ELASTODYN INPUT FILE ----------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sed['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['IntMethod'], 'IntMethod', '- Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['DT'], 'DT', '- Integration time step (s)\n'))
            f.write('---------------------- DEGREES OF FREEDOM --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sed['GenDOF'], 'GenDOF', '- Generator DOF (flag)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sed['YawDOF'], 'YawDOF', '- Yaw degree of freedom -- controlled by controller (flag)\n'))
            f.write('---------------------- INITIAL CONDITIONS --------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sed['Azimuth'], 'Azimuth', '- Initial azimuth angle for blades (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['BlPitch'], 'BlPitch', '- Blades initial pitch (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['RotSpeed'], 'RotSpeed', '- Initial or fixed rotor speed (rpm)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['NacYaw'], 'NacYaw', '- Initial or fixed nacelle-yaw angle (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['PtfmPitch'], 'PtfmPitch', '- Fixed pitch tilt rotational displacement of platform (degrees)\n'))
            f.write('---------------------- TURBINE CONFIGURATION -----------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sed['NumBl'], 'NumBl', '- Number of blades (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['TipRad'], 'TipRad', '- The distance from the rotor apex to the blade tip (meters)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['HubRad'], 'HubRad', '- The distance from the rotor apex to the blade root (meters)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['PreCone'], 'PreCone', '- Blades cone angle (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['OverHang'], 'OverHang', '- Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['ShftTilt'], 'ShftTilt', '- Rotor shaft tilt angle (degrees)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['Twr2Shft'], 'Twr2Shft', '- Vertical distance from the tower-top to the rotor shaft (meters)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['TowerHt'], 'TowerHt', '- Height of tower above ground level [onshore] or MSL [offshore] (meters)\n'))
            f.write('---------------------- MASS AND INERTIA ----------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sed['RotIner'], 'RotIner', '- Rot inertia about rotor axis [blades + hub] (kg m^2)\n'))
            f.write('{:<22} {:<11} {:}'.format(sed['GenIner'], 'GenIner', '- Generator inertia about HSS (kg m^2)\n'))
            f.write('---------------------- DRIVETRAIN ----------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(sed['GBoxRatio'], 'GBoxRatio', '- Gearbox ratio (-)\n'))
            f.write('---------------------- OUTPUT --------------------------------------------------\n')
            f.write('                   OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')
            out_channels = sed.get('_outlist', {})
            if outlist and 'SimpleElastoDyn' in outlist:
                out_channels = outlist['SimpleElastoDyn']
            for ch in out_channels:
                f.write('"' + ch + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
            f.write('---------------------------------------------------------------------------------------\n')
