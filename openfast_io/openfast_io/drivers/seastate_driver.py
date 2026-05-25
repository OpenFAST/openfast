"""SeaState standalone driver — reads/writes SeaState driver input files (.inp).

The simplest of the standalone drivers — just environmental conditions,
SeaState file reference, wave output options, and simulation timing.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict

from ..io.seastate import SeaStateIO
from ..parsing import bool_read, float_read, int_read, quoted_read


class SeaStateStandaloneDriver:
    """Reads and writes SeaState standalone driver input files (.inp)."""

    def __init__(self) -> None:
        self._seastate = SeaStateIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, inp_path: Path) -> dict:
        """Read a SeaState driver file and the referenced SeaState input.

        Returns a dict with keys:
          ``'SeaStateDriver'`` — driver-level parameters
          ``'SeaState'`` — from SeaStateIO
        """
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr: Dict[str, Any] = {}

        with open(inp_path) as f:
            # --- Header ---
            f.readline()  # title line
            f.readline()  # compatibility line
            dvr['Echo'] = bool_read(f.readline().split()[0])

            # --- Environmental Conditions ---
            f.readline()  # section header
            dvr['Gravity'] = float_read(f.readline().split()[0])
            dvr['WtrDens'] = float_read(f.readline().split()[0])
            dvr['WtrDpth'] = float_read(f.readline().split()[0])
            dvr['MSL2SWL'] = float_read(f.readline().split()[0])

            # --- SeaState ---
            f.readline()  # section header
            dvr['SeaStateInputFile'] = quoted_read(f.readline().split()[0])
            dvr['OutRootName'] = quoted_read(f.readline().split()[0])
            dvr['WrWvKinMod'] = int_read(f.readline().split()[0])
            dvr['NSteps'] = int_read(f.readline().split()[0])
            dvr['TimeInterval'] = float_read(f.readline().split()[0])

            # --- Waves multipoint elevation output ---
            f.readline()  # section header
            dvr['WaveElevSeriesFlag'] = bool_read(f.readline().split()[0])

        result: Dict[str, Any] = {'SeaStateDriver': dvr}

        # --- Delegate to SeaStateIO ---
        ss_file = dvr.get('SeaStateInputFile', '')
        ss_path = os.path.normpath(os.path.join(str(base_dir), ss_file))
        if ss_file and os.path.isfile(ss_path):
            ss_data = self._seastate.read(Path(ss_path), base_dir)
            result.update(ss_data)

        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, inp_path: Path) -> None:
        """Write a SeaState driver file and all referenced module files."""
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr = data['SeaStateDriver']

        with open(inp_path, 'w') as f:
            f.write('Seastate driver file for stand-alone applications.\n')
            f.write('Compatible with SeaState v1.00\n')
            f.write('{!s:<17}{:<19}{:}\n'.format(dvr['Echo'], 'Echo', '- Echo the input file data (flag)'))
            f.write('---------------------- ENVIRONMENTAL CONDITIONS -------------------------------\n')
            f.write('{:<17}{:<19}{:}\n'.format(dvr['Gravity'], 'Gravity', '- Gravity (m/s^2)'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['WtrDens'], 'WtrDens', '- Water density (kg/m^3)'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['WtrDpth'], 'WtrDpth', '- Water depth (m)'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['MSL2SWL'], 'MSL2SWL', '- Offset between still-water level and mean sea level (m) [positive upward]'))
            f.write('---------------------- SEASTATE -----------------------------------------------\n')

            ss_name = dvr.get('SeaStateInputFile', 'SeaState.dat')
            f.write('{:<17}{:<19}{:}\n'.format('"' + ss_name + '"', 'SeaStateInputFile', '- Primary SeaState input file name (quoted string)'))
            f.write('{:<17}{:<19}{:}\n'.format('"' + dvr.get('OutRootName', './seastate') + '"', 'OutRootName', '- The name which prefixes all SeaState generated files (quoted string)'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['WrWvKinMod'], 'WrWvKinMod', '- Write Wave Kinematics? [0: none, 1: (0,0) elevations, 2: complete]'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['NSteps'], 'NSteps', '- Number of time steps in the simulations (-)'))
            f.write('{:<17}{:<19}{:}\n'.format(dvr['TimeInterval'], 'TimeInterval', '- TimeInterval for the simulation (sec)'))
            f.write('---------------------- Waves multipoint elevation output ----------------------\n')
            f.write('{!s:<17}{:<19}{:}\n'.format(dvr['WaveElevSeriesFlag'], 'WaveElevSeriesFlag', '- T/F flag to output the wave elevation field (for movies)'))
            f.write('END of driver input file\n')

        # --- Delegate to SeaStateIO ---
        if 'SeaState' in data:
            ss_path = os.path.normpath(os.path.join(str(base_dir), ss_name))
            self._seastate.write({'SeaState': data['SeaState']}, ss_path, str(base_dir))
