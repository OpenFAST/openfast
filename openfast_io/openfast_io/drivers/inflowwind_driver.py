"""InflowWind standalone driver — reads/writes InflowWind driver input files (.inp).

Uses ``=====`` section separators and ``--`` comment syntax (different from most
OpenFAST drivers). Handles file conversion options, interpolation tests, points
file input, output grid, and VTK slice output.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict

from ..io.inflowwind import InflowWindIO
from ..parsing import bool_read, float_read, int_read, quoted_read


class InflowWindStandaloneDriver:
    """Reads and writes InflowWind standalone driver input files (.inp)."""

    def __init__(self) -> None:
        self._inflowwind = InflowWindIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, inp_path: Path) -> dict:
        """Read an InflowWind driver file and the referenced primary input.

        Returns a dict with keys:
          ``'InflowWindDriver'`` — driver-level parameters
          ``'InflowWind'`` — from InflowWindIO
        """
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr: Dict[str, Any] = {}

        with open(inp_path) as f:
            # --- Header ---
            f.readline()  # title line 1
            f.readline()  # title line 2
            dvr['Echo'] = bool_read(f.readline().split()[0])

            # --- InflowWind input filename ---
            f.readline()  # separator
            dvr['IfWFileName'] = quoted_read(f.readline().split()[0])

            # --- File Conversion Options ---
            f.readline()  # separator
            dvr['WrHAWC'] = bool_read(f.readline().split()[0])
            dvr['WrBladed'] = bool_read(f.readline().split()[0])
            dvr['WrVTK'] = bool_read(f.readline().split()[0])
            dvr['WrUniform'] = bool_read(f.readline().split()[0])

            # --- Interpolation Options ---
            f.readline()  # separator
            dvr['NumTSteps'] = int_read(f.readline().split()[0])
            dvr['TStart'] = float_read(f.readline().split()[0])
            dvr['DT'] = float_read(f.readline().split()[0])
            dvr['Summary'] = bool_read(f.readline().split()[0])
            dvr['SummaryFile'] = bool_read(f.readline().split()[0])
            dvr['BoxExceedAllow'] = bool_read(f.readline().split()[0])

            # --- Points file input ---
            f.readline()  # separator
            dvr['PointsFlag'] = bool_read(f.readline().split()[0])
            dvr['PointsFileName'] = quoted_read(f.readline().split()[0])
            dvr['CalcAccel'] = bool_read(f.readline().split()[0])

            # --- Output grid ---
            f.readline()  # separator
            dvr['WindGrid'] = bool_read(f.readline().split()[0])
            dvr['GridCtrCoord'] = _read_csv_vec(f.readline())
            dvr['GridDXYZ'] = _read_csv_vec(f.readline())
            dvr['GridNXYZ'] = _read_csv_vec(f.readline())

            # --- Output VTK slices ---
            f.readline()  # separator
            dvr['NOutWindXY'] = int_read(f.readline().split()[0])
            dvr['OutWindZ'] = float_read(f.readline().split()[0])

        result: Dict[str, Any] = {'InflowWindDriver': dvr}

        # --- Delegate to InflowWindIO ---
        ifw_file = dvr.get('IfWFileName', '')
        ifw_path = os.path.normpath(os.path.join(str(base_dir), ifw_file))
        if ifw_file and os.path.isfile(ifw_path):
            ifw_data = self._inflowwind.read(Path(ifw_path), base_dir)
            result.update(ifw_data)

        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, inp_path: Path) -> None:
        """Write an InflowWind driver file and all referenced module files."""
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr = data['InflowWindDriver']

        with open(inp_path, 'w') as f:
            f.write('InflowWind driver input file\n')
            f.write('InflowWind driver input file. V1.00\n')
            f.write('{!s:<8}{:}\n'.format(dvr['Echo'], 'echo (flag)'))
            f.write('===============================================================================\n')

            ifw_name = dvr.get('IfWFileName', 'InflowWind.dat')
            f.write('{:<16}{:}\n'.format('"' + ifw_name + '"', 'IfWFileName -- IfW input filename (-)'))
            f.write('===================== File Conversion Options =================================\n')
            f.write(' {!s:<14}{:<15}{:}\n'.format(dvr['WrHAWC'], 'WrHAWC', '-- Convert all data to HAWC2 format? (flag)'))
            f.write(' {!s:<14}{:<15}{:}\n'.format(dvr['WrBladed'], 'WrBladed', '-- Convert all data to Bladed format? (flag)'))
            f.write(' {!s:<14}{:<15}{:}\n'.format(dvr['WrVTK'], 'WrVTK', '-- Convert all data to VTK format? (flag)'))
            f.write(' {!s:<14}{:<15}{:}\n'.format(dvr.get('WrUniform', False), 'WrUniform', '-- Convert data to Uniform wind format? (flag)'))
            f.write('=====================  Tests of Interpolation Options =========================\n')
            f.write('{:<15} {:<15}{:}\n'.format(dvr['NumTSteps'], 'NumTSteps', '-- number of timesteps to run (DEFAULT for all) (-)'))
            f.write('{:<15} {:<15}{:}\n'.format(dvr['TStart'], 'TStart', '-- Start time (s)'))
            f.write('{:<15} {:<15}{:}\n'.format(dvr['DT'], 'DT', '-- timestep size for driver to take (s, or DEFAULT for what the file contains)'))
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr['Summary'], 'Summary', '-- Summarize the data extents in the windfile (flag)'))
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr['SummaryFile'], 'SummaryFile', '-- Write summary to file .dvr.sum (flag)'))
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr.get('BoxExceedAllow', False), 'BoxExceedAllow', '-- Allow point sampling outside grid'))
            f.write('---- Points file input (output given as PointsFileName.Velocity.dat) --------\n')
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr['PointsFlag'], 'PointsFileName', '-- read in a list of points from a file (flag)'))
            f.write('{:<15} {:<15}{:}\n'.format('"' + dvr['PointsFileName'] + '"', 'PointsFileName', '-- name of points file (-)'))
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr.get('CalcAccel', False), 'CalcAccel', '-- calculate and output acceleration at points'))
            f.write('---- Output grid (Points below ground will simply be ignored) ---------------\n')
            f.write('{!s:<15} {:<15}{:}\n'.format(dvr['WindGrid'], 'WindGrid', '-- report wind data at set of Y,Z coordinates (flag)'))
            f.write('{:<15} {:<15}{:}\n'.format(','.join(str(v) for v in dvr['GridCtrCoord']), 'GridCtrCoord', '-- coordinates of center of grid (m)'))
            f.write('{:<15} {:<15}{:}\n'.format(','.join(str(v) for v in dvr['GridDXYZ']), 'GridDX,GridDY,GridDZ', '-- Stepsize of grid (m)'))
            f.write('{:<15} {:<15}{:}\n'.format(','.join(str(v) for v in dvr['GridNXYZ']), 'GridNX,GridNY,GridNZ', '-- number of grid points in X, Y and Z directions (-)'))
            f.write('----  Output VTK slices  ------------------------------------------------------\n')
            f.write('{:>4}            {:<14}{:}\n'.format(dvr['NOutWindXY'], 'NOutWindXY', '-- Number of XY planes for output (-) [0 to 9]'))
            f.write('{:>4}           {:<14}{:}\n'.format(dvr['OutWindZ'], 'OutWindZ', '-- Z coordinates of XY planes for output (m)'))
            f.write('END of driver input file\n')

        # --- Copy referenced data files (Points.inp, wind files, etc.) ---
        # These are external data files the driver references but doesn't generate.

        # --- Delegate to InflowWindIO ---
        if 'InflowWind' in data:
            ifw_path = os.path.normpath(os.path.join(str(base_dir), ifw_name))
            self._inflowwind.write({'InflowWind': data['InflowWind']}, Path(ifw_path), base_dir)


def _read_csv_vec(line: str) -> list:
    """Parse a comma-separated vector from a line like '0,0,150  name ...'."""
    raw = line.split()[0]
    return [float_read(p.strip()) for p in raw.split(',')]
