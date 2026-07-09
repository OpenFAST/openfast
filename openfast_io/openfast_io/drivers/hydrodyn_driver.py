"""HydroDyn standalone driver — reads/writes HydroDyn driver input files (.inp).

Handles environmental conditions, primary HD + SeaState file references,
PRP (Platform Reference Point) inputs, and steady-state 6-DOF vectors.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict

from ..io.hydrodyn import HydroDynIO
from ..io.seastate import SeaStateIO
from ..outlist import capture_outlist
from ..parsing import bool_read, float_read, int_read, quoted_read

try:
    from ..FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}


class HydroDynStandaloneDriver:
    """Reads and writes HydroDyn standalone driver input files (.inp)."""

    def __init__(self) -> None:
        self._hydrodyn = HydroDynIO()
        self._seastate = SeaStateIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, inp_path: Path) -> dict:
        """Read a HydroDyn driver file and all referenced module files.

        Returns a dict with keys:
          ``'HydroDynDriver'`` — driver-level parameters
          ``'HydroDyn'`` — from HydroDynIO
          ``'SeaState'`` — from SeaStateIO
        """
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr: Dict[str, Any] = {}

        with open(inp_path) as f:
            # --- Header ---
            f.readline()  # title line
            dvr['description'] = f.readline().rstrip()

            # --- Environmental Conditions ---
            dvr['Echo'] = bool_read(f.readline().split()[0])
            f.readline()  # section header
            dvr['Gravity'] = float_read(f.readline().split()[0])
            dvr['WtrDens'] = float_read(f.readline().split()[0])
            dvr['WtrDpth'] = float_read(f.readline().split()[0])
            dvr['MSL2SWL'] = float_read(f.readline().split()[0])

            # --- HydroDyn ---
            f.readline()  # section header
            dvr['HDInputFile'] = quoted_read(f.readline().split()[0])
            dvr['SeaStateInputFile'] = quoted_read(f.readline().split()[0])
            dvr['OutRootName'] = quoted_read(f.readline().split()[0])
            dvr['Linearize'] = bool_read(f.readline().split()[0])
            dvr['NSteps'] = int_read(f.readline().split()[0])
            dvr['TimeInterval'] = float_read(f.readline().split()[0])

            # --- PRP Inputs ---
            f.readline()  # section header
            dvr['PRPInputsMod'] = int_read(f.readline().split()[0])
            dvr['NAddDOF'] = int_read(f.readline().split()[0])
            dvr['PtfmRefzt'] = float_read(f.readline().split()[0])
            dvr['PRPInputsFile'] = quoted_read(f.readline().split()[0])

            # --- PRP Steady State Inputs ---
            f.readline()  # section header
            dvr['uPRPInSteady'] = _read_6dof_vec(f.readline())
            dvr['uDotPRPInSteady'] = _read_6dof_vec(f.readline())
            dvr['uDotDotPRPInSteady'] = _read_6dof_vec(f.readline())

        result: Dict[str, Any] = {'HydroDynDriver': dvr}

        # Shared OutList registry (mirrors OpenFASTDriver.read) so both the
        # HydroDyn and SeaState output channel sections survive a standalone
        # read → write roundtrip.
        outlist: Dict[str, Any] = copy.deepcopy(FstOutput) if FstOutput else {}

        def _cap(f, module, freeform=False):
            return capture_outlist(f, outlist, module, freeform=freeform)

        # --- Delegate to HydroDynIO ---
        hd_file = dvr.get('HDInputFile', '')
        hd_path = os.path.normpath(os.path.join(str(base_dir), hd_file))
        if hd_file and os.path.isfile(hd_path):
            hd_data = self._hydrodyn.read(Path(hd_path), base_dir,
                                           outlist=outlist, read_outlist_fn=_cap)
            result.update(hd_data)

        # --- Delegate to SeaStateIO ---
        ss_file = dvr.get('SeaStateInputFile', '')
        ss_path = os.path.normpath(os.path.join(str(base_dir), ss_file))
        if ss_file and os.path.isfile(ss_path):
            ss_data = self._seastate.read(
                Path(ss_path), base_dir,
                outlist=outlist,
                read_outlist_fn=lambda f, module: _cap(f, module, freeform=True),
            )
            result.update(ss_data)

        result['outlist'] = outlist
        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, inp_path: Path) -> None:
        """Write a HydroDyn driver file and all referenced module files."""
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr = data['HydroDynDriver']

        with open(inp_path, 'w') as f:
            f.write('------- HydroDyn Driver Input File --------------------------------------------\n')
            f.write('Generated by OpenFAST_IO\n')
            f.write('{!s:<14} {:<20} {:}\n'.format(dvr['Echo'], 'Echo', '- Echo the input file data (flag)'))
            f.write('---------------------- ENVIRONMENTAL CONDITIONS -------------------------------\n')
            f.write('{:<14} {:<20} {:}\n'.format(dvr['Gravity'], 'Gravity', '- Gravity (m/s^2)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['WtrDens'], 'WtrDens', '- Water density (kg/m^3)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['WtrDpth'], 'WtrDpth', '- Water depth (m)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['MSL2SWL'], 'MSL2SWL', '- Offset between still-water level and mean sea level (m) [positive upward]'))
            f.write('---------------------- HYDRODYN -----------------------------------------------\n')

            hd_name = dvr.get('HDInputFile', 'HydroDyn.dat')
            ss_name = dvr.get('SeaStateInputFile', 'SeaState.dat')
            f.write('{:<14} {:<20} {:}\n'.format('"' + hd_name + '"', 'HDInputFile', '- Primary HydroDyn input file name (quoted string)'))
            f.write('{:<14} {:<20} {:}\n'.format('"' + ss_name + '"', 'SeaStateInputFile', '- Primary SeaState input file name (quoted string)'))
            f.write('{:<14} {:<20} {:}\n'.format('"' + dvr.get('OutRootName', './driver') + '"', 'OutRootName', '- The name which prefixes all HydroDyn generated files (quoted string)'))
            f.write('{!s:<14} {:<20} {:}\n'.format(dvr['Linearize'], 'Linearize', '- Flag to enable linearization'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['NSteps'], 'NSteps', '- Number of time steps in the simulations (-)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['TimeInterval'], 'TimeInterval', '- TimeInterval for the simulation (sec)'))
            f.write('---------------------- PRP INPUTS (Platform Reference Point) ------------------\n')
            f.write('{:<14} {:<20} {:}\n'.format(dvr['PRPInputsMod'], 'PRPInputsMod', '- Model for the PRP inputs (switch)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['NAddDOF'], 'NAddDOF', '- Number of additional generalized DOF (-)'))
            f.write('{:<14} {:<20} {:}\n'.format(dvr['PtfmRefzt'], 'PtfmRefzt', '- Vertical distance from ground to platform reference point (m)'))
            f.write('{:<14} {:<20} {:}\n'.format('"' + dvr.get('PRPInputsFile', '') + '"', 'PRPInputsFile', '- Filename for the PRP inputs'))
            f.write('---------------------- PRP STEADY STATE INPUTS  -------------------------------\n')
            f.write('{:>10},{:>10},{:>10},{:>10},{:>10},{:>10}    {:<25} {:}\n'.format(
                *dvr['uPRPInSteady'], 'uPRPInSteady', '- PRP Steady-state displacements and rotations'))
            f.write('{:>10},{:>10},{:>10},{:>10},{:>10},{:>10}    {:<25} {:}\n'.format(
                *dvr['uDotPRPInSteady'], 'uDotPRPInSteady', '- PRP Steady-state velocities'))
            f.write('{:>10},{:>10},{:>10},{:>10},{:>10},{:>10}    {:<25} {:}\n'.format(
                *dvr['uDotDotPRPInSteady'], 'uDotDotPRPInSteady', '- PRP Steady-state accelerations'))

        outlist = data.get('outlist') or (copy.deepcopy(FstOutput) if FstOutput else {})

        # --- Delegate to HydroDynIO ---
        if 'HydroDyn' in data:
            hd_path = os.path.normpath(os.path.join(str(base_dir), hd_name))
            self._hydrodyn.write({'HydroDyn': data['HydroDyn']}, hd_path, str(base_dir), outlist=outlist)

        # --- Delegate to SeaStateIO ---
        if 'SeaState' in data:
            ss_path = os.path.normpath(os.path.join(str(base_dir), ss_name))
            self._seastate.write({'SeaState': data['SeaState']}, ss_path, str(base_dir), outlist=outlist)


def _read_6dof_vec(line: str) -> list:
    """Parse a 6-value comma-separated vector from a line like '1, 2, -3, 0, 0, 0  name'."""
    # Split on whitespace, take tokens that look numeric (before the parameter name)
    parts = line.replace(',', ' ').split()
    values = []
    for p in parts:
        try:
            values.append(float_read(p))
        except (ValueError, IndexError):
            break
        if len(values) == 6:
            break
    return values
