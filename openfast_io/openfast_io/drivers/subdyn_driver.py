"""SubDyn standalone driver — reads/writes SubDyn driver input files (.dvr).

Handles environmental conditions, SubDyn input file reference, TP reference
points, input modes, steady-state 6-DOF vectors, and applied loads table.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict, List

from ..io.subdyn import SubDynIO
from ..outlist import capture_outlist
from ..parsing import bool_read, float_read, int_read, quoted_read

try:
    from ..FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}


class SubDynStandaloneDriver:
    """Reads and writes SubDyn standalone driver input files (.dvr)."""

    def __init__(self) -> None:
        self._subdyn = SubDynIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, dvr_path: Path) -> dict:
        """Read a SubDyn driver file and the referenced SubDyn input.

        Returns a dict with keys:
          ``'SubDynDriver'`` — driver-level parameters
          ``'SubDyn'`` — from SubDynIO
        """
        dvr_path = Path(dvr_path)
        base_dir = dvr_path.parent
        dvr: Dict[str, Any] = {}

        with open(dvr_path) as f:
            # --- Header ---
            f.readline()  # title line
            f.readline()  # compatibility line
            dvr['Echo'] = bool_read(f.readline().split()[0])

            # --- Environmental Conditions ---
            f.readline()  # section header
            dvr['Gravity'] = float_read(f.readline().split()[0])
            dvr['WtrDpth'] = float_read(f.readline().split()[0])

            # --- SubDyn ---
            f.readline()  # section header
            dvr['SDInputFile'] = quoted_read(f.readline().split()[0])
            dvr['OutRootName'] = quoted_read(f.readline().split()[0])
            dvr['NSteps'] = int_read(f.readline().split()[0])
            dvr['TimeInterval'] = float_read(f.readline().split()[0])
            dvr['NTPs'] = int_read(f.readline().split()[0])
            dvr['TP_RefPoint_X'] = float_read(f.readline().split()[0])
            dvr['TP_RefPoint_Y'] = float_read(f.readline().split()[0])
            dvr['TP_RefPoint_Z'] = float_read(f.readline().split()[0])
            dvr['SubRotateZ'] = float_read(f.readline().split()[0])

            # --- Inputs ---
            f.readline()  # section header
            dvr['InputsMod'] = int_read(f.readline().split()[0])
            dvr['InputsFile'] = quoted_read(f.readline().split()[0])

            # --- Steady Inputs ---
            f.readline()  # section header
            dvr['uTPInSteady'] = _read_6dof_vec(f.readline())
            dvr['uDotTPInSteady'] = _read_6dof_vec(f.readline())
            dvr['uDotDotTPInSteady'] = _read_6dof_vec(f.readline())

            # --- Loads ---
            f.readline()  # section header
            dvr['nAppliedLoads'] = int_read(f.readline().split()[0])
            f.readline()  # column headers
            f.readline()  # column units
            loads: List[Dict[str, Any]] = []
            for _ in range(dvr['nAppliedLoads']):
                ln = f.readline().split()
                loads.append({
                    'ALJointID': int_read(ln[0]),
                    'Fx': float_read(ln[1]),
                    'Fy': float_read(ln[2]),
                    'Fz': float_read(ln[3]),
                    'Mx': float_read(ln[4]),
                    'My': float_read(ln[5]),
                    'Mz': float_read(ln[6]),
                    'UnsteadyFile': ln[7] if len(ln) > 7 else '',
                })
            dvr['AppliedLoads'] = loads

        result: Dict[str, Any] = {'SubDynDriver': dvr}

        # Shared OutList registry (mirrors OpenFASTDriver.read) so the SubDyn
        # SSOutList section survives a standalone read → write roundtrip.
        outlist: Dict[str, Any] = copy.deepcopy(FstOutput) if FstOutput else {}

        def _cap(f, module, freeform=False):
            return capture_outlist(f, outlist, module, freeform=freeform)

        # --- Delegate to SubDynIO ---
        sd_file = dvr.get('SDInputFile', '')
        sd_path = os.path.normpath(os.path.join(str(base_dir), sd_file))
        if sd_file and os.path.isfile(sd_path):
            sd_data = self._subdyn.read(
                Path(sd_path), base_dir,
                outlist=outlist,
                read_outlist_fn=lambda f, module: _cap(f, module, freeform=True),
            )
            result.update(sd_data)

        result['outlist'] = outlist
        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, dvr_path: Path) -> None:
        """Write a SubDyn driver file and the referenced SubDyn input."""
        dvr_path = Path(dvr_path)
        base_dir = dvr_path.parent
        dvr = data['SubDynDriver']

        with open(dvr_path, 'w') as f:
            f.write('SubDyn Driver file for stand-alone applications\n')
            f.write('Compatible with SubDyn v1.xx.x\n')
            f.write('{!s:<19} {:<15} {:}\n'.format(dvr['Echo'], 'Echo', '- Echo the input file data (flag).'))
            f.write('---------------------- ENVIRONMENTAL CONDITIONS -------------------------------------------------\n')
            f.write('{:<19} {:<15} {:}\n'.format(dvr['Gravity'], 'Gravity', '- Gravity (m/s^2).'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['WtrDpth'], 'WtrDpth', '- Water Depth (m) positive value.'))
            f.write('---------------------- SubDyn -------------------------------------------------------------------\n')

            sd_name = dvr.get('SDInputFile', 'SubDyn.dat')
            f.write('{:<19} {:<15} {:}\n'.format('"' + sd_name + '"', 'SDInputFile', '- SubDyn input file.'))
            f.write('{:<19} {:<15} {:}\n'.format('"' + dvr.get('OutRootName', 'SubDyn') + '"', 'OutRootName', '- All the output files will have this name.'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['NSteps'], 'NSteps', '- Number of time steps in the simulations (-).'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['TimeInterval'], 'TimeInterval', '- TimeInterval for the simulation (sec).'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['NTPs'], 'NTPs', '- Number of transition pieces'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['TP_RefPoint_X'], 'TP_RefPoint_X', '- X location of the TP reference points in global coordinates (m)'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['TP_RefPoint_Y'], 'TP_RefPoint_Y', '- Y location of the TP reference points in global coordinates (m)'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['TP_RefPoint_Z'], 'TP_RefPoint_Z', '- Z location of the TP reference points in global coordinates (m)'))
            f.write('{:<19} {:<15} {:}\n'.format(dvr['SubRotateZ'], 'SubRotateZ', '- Rotation angle of the structure geometry in [deg] about the global Z axis.'))
            f.write('---------------------- INPUTS -------------------------------------------------------------------\n')
            f.write('{:<19} {:<15} {:}\n'.format(dvr['InputsMod'], 'InputsMod', '- Inputs model {0: zero, 1: steady state, 2: from file} (switch)'))
            f.write('{:<19} {:<15} {:}\n'.format('"' + dvr.get('InputsFile', '') + '"', 'InputsFile', '- Name of the inputs file if InputsMod = 2.'))
            f.write('---------------------- STEADY INPUTS (for InputsMod = 1) ----------------------------------------\n')
            f.write('{}   {}\n'.format('   '.join('{}'.format(v) for v in dvr['uTPInSteady']), 'uTPInSteady     - input displacements and rotations ([m], [rad])'))
            f.write('{}   {}\n'.format('   '.join('{}'.format(v) for v in dvr['uDotTPInSteady']), 'uDotTPInSteady  - input translational and rotational velocities ([m/s], [rad/s])'))
            f.write('{}   {}\n'.format('   '.join('{}'.format(v) for v in dvr['uDotDotTPInSteady']), 'uDotTPInSteady  - input translational and rotational accelerations([m/s^2], [rad/s^2])'))
            f.write('---------------------- LOADS --------------------------------------------------------------------\n')
            f.write('{:<5} {:<15} {:}\n'.format(dvr['nAppliedLoads'], 'nAppliedLoads', '- Number of applied loads at given nodes'))
            f.write('ALJointID    Fx     Fy    Fz     Mx     My     Mz   UnsteadyFile\n')
            f.write('   (-)       (N)    (N)   (N)   (Nm)   (Nm)   (Nm)     (-)\n')
            for load in dvr.get('AppliedLoads', []):
                f.write('{:>5}  {:>8}  {:>8}  {:>8}  {:>8}  {:>8}  {:>8}  {}\n'.format(
                    load['ALJointID'], load['Fx'], load['Fy'], load['Fz'],
                    load['Mx'], load['My'], load['Mz'], load.get('UnsteadyFile', '')))
            f.write('END of driver input file\n')

        # --- Delegate to SubDynIO ---
        if 'SubDyn' in data:
            sd_path = os.path.normpath(os.path.join(str(base_dir), sd_name))
            outlist = data.get('outlist') or (copy.deepcopy(FstOutput) if FstOutput else {})
            self._subdyn.write({'SubDyn': data['SubDyn']}, sd_path, str(base_dir), outlist=outlist)


def _read_6dof_vec(line: str) -> list:
    """Parse a 6-value space-separated vector from a line like '0.5 0.0 0.0 0.0 0.0 0.0  name'."""
    parts = line.split()
    return [float_read(parts[i]) for i in range(min(6, len(parts)))]
