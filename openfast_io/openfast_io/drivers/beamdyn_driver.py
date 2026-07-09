"""BeamDyn standalone driver — reads/writes BeamDyn driver input files (.inp).

Handles simulation control, gravity, frame parameters (including 3x3 DCM),
root velocity, applied forces, point loads, and primary input file reference.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict, List

from ..io.beamdyn import BeamDynIO
from ..outlist import capture_outlist
from ..parsing import bool_read, float_read, fmt_field, int_read, quoted_read

try:
    from ..FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}


def _fw(val, width: int = 16) -> str:
    """Format a value field with guaranteed 2-space separator."""
    return fmt_field(val, min_width=width)


class BeamDynStandaloneDriver:
    """Reads and writes BeamDyn standalone driver input files (.inp)."""

    def __init__(self) -> None:
        self._beamdyn = BeamDynIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, inp_path: Path) -> dict:
        """Read a BeamDyn driver file and the referenced primary input file.

        Returns a dict with keys:
          ``'BeamDynDriver'`` — driver-level parameters
          ``'BeamDyn'``, ``'BeamDynBlade'`` — from BeamDynIO
        """
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr: Dict[str, Any] = {}

        with open(inp_path) as f:
            # --- Header ---
            f.readline()  # title line
            dvr['description'] = f.readline().rstrip()

            # --- Simulation Control ---
            f.readline()  # section header
            dvr['DynamicSolve'] = bool_read(f.readline().split()[0])
            dvr['t_initial'] = float_read(f.readline().split()[0])
            dvr['t_final'] = float_read(f.readline().split()[0])
            dvr['dt'] = float_read(f.readline().split()[0])

            # --- Gravity Parameter ---
            f.readline()  # section header
            dvr['Gx'] = float_read(f.readline().split()[0])
            dvr['Gy'] = float_read(f.readline().split()[0])
            dvr['Gz'] = float_read(f.readline().split()[0])

            # --- Frame Parameter ---
            f.readline()  # section header
            dvr['GlbPos'] = [
                float_read(f.readline().split()[0]),
                float_read(f.readline().split()[0]),
                float_read(f.readline().split()[0]),
            ]

            # 3x3 direction cosine matrix
            f.readline()  # DCM comment line 1
            f.readline()  # DCM comment line 2
            dcm: List[List[float]] = []
            for _ in range(3):
                row = [float_read(v) for v in f.readline().split()]
                dcm.append(row)
            dvr['GlbDCM'] = dcm

            dvr['GlbRotBladeT0'] = bool_read(f.readline().split()[0])

            # --- Root Velocity Parameter ---
            f.readline()  # section header
            dvr['RootVel'] = [
                float_read(f.readline().split()[0]),
                float_read(f.readline().split()[0]),
                float_read(f.readline().split()[0]),
            ]

            # --- Applied Force ---
            f.readline()  # section header
            dvr['DistrLoad'] = [float_read(f.readline().split()[0]) for _ in range(6)]
            dvr['TipLoad'] = [float_read(f.readline().split()[0]) for _ in range(6)]
            dvr['NumPointLoads'] = int_read(f.readline().split()[0])

            # Point loads table
            f.readline()  # column headers
            f.readline()  # column units
            dvr['PointLoads'] = []
            for _ in range(dvr['NumPointLoads']):
                ln = f.readline().split()
                dvr['PointLoads'].append({
                    'eta': float_read(ln[0]),
                    'Fx': float_read(ln[1]),
                    'Fy': float_read(ln[2]),
                    'Fz': float_read(ln[3]),
                    'Mx': float_read(ln[4]),
                    'My': float_read(ln[5]),
                    'Mz': float_read(ln[6]),
                })

            # --- Primary Input File ---
            f.readline()  # section header
            dvr['InputFile'] = quoted_read(f.readline().split()[0])

            # --- Output Settings ---
            f.readline()  # section header
            dvr['WrVTK'] = int_read(f.readline().split()[0])
            try:
                dvr['VTK_fps'] = int_read(f.readline().split()[0])
            except (IndexError, ValueError):
                dvr['VTK_fps'] = 15

        result: Dict[str, Any] = {'BeamDynDriver': dvr}

        # Shared OutList registry (mirrors OpenFASTDriver.read) so the BeamDyn
        # OutList section survives a standalone read → write roundtrip.
        outlist: Dict[str, Any] = copy.deepcopy(FstOutput) if FstOutput else {}

        def _cap(f, module, freeform=False):
            return capture_outlist(f, outlist, module, freeform=freeform)

        # --- Delegate to BeamDynIO ---
        bd_file = dvr.get('InputFile', '')
        bd_path = os.path.normpath(os.path.join(str(base_dir), bd_file))
        if bd_file and os.path.isfile(bd_path):
            bd_data = self._beamdyn.read(Path(bd_path), base_dir,
                                          outlist=outlist, read_outlist_fn=_cap)
            result.update(bd_data)

        result['outlist'] = outlist
        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, inp_path: Path) -> None:
        """Write a BeamDyn driver file and the referenced primary input."""
        inp_path = Path(inp_path)
        base_dir = inp_path.parent
        dvr = data['BeamDynDriver']

        with open(inp_path, 'w') as f:
            f.write('------- BEAMDYN Driver with OpenFAST INPUT FILE --------------------------------\n')
            f.write(dvr.get('description', 'Generated by OpenFAST_IO') + '\n')
            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write(_fw(dvr['DynamicSolve']) + '{:<14} {:}\n'.format('DynamicSolve', '- Dynamic solve (false for static solve) (-)'))
            f.write(_fw(dvr['t_initial']) + '{:<14} {:}\n'.format('t_initial', '- Starting time of simulation (s)'))
            f.write(_fw(dvr['t_final']) + '{:<14} {:}\n'.format('t_final', '- Ending time of simulation (s)'))
            f.write(_fw(dvr['dt']) + '{:<14} {:}\n'.format('dt', '- Time increment size (s)'))
            f.write('---------------------- GRAVITY PARAMETER --------------------------------------\n')
            f.write(_fw(dvr['Gx']) + '{:<14} {:}\n'.format('Gx', '- Component of gravity vector along X direction (m/s^2)'))
            f.write(_fw(dvr['Gy']) + '{:<14} {:}\n'.format('Gy', '- Component of gravity vector along Y direction (m/s^2)'))
            f.write(_fw(dvr['Gz']) + '{:<14} {:}\n'.format('Gz', '- Component of gravity vector along Z direction (m/s^2)'))
            f.write('---------------------- FRAME PARAMETER --------------------------------------\n')
            f.write(_fw(dvr['GlbPos'][0]) + '{:<14} {:}\n'.format('GlbPos(1)', '- Component of position vector along X direction (m)'))
            f.write(_fw(dvr['GlbPos'][1]) + '{:<14} {:}\n'.format('GlbPos(2)', '- Component of position vector along Y direction (m)'))
            f.write(_fw(dvr['GlbPos'][2]) + '{:<14} {:}\n'.format('GlbPos(3)', '- Component of position vector along Z direction (m)'))
            f.write('---The following 3 by 3 matrix is the direction cosine matirx ,GlbDCM(3,3),\n')
            f.write('---relates global frame to the initial blade root frame\n')
            for row in dvr['GlbDCM']:
                f.write('  '.join('{:.7E}'.format(v) for v in row) + '\n')
            f.write('{!s:<14} {:<14} {:}\n'.format(dvr['GlbRotBladeT0'], 'GlbRotBladeT0', '- Reference orientation aligned with initial blade root?'))
            f.write('---------------------- ROOT VELOCITY PARAMETER ----------------------------------\n')
            f.write(_fw(dvr['RootVel'][0]) + '{:<14} {:}\n'.format('RootVel(4)', '- Component of angular velocity about X axis (rad/s)'))
            f.write(_fw(dvr['RootVel'][1]) + '{:<14} {:}\n'.format('RootVel(5)', '- Component of angular velocity about Y axis (rad/s)'))
            f.write(_fw(dvr['RootVel'][2]) + '{:<14} {:}\n'.format('RootVel(6)', '- Component of angular velocity about Z axis (rad/s)'))
            f.write('---------------------- APPLIED FORCE ----------------------------------\n')
            labels = ['DistrLoad(1)', 'DistrLoad(2)', 'DistrLoad(3)', 'DistrLoad(4)', 'DistrLoad(5)', 'DistrLoad(6)']
            descs = [
                '- Component of distributed force vector along X direction (N/m)',
                '- Component of distributed force vector along Y direction (N/m)',
                '- Component of distributed force vector along Z direction (N/m)',
                '- Component of distributed moment vector along X direction (N-m/m)',
                '- Component of distributed moment vector along Y direction (N-m/m)',
                '- Component of distributed moment vector along Z direction (N-m/m)',
            ]
            for i in range(6):
                f.write(_fw(dvr['DistrLoad'][i]) + '{:<14} {:}\n'.format(labels[i], descs[i]))
            tip_labels = ['TipLoad(1)', 'TipLoad(2)', 'TipLoad(3)', 'TipLoad(4)', 'TipLoad(5)', 'TipLoad(6)']
            tip_descs = [
                '- Component of concentrated force at blade tip along X direction (N)',
                '- Component of concentrated force at blade tip along Y direction (N)',
                '- Component of concentrated force at blade tip along Z direction (N)',
                '- Component of concentrated moment at blade tip along X direction (N-m)',
                '- Component of concentrated moment at blade tip along Y direction (N-m)',
                '- Component of concentrated moment at blade tip along Z direction (N-m)',
            ]
            for i in range(6):
                f.write(_fw(dvr['TipLoad'][i]) + '{:<14} {:}\n'.format(tip_labels[i], tip_descs[i]))
            f.write(_fw(dvr['NumPointLoads']) + '{:<14} {:}\n'.format('NumPointLoads', '- Number of point loads along blade'))
            f.write('Non-dim blade-span eta   Fx          Fy            Fz           Mx           My           Mz\n')
            f.write('(-)                      (N)         (N)           (N)          (N-m)        (N-m)        (N-m)\n')
            for pl in dvr.get('PointLoads', []):
                f.write('{:<10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}\n'.format(
                    pl['eta'], pl['Fx'], pl['Fy'], pl['Fz'], pl['Mx'], pl['My'], pl['Mz']))

            f.write('---------------------- PRIMARY INPUT FILE --------------------------------------\n')
            bd_name = dvr.get('InputFile', 'bd_primary.inp')
            f.write('{:<14} {:<14} {:}\n'.format('"' + bd_name + '"', 'InputFile', '- Name of the primary BeamDyn input file'))

            f.write('----- Output Settings -------------------------------------------------------------------\n')
            f.write(_fw(dvr.get('WrVTK', 0)) + '{:<14} {:}\n'.format('WrVTK', '- VTK visualization data output (switch)'))
            f.write(_fw(dvr.get('VTK_fps', 15)) + '{:<14} {:}\n'.format('VTK_fps', '- Frame rate for VTK output (frames per second)'))

        # --- Delegate to BeamDynIO ---
        if 'BeamDyn' in data:
            bd_path = os.path.normpath(os.path.join(str(base_dir), bd_name))
            bd_data = {'BeamDyn': data['BeamDyn']}
            if 'BeamDynBlade' in data:
                bd_data['BeamDynBlade'] = data['BeamDynBlade']
            outlist = data.get('outlist') or (copy.deepcopy(FstOutput) if FstOutput else {})
            self._beamdyn.write(bd_data, Path(bd_path), base_dir, outlist=outlist)
