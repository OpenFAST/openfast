"""FASTFarmDriver — reads/writes FAST.Farm (.fstf) simulation decks.

Composes ``OpenFASTDriver`` per turbine, plus farm-level parameters.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List

from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
)
from .openfast import OpenFASTDriver


class FASTFarmDriver:
    """Reads and writes a FAST.Farm simulation deck.

    The result dict has:
      - ``'FASTFarm'``: farm-level scalars (sim control, ambient wind, wake, output, …)
      - ``'Turbines'``: list of per-turbine dicts, each containing a full ``fst_vt``
    """

    def __init__(self):
        self._of_driver = OpenFASTDriver()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, fstf_path: Path) -> dict:
        fstf_path = Path(fstf_path)
        base_dir = fstf_path.parent
        ff: Dict[str, Any] = {}

        f = open(fstf_path)
        f.readline()  # header 1
        ff['description'] = f.readline().rstrip()

        # --- SIMULATION CONTROL ---
        f.readline()
        ff['Echo'] = bool_read(f.readline().split()[0])
        ff['AbortLevel'] = f.readline().split()[0]
        ff['TMax'] = float_read(f.readline().split()[0])
        ff['Mod_AmbWind'] = int_read(f.readline().split()[0])
        ff['Mod_WaveField'] = int_read(f.readline().split()[0])
        ff['Mod_SharedMooring'] = int_read(f.readline().split()[0])

        # --- SHARED MOORING ---
        f.readline()
        ff['SharedMoorFile'] = quoted_read(f.readline().split()[0])
        ff['DT_Mooring'] = float_read(f.readline().split()[0])
        ff['WrMooringVis'] = bool_read(f.readline().split()[0])

        # --- AMBIENT WIND: VTK ---
        f.readline()
        ff['DT_Low-VTK'] = float_read(f.readline().split()[0])
        ff['DT_High-VTK'] = float_read(f.readline().split()[0])
        ff['WindFilePath'] = quoted_read(f.readline().split()[0])
        ff['ChkWndFiles'] = bool_read(f.readline().split()[0])

        # --- AMBIENT WIND: INFLOWWIND ---
        f.readline()
        ff['DT_Low'] = float_read(f.readline().split()[0])
        ff['DT_High'] = float_read(f.readline().split()[0])
        ff['NX_Low'] = int_read(f.readline().split()[0])
        ff['NY_Low'] = int_read(f.readline().split()[0])
        ff['NZ_Low'] = int_read(f.readline().split()[0])
        ff['X0_Low'] = float_read(f.readline().split()[0])
        ff['Y0_Low'] = float_read(f.readline().split()[0])
        ff['Z0_Low'] = float_read(f.readline().split()[0])
        ff['dX_Low'] = float_read(f.readline().split()[0])
        ff['dY_Low'] = float_read(f.readline().split()[0])
        ff['dZ_Low'] = float_read(f.readline().split()[0])
        ff['NX_High'] = int_read(f.readline().split()[0])
        ff['NY_High'] = int_read(f.readline().split()[0])
        ff['NZ_High'] = int_read(f.readline().split()[0])
        ff['InflowFile'] = quoted_read(f.readline().split()[0])

        # --- AMBIENT WIND: AMReX ---
        f.readline()
        ff['WindDirPrefix'] = quoted_read(f.readline().split()[0])
        ff['DirStartIndex'] = int_read(f.readline().split()[0])
        ff['DT_Low-AMReX'] = float_read(f.readline().split()[0])
        ff['DT_High-AMReX'] = float_read(f.readline().split()[0])

        # --- WIND TURBINES ---
        f.readline()
        ff['NumTurbines'] = int_read(f.readline().split()[0])
        f.readline()  # column headers (WT_X  WT_Y  WT_Z  WT_FASTInFile …)
        f.readline()  # column units   (m)    (m)   (m)   (string)      …

        turbine_rows = []
        for _ in range(ff['NumTurbines']):
            ln = f.readline().split()
            row = {
                'WT_X': float_read(ln[0]),
                'WT_Y': float_read(ln[1]),
                'WT_Z': float_read(ln[2]),
                'WT_FASTInFile': ln[3].strip('"'),
                'X0_High': float_read(ln[4]) if len(ln) > 4 else 0.0,
                'Y0_High': float_read(ln[5]) if len(ln) > 5 else 0.0,
                'Z0_High': float_read(ln[6]) if len(ln) > 6 else 0.0,
                'dX_High': float_read(ln[7]) if len(ln) > 7 else 0.0,
                'dY_High': float_read(ln[8]) if len(ln) > 8 else 0.0,
                'dZ_High': float_read(ln[9]) if len(ln) > 9 else 0.0,
            }
            turbine_rows.append(row)
        ff['TurbineRows'] = turbine_rows

        # Read remaining lines as raw text (wake dynamics, curl, WAT, viz, output)
        remaining = f.read()
        ff['_remaining'] = remaining
        f.close()

        # --- Read per-turbine .fst files ---
        turbines: List[dict] = []
        for row in turbine_rows:
            fst_rel = row['WT_FASTInFile']
            fst_file = os.path.normpath(os.path.join(str(base_dir), fst_rel))
            if os.path.isfile(fst_file):
                turb_vt = self._of_driver.read(Path(fst_file))
            else:
                turb_vt = {}
            turb_vt['_farm_position'] = {
                'WT_X': row['WT_X'],
                'WT_Y': row['WT_Y'],
                'WT_Z': row['WT_Z'],
            }
            turbines.append(turb_vt)

        return {'FASTFarm': ff, 'Turbines': turbines}

    # ------------------------------------------------------------------
    # WRITE (stub — farm-level write is complex, for now just the .fstf)
    # ------------------------------------------------------------------
    def write(
        self,
        data: dict,
        output_dir: Path,
        case_name: str,
    ) -> list:
        """Write the .fstf file and per-turbine files.

        Currently writes only the .fstf file (no per-turbine write delegation yet).
        Returns list of written paths.
        """
        raise NotImplementedError(
            'FASTFarmDriver.write() is not yet implemented. '
            'Use OpenFASTDriver.write() for individual turbine files.'
        )
