"""UnsteadyAero standalone driver — reads/writes UA driver input files (.dvr).

The most complex standalone driver format. Supports three simulation modes:
  1 = reduced frequency
  2 = prescribed-aero time series
  3 = elastic cross section (with 3x3 matrices)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from ..parsing import bool_read, float_read, fmt_field, int_read, quoted_read


def _fw(val, width: int = 28) -> str:
    """Format a value field with guaranteed 2-space separator."""
    return fmt_field(val, min_width=width)


class UnsteadyAeroStandaloneDriver:
    """Reads and writes UnsteadyAero standalone driver input files (.dvr)."""

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, dvr_path: Path) -> dict:
        """Read an UnsteadyAero driver file.

        No sub-module delegation — the airfoil table is referenced by filename
        but is a simple polar table, not a full ModuleIO file.

        Returns a dict with key ``'UnsteadyAeroDriver'``.
        """
        dvr_path = Path(dvr_path)
        dvr: Dict[str, Any] = {}

        with open(dvr_path) as f:
            # --- Header ---
            f.readline()  # title line
            f.readline()  # separator / comment
            dvr['Echo'] = bool_read(f.readline().split()[0])

            # --- Environmental Conditions ---
            f.readline()  # section header
            dvr['FldDens'] = float_read(f.readline().split()[0])
            dvr['KinVisc'] = float_read(f.readline().split()[0])
            dvr['SpdSound'] = float_read(f.readline().split()[0])

            # --- UnsteadyAero ---
            f.readline()  # section header
            dvr['UAMod'] = int_read(f.readline().split()[0])
            dvr['FLookup'] = bool_read(f.readline().split()[0])

            # --- Airfoil Properties ---
            f.readline()  # section header
            dvr['AirFoil'] = quoted_read(f.readline().split()[0])
            dvr['Chord'] = float_read(f.readline().split()[0])
            dvr['Vec_AQ'] = _read_csv_vec(f.readline())
            dvr['Vec_AT'] = _read_csv_vec(f.readline())
            dvr['UseCm'] = bool_read(f.readline().split()[0])

            # --- Simulation Control ---
            f.readline()  # section header
            dvr['SimMod'] = int_read(f.readline().split()[0])

            # --- Reduced-Frequency Simulation (SimMod=1) ---
            f.readline()  # section header
            dvr['InflowVel'] = float_read(f.readline().split()[0])
            dvr['NCycles'] = int_read(f.readline().split()[0])
            dvr['StepsPerCycle'] = int_read(f.readline().split()[0])
            dvr['Frequency'] = float_read(f.readline().split()[0])
            dvr['Amplitude'] = float_read(f.readline().split()[0])
            dvr['Mean'] = float_read(f.readline().split()[0])
            dvr['Phase'] = float_read(f.readline().split()[0])

            # --- Prescribed-Aero Simulation (SimMod=2) ---
            f.readline()  # section header
            dvr['TMax_PA'] = float_read(f.readline().split()[0])
            dvr['DT_PA'] = float_read(f.readline().split()[0])
            dvr['AeroTSFile'] = quoted_read(f.readline().split()[0])

            # --- Aero-Elastic Simulation (SimMod=3) ---
            f.readline()  # section header
            dvr['TMax'] = float_read(f.readline().split()[0])
            dvr['DT'] = float_read(f.readline().split()[0])
            dvr['ActiveDOF'] = _read_csv_bools(f.readline())
            dvr['InitPos'] = _read_csv_vec(f.readline())
            dvr['InitVel'] = _read_csv_vec(f.readline())

            # 3x3 matrices: GFScaling, Mass, Damp, Stif
            dvr['GFScaling'] = _read_3x3(f)
            dvr['MassMatrix'] = _read_3x3(f)
            dvr['DampMatrix'] = _read_3x3(f)
            dvr['StifMatrix'] = _read_3x3(f)

            dvr['Twist'] = float_read(f.readline().split()[0])
            dvr['InflowMod'] = int_read(f.readline().split()[0])
            dvr['Inflow'] = _read_csv_vec(f.readline())
            dvr['InflowTSFile'] = quoted_read(f.readline().split()[0])
            dvr['MotionMod'] = int_read(f.readline().split()[0])
            dvr['MotionTSFile'] = quoted_read(f.readline().split()[0])

            # --- Output Control ---
            f.readline()  # section header
            dvr['SumPrint'] = bool_read(f.readline().split()[0])
            dvr['WrAFITables'] = bool_read(f.readline().split()[0])

        return {'UnsteadyAeroDriver': dvr}

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, dvr_path: Path) -> None:
        """Write an UnsteadyAero driver file.

        No sub-module delegation — all data lives in 'UnsteadyAeroDriver'.
        """
        dvr_path = Path(dvr_path)
        dvr = data['UnsteadyAeroDriver']

        with open(dvr_path, 'w') as f:
            f.write('----------- Unsteady Aerodynamics standalone driver ---------------------------\n')
            f.write('-------------------------------------------------------------------------------\n')
            f.write('{!s:<14} {:<20} {:}\n'.format(dvr['Echo'], 'Echo', '- Echo the input file data (flag)'))
            f.write('---------------------- ENVIRONMENTAL CONDITIONS -------------------------------\n')
            f.write('{:<13} {:<20} {:}\n'.format(dvr['FldDens'], 'FldDens', '- Density of working fluid (kg/m^3)'))
            f.write('{:<13} {:<20} {:}\n'.format(dvr['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)'))
            f.write('{:<13} {:<20} {:}\n'.format(dvr['SpdSound'], 'SpdSound', '- Speed of sound of working fluid (m/s)'))
            f.write('---------------------- UNSTEADYAERO -------------------------------------------\n')
            f.write('{:<13} {:<20} {:}\n'.format(dvr['UAMod'], 'UAMod', '- Unsteady Aero Model Switch (switch)'))
            f.write('{!s:<13} {:<20} {:}\n'.format(dvr['FLookup'], 'FLookup', '- Flag for f\' lookup table (flag)'))
            f.write('------------------- AIRFOIL PROPERTIES ----------------------------------------\n')
            f.write(_fw('"' + dvr['AirFoil'] + '"') + '{:<22} {:}\n'.format('AirFoil', '- Airfoil table file'))
            f.write(_fw(dvr['Chord']) + '{:<22} {:}\n'.format('Chord', '- Chord length (m)'))
            f.write(_fw(', '.join(str(v) for v in dvr['Vec_AQ'])) + '{:<22} {:}\n'.format('Vec_AQ', '- Vector from reference point A to aerodynamic center Q'))
            f.write(_fw(', '.join(str(v) for v in dvr['Vec_AT'])) + '{:<22} {:}\n'.format('Vec_AT', '- Vector from reference point A to three-quarter chord point T'))
            f.write('{!s:<14} {:<22} {:}\n'.format(dvr['UseCm'], 'UseCm', '- Use Cm (moment coefficient) data (flag)'))
            f.write('------------------- SIMULATION CONTROL ----------------------------------------\n')
            f.write(_fw(dvr['SimMod']) + '{:<22} {:}\n'.format('SimMod', '- Simulation model {1=reduced frequency, 2=prescribed-aero, 3=elastic}'))
            f.write('---------- REDUCED-FREQUENCY SIMULATION [used only when SimMod=1] -------------\n')
            f.write(_fw(dvr['InflowVel']) + '{:<22} {:}\n'.format('InflowVel', '- Inflow velocity (m/s)'))
            f.write(_fw(dvr['NCycles']) + '{:<22} {:}\n'.format('NCycles', '- Number of angle-of-attack oscillations (-)'))
            f.write(_fw(dvr['StepsPerCycle']) + '{:<22} {:}\n'.format('StepsPerCycle', '- Number of timesteps per cycle (-)'))
            f.write(_fw(dvr['Frequency']) + '{:<22} {:}\n'.format('Frequency', '- Frequency for the airfoil oscillations (Hz)'))
            f.write(_fw(dvr['Amplitude']) + '{:<22} {:}\n'.format('Amplitude', '- Amplitude of the angle of attack oscillations (deg)'))
            f.write(_fw(dvr['Mean']) + '{:<22} {:}\n'.format('Mean', '- Cycle mean (deg)'))
            f.write(_fw(dvr['Phase']) + '{:<22} {:}\n'.format('Phase', '- Initial phase (num steps).'))
            f.write('---------- PRESCRIBED-AERO INPUTS SIMULATION [used only when SimMod=2] --------\n')
            f.write(_fw(dvr['TMax_PA']) + '{:<22} {:}\n'.format('TMax_PA', '- Total run time (s)'))
            f.write(_fw(dvr['DT_PA']) + '{:<22} {:}\n'.format('DT_PA', '- Recommended module time step (s)'))
            f.write(_fw('"' + dvr['AeroTSFile'] + '"') + '{:<22} {:}\n'.format('AeroTSFile', '- Time series data input file'))
            f.write('---------- AERO-ELASTIC SIMULATION [used only when SimMod=3] ------------------\n')
            f.write(_fw(dvr['TMax']) + '{:<22} {:}\n'.format('TMax', '- Total run time (s)'))
            f.write(_fw(dvr['DT']) + '{:<22} {:}\n'.format('DT', '- Recommended module time step (s)'))
            f.write(_fw(', '.join('T' if v else 'F' for v in dvr['ActiveDOF'])) + '{:<22} {:}\n'.format('ActiveDOF', '- List of active degrees of freedom (true or false)'))
            f.write(_fw(', '.join(str(v) for v in dvr['InitPos'])) + '{:<22} {:}\n'.format('InitPos', '- List of initial positions for elastic DOFs'))
            f.write(_fw(', '.join(str(v) for v in dvr['InitVel'])) + '{:<22} {:}\n'.format('InitVel', '- List of initial velocities for elastic DOFs'))
            for i, row in enumerate(dvr['GFScaling']):
                f.write(_fw(', '.join(str(v) for v in row)) + '{:<22} {:}\n'.format('GFScalingL{}'.format(i+1), '- Generalized force scaling factors'))
            for i, row in enumerate(dvr['MassMatrix']):
                f.write(_fw(', '.join(str(v) for v in row)) + '{:<22} {:}\n'.format('MassMatrixL{}'.format(i+1), '- Mass matrix'))
            for i, row in enumerate(dvr['DampMatrix']):
                f.write(_fw(', '.join(str(v) for v in row)) + '{:<22} {:}\n'.format('DampMatrixL{}'.format(i+1), '- Damping matrix'))
            for i, row in enumerate(dvr['StifMatrix']):
                f.write(_fw(', '.join(str(v) for v in row)) + '{:<22} {:}\n'.format('StifMatrixL{}'.format(i+1), '- Stiffness matrix'))
            f.write(_fw(dvr['Twist']) + '{:<22} {:}\n'.format('Twist', '- Fixed twist of the section (deg)'))
            f.write(_fw(dvr['InflowMod']) + '{:<22} {:}\n'.format('InflowMod', '- Model for the inflow velocity'))
            f.write(_fw(', '.join(str(v) for v in dvr['Inflow'])) + '{:<22} {:}\n'.format('Inflow', '- Inflow velocity in x and y direction'))
            f.write(_fw('"' + dvr['InflowTSFile'] + '"') + '{:<22} {:}\n'.format('InflowTSFile', '- Input file for inflow velocity time series'))
            f.write(_fw(dvr['MotionMod']) + '{:<22} {:}\n'.format('MotionMod', '- Model for the motion of the degrees of freedom'))
            f.write(_fw('"' + dvr['MotionTSFile'] + '"') + '{:<22} {:}\n'.format('MotionTSFile', '- Input file for prescribed motion'))
            f.write('------------------- OUTPUT CONTROL --------------------------------------------\n')
            f.write('{!s:<14} {:<22} {:}\n'.format(dvr['SumPrint'], 'SumPrint', '- Write unsteady aerodynamics summary file (flag)'))
            f.write('{!s:<14} {:<22} {:}\n'.format(dvr['WrAFITables'], 'WrAFITables', '- Write back the aerodynamic coefficients used internally (flag)'))
            f.write('END of driver input file\n')


def _read_csv_vec(line: str) -> list:
    """Parse comma-separated numeric values from the start of a line."""
    raw = line.split()[0]
    # Handle case where values span multiple tokens (e.g. "1, 10")
    # by re-joining and splitting on comma
    tokens = line.strip().split()
    # Collect numeric comma-separated tokens until we hit a non-numeric word
    numeric_part = []
    for t in tokens:
        clean = t.rstrip(',').lstrip(',')
        try:
            float(clean)
            numeric_part.append(t)
        except ValueError:
            break
    joined = ' '.join(numeric_part).replace(',', ' ')
    return [float_read(p) for p in joined.split() if p]


def _read_csv_bools(line: str) -> list:
    """Parse comma-separated bool values like 'T, T, T  ActiveDOF ...'."""
    tokens = line.strip().split()
    bools = []
    for t in tokens:
        clean = t.rstrip(',').lstrip(',')
        if clean.upper() in ('T', 'TRUE'):
            bools.append(True)
        elif clean.upper() in ('F', 'FALSE'):
            bools.append(False)
        else:
            break
    return bools


def _read_3x3(f) -> List[List[float]]:
    """Read 3 rows of 3 comma-separated floats each."""
    mat: List[List[float]] = []
    for _ in range(3):
        line = f.readline()
        tokens = line.strip().split()
        numeric_part = []
        for t in tokens:
            clean = t.rstrip(',').lstrip(',')
            try:
                float(clean)
                numeric_part.append(clean)
            except ValueError:
                break
        mat.append([float_read(p) for p in numeric_part])
    return mat
