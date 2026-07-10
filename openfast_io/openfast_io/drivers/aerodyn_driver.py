"""AeroDyn standalone driver — reads/writes AeroDyn driver input files (.dvr).

Supports all three analysis types:
  1 = multiple turbines, one simulation
  2 = one turbine, time-dependent simulation
  3 = one turbine, combined cases

Both BasicHAWTFormat and advanced (generic) turbine geometry are supported.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..io.aerodyn import AeroDynIO
from ..io.inflowwind import InflowWindIO
from ..io.seastate import SeaStateIO
from ..outlist import capture_outlist
from ..parsing import bool_read, float_read, fmt_field, int_read, quoted_read

try:
    from ..FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}


def _fw(val, width: int = 22) -> str:
    """Format a value field with guaranteed 2-space separator."""
    return fmt_field(val, min_width=width)


class AeroDynStandaloneDriver:
    """Reads and writes AeroDyn standalone driver input files (.dvr)."""

    def __init__(self) -> None:
        self._aerodyn = AeroDynIO()
        self._inflowwind = InflowWindIO()
        self._seastate = SeaStateIO()

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------
    def read(self, dvr_path: Path) -> dict:
        """Read an AeroDyn driver file and all referenced module files.

        Returns a dict with keys:
          ``'AeroDynDriver'`` — driver-level parameters
          ``'AeroDyn'``, ``'AeroDynBlade'``, etc. — from AeroDynIO
          ``'InflowWind'`` — from InflowWindIO (if referenced)
        """
        dvr_path = Path(dvr_path)
        base_dir = dvr_path.parent
        dvr: Dict[str, Any] = {}

        with open(dvr_path) as f:
            # --- Header ---
            f.readline()  # title line
            dvr['description'] = f.readline().rstrip()

            # --- Input Configuration ---
            f.readline()  # section header
            dvr['Echo'] = bool_read(f.readline().split()[0])
            dvr['MHK'] = int_read(f.readline().split()[0])
            dvr['AnalysisType'] = int_read(f.readline().split()[0])
            dvr['TMax'] = float_read(f.readline().split()[0])
            dvr['DT'] = float_read(f.readline().split()[0])
            dvr['AeroFile'] = quoted_read(f.readline().split()[0])

            # --- Environmental Conditions ---
            f.readline()  # section header
            dvr['FldDens'] = float_read(f.readline().split()[0])
            dvr['KinVisc'] = float_read(f.readline().split()[0])
            dvr['SpdSound'] = float_read(f.readline().split()[0])
            dvr['Patm'] = float_read(f.readline().split()[0])
            dvr['Pvap'] = float_read(f.readline().split()[0])
            dvr['WtrDpth'] = float_read(f.readline().split()[0])

            # --- Inflow Data ---
            f.readline()  # section header
            dvr['CompInflow'] = int_read(f.readline().split()[0])
            dvr['InflowFile'] = quoted_read(f.readline().split()[0])
            dvr['HWindSpeed'] = float_read(f.readline().split()[0])
            dvr['RefHt'] = float_read(f.readline().split()[0])
            dvr['PLExp'] = float_read(f.readline().split()[0])

            # --- SeaState Data ---
            f.readline()  # section header
            dvr['CompSeaSt'] = int_read(f.readline().split()[0])
            dvr['SeaStFile'] = quoted_read(f.readline().split()[0])

            # --- Turbine Data ---
            f.readline()  # section header
            dvr['NumTurbines'] = int_read(f.readline().split()[0])

            # --- Per-turbine data ---
            turbines: List[Dict[str, Any]] = []
            for i_turb in range(dvr['NumTurbines']):
                turb: Dict[str, Any] = {}
                f.readline()  # turbine section header

                ln = f.readline().split()
                turb['BasicHAWTFormat'] = bool_read(ln[0])

                # Read origin common to both formats
                turb['BaseOriginInit'] = _read_csv_vec(f.readline())

                if turb['BasicHAWTFormat']:
                    # Basic HAWT format: 8 more lines
                    turb['NumBlades'] = int_read(f.readline().split()[0])
                    turb['HubRad'] = float_read(f.readline().split()[0])
                    turb['HubHt'] = float_read(f.readline().split()[0])
                    turb['Overhang'] = float_read(f.readline().split()[0])
                    turb['ShftTilt'] = float_read(f.readline().split()[0])
                    turb['Precone'] = float_read(f.readline().split()[0])
                    turb['Twr2Shft'] = float_read(f.readline().split()[0])
                else:
                    # Advanced (generic) format
                    turb['BaseOrientationInit'] = _read_csv_vec(f.readline())
                    turb['HasTower'] = bool_read(f.readline().split()[0])
                    turb['HAWTprojection'] = bool_read(f.readline().split()[0])
                    turb['TwrOrigin_t'] = _read_csv_vec(f.readline())
                    turb['NacOrigin_t'] = _read_csv_vec(f.readline())
                    turb['HubOrigin_n'] = _read_csv_vec(f.readline())
                    turb['HubOrientation_n'] = _read_csv_vec(f.readline())

                    # Blades section
                    f.readline()  # blade section header
                    turb['NumBlades'] = int_read(f.readline().split()[0])
                    n_bld = turb['NumBlades']
                    turb['BldOrigin_h'] = [_read_csv_vec(f.readline()) for _ in range(n_bld)]
                    turb['BldOrientation_h'] = [_read_csv_vec(f.readline()) for _ in range(n_bld)]
                    turb['BldHubRad_bl'] = [float_read(f.readline().split()[0]) for _ in range(n_bld)]

                # --- Motion (AnalysisType=1) ---
                f.readline()  # motion section header
                if turb['BasicHAWTFormat']:
                    turb['BaseMotionType'] = int_read(f.readline().split()[0])
                    turb['DegreeOfFreedom'] = int_read(f.readline().split()[0])
                    turb['Amplitude'] = float_read(f.readline().split()[0])
                    turb['Frequency'] = float_read(f.readline().split()[0])
                    turb['BaseMotionFileName'] = quoted_read(f.readline().split()[0])
                    turb['NacYaw'] = float_read(f.readline().split()[0])
                    turb['RotSpeed'] = float_read(f.readline().split()[0])
                    turb['BldPitch'] = float_read(f.readline().split()[0])
                else:
                    turb['BaseMotionType'] = int_read(f.readline().split()[0])
                    turb['DegreeOfFreedom'] = int_read(f.readline().split()[0])
                    turb['Amplitude'] = float_read(f.readline().split()[0])
                    turb['Frequency'] = float_read(f.readline().split()[0])
                    turb['BaseMotionFileName'] = quoted_read(f.readline().split()[0])
                    turb['NacMotionType'] = int_read(f.readline().split()[0])
                    turb['NacYaw'] = float_read(f.readline().split()[0])
                    turb['NacMotionFileName'] = quoted_read(f.readline().split()[0])
                    turb['RotMotionType'] = int_read(f.readline().split()[0])
                    turb['RotSpeed'] = float_read(f.readline().split()[0])
                    turb['RotMotionFileName'] = quoted_read(f.readline().split()[0])
                    turb['BldMotionType'] = int_read(f.readline().split()[0])
                    n_bld = turb['NumBlades']
                    turb['BldPitch'] = [float_read(f.readline().split()[0]) for _ in range(n_bld)]
                    turb['BldMotionFileName'] = [quoted_read(f.readline().split()[0]) for _ in range(n_bld)]

                turbines.append(turb)

            dvr['Turbines'] = turbines

            # --- Time-dependent Analysis ---
            f.readline()  # section header
            dvr['TimeAnalysisFileName'] = quoted_read(f.readline().split()[0])

            # --- Combined-Case Analysis ---
            f.readline()  # section header
            dvr['NumCases'] = int_read(f.readline().split()[0])
            f.readline()  # column headers
            f.readline()  # column units
            cases: List[Dict[str, float]] = []
            for _ in range(dvr['NumCases']):
                ln = f.readline().split()
                cases.append({
                    'HWndSpeed': float_read(ln[0]),
                    'PLExp': float_read(ln[1]),
                    'RotSpd': float_read(ln[2]),
                    'Pitch': float_read(ln[3]),
                    'Yaw': float_read(ln[4]),
                    'dT': float_read(ln[5]),
                    'Tmax': float_read(ln[6]),
                    'DOF': int_read(ln[7]),
                    'Amplitude': float_read(ln[8]),
                    'Frequency': float_read(ln[9]),
                })
            dvr['Cases'] = cases

            # --- Output Settings ---
            f.readline()  # section header
            dvr['OutFmt'] = quoted_read(f.readline().split()[0])
            dvr['OutFileFmt'] = int_read(f.readline().split()[0])
            dvr['WrVTK'] = int_read(f.readline().split()[0])
            dvr['WrVTK_Type'] = int_read(f.readline().split()[0])
            dvr['VTKHubRad'] = float_read(f.readline().split()[0])
            dvr['VTKNacDim'] = _read_csv_vec(f.readline())

        result: Dict[str, Any] = {'AeroDynDriver': dvr}

        # Shared OutList registry (mirrors OpenFASTDriver.read) so the AeroDyn
        # and InflowWind OutList sections survive a standalone read → write
        # roundtrip.
        outlist: Dict[str, Any] = copy.deepcopy(FstOutput) if FstOutput else {}

        def _cap(f, module, freeform=False):
            return capture_outlist(f, outlist, module, freeform=freeform)

        # --- Delegate to AeroDynIO for the primary AeroDyn file ---
        aero_file = dvr.get('AeroFile', '')
        aero_path = os.path.normpath(os.path.join(str(base_dir), aero_file))
        if aero_file and os.path.isfile(aero_path):
            n_bld = turbines[0]['NumBlades'] if turbines else 3
            try:
                ad_data = self._aerodyn.read(Path(aero_path), base_dir, num_blades=n_bld,
                                              outlist=outlist, read_outlist_fn=_cap)
                result.update(ad_data)
            except Exception:
                pass  # module file parse failure — driver data still valid

        # --- Delegate to InflowWindIO ---
        if dvr.get('CompInflow', 0) == 1:
            ifw_file = dvr.get('InflowFile', '')
            ifw_path = os.path.normpath(os.path.join(str(base_dir), ifw_file))
            if ifw_file and os.path.isfile(ifw_path):
                try:
                    ifw_data = self._inflowwind.read(Path(ifw_path), base_dir,
                                                      outlist=outlist, read_outlist_fn=_cap)
                    result.update(ifw_data)
                except Exception:
                    pass

        result['outlist'] = outlist
        return result

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------
    def write(self, data: dict, dvr_path: Path) -> None:
        """Write an AeroDyn driver file and all referenced module files."""
        dvr_path = Path(dvr_path)
        base_dir = dvr_path.parent
        dvr = data['AeroDynDriver']

        with open(dvr_path, 'w') as f:
            f.write('----- AeroDyn Driver Input File ---------------------------------------------------------\n')
            f.write(dvr.get('description', 'Generated by OpenFAST_IO') + '\n')
            f.write('----- Input Configuration ---------------------------------------------------------------\n')
            f.write('{!s:<16} {:<12} {:}\n'.format(dvr['Echo'], 'Echo', '- Echo input parameters to "<rootname>.ech"?'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['MHK'], 'MHK', '- MHK turbine type (switch)'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['AnalysisType'], 'AnalysisType', '- {1: multiple turbines, 2: time-dependent, 3: combined cases}'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['TMax'], 'TMax', '- Total run time (s)'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['DT'], 'DT', '- Simulation time step (s)'))

            ad_file = dvr.get('AeroFile', 'AeroDyn.dat')
            f.write('{:<16} {:<12} {:}\n'.format('"' + ad_file + '"', 'AeroFile', '- Name of the primary AeroDyn input file'))
            f.write('----- Environmental Conditions ----------------------------------------------------------\n')
            f.write('{:<26} {:<13} {:}\n'.format(dvr['FldDens'], 'FldDens', '- Density of working fluid (kg/m^3)'))
            f.write('{:<26} {:<13} {:}\n'.format(dvr['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)'))
            f.write('{:<26} {:<13} {:}\n'.format(dvr['SpdSound'], 'SpdSound', '- Speed of sound in working fluid (m/s)'))
            f.write('{:<26} {:<13} {:}\n'.format(dvr['Patm'], 'Patm', '- Atmospheric pressure (Pa)'))
            f.write('{:<26} {:<13} {:}\n'.format(dvr['Pvap'], 'Pvap', '- Vapour pressure of working fluid (Pa)'))
            f.write('{:<26} {:<13} {:}\n'.format(dvr['WtrDpth'], 'WtrDpth', '- Water depth (m)'))
            f.write('----- Inflow Data -----------------------------------------------------------------------\n')
            f.write('{:<16} {:<12} {:}\n'.format(dvr['CompInflow'], 'CompInflow', '- Compute inflow wind velocities (switch)'))
            ifw_file = dvr.get('InflowFile', 'unused')
            f.write('{:<16} {:<12} {:}\n'.format('"' + ifw_file + '"', 'InflowFile', '- Name of the InflowWind input file'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['HWindSpeed'], 'HWindSpeed', '- Horizontal wind speed (m/s)'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['RefHt'], 'RefHt', '- Reference height for horizontal wind speed (m)'))
            f.write('{:<16} {:<12} {:}\n'.format(dvr['PLExp'], 'PLExp', '- Power law exponent (-)'))
            f.write('----- SeaState Data [used only when MHK = 1 or 2] ---------------------------------------\n')
            f.write('{:<16} {:<27} {:}\n'.format(dvr['CompSeaSt'], 'CompSeaSt', '- Compute wave velocities (switch)'))
            ss_file = dvr.get('SeaStFile', 'unused')
            f.write('{:<16} {:<27} {:}\n'.format('"' + ss_file + '"', 'SeaStFile', '- Name of the SeaState input file'))
            f.write('----- Turbine Data ----------------------------------------------------------------------\n')
            f.write('{:<16} {:<12} {:}\n'.format(dvr['NumTurbines'], 'NumTurbines', '- Number of turbines'))

            for i_turb, turb in enumerate(dvr.get('Turbines', [])):
                f.write('----- Turbine({}) ------------------------------------------------------------------------\n'.format(i_turb + 1))
                f.write('{!s:<16} {:<22} {:}\n'.format(turb['BasicHAWTFormat'], 'BasicHAWTFormat({})'.format(i_turb + 1), '- Flag to switch between basic or generic input format'))
                f.write(_fw(','.join(str(v) for v in turb['BaseOriginInit'])) + '{:<22} {:}\n'.format('BaseOriginInit({})'.format(i_turb + 1), '- x,y,z coordinates of base origin (m)'))

                if turb['BasicHAWTFormat']:
                    f.write(_fw(turb['NumBlades']) + '{:<22} {:}\n'.format('NumBlades({})'.format(i_turb + 1), '- Number of blades'))
                    f.write(_fw(turb['HubRad']) + '{:<22} {:}\n'.format('HubRad({})'.format(i_turb + 1), '- Hub radius (m)'))
                    f.write(_fw(turb['HubHt']) + '{:<22} {:}\n'.format('HubHt({})'.format(i_turb + 1), '- Hub height (m)'))
                    f.write(_fw(turb['Overhang']) + '{:<22} {:}\n'.format('Overhang({})'.format(i_turb + 1), '- Overhang (m)'))
                    f.write(_fw(turb['ShftTilt']) + '{:<22} {:}\n'.format('ShftTilt({})'.format(i_turb + 1), '- Shaft tilt (deg)'))
                    f.write(_fw(turb['Precone']) + '{:<22} {:}\n'.format('Precone({})'.format(i_turb + 1), '- Precone (deg)'))
                    f.write(_fw(turb['Twr2Shft']) + '{:<22} {:}\n'.format('Twr2Shft({})'.format(i_turb + 1), '- Twr2Shft (m)'))
                else:
                    f.write(_fw(','.join(str(v) for v in turb['BaseOrientationInit'])) + '{:<22} {:}\n'.format('BaseOrientationInit({})'.format(i_turb + 1), '- successive rotations defining initial orientation of the base frame (deg)'))
                    f.write('{!s:<16} {:<22} {:}\n'.format(turb['HasTower'], 'HasTower({})'.format(i_turb + 1), '- True if turbine has a tower (flag)'))
                    f.write('{!s:<16} {:<22} {:}\n'.format(turb['HAWTprojection'], 'HAWTprojection({})'.format(i_turb + 1), '- True if turbine is a horizontal axis turbine (flag)'))
                    f.write(_fw(','.join(str(v) for v in turb['TwrOrigin_t'])) + '{:<22} {:}\n'.format('TwrOrigin_t({})'.format(i_turb + 1), '- Coordinate of tower base in base coordinates (m)'))
                    f.write(_fw(','.join(str(v) for v in turb['NacOrigin_t'])) + '{:<22} {:}\n'.format('NacOrigin_t({})'.format(i_turb + 1), '- x,y,z coordinates of nacelle origin from base (m)'))
                    f.write(_fw(','.join(str(v) for v in turb['HubOrigin_n'])) + '{:<22} {:}\n'.format('HubOrigin_n({})'.format(i_turb + 1), '- x,y,z coordinates of hub origin from nacelle (m)'))
                    f.write(_fw(','.join(str(v) for v in turb['HubOrientation_n'])) + '{:<22} {:}\n'.format('HubOrientation_n({})'.format(i_turb + 1), '- successive rotations defining hub frame from nacelle frame (deg)'))

                    f.write('----- Turbine({}) Blades -----------------------------------------------------------------\n'.format(i_turb + 1))
                    f.write(_fw(turb['NumBlades']) + '{:<22} {:}\n'.format('NumBlades({})'.format(i_turb + 1), '- Number of blades for current rotor (-)'))
                    n_bld = turb['NumBlades']
                    for j in range(n_bld):
                        f.write(_fw(','.join(str(v) for v in turb['BldOrigin_h'][j])) + '{:<22} {:}\n'.format('BldOrigin_h({0}_{1})'.format(i_turb + 1, j + 1), '- Orign of blade {:d} wrt. hub origin in hub coordinates (m)'.format(j + 1)))
                    for j in range(n_bld):
                        f.write(_fw(','.join(str(v) for v in turb['BldOrientation_h'][j])) + '{:<22} {:}\n'.format('BldOrientation_h({0}_{1})'.format(i_turb + 1, j + 1), '- successive rotations defining blade {:d} frame from hub frame (deg)'.format(j + 1)))
                    for j in range(n_bld):
                        f.write(_fw(turb['BldHubRad_bl'][j]) + '{:<22} {:}\n'.format('BldHubRad_bl({0}_{1})'.format(i_turb + 1, j + 1), '- z-offset where radial input data start for blade {:d} (m)'.format(j + 1)))

                f.write('----- Turbine({}) Motion [used only when AnalysisType=1] ---------------------------------\n'.format(i_turb + 1))
                if turb['BasicHAWTFormat']:
                    f.write('{:<16} {:<28} {:}\n'.format(turb['BaseMotionType'], 'BaseMotionType({})'.format(i_turb + 1), '- Type of motion prescribed for this base (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['DegreeOfFreedom'], 'DegreeOfFreedom({})'.format(i_turb + 1), '- Degree of freedom (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['Amplitude'], 'Amplitude({})'.format(i_turb + 1), '- Amplitude of sinusoidal motion (m or rad)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['Frequency'], 'Frequency({})'.format(i_turb + 1), '- Frequency of sinusoidal motion (Hz)'))
                    f.write('{:<16} {:<28} {:}\n'.format('"' + turb['BaseMotionFileName'] + '"', 'BaseMotionFileName({})'.format(i_turb + 1), '- Filename for arbitrary base motion'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['NacYaw'], 'NacYaw({})'.format(i_turb + 1), '- Yaw angle of the nacelle (deg)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['RotSpeed'], 'RotSpeed({})'.format(i_turb + 1), '- Rotational speed (rpm)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['BldPitch'], 'BldPitch({})'.format(i_turb + 1), '- Blade pitch (deg)'))
                else:
                    f.write('{:<16} {:<28} {:}\n'.format(turb['BaseMotionType'], 'BaseMotionType({})'.format(i_turb + 1), '- Type of motion prescribed for this base (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb.get('DegreeOfFreedom', 1), 'DegreeOfFreedom({})'.format(i_turb + 1), '- {1:xt, 2:yt, 3:zt, 4:theta_xt, 5:theta_yt, 6:theta_zt} (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['Amplitude'], 'Amplitude({})'.format(i_turb + 1), '- Amplitude of sinusoidal motion (m or rad)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['Frequency'], 'Frequency({})'.format(i_turb + 1), '- Frequency of sinusoidal motion (Hz)'))
                    f.write('{:<16} {:<28} {:}\n'.format('"' + turb['BaseMotionFileName'] + '"', 'BaseMotionFileName({})'.format(i_turb + 1), '- Filename for arbitrary base motion'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['NacMotionType'], 'NacMotionType({})'.format(i_turb + 1), '- Type of motion prescribed for the nacelle (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['NacYaw'], 'NacYaw({})'.format(i_turb + 1), '- Yaw angle of the nacelle (deg)'))
                    f.write('{:<16} {:<28} {:}\n'.format('"' + turb['NacMotionFileName'] + '"', 'NacMotionFileName({})'.format(i_turb + 1), '- Filename for yaw motion'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['RotMotionType'], 'RotMotionType({})'.format(i_turb + 1), '- Type of motion prescribed for this rotor (flag)'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['RotSpeed'], 'RotSpeed({})'.format(i_turb + 1), '- Rotational speed (rpm)'))
                    f.write('{:<16} {:<28} {:}\n'.format('"' + turb['RotMotionFileName'] + '"', 'RotMotionFileName({})'.format(i_turb + 1), '- Filename for rotor motion'))
                    f.write('{:<16} {:<28} {:}\n'.format(turb['BldMotionType'], 'BldMotionType({})'.format(i_turb + 1), '- Type of pitch motion prescribed for the blades (flag)'))
                    n_bld = turb['NumBlades']
                    for j in range(n_bld):
                        f.write('{:<16} {:<28} {:}\n'.format(turb['BldPitch'][j], 'BldPitch({0}_{1})'.format(i_turb + 1, j + 1), '- Blade {:d} pitch (deg)'.format(j + 1)))
                    for j in range(n_bld):
                        f.write('{:<16} {:<28} {:}\n'.format('"' + turb['BldMotionFileName'][j] + '"', 'BldMotionFileName({0}_{1})'.format(i_turb + 1, j + 1), '- Filename containing blade pitch motion'))

            f.write('----- Time-dependent Analysis [used only when AnalysisType=2, numTurbines=1] ------------\n')
            f.write('{:<17} {:<22} {:}\n'.format('"' + dvr.get('TimeAnalysisFileName', 'unused') + '"', 'TimeAnalysisFileName', '- Filename containing time series'))
            f.write('----- Combined-Case Analysis [used only when AnalysisType=3, numTurbines=1] -------------\n')
            f.write('{:<5} {:<13} {:}\n'.format(dvr.get('NumCases', 0), 'NumCases', '- Number of cases to run'))
            f.write('HWndSpeed    PLExp        RotSpd       Pitch        Yaw     dT      Tmax  DOF  Amplitude Frequency\n')
            f.write('(m/s)        (-)          (rpm)        (deg)        (deg)   (s)     (s)    (-)  (m or rad) (Hz)\n')
            for case in dvr.get('Cases', []):
                f.write('{:<13} {:<13} {:<13} {:<13} {:<8} {:<8} {:<6} {:<5} {:<10} {}\n'.format(
                    case['HWndSpeed'], case['PLExp'], case['RotSpd'], case['Pitch'],
                    case['Yaw'], case['dT'], case['Tmax'], case['DOF'],
                    case['Amplitude'], case['Frequency']))
            f.write('----- Output Settings -------------------------------------------------------------------\n')
            f.write('{:<12} {:<12} {:}\n'.format('"' + dvr.get('OutFmt', 'ES15.8E2') + '"', 'OutFmt', '- Format used for text tabular output'))
            f.write('{:<12} {:<12} {:}\n'.format(dvr['OutFileFmt'], 'OutFileFmt', '- Format for tabular output file'))
            f.write('{:<12} {:<12} {:}\n'.format(dvr['WrVTK'], 'WrVTK', '- VTK visualization data output'))
            f.write('{:<12} {:<12} {:}\n'.format(dvr['WrVTK_Type'], 'WrVTK_Type', '- VTK visualization data type'))
            f.write('{:<12} {:<12} {:}\n'.format(dvr['VTKHubRad'], 'VTKHubRad', '- HubRadius for VTK visualization (m)'))
            f.write(_fw(','.join(str(v) for v in dvr['VTKNacDim'])) + '{:<12} {:}\n'.format('VTKNacDim', '- Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)'))

        outlist = data.get('outlist') or (copy.deepcopy(FstOutput) if FstOutput else {})

        # --- Delegate to AeroDynIO ---
        if 'AeroDyn' in data:
            n_bld = dvr['Turbines'][0]['NumBlades'] if dvr.get('Turbines') else 3
            ad_path = os.path.normpath(os.path.join(str(base_dir), ad_file))
            ad_write_data = {k: data[k] for k in ('AeroDyn', 'AeroDynBlade', 'af_data', 'af_coord', 'ac') if k in data}
            try:
                self._aerodyn.write(ad_write_data, Path(ad_path), base_dir,
                                    naming_out=dvr_path.stem, outlist=outlist)
            except Exception:
                pass

        # --- Delegate to InflowWindIO ---
        if dvr.get('CompInflow', 0) == 1 and 'InflowWind' in data:
            ifw_path = os.path.normpath(os.path.join(str(base_dir), ifw_file))
            try:
                self._inflowwind.write({'InflowWind': data['InflowWind']}, Path(ifw_path), base_dir, outlist=outlist)
            except Exception:
                pass


def _read_csv_vec(line: str) -> list:
    """Parse a comma- or space-separated vector from a line like '0,0,0  name ...'."""
    # Take everything before the first alphabetic token (the param name)
    raw = line.split()[0]
    parts = raw.split(',')
    return [float_read(p.strip()) for p in parts]
