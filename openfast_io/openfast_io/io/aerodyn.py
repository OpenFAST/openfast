"""AeroDyn module IO — reads and writes AeroDyn15 input files.

Extracted from FAST_reader.py and FAST_writer.py.  Covers:
  - Main AeroDyn15 input (.dat)
  - AeroDyn blade files
  - Airfoil polar files
  - Airfoil coordinate files
  - OLAF input file

All parsing logic is identical to the original; only the data target changes
(local dict vs self.fst_vt).
"""
from __future__ import annotations

import copy
import os
import random
import time
from pathlib import Path

import numpy as np

from .base import ModuleIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    fix_path,
    readline_filterComments,
)


# ---------------------------------------------------------------------------
# Writer helpers (copied from FAST_writer.py to avoid coupling)
# ---------------------------------------------------------------------------

def _float_default_out(val, trim=False):
    if isinstance(val, float):
        return '{:.4f}'.format(val) if trim else '{: 22f}'.format(val)
    else:
        return '{:}'.format(val) if trim else '{:<22}'.format(val)


def _get_outlist(outlist_dict, channel_list):
    """Extract True channels from an outlist dict for the given modules."""
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


class AeroDynIO(ModuleIO):
    """Reads and writes AeroDyn v15.03+ input files.

    read() returns::

        {
            'AeroDyn': { ... main params, tower arrays, af_data, af_coord, OLAF ... },
            'AeroDynBlade': [ {blade0}, {blade1}, {blade2} ]  or  {single_blade},
        }

    write() accepts data with keys 'AeroDyn' and 'AeroDynBlade', plus outlist
    and Fst context.  Writes the main AD file, blade files, polars, coords, OLAF.
    """

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------

    def read(self, file_path: Path, base_dir: Path, *, num_blades: int = 3,
             aero_file_path: str = '', outlist: dict | None = None,
             read_outlist_fn=None) -> dict:
        """Read AeroDyn main file and all sub-files.

        Parameters
        ----------
        file_path : path to AeroDyn main input file
        base_dir : root directory for resolving relative paths
        num_blades : number of blades (from ElastoDyn)
        aero_file_path : relative directory from base_dir where sub-files live
        outlist : optional outlist dict to populate via read_outlist_fn
        read_outlist_fn : callback(f, module_name) to read outlist sections
        """
        ad = {}
        ad_blade = [{}, {}, {}]
        file_path = str(file_path)
        base_dir = str(base_dir) if base_dir else ''

        f = open(file_path)

        # General Options
        f.readline(); f.readline(); f.readline()
        ad['Echo']          = bool_read(f.readline().split()[0])
        ad['DTAero']        = float_read(f.readline().split()[0])
        ad['Wake_Mod']      = int(f.readline().split()[0])
        ad['TwrPotent']     = int(f.readline().split()[0])
        ad['TwrShadow']     = int(f.readline().split()[0])
        ad['TwrAero']       = bool_read(f.readline().split()[0])
        ad['CavitCheck']    = bool_read(f.readline().split()[0])
        ad['NacelleDrag']   = bool_read(f.readline().split()[0])
        ad['CompAA']        = bool_read(f.readline().split()[0])
        ad['AA_InputFile']  = f.readline().split()[0]

        # Environmental Conditions
        f.readline()
        ad['AirDens']       = float_read(f.readline().split()[0])
        ad['KinVisc']       = float_read(f.readline().split()[0])
        ad['SpdSound']      = float_read(f.readline().split()[0])
        ad['Patm']          = float_read(f.readline().split()[0])
        ad['Pvap']          = float_read(f.readline().split()[0])

        f.readline()
        ad['BEM_Mod']       = int(f.readline().split()[0])

        # BEM Options
        f.readline()
        ad['Skew_Mod']              = int_read(f.readline().split()[0])
        ad['SkewMomCorr']           = bool_read(f.readline().split()[0])
        ad['SkewRedistr_Mod']       = int_read(f.readline().split()[0])
        ad['SkewRedistrFactor']     = float_read(f.readline().split()[0])
        f.readline()
        ad['TipLoss']               = bool_read(f.readline().split()[0])
        ad['HubLoss']               = bool_read(f.readline().split()[0])
        ad['TanInd']                = bool_read(f.readline().split()[0])
        ad['AIDrag']                = bool_read(f.readline().split()[0])
        ad['TIDrag']                = bool_read(f.readline().split()[0])
        ad['IndToler']              = float_read(f.readline().split()[0])
        ad['MaxIter']               = int(f.readline().split()[0])
        f.readline()
        ad['SectAvg']               = bool_read(f.readline().split()[0])
        ad['SectAvgWeighting']      = int_read(f.readline().split()[0])
        ad['SectAvgNPoints']        = int_read(f.readline().split()[0])
        ad['SectAvgPsiBwd']         = float_read(f.readline().split()[0])
        ad['SectAvgPsiFwd']         = float_read(f.readline().split()[0])

        # Dynamic BEM
        f.readline()
        ad['DBEMT_Mod']     = int(f.readline().split()[0])
        ad['tau1_const']    = float_read(f.readline().split()[0])

        # OLAF
        f.readline()
        ad['OLAFInputFileName'] = quoted_read(f.readline().split()[0])

        # Unsteady Airfoil Aero
        f.readline()
        ad['AoA34']             = bool_read(f.readline().split()[0])
        ad['UA_Mod']            = int(f.readline().split()[0])
        ad['FLookup']           = bool_read(f.readline().split()[0])
        ad['IntegrationMethod'] = int(f.readline().split()[0])

        file_pos = f.tell()
        line = f.readline()
        if 'UAStartRad' in line:
            ad['UAStartRad'] = float_read(line.split()[0])
        else:
            f.seek(file_pos)

        file_pos = f.tell()
        line = f.readline()
        if 'UAEndRad' in line:
            ad['UAEndRad'] = float_read(line.split()[0])
        else:
            f.seek(file_pos)

        # Airfoil Information
        f.readline()
        ad['AFTabMod']      = int(f.readline().split()[0])
        ad['InCol_Alfa']    = int(f.readline().split()[0])
        ad['InCol_Cl']      = int(f.readline().split()[0])
        ad['InCol_Cd']      = int(f.readline().split()[0])
        ad['InCol_Cm']      = int(f.readline().split()[0])
        ad['InCol_Cpmin']   = int(f.readline().split()[0])
        ad['NumAFfiles']    = int(f.readline().split()[0])
        ad['AFNames']       = [None] * ad['NumAFfiles']
        for i in range(ad['NumAFfiles']):
            af_filename = fix_path(f.readline().split()[0])[1:-1]
            ad['AFNames'][i] = os.path.abspath(os.path.join(base_dir, aero_file_path, af_filename))

        # Rotor/Blade Properties
        f.readline()
        ad['UseBlCm']  = bool_read(f.readline().split()[0])
        ad['ADBlFile1'] = quoted_read(f.readline().split()[0])
        ad['ADBlFile2'] = quoted_read(f.readline().split()[0])
        ad['ADBlFile3'] = quoted_read(f.readline().split()[0])

        # Hub, nacelle, tail fin
        f.readline()
        ad['VolHub']    = float_read(f.readline().split()[0])
        ad['HubCenBx']  = float_read(f.readline().split()[0])
        f.readline()
        ad['VolNac']    = float_read(f.readline().split()[0])
        ad['NacCenB']   = [idx.strip() for idx in f.readline().split('NacCenB')[0].split(',')]
        ad['NacArea']   = [idx.strip() for idx in f.readline().split('NacArea')[0].split(',')]
        ad['NacCd']     = [idx.strip() for idx in f.readline().split('NacCd')[0].split(',')]
        ad['NacDragAC'] = [idx.strip() for idx in f.readline().split('NacDragAC')[0].split(',')]
        f.readline()
        ad['TFinAero']  = bool_read(f.readline().split()[0])
        tfa_filename    = fix_path(f.readline().split()[0])[1:-1]
        ad['TFinFile']  = os.path.abspath(os.path.join(base_dir, tfa_filename))

        # Tower Influence and Aerodynamics
        f.readline()
        ad['NumTwrNds'] = int(f.readline().split()[0])
        f.readline(); f.readline()
        for k in ('TwrElev', 'TwrDiam', 'TwrCd', 'TwrTI', 'TwrCb', 'TwrCp', 'TwrCa'):
            ad[k] = [None] * ad['NumTwrNds']
        for i in range(ad['NumTwrNds']):
            data = [float(val) for val in f.readline().split()]
            ad['TwrElev'][i] = data[0]
            ad['TwrDiam'][i] = data[1]
            ad['TwrCd'][i]   = data[2]
            ad['TwrTI'][i]   = data[3]
            ad['TwrCb'][i]   = data[4]
            ad['TwrCp'][i]   = data[5]
            ad['TwrCa'][i]   = data[6]

        # Outputs
        f.readline()
        ad['SumPrint']  = bool_read(f.readline().split()[0])
        ad['NBlOuts']   = int(f.readline().split()[0])
        ad['BlOutNd']   = [idx.strip() for idx in f.readline().split('BlOutNd')[0].split(',')]
        ad['NTwOuts']   = int(f.readline().split()[0])
        ad['TwOutNd']   = [idx.strip() for idx in f.readline().split('TwOutNd')[0].split(',')]

        # AeroDyn OutList
        f.readline()
        if read_outlist_fn is not None and outlist is not None:
            read_outlist_fn(f, 'AeroDyn')
        else:
            # skip outlist section
            line = f.readline()
            while line and 'END' not in line.split('!')[0].upper()[:3]:
                line = f.readline()

        # Optional nodal output
        try:
            f.readline()
            ad['BldNd_BladesOut'] = int(f.readline().split()[0])
            ad['BldNd_BlOutNd']   = f.readline().split()[0]
            f.readline()
            if read_outlist_fn is not None and outlist is not None:
                read_outlist_fn(f, 'AeroDyn_Nodes')
            else:
                line = f.readline()
                while line and 'END' not in line.split('!')[0].upper()[:3]:
                    line = f.readline()
        except Exception:
            pass

        f.close()

        # ── Read blade files ──
        ad_bld_file1 = os.path.join(base_dir, aero_file_path, ad['ADBlFile1'])
        ad_bld_file2 = os.path.join(base_dir, aero_file_path, ad['ADBlFile2'])
        ad_bld_file3 = os.path.join(base_dir, aero_file_path, ad['ADBlFile3'])

        if ad_bld_file1 == ad_bld_file2 and ad_bld_file1 == ad_bld_file3:
            self._read_blade(ad_blade, ad_bld_file1, 0)
            ad_blade = ad_blade[0]
        elif num_blades == 2 and ad_bld_file1 == ad_bld_file2:
            self._read_blade(ad_blade, ad_bld_file1, 0)
            ad_blade = ad_blade[0]
        else:
            self._read_blade(ad_blade, ad_bld_file1, 0)
            if num_blades > 1:
                self._read_blade(ad_blade, ad_bld_file2, 1)
            if num_blades > 2:
                self._read_blade(ad_blade, ad_bld_file3, 2)
            else:
                ad_blade = ad_blade[0]

        # ── Read polars ──
        self._read_polars(ad)

        # ── Read coords ──
        self._read_coords(ad)

        # ── Read OLAF if present ──
        olaf_filename = os.path.join(base_dir, ad['OLAFInputFileName'])
        if os.path.isfile(olaf_filename):
            self._read_olaf(ad, olaf_filename, base_dir)

        return {'AeroDyn': ad, 'AeroDynBlade': ad_blade}

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------

    def write(self, data: dict, file_path: Path, base_dir: Path, *,
              naming_out: str = 'openfast', outlist: dict | None = None) -> None:
        """Write AeroDyn main file and all sub-files.

        Parameters
        ----------
        data : dict with 'AeroDyn', 'AeroDynBlade' keys (and optionally 'outlist')
        file_path : output path for main AeroDyn .dat file
        base_dir : directory where sub-files will be written
        naming_out : base naming for generated files
        outlist : optional outlist dict for channel output sections
        """
        ad = data['AeroDyn']
        ad_blade = data['AeroDynBlade']
        run_dir = str(base_dir)

        # ── Write blade files ──
        if isinstance(ad_blade, list):
            for i_bld, _ in enumerate(ad_blade):
                ad['ADBlFile%d' % (i_bld + 1)] = naming_out + '_AeroDyn_blade_%d.dat' % (i_bld + 1)
                self._write_blade(ad, ad_blade, run_dir, bld_ind=i_bld)
        elif isinstance(ad_blade, dict):
            ad['ADBlFile1'] = naming_out + '_AeroDyn_blade.dat'
            ad['ADBlFile2'] = ad['ADBlFile1']
            ad['ADBlFile3'] = ad['ADBlFile1']
            self._write_blade(ad, ad_blade, run_dir)

        # ── Write polars ──
        self._write_polars(ad, run_dir, naming_out)

        # ── Write coords ──
        if any(ad['af_data'][i][0]['NumCoords'] != '0' for i in range(len(ad['af_data']))):
            af_coords = [i for i in range(len(ad['af_data'])) if ad['af_data'][i][0]['NumCoords'] != '0']
            self._write_coords(ad, af_coords, run_dir, naming_out)

        # ── Write OLAF ──
        if ad['Wake_Mod'] == 3:
            self._write_olaf(ad, run_dir, naming_out)

        # ── Write main AeroDyn file ──
        f = open(str(file_path), 'w')

        f.write('------- AERODYN15 INPUT FILE ------------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('======  General Options  ============================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['Echo'], 'Echo', '- Echo the input to "<rootname>.AD.ech"?  (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['DTAero'], 'DTAero', '- Time interval for aerodynamic calculations {or "default"} (s)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['Wake_Mod'], 'Wake_Mod', '- Wake/induction model (switch) {0=none, 1=BEMT, 3=OLAF} [Wake_Mod cannot be 2 or 3 when linearizing]\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['TwrPotent'], 'TwrPotent', '- Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['TwrShadow'], 'TwrShadow', '- Calculate tower influence on wind based on downstream tower shadow (switch) {0=none, 1=Powles model, 2=Eames model}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['TwrAero'], 'TwrAero', '- Calculate tower aerodynamic loads? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['CavitCheck'], 'CavitCheck', '- Perform cavitation check? (flag) [UA_Mod must be 0 when CavitCheck=true]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['NacelleDrag'], 'NacelleDrag', '- Include Nacelle Drag effects? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['CompAA'], 'CompAA', '- Flag to compute AeroAcoustics calculation [used only when Wake_Mod = 1 or 2]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['AA_InputFile'], 'AA_InputFile', '- AeroAcoustics input file [used only when CompAA=true]\n'))
        f.write('======  Environmental Conditions  ===================================================================\n')
        f.write('{:<22} {:<11} {:}'.format(ad['AirDens'], 'AirDens', '- Air density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['SpdSound'], 'SpdSound', '- Speed of sound in working fluid (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['Patm'], 'Patm', '- Atmospheric pressure (Pa) [used only when CavitCheck=True]\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['Pvap'], 'Pvap', '- Vapour pressure of working fluid (Pa) [used only when CavitCheck=True]\n'))
        f.write('======  Blade-Element/Momentum Theory Options  ====================================================== [unused when Wake_Mod=0 or 3, except for BEM_Mod]\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['BEM_Mod'], 'BEM_Mod', '- BEM model {1=legacy NoSweepPitchTwist, 2=polar} (switch) [used for all Wake_Mod to determine output coordinate system]\n'))
        f.write('--- Skew correction\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['Skew_Mod'], 'Skew_Mod', '- Skew model {0=No skew model, -1=Remove non-normal component for linearization, 1=skew model active}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SkewMomCorr'], 'SkewMomCorr', '- Turn the skew momentum correction on or off [used only when Skew_Mod=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SkewRedistr_Mod'], 'SkewRedistr_Mod', '- Type of skewed-wake correction model (switch) {0=no redistribution, 1=Glauert/Pitt/Peters, default=1} [used only when Skew_Mod=1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SkewRedistrFactor'], 'SkewRedistrFactor', '- Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when Skew_Mod=1 and SkewRedistr_Mod=1]\n'))
        f.write('--- BEM algorithm\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['TipLoss'], 'TipLoss', '- Use the Prandtl tip-loss model? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['HubLoss'], 'HubLoss', '- Use the Prandtl hub-loss model? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['TanInd'], 'TanInd', '- Include tangential induction in BEMT calculations? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['AIDrag'], 'AIDrag', '- Include the drag term in the axial-induction calculation? (flag) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['TIDrag'], 'TIDrag', '- Include the drag term in the tangential-induction calculation? (flag) [unused when Wake_Mod=0,3 or TanInd=FALSE]\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['IndToler'], 'IndToler', '- Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [unused when Wake_Mod=0 or 3]\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['MaxIter'], 'MaxIter', '- Maximum number of iteration steps (-) [unused when Wake_Mod=0]\n'))
        f.write('--- Shear correction\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['SectAvg'], 'SectAvg', '- Use sector averaging (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SectAvgWeighting'], 'SectAvgWeighting', '- Weighting function for sector average {1=Uniform, default=1} within a sector centered on the blade (switch) [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SectAvgNPoints'], 'SectAvgNPoints', '- Number of points per sectors (-) {default=5} [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SectAvgPsiBwd'], 'SectAvgPsiBwd', '- Backward azimuth relative to blade where the sector starts (<=0) {default=-60} (deg) [used only when SectAvg=True]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['SectAvgPsiFwd'], 'SectAvgPsiFwd', '- Forward azimuth relative to blade where the sector ends (>=0) {default=60} (deg) [used only when SectAvg=True]\n'))
        f.write('--- Dynamic wake/inflow\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['DBEMT_Mod'], 'DBEMT_Mod', '- Type of dynamic BEMT (DBEMT) model {0=No Dynamic Wake, -1=Frozen Wake for linearization, 1:constant tau1, 2=time-dependent tau1, 3=constant tau1 with continuous formulation} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['tau1_const'], 'tau1_const', '- Time constant for DBEMT (s) [used only when DBEMT_Mod=1 or 3]\n'))
        f.write('======  OLAF -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when Wake_Mod=3]\n')
        olaf_file = naming_out + '_OLAF.dat'
        f.write('{!s:<22} {:<11} {:}'.format(olaf_file, 'OLAFInputFileName', '- Input file for OLAF [used only when Wake_Mod=3]\n'))
        f.write('======  Unsteady Airfoil Aerodynamics Options  ===================================== \n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['AoA34'], 'AoA34', "- Sample the angle of attack (AoA) at the 3/4 chord or the AC point {default=True} [always used]\n"))
        f.write('{:<22d} {:<11} {:}'.format(ad['UA_Mod'], 'UA_Mod', "- Unsteady Aero Model Switch (switch) {0=Quasi-steady (no UA), 2=B-L Gonzalez, 3=B-L Minnema/Pierce, 4=B-L HGM 4-states, 5=B-L HGM+vortex 5 states, 6=Oye, 7=Boeing-Vertol}\n"))
        f.write('{!s:<22} {:<11} {:}'.format(ad['FLookup'], 'FLookup', "- Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when UA_Mod>0]\n"))
        f.write('{!s:<22} {:<11} {:}'.format(ad['IntegrationMethod'], 'IntegrationMethod', "- Switch to indicate which integration method UA uses (1=RK4, 2=AB4, 3=ABM4, 4=BDF2)\n"))
        if 'UAStartRad' in ad and 'UAEndRad' in ad:
            f.write('{:<22} {:<11} {:}'.format(ad['UAStartRad'], 'UAStartRad', '- Starting radius for dynamic stall (fraction of rotor radius [0.0,1.0]) [used only when UA_Mod>0; if line is missing UAStartRad=0]\n'))
            f.write('{:<22} {:<11} {:}'.format(ad['UAEndRad'], 'UAEndRad', '- Ending radius for dynamic stall (fraction of rotor radius [0.0,1.0]) [used only when UA_Mod>0; if line is missing UAEndRad=1]\n'))
        f.write('======  Airfoil Information =========================================================================\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['AFTabMod'], 'AFTabMod', '- Interpolation method for multiple airfoil tables {1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp} (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['InCol_Alfa'], 'InCol_Alfa', '- The column in the airfoil tables that contains the angle of attack (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['InCol_Cl'], 'InCol_Cl', '- The column in the airfoil tables that contains the lift coefficient (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['InCol_Cd'], 'InCol_Cd', '- The column in the airfoil tables that contains the drag coefficient (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['InCol_Cm'], 'InCol_Cm', '- The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['InCol_Cpmin'], 'InCol_Cpmin', '- The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['NumAFfiles'], 'NumAFfiles', '- Number of airfoil files used (-)\n'))
        for i in range(ad['NumAFfiles']):
            if i == 0:
                f.write('"' + ad['AFNames'][i] + '"    AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)\n')
            else:
                f.write('"' + ad['AFNames'][i] + '"\n')
        f.write('======  Rotor/Blade Properties  =====================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['UseBlCm'], 'UseBlCm', '- Include aerodynamic pitching moment in calculations?  (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ad['ADBlFile1'] + '"', 'ADBlFile(1)', '- Name of file containing distributed aerodynamic properties for Blade #1 (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ad['ADBlFile2'] + '"', 'ADBlFile(2)', '- Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ad['ADBlFile3'] + '"', 'ADBlFile(3)', '- Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]\n'))
        f.write('======  Hub Properties ============================================================================== [used only when MHK=1 or 2]\n')
        f.write('{:<22} {:<11} {:}'.format(ad['VolHub'], 'VolHub', '- Hub volume (m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['HubCenBx'], 'HubCenBx', '- Hub center of buoyancy x direction offset (m)\n'))
        f.write('======  Nacelle Properties ========================================================================== [used only when MHK=1 or 2 or when NacelleDrag=True]\n')
        f.write('{:<22} {:<11} {:}'.format(ad['VolNac'], 'VolNac', '- Nacelle volume (m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ad['NacCenB'], dtype=str)), 'NacCenB', '- Position of nacelle center of buoyancy from yaw bearing in nacelle coordinates (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ad['NacArea'], dtype=str)), 'NacArea', '- Projected area of the nacelle in X, Y, Z in the nacelle coordinate system (m^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ad['NacCd'], dtype=str)), 'NacCd', '- Drag coefficient for the nacelle areas defined above (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ad['NacDragAC'], dtype=str)), 'NacDragAC', '- Position of aerodynamic center of nacelle drag in nacelle coordinates (m)\n'))
        f.write('======  Tail Fin Aerodynamics ========================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['TFinAero'], 'TFinAero', '- Calculate tail fin aerodynamics model (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ad['TFinFile'] + '"', 'TFinFile', '- Input file for tail fin aerodynamics [used only when TFinAero=True]\n'))
        f.write('======  Tower Influence and Aerodynamics ============================================================ [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or MHK=1 or 2]\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['NumTwrNds'], 'NumTwrNds', '- Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or MHK=1 or 2]\n'))
        f.write('TwrElev        TwrDiam        TwrCd          TwrTI          TwrCb          TwrCp          TwrCa !TwrTI used only with TwrShadow=2, TwrCb/TwrCp/TwrCa used only with MHK=1 or 2\n')
        f.write('(m)            (m)            (-)            (-)            (-)            (-)            (-)\n')
        for TwrElev, TwrDiam, TwrCd, TwrTI, TwrCb, TwrCp, TwrCa in zip(
                ad['TwrElev'], ad['TwrDiam'], ad['TwrCd'], ad['TwrTI'],
                ad['TwrCb'], ad['TwrCp'], ad['TwrCa']):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} \n'.format(
                TwrElev, TwrDiam, TwrCd, TwrTI, TwrCb, TwrCp, TwrCa))
        f.write('======  Outputs  ====================================================================================\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['SumPrint'], 'SumPrint', '- Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['NBlOuts'], 'NBlOuts', '- Number of blade node outputs [0 - 9] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(ad['BlOutNd']), 'BlOutNd', '- Blade nodes whose values will be output  (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['NTwOuts'], 'NTwOuts', '- Number of tower node outputs [0 - 9]  (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(np.array(ad['TwOutNd'], dtype=str)), 'TwOutNd', '- Tower nodes whose values will be output  (-)\n'))
        f.write('                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        if outlist is not None:
            ol = _get_outlist(outlist, ['AeroDyn'])
            for channel_list in ol:
                for ch in channel_list:
                    f.write('"' + ch + '"\n')
        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')

        # Optional nodal output
        if 'BldNd_BladesOut' in ad:
            f.write('====== Outputs for all blade stations (same ending as above for B1N1.... =========================== [optional section]\n')
            f.write('{:<22d} {:<11} {:}'.format(ad['BldNd_BladesOut'], 'BldNd_BladesOut', '- Number of blades to output all node information at.  Up to number of blades on turbine. (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ad['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)\n'))
            f.write('                   OutList_Nodal     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, AeroDyn_Nodes tab for a listing of available output channels, (-)\n')
            if outlist is not None:
                opt_ol = _get_outlist(outlist, ['AeroDyn_Nodes'])
                for opt_channel_list in opt_ol:
                    for ch in opt_channel_list:
                        f.write('"' + ch + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')

        f.write('---------------------------------------------------------------------------------------\n')
        f.flush()
        os.fsync(f)
        f.close()

    # ------------------------------------------------------------------
    # PRIVATE: Blade read/write
    # ------------------------------------------------------------------

    @staticmethod
    def _read_blade(ad_blade_list, blade_file, blade_number):
        bld = ad_blade_list[blade_number]
        f = open(blade_file)
        f.readline(); f.readline(); f.readline()

        bld['NumBlNds'] = int(f.readline().split()[0])
        f.readline(); f.readline()

        for k in ('BlSpn', 't_c', 'BlCrvAC', 'BlSwpAC', 'BlCrvAng', 'BlTwist',
                   'BlChord', 'BlAFID', 'BlCb', 'BlCenBn', 'BlCenBt',
                   'BlCpn', 'BlCpt', 'BlCan', 'BlCat', 'BlCam'):
            bld[k] = [None] * bld['NumBlNds']

        for i in range(bld['NumBlNds']):
            data = [float(val) for val in f.readline().split()]
            bld['BlSpn'][i]    = data[0]
            bld['BlCrvAC'][i]  = data[1]
            bld['BlSwpAC'][i]  = data[2]
            bld['BlCrvAng'][i] = data[3]
            bld['BlTwist'][i]  = data[4]
            bld['BlChord'][i]  = data[5]
            bld['BlAFID'][i]   = data[6]
            if len(data) == 16:
                bld['t_c'][i]     = data[7]
                bld['BlCb'][i]    = data[8]
                bld['BlCenBn'][i] = data[9]
                bld['BlCenBt'][i] = data[10]
                bld['BlCpn'][i]   = data[11]
                bld['BlCpt'][i]   = data[12]
                bld['BlCan'][i]   = data[13]
                bld['BlCat'][i]   = data[14]
                bld['BlCam'][i]   = data[15]
            else:
                bld['t_c'][i] = 0.0
                bld['BlCb'][i] = 0.0
                bld['BlCenBn'][i] = 0.0
                bld['BlCenBt'][i] = 0.0
                bld['BlCpn'][i] = 0.0
                bld['BlCpt'][i] = 0.0
                bld['BlCan'][i] = 0.0
                bld['BlCat'][i] = 0.0
                bld['BlCam'][i] = 0.0

        f.close()

    @staticmethod
    def _write_blade(ad, ad_blade, run_dir, bld_ind=None):
        if bld_ind is None:
            filename = os.path.join(run_dir, ad['ADBlFile1'])
            bld_dict = ad_blade
        else:
            filename = os.path.join(run_dir, ad['ADBlFile%d' % (bld_ind + 1)])
            bld_dict = ad_blade[bld_ind]

        f = open(filename, 'w')
        f.write('------- AERODYN15 BLADE DEFINITION INPUT FILE -------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('======  Blade Properties =================================================================\n')
        f.write('{:<11d} {:<11} {:}'.format(bld_dict['NumBlNds'], 'NumBlNds', '- Number of blade nodes used in the analysis (-)\n'))
        f.write('    BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID       t_c       BlCb        BlCenBn      BlCenBt       BlCpn     BlCpt     BlCan     BlCat     BlCam\n')
        f.write('     (m)           (m)            (m)            (deg)         (deg)           (m)              (-)        (-)       (-)         (m)             (m)         (-)       (-)       (-)       (-)       (-)\n')

        for Spn, CrvAC, SwpAC, CrvAng, Twist, Chord, AFID, tc, Cb, CenBn, CenBt, Cpn, Cpt, Can, Cat, Cam in zip(
                bld_dict['BlSpn'], bld_dict['BlCrvAC'], bld_dict['BlSwpAC'],
                bld_dict['BlCrvAng'], bld_dict['BlTwist'], bld_dict['BlChord'],
                bld_dict['BlAFID'], bld_dict['t_c'], bld_dict['BlCb'],
                bld_dict['BlCenBn'], bld_dict['BlCenBt'],
                bld_dict['BlCpn'], bld_dict['BlCpt'],
                bld_dict['BlCan'], bld_dict['BlCat'], bld_dict['BlCam']):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 8d} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(
                Spn, CrvAC, SwpAC, CrvAng, Twist, Chord, int(AFID),
                tc, Cb, CenBn, CenBt, Cpn, Cpt, Can, Cat, Cam))

        f.flush()
        os.fsync(f)
        f.close()

    # ------------------------------------------------------------------
    # PRIVATE: Polar read/write
    # ------------------------------------------------------------------

    @staticmethod
    def _read_polars(ad):
        ad['af_data'] = [None] * ad['NumAFfiles']

        for afi, af_filename in enumerate(ad['AFNames']):
            f = open(af_filename)
            polar = {}

            polar['InterpOrd']  = int_read(readline_filterComments(f).split()[0])
            temp = readline_filterComments(f).split()
            if temp[1] == "RelThickness":
                polar['RelThickness'] = float_read(temp[0])
                polar['NonDimArea'] = float_read(readline_filterComments(f).split()[0])
            else:
                polar['NonDimArea'] = float_read(temp[0])
            polar['NumCoords']  = readline_filterComments(f).split()[0]
            polar['BL_file']    = readline_filterComments(f).split()[0]
            polar['NumTabs']    = int_read(readline_filterComments(f).split()[0])
            ad['af_data'][afi]  = [None] * polar['NumTabs']

            for tab in range(polar['NumTabs']):
                polar['Re']         = float_read(readline_filterComments(f).split()[0]) * 1.e+6
                polar['UserProp']   = int_read(readline_filterComments(f).split()[0])
                polar['InclUAdata'] = bool_read(readline_filterComments(f).split()[0])

                if polar['InclUAdata']:
                    polar['alpha0']     = float_read(readline_filterComments(f).split()[0])
                    polar['alpha1']     = float_read(readline_filterComments(f).split()[0])
                    polar['alpha2']     = float_read(readline_filterComments(f).split()[0])
                    polar['eta_e']      = float_read(readline_filterComments(f).split()[0])
                    polar['C_nalpha']   = float_read(readline_filterComments(f).split()[0])
                    polar['T_f0']       = float_read(readline_filterComments(f).split()[0])
                    polar['T_V0']       = float_read(readline_filterComments(f).split()[0])
                    polar['T_p']        = float_read(readline_filterComments(f).split()[0])
                    polar['T_VL']       = float_read(readline_filterComments(f).split()[0])
                    polar['b1']         = float_read(readline_filterComments(f).split()[0])
                    polar['b2']         = float_read(readline_filterComments(f).split()[0])
                    polar['b5']         = float_read(readline_filterComments(f).split()[0])
                    polar['A1']         = float_read(readline_filterComments(f).split()[0])
                    polar['A2']         = float_read(readline_filterComments(f).split()[0])
                    polar['A5']         = float_read(readline_filterComments(f).split()[0])
                    polar['S1']         = float_read(readline_filterComments(f).split()[0])
                    polar['S2']         = float_read(readline_filterComments(f).split()[0])
                    polar['S3']         = float_read(readline_filterComments(f).split()[0])
                    polar['S4']         = float_read(readline_filterComments(f).split()[0])
                    polar['Cn1']        = float_read(readline_filterComments(f).split()[0])
                    polar['Cn2']        = float_read(readline_filterComments(f).split()[0])
                    polar['St_sh']      = float_read(readline_filterComments(f).split()[0])
                    polar['Cd0']        = float_read(readline_filterComments(f).split()[0])
                    polar['Cm0']        = float_read(readline_filterComments(f).split()[0])
                    polar['k0']         = float_read(readline_filterComments(f).split()[0])
                    polar['k1']         = float_read(readline_filterComments(f).split()[0])
                    polar['k2']         = float_read(readline_filterComments(f).split()[0])
                    polar['k3']         = float_read(readline_filterComments(f).split()[0])
                    polar['k1_hat']     = float_read(readline_filterComments(f).split()[0])
                    polar['x_cp_bar']   = float_read(readline_filterComments(f).split()[0])
                    polar['UACutout']   = float_read(readline_filterComments(f).split()[0])
                    polar['filtCutOff'] = float_read(readline_filterComments(f).split()[0])

                polar['NumAlf'] = int_read(readline_filterComments(f).split()[0])
                polar['Alpha']  = [None] * polar['NumAlf']
                polar['Cl']     = [None] * polar['NumAlf']
                polar['Cd']     = [None] * polar['NumAlf']
                polar['Cm']     = [None] * polar['NumAlf']
                polar['Cpmin']  = [None] * polar['NumAlf']
                for i in range(polar['NumAlf']):
                    data = [float(val) for val in readline_filterComments(f).split()]
                    if ad['InCol_Alfa'] > 0:
                        polar['Alpha'][i] = data[ad['InCol_Alfa'] - 1]
                    if ad['InCol_Cl'] > 0:
                        polar['Cl'][i]    = data[ad['InCol_Cl'] - 1]
                    if ad['InCol_Cd'] > 0:
                        polar['Cd'][i]    = data[ad['InCol_Cd'] - 1]
                    if ad['InCol_Cm'] > 0:
                        polar['Cm'][i]    = data[ad['InCol_Cm'] - 1]
                    if ad['InCol_Cpmin'] > 0:
                        polar['Cpmin'][i] = data[ad['InCol_Cpmin'] - 1]

                ad['af_data'][afi][tab] = copy.copy(polar)

            f.close()

    @staticmethod
    def _write_polars(ad, run_dir, naming_out):
        airfoils_dir = os.path.join(run_dir, 'Airfoils')
        if not os.path.isdir(airfoils_dir):
            try:
                os.makedirs(airfoils_dir)
            except Exception:
                try:
                    time.sleep(random.random())
                    if not os.path.isdir(airfoils_dir):
                        os.makedirs(airfoils_dir)
                except Exception:
                    print("Error trying to make '%s'!" % airfoils_dir)

        ad['NumAFfiles'] = len(ad['af_data'])
        ad['AFNames'] = [''] * ad['NumAFfiles']

        for afi in range(int(ad['NumAFfiles'])):
            ad['AFNames'][afi] = os.path.join('Airfoils', naming_out + '_AeroDyn_Polar_%02d.dat' % afi)
            af_file = os.path.join(run_dir, ad['AFNames'][afi])
            f = open(af_file, 'w')

            f.write('! ------------ AirfoilInfo Input File ----------------------------------\n')
            f.write('! Generated with OpenFAST_IO\n')
            f.write('! line\n')
            f.write('! line\n')
            f.write('! ------------------------------------------------------------------------------\n')
            f.write('{:<22}   {:<11} {:}'.format(ad['af_data'][afi][0]['InterpOrd'], 'InterpOrd', '! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]\n'))
            if 'RelThickness' in ad['af_data'][afi][0]:
                f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][0]['RelThickness'], 'RelThickness', '! The non-dimensional thickness of the airfoil (thickness/chord) [only used if UAMod=7] [default=0.2] (-)\n'))
            f.write('{:<22}   {:<11} {:}'.format(ad['af_data'][afi][0]['NonDimArea'], 'NonDimArea', '! The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)\n'))
            if ad['af_data'][afi][0]['NumCoords'] != '0':
                f.write('@"{:}_AF{:02d}_Coords.txt"       {:<11} {:}'.format(naming_out, afi, 'NumCoords', '! The number of coordinates in the airfoil shape file. Set to zero if coordinates not included.\n'))
            else:
                f.write('{:<22d}       {:<11} {:}'.format(0, 'NumCoords', '! The number of coordinates in the airfoil shape file. Set to zero if coordinates not included.\n'))
            f.write('AF{:02d}_BL.txt              {:<11} {:}'.format(afi, 'BL_file', '! The file name including the boundary layer characteristics of the profile. Ignored if the aeroacoustic module is not called.\n'))

            # Determine number of tabs to write
            if ad['AFTabMod'] == 2:
                num_tab = len(ad['af_data'][afi])
            elif ad['AFTabMod'] == 3:
                if len(ad['af_data'][afi]) == 1 or \
                   ad['af_data'][afi][0]['UserProp'] == ad['af_data'][afi][1]['UserProp']:
                    num_tab = 1
                else:
                    num_tab = ad['af_data'][afi][0]['NumTabs']
            else:
                num_tab = 1

            f.write('{:<22d}   {:<11} {:}'.format(num_tab, 'NumTabs', '! Number of airfoil tables in this file.  Each table must have lines for Re and UserProp.\n'))

            for tab in range(num_tab):
                f.write('! ------------------------------------------------------------------------------\n')
                f.write("! data for table %i \n" % (tab + 1))
                f.write('! ------------------------------------------------------------------------------\n')
                f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['Re'] * 1.e-6, 'Re', '! Reynolds number in millions\n'))
                f.write('{:<22d}   {:<11} {:}'.format(int(ad['af_data'][afi][tab]['UserProp']), 'UserProp', '! User property (control) setting\n'))
                f.write('{!s:<22}   {:<11} {:}'.format(ad['af_data'][afi][tab]['InclUAdata'], 'InclUAdata', '! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n'))
                f.write('!........................................\n')
                if ad['af_data'][afi][tab]['InclUAdata']:
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['alpha0'], 'alpha0', '! 0-lift angle of attack, depends on airfoil.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['alpha1'], 'alpha1', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['alpha2'], 'alpha2', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['eta_e'], 'eta_e', '! Recovery factor in the range [0.85 - 0.95] used only for UA_Mod=1, it is set to 1 in the code when flookup=True. (-)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['C_nalpha'], 'C_nalpha', '! Slope of the 2D normal force coefficient curve. (1/rad)\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['T_f0']) + '   {:<11} {:}'.format('T_f0', '! Initial value of the time constant associated with Df in the expression of Df and f\'\'. [default = 3]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['T_V0']) + '   {:<11} {:}'.format('T_V0', '! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['T_p']) + '   {:<11} {:}'.format('T_p', '! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['T_VL']) + '   {:<11} {:}'.format('T_VL', '! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['b1']) + '   {:<11} {:}'.format('b1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['b2']) + '   {:<11} {:}'.format('b2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['b5']) + '   {:<11} {:}'.format('b5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n"))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['A1']) + '   {:<11} {:}'.format('A1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['A2']) + '   {:<11} {:}'.format('A2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['A5']) + '   {:<11} {:}'.format('A5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['S1'], 'S1', '! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['S2'], 'S2', '! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['S3'], 'S3', '! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['S4'], 'S4', '! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['Cn1'], 'Cn1', '! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['Cn2'], 'Cn2', '! As Cn1 for negative AOAs.\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['St_sh']) + '   {:<11} {:}'.format('St_sh', "! Strouhal's shedding frequency constant.  [default = 0.19]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['Cd0'], 'Cd0', '! 2D drag coefficient value at 0-lift.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['Cm0'], 'Cm0', '! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['k0'], 'k0', '! Constant in the \\hat(x)_cp curve best-fit; = (\\hat(x)_AC-0.25).  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['k1'], 'k1', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['k2'], 'k2', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['k3'], 'k3', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UA_Mod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(ad['af_data'][afi][tab]['k1_hat'], 'k1_hat', '! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UA_Mod<>1]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['x_cp_bar']) + '   {:<11} {:}'.format('x_cp_bar', '! Constant in the expression of \\hat(x)_cp^v. [ignored if UA_Mod<>1, default = 0.2]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['UACutout']) + '   {:<11} {:}'.format('UACutout', '! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n'))
                    f.write(_float_default_out(ad['af_data'][afi][tab]['filtCutOff']) + '   {:<11} {:}'.format('filtCutOff', '! Reduced frequency cut-off for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (-) [default = 0.5]\n'))

                f.write('!........................................\n')
                f.write('! Table of aerodynamics coefficients\n')
                f.write('{:<22d}   {:<11} {:}'.format(ad['af_data'][afi][tab]['NumAlf'], 'NumAlf', '! Number of data lines in the following table\n'))
                f.write('!    Alpha      Cl      Cd        Cm\n')
                f.write('!    (deg)      (-)     (-)       (-)\n')

                polar_map = [ad['InCol_Alfa'], ad['InCol_Cl'], ad['InCol_Cd'], ad['InCol_Cm'], ad['InCol_Cpmin']]
                polar_map.remove(0)
                polar_map = [i - 1 for i in polar_map]

                alpha = np.asarray(ad['af_data'][afi][tab]['Alpha'])
                cl = np.asarray(ad['af_data'][afi][tab]['Cl'])
                cd = np.asarray(ad['af_data'][afi][tab]['Cd'])
                cm = np.asarray(ad['af_data'][afi][tab]['Cm'])
                cpmin = np.asarray(ad['af_data'][afi][tab]['Cpmin'])

                if alpha[0] != -180.:
                    alpha[0] = -180.
                if alpha[-1] != 180.:
                    alpha[-1] = 180.
                if cl[0] != cl[-1]:
                    cl[0] = cl[-1]
                if cd[0] != cd[-1]:
                    cd[0] = cd[-1]
                if cm[0] != cm[-1]:
                    cm[0] = cm[-1]

                if ad['InCol_Cm'] == 0:
                    cm = np.zeros_like(cl)
                if ad['InCol_Cpmin'] == 0:
                    cpmin = np.zeros_like(cl)
                polar = np.column_stack((alpha, cl, cd, cm, cpmin))
                polar = polar[:, polar_map]

                for row in polar:
                    f.write(' '.join(['{: 2.14e}'.format(val) for val in row]) + '\n')

            f.flush()
            os.fsync(f)
            f.close()

    # ------------------------------------------------------------------
    # PRIVATE: Coord read/write
    # ------------------------------------------------------------------

    @staticmethod
    def _read_coords(ad):
        ad['af_coord'] = []
        ad['ac'] = np.zeros(len(ad['AFNames']))

        for afi, af_filename in enumerate(ad['AFNames']):
            ad['af_coord'].append({})
            if not (ad['af_data'][afi][0]['NumCoords'] == 0 or ad['af_data'][afi][0]['NumCoords'] == '0'):
                coord_filename = af_filename[:af_filename.rfind(os.sep)] + os.sep + ad['af_data'][afi][0]['NumCoords'][2:-1]

                f = open(coord_filename)
                lines = f.readlines()
                f.close()
                lines = [line for line in lines if not line.strip().startswith('!')]
                n_coords = int(lines[0].split()[0])

                x = np.zeros(n_coords - 1)
                y = np.zeros(n_coords - 1)
                ad['ac'][afi] = float(lines[1].split()[0])

                for j in range(2, n_coords + 1):
                    x[j - 2], y[j - 2] = map(float, lines[j].split())

                ad['af_coord'][afi]['x'] = x
                ad['af_coord'][afi]['y'] = y

    @staticmethod
    def _write_coords(ad, af_coord_indices, run_dir, naming_out):
        coord_names = [''] * ad['NumAFfiles']

        for afi in af_coord_indices:
            coord_names[afi] = os.path.join('Airfoils', naming_out + '_AF%02d_Coords.txt' % afi)

            x = ad['af_coord'][afi]['x']
            y = ad['af_coord'][afi]['y']
            coord = np.vstack((x, y)).T

            af_file = os.path.join(run_dir, coord_names[afi])
            f = open(af_file, 'w')

            f.write('{: 22d}   {:<11} {:}'.format(len(x) + 1, 'NumCoords', '! The number of coordinates in the airfoil shape file (including an extra coordinate for airfoil reference).  Set to zero if coordinates not included.\n'))
            f.write('! ......... x-y coordinates are next if NumCoords > 0 .............\n')
            f.write('! x-y coordinate of airfoil reference\n')
            f.write('!  x/c        y/c\n')
            f.write('{: 5f}       0\n'.format(ad['ac'][afi]))
            f.write('! coordinates of airfoil shape\n')
            f.write('! interpolation to 200 points\n')
            f.write('!  x/c        y/c\n')
            for row in coord:
                f.write(' '.join(['{: 2.14e}'.format(val) for val in row]) + '\n')

            f.flush()
            os.fsync(f)
            f.close()

    # ------------------------------------------------------------------
    # PRIVATE: OLAF read/write
    # ------------------------------------------------------------------

    @staticmethod
    def _read_olaf(ad, olaf_filename, base_dir):
        ad['OLAF'] = {}
        f = open(olaf_filename)
        f.readline(); f.readline(); f.readline()
        ad['OLAF']['IntMethod']       = int_read(f.readline().split()[0])
        ad['OLAF']['DTfvw']           = float_read(f.readline().split()[0])
        ad['OLAF']['FreeWakeStart']   = float_read(f.readline().split()[0])
        ad['OLAF']['FullCircStart']   = float_read(f.readline().split()[0])
        f.readline()
        ad['OLAF']['CircSolvMethod']       = int_read(f.readline().split()[0])
        ad['OLAF']['CircSolvConvCrit']     = float_read(f.readline().split()[0])
        ad['OLAF']['CircSolvRelaxation']   = float_read(f.readline().split()[0])
        ad['OLAF']['CircSolvMaxIter']      = int_read(f.readline().split()[0])
        ad['OLAF']['PrescribedCircFile']   = os.path.join(str(base_dir), quoted_read(f.readline().split()[0]))
        f.readline(); f.readline(); f.readline()
        ad['OLAF']['nNWPanels']       = int_read(f.readline().split()[0])
        ad['OLAF']['nNWPanelsFree']   = int_read(f.readline().split()[0])
        ad['OLAF']['nFWPanels']       = int_read(f.readline().split()[0])
        ad['OLAF']['nFWPanelsFree']   = int_read(f.readline().split()[0])
        ad['OLAF']['FWShedVorticity'] = bool_read(f.readline().split()[0])
        f.readline()
        ad['OLAF']['DiffusionMethod'] = int_read(f.readline().split()[0])
        ad['OLAF']['RegDeterMethod']  = int_read(f.readline().split()[0])
        ad['OLAF']['RegFunction']     = int_read(f.readline().split()[0])
        ad['OLAF']['WakeRegMethod']   = int_read(f.readline().split()[0])
        ad['OLAF']['WakeRegFactor']   = float(f.readline().split()[0])
        ad['OLAF']['WingRegFactor']   = float(f.readline().split()[0])
        ad['OLAF']['CoreSpreadEddyVisc'] = int(f.readline().split()[0])
        f.readline()
        ad['OLAF']['TwrShadowOnWake'] = bool_read(f.readline().split()[0])
        ad['OLAF']['ShearModel']      = int_read(f.readline().split()[0])
        f.readline()
        ad['OLAF']['VelocityMethod']  = int_read(f.readline().split()[0])
        ad['OLAF']['TreeBranchFactor'] = float_read(f.readline().split()[0])
        ad['OLAF']['PartPerSegment']  = int_read(f.readline().split()[0])
        f.readline(); f.readline()
        ad['OLAF']['WrVTk']       = int_read(f.readline().split()[0])
        ad['OLAF']['nVTKBlades']  = int_read(f.readline().split()[0])
        ad['OLAF']['VTKCoord']    = int_read(f.readline().split()[0])
        ad['OLAF']['VTK_fps']     = float_read(f.readline().split()[0])
        ad['OLAF']['nGridOut']    = int_read(f.readline().split()[0])
        f.readline()
        f.close()

    @staticmethod
    def _write_olaf(ad, run_dir, naming_out):
        olaf_file = os.path.join(run_dir, naming_out + '_OLAF.dat')
        f = open(olaf_file, 'w')

        f.write('--------------------------- OLAF (cOnvecting LAgrangian Filaments) INPUT FILE -----------------\n')
        f.write('Generated by OpenFAST_IO\n')
        f.write('--------------------------- GENERAL OPTIONS ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['IntMethod'], 'IntMethod', '- Integration method {1: RK4, 5: Forward Euler 1st order, default: 5} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['DTfvw'], 'DTfvw', '- Time interval for wake propagation. {default: dtaero} (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['FreeWakeStart'], 'FreeWakeStart', '- Time when wake is free. (-) value = always free. {default: 0.0} (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['FullCircStart'], 'FullCircStart', '- Time at which full circulation is reached. {default: 0.0} (s)\n'))
        f.write('--------------------------- CIRCULATION SPECIFICATIONS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['CircSolvMethod'], 'CircSolvingMethod', '- Circulation solving method {1: Cl-Based, 2: No-Flow Through, 3: Prescribed, default: 1 }(switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['CircSolvConvCrit'], 'CircSolvConvCrit', ' - Convergence criteria {default: 0.001} [only if CircSolvMethod=1] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['CircSolvRelaxation'], 'CircSolvRelaxation', '- Relaxation factor {default: 0.1} [only if CircSolvMethod=1] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['CircSolvMaxIter'], 'CircSolvMaxIter', ' - Maximum number of iterations for circulation solving {default: 30} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + ad['OLAF']['PrescribedCircFile'] + '"', 'PrescribedCircFile', '- File containing prescribed circulation [only if CircSolvMethod=3] (quoted string)\n'))
        f.write('===============================================================================================\n')
        f.write('--------------------------- WAKE OPTIONS ------------------------------------------------------\n')
        f.write('------------------- WAKE EXTENT AND DISCRETIZATION --------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['OLAF']['nNWPanels'], 'nNWPanels', '- Number of near-wake panels (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['nNWPanelsFree'], 'nNWPanelsFree', '- Number of free near-wake panels (-) {default: nNWPanels}\n'))
        f.write('{:<22d} {:<11} {:}'.format(ad['OLAF']['nFWPanels'], 'nFWPanels', '- Number of far-wake panels (-) {default: 0}\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['nFWPanelsFree'], 'nFWPanelsFree', '- Number of free far-wake panels (-) {default: nFWPanels}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ad['OLAF']['FWShedVorticity'], 'FWShedVorticity', '- Include shed vorticity in the far wake {default: False}\n'))
        f.write('------------------- WAKE REGULARIZATIONS AND DIFFUSION -----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['DiffusionMethod'], 'DiffusionMethod', '- Diffusion method to account for viscous effects {0: None, 1: Core Spreading, "default": 0}\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['RegDeterMethod'], 'RegDeterMethod', '- Method to determine the regularization parameters {0:  Manual, 1: Optimized, 2: Chord, 3: Span, default: 0 }\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['RegFunction'], 'RegFunction', '- Viscous diffusion function {0: None, 1: Rankine, 2: LambOseen, 3: Vatistas, 4: Denominator, "default": 3} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['WakeRegMethod'], 'WakeRegMethod', '- Wake regularization method {1: Constant, 2: Stretching, 3: Age, default: 3} (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['WakeRegFactor'], 'WakeRegFactor', '- Wake regularization factor (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['WingRegFactor'], 'WingRegFactor', '- Wing regularization factor (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['CoreSpreadEddyVisc'], 'CoreSpreadEddyVisc', '- Eddy viscosity in core spreading methods, typical values 1-1000\n'))
        f.write('------------------- WAKE TREATMENT OPTIONS ---------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(ad['OLAF']['TwrShadowOnWake'], 'TwrShadowOnWake', '- Include tower flow disturbance effects on wake convection {default:false} [only if TwrPotent or TwrShadow]\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['ShearModel'], 'ShearModel', '- Shear Model {0: No treatment, 1: Mirrored vorticity, default: 0}\n'))
        f.write('------------------- SPEEDUP OPTIONS -----------------------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(ad['OLAF']['VelocityMethod'], 'VelocityMethod', '- Method to determine the velocity {1:Segment N^2, 2:Particle tree, 3:Particle N^2, 4:Segment tree, default: 2}\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['TreeBranchFactor'], 'TreeBranchFactor', '- Branch radius fraction above which a multipole calculation is used {default: 1.5} [only if VelocityMethod=2,4]\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['PartPerSegment'], 'PartPerSegment', '- Number of particles per segment {default: 1} [only if VelocityMethod=2,3]\n'))
        f.write('===============================================================================================\n')
        f.write('--------------------------- OUTPUT OPTIONS  ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['WrVTk'], 'WrVTk', '- Outputs Visualization Toolkit (VTK) (independent of .fst option) {0: NoVTK, 1: Write VTK at VTK_fps, 2: Write VTK at init and final, default: 0} (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['nVTKBlades'], 'nVTKBlades', '- Number of blades for which VTK files are exported {0: No VTK per blade, n: VTK for blade 1 to n, default: 0} (-) \n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['VTKCoord'], 'VTKCoord', '- Coordinate system used for VTK export. {1: Global, 2: Hub, 3: Both, default: 1} \n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['VTK_fps'], 'VTK_fps', '- Frame rate for VTK output (frames per second) {"all" for all glue code timesteps, "default" for all OLAF timesteps} [only if WrVTK=1]\n'))
        f.write('{:<22} {:<11} {:}'.format(ad['OLAF']['nGridOut'], 'nGridOut', '- Number of grid outputs\n'))
        f.write('GridName  GridType  TStart   TEnd    DTGrid   XStart  XEnd   nX   YStart   YEnd   nY   ZStart   ZEnd   nZ\n')
        f.write('(-)         (-)      (s)     (s)      (s)       (m)    (m)   (-)   (m)     (m)    (-)   (m)     (m)    (-)\n')
        f.write('===============================================================================================\n')
        f.write('--------------------------- ADVANCED OPTIONS --------------------------------------------------\n')
        f.write('===============================================================================================\n')

        f.flush()
        os.fsync(f)
        f.close()
