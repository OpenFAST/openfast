"""ElastoDyn module IO — reads and writes ElastoDyn, blade, and tower input files.

Extracted from FAST_reader.py and FAST_writer.py. All parsing logic is identical
to the original; only the data target changes (local dict vs self.fst_vt).
"""
import os
from pathlib import Path

from .base import ModuleIO
from ..outlist import capture_outlist, emit_outlist
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
    fix_path,
)


class ElastoDynIO(ModuleIO):
    """Reads and writes ElastoDyn input files.

    read() returns::

        {
            'ElastoDyn': { ... },
            'ElastoDynBlade': [ {blade0}, {blade1}, {blade2} ],
            'ElastoDynTower': { ... },
        }

    write() accepts data with key 'ElastoDyn' (dict), optional 'ElastoDynBlade'
    and 'ElastoDynTower'. Writes the main ED file plus blade/tower sub-files.
    """

    # ------------------------------------------------------------------
    # READ
    # ------------------------------------------------------------------

    def read(self, file_path: Path, base_dir: Path, outlist: dict = None) -> dict:
        ed = {}
        file_path = Path(file_path)
        base_dir = Path(base_dir)

        f = open(file_path)

        f.readline()
        f.readline()

        # Simulation Control
        f.readline()
        ed['Echo'] = bool_read(f.readline().split()[0])
        ed['Method'] = int(f.readline().split()[0])
        ed['DT'] = float_read(f.readline().split()[0])

        # Degrees of Freedom
        f.readline()
        ed['FlapDOF1'] = bool_read(f.readline().split()[0])
        ed['FlapDOF2'] = bool_read(f.readline().split()[0])
        ed['EdgeDOF'] = bool_read(f.readline().split()[0])
        ed['PitchDOF'] = bool_read(f.readline().split()[0])
        ed['TeetDOF'] = bool_read(f.readline().split()[0])
        ed['DrTrDOF'] = bool_read(f.readline().split()[0])
        ed['GenDOF'] = bool_read(f.readline().split()[0])
        ed['YawDOF'] = bool_read(f.readline().split()[0])
        ed['TwFADOF1'] = bool_read(f.readline().split()[0])
        ed['TwFADOF2'] = bool_read(f.readline().split()[0])
        ed['TwSSDOF1'] = bool_read(f.readline().split()[0])
        ed['TwSSDOF2'] = bool_read(f.readline().split()[0])
        ed['PtfmSgDOF'] = bool_read(f.readline().split()[0])
        ed['PtfmSwDOF'] = bool_read(f.readline().split()[0])
        ed['PtfmHvDOF'] = bool_read(f.readline().split()[0])
        ed['PtfmRDOF'] = bool_read(f.readline().split()[0])
        ed['PtfmPDOF'] = bool_read(f.readline().split()[0])
        ed['PtfmYDOF'] = bool_read(f.readline().split()[0])

        # Initial Conditions
        f.readline()
        ed['OoPDefl'] = float_read(f.readline().split()[0])
        ed['IPDefl'] = float_read(f.readline().split()[0])
        ed['BlPitch1'] = float_read(f.readline().split()[0])
        ed['BlPitch2'] = float_read(f.readline().split()[0])
        ed['BlPitch3'] = float_read(f.readline().split()[0])
        ed['TeetDefl'] = float_read(f.readline().split()[0])
        ed['Azimuth'] = float_read(f.readline().split()[0])
        ed['RotSpeed'] = float_read(f.readline().split()[0])
        ed['NacYaw'] = float_read(f.readline().split()[0])
        ed['TTDspFA'] = float_read(f.readline().split()[0])
        ed['TTDspSS'] = float_read(f.readline().split()[0])
        ed['PtfmSurge'] = float_read(f.readline().split()[0])
        ed['PtfmSway'] = float_read(f.readline().split()[0])
        ed['PtfmHeave'] = float_read(f.readline().split()[0])
        ed['PtfmRoll'] = float_read(f.readline().split()[0])
        ed['PtfmPitch'] = float_read(f.readline().split()[0])
        ed['PtfmYaw'] = float_read(f.readline().split()[0])

        # Turbine Configuration
        f.readline()
        ed['NumBl'] = int(f.readline().split()[0])
        ed['TipRad'] = float_read(f.readline().split()[0])
        ed['HubRad'] = float_read(f.readline().split()[0])
        ed['PreCone(1)'] = float_read(f.readline().split()[0])
        ed['PreCone(2)'] = float_read(f.readline().split()[0])
        ed['PreCone(3)'] = float_read(f.readline().split()[0])
        ed['HubCM'] = float_read(f.readline().split()[0])
        ed['UndSling'] = float_read(f.readline().split()[0])
        ed['Delta3'] = float_read(f.readline().split()[0])
        ed['AzimB1Up'] = float_read(f.readline().split()[0])
        ed['OverHang'] = float_read(f.readline().split()[0])
        ed['ShftGagL'] = float_read(f.readline().split()[0])
        ed['ShftTilt'] = float_read(f.readline().split()[0])
        ed['NacCMxn'] = float_read(f.readline().split()[0])
        ed['NacCMyn'] = float_read(f.readline().split()[0])
        ed['NacCMzn'] = float_read(f.readline().split()[0])
        ed['NcIMUxn'] = float_read(f.readline().split()[0])
        ed['NcIMUyn'] = float_read(f.readline().split()[0])
        ed['NcIMUzn'] = float_read(f.readline().split()[0])
        ed['Twr2Shft'] = float_read(f.readline().split()[0])
        ed['TowerHt'] = float_read(f.readline().split()[0])
        ed['TowerBsHt'] = float_read(f.readline().split()[0])
        ed['PtfmCMxt'] = float_read(f.readline().split()[0])
        ed['PtfmCMyt'] = float_read(f.readline().split()[0])
        ed['PtfmCMzt'] = float_read(f.readline().split()[0])
        ed['PtfmRefxt'] = float_read(f.readline().split()[0])
        ed['PtfmRefyt'] = float_read(f.readline().split()[0])
        ed['PtfmRefzt'] = float_read(f.readline().split()[0])

        # Mass and Inertia
        f.readline()
        ed['TipMass(1)'] = float_read(f.readline().split()[0])
        ed['TipMass(2)'] = float_read(f.readline().split()[0])
        ed['TipMass(3)'] = float_read(f.readline().split()[0])
        ed['PBrIner(1)'] = float_read(f.readline().split()[0])
        ed['PBrIner(2)'] = float_read(f.readline().split()[0])
        ed['PBrIner(3)'] = float_read(f.readline().split()[0])
        ed['BlPIner(1)'] = float_read(f.readline().split()[0])
        ed['BlPIner(2)'] = float_read(f.readline().split()[0])
        ed['BlPIner(3)'] = float_read(f.readline().split()[0])
        ed['HubMass'] = float_read(f.readline().split()[0])
        ed['HubIner'] = float_read(f.readline().split()[0])
        ed['HubIner_Teeter'] = float_read(f.readline().split()[0])
        ed['GenIner'] = float_read(f.readline().split()[0])
        ed['NacMass'] = float_read(f.readline().split()[0])
        ed['NacYIner'] = float_read(f.readline().split()[0])
        ed['YawBrMass'] = float_read(f.readline().split()[0])
        ed['PtfmMass'] = float_read(f.readline().split()[0])
        ed['PtfmRIner'] = float_read(f.readline().split()[0])
        ed['PtfmPIner'] = float_read(f.readline().split()[0])
        ed['PtfmYIner'] = float_read(f.readline().split()[0])
        ed['PtfmXYIner'] = float_read(f.readline().split()[0])
        ed['PtfmYZIner'] = float_read(f.readline().split()[0])
        ed['PtfmXZIner'] = float_read(f.readline().split()[0])

        # Blade
        f.readline()
        ed['BldNodes'] = int(f.readline().split()[0])
        ed['BldFile1'] = quoted_read(f.readline().split()[0])
        ed['BldFile2'] = quoted_read(f.readline().split()[0])
        ed['BldFile3'] = quoted_read(f.readline().split()[0])

        # Rotor-Teeter
        f.readline()
        ed['TeetMod'] = int(f.readline().split()[0])
        ed['TeetDmpP'] = float_read(f.readline().split()[0])
        ed['TeetDmp'] = float_read(f.readline().split()[0])
        ed['TeetCDmp'] = float_read(f.readline().split()[0])
        ed['TeetSStP'] = float_read(f.readline().split()[0])
        ed['TeetHStP'] = float_read(f.readline().split()[0])
        ed['TeetSSSp'] = float_read(f.readline().split()[0])
        ed['TeetHSSp'] = float_read(f.readline().split()[0])

        # Yaw Friction
        f.readline()
        ed['YawFrctMod'] = int(f.readline().split()[0])
        ed['M_CSmax'] = float_read(f.readline().split()[0])
        ed['M_FCSmax'] = float_read(f.readline().split()[0])
        ed['M_MCSmax'] = float_read(f.readline().split()[0])
        ed['M_CD'] = float_read(f.readline().split()[0])
        ed['M_FCD'] = float_read(f.readline().split()[0])
        ed['M_MCD'] = float_read(f.readline().split()[0])
        ed['sig_v'] = float_read(f.readline().split()[0])
        ed['sig_v2'] = float_read(f.readline().split()[0])
        ed['OmgCut'] = float_read(f.readline().split()[0])

        # Drivetrain
        f.readline()
        ed['GBoxEff'] = float_read(f.readline().split()[0])
        ed['GBRatio'] = float_read(f.readline().split()[0])
        ed['DTTorSpr'] = float_read(f.readline().split()[0])
        ed['DTTorDmp'] = float_read(f.readline().split()[0])

        # Furling
        f.readline()
        ed['Furling'] = bool_read(f.readline().split()[0])
        ed['FurlFile'] = os.path.join(str(base_dir), quoted_read(f.readline().split()[0]))

        # Tower
        f.readline()
        ed['TwrNodes'] = int(f.readline().split()[0])
        ed['TwrFile'] = quoted_read(f.readline().split()[0])

        # Output Parameters
        f.readline()
        ed['SumPrint'] = bool_read(f.readline().split()[0])
        ed['OutFile'] = int(f.readline().split()[0])
        ed['TabDelim'] = bool_read(f.readline().split()[0])
        ed['OutFmt'] = quoted_read(f.readline().split()[0])
        ed['TStart'] = float_read(f.readline().split()[0])
        ed['DecFact'] = int(f.readline().split()[0])
        ed['NTwGages'] = int(f.readline().split()[0])
        if ed['NTwGages'] != 0:
            ed['TwrGagNd'] = read_array(f, ed['NTwGages'], array_type=int)
        else:
            ed['TwrGagNd'] = 0
            f.readline()
        ed['NBlGages'] = int(f.readline().split()[0])
        if ed['NBlGages'] != 0:
            ed['BldGagNd'] = read_array(f, ed['NBlGages'], array_type=int)
        else:
            ed['BldGagNd'] = 0

        # OutList — capture into the shared registry (mirrors baseline read_outlist)
        f.readline()
        if outlist is not None:
            capture_outlist(f, outlist, 'ElastoDyn')
        else:
            self._read_outlist(f)

        # Optional nodal output
        try:
            f.readline()
            ed['BldNd_BladesOut'] = int(f.readline().split()[0])
            ed['BldNd_BlOutNd'] = f.readline().split()[0]
            f.readline()
            if outlist is not None:
                capture_outlist(f, outlist, 'ElastoDyn')
            else:
                self._read_outlist(f)
        except:
            None

        f.close()

        # ── Read blade files ──
        blades = [{}, {}, {}]
        bld_files = [ed.get('BldFile1', ''), ed.get('BldFile2', ''), ed.get('BldFile3', '')]
        num_bl = ed.get('NumBl', 3)
        for i in range(num_bl):
            bld_path = base_dir / fix_path(bld_files[i])
            if bld_path.exists():
                blades[i] = self._read_blade(bld_path)

        # ── Read tower file ──
        tower = {}
        twr_file = ed.get('TwrFile', '')
        if twr_file:
            twr_path = base_dir / fix_path(twr_file)
            if twr_path.exists():
                tower = self._read_tower(twr_path)

        return {
            'ElastoDyn': ed,
            'ElastoDynBlade': blades,
            'ElastoDynTower': tower,
        }

    @staticmethod
    def _read_outlist(f):
        """Skip over an OutList section (consumed but not stored here — driver handles it)."""
        data = f.readline()
        while data.strip() == '':
            data = f.readline()
        while data and not data.strip().startswith('END'):
            data = f.readline()
            while data.strip() == '':
                data = f.readline()

    @staticmethod
    def _read_blade(blade_path: Path) -> dict:
        blade = {}
        f = open(blade_path)

        f.readline()
        f.readline()
        f.readline()

        # Blade Parameters
        blade['NBlInpSt'] = int(f.readline().split()[0])
        blade['BldFlDmp1'] = float_read(f.readline().split()[0])
        blade['BldFlDmp2'] = float_read(f.readline().split()[0])
        blade['BldEdDmp1'] = float_read(f.readline().split()[0])

        # Blade Adjustment Factors
        f.readline()
        blade['FlStTunr1'] = float_read(f.readline().split()[0])
        blade['FlStTunr2'] = float_read(f.readline().split()[0])
        blade['AdjBlMs'] = float_read(f.readline().split()[0])
        blade['AdjFlSt'] = float_read(f.readline().split()[0])
        blade['AdjEdSt'] = float_read(f.readline().split()[0])

        # Distributed Blade Properties
        f.readline()
        f.readline()
        f.readline()
        n = blade['NBlInpSt']
        blade['BlFract'] = [None] * n
        blade['StrcTwst'] = [None] * n
        blade['BMassDen'] = [None] * n
        blade['FlpStff'] = [None] * n
        blade['EdgStff'] = [None] * n

        for i in range(n):
            data = f.readline().split()
            blade['BlFract'][i] = float_read(data[0])
            blade['StrcTwst'][i] = float_read(data[1])
            blade['BMassDen'][i] = float_read(data[2])
            blade['FlpStff'][i] = float_read(data[3])
            blade['EdgStff'][i] = float_read(data[4])

        f.readline()
        blade['BldFl1Sh'] = [None] * 5
        blade['BldFl2Sh'] = [None] * 5
        blade['BldEdgSh'] = [None] * 5
        for i in range(5):
            blade['BldFl1Sh'][i] = float_read(f.readline().split()[0])
        for i in range(5):
            blade['BldFl2Sh'][i] = float_read(f.readline().split()[0])
        for i in range(5):
            blade['BldEdgSh'][i] = float_read(f.readline().split()[0])

        f.close()
        return blade

    @staticmethod
    def _read_tower(tower_path: Path) -> dict:
        tower = {}
        f = open(tower_path)

        f.readline()
        f.readline()

        # General Tower Parameters
        f.readline()
        tower['NTwInpSt'] = int(f.readline().split()[0])
        tower['TwrFADmp1'] = float_read(f.readline().split()[0])
        tower['TwrFADmp2'] = float_read(f.readline().split()[0])
        tower['TwrSSDmp1'] = float_read(f.readline().split()[0])
        tower['TwrSSDmp2'] = float_read(f.readline().split()[0])

        # Tower Adjustment Factors
        f.readline()
        tower['FAStTunr1'] = float_read(f.readline().split()[0])
        tower['FAStTunr2'] = float_read(f.readline().split()[0])
        tower['SSStTunr1'] = float_read(f.readline().split()[0])
        tower['SSStTunr2'] = float_read(f.readline().split()[0])
        tower['AdjTwMa'] = float_read(f.readline().split()[0])
        tower['AdjFASt'] = float_read(f.readline().split()[0])
        tower['AdjSSSt'] = float_read(f.readline().split()[0])

        # Distributed Tower Properties
        f.readline()
        f.readline()
        f.readline()
        n = tower['NTwInpSt']
        tower['HtFract'] = [None] * n
        tower['TMassDen'] = [None] * n
        tower['TwFAStif'] = [None] * n
        tower['TwSSStif'] = [None] * n

        for i in range(n):
            data = f.readline().split()
            tower['HtFract'][i] = float_read(data[0])
            tower['TMassDen'][i] = float_read(data[1])
            tower['TwFAStif'][i] = float_read(data[2])
            tower['TwSSStif'][i] = float_read(data[3])

        # Tower Mode Shapes
        f.readline()
        tower['TwFAM1Sh'] = [None] * 5
        tower['TwFAM2Sh'] = [None] * 5
        for i in range(5):
            tower['TwFAM1Sh'][i] = float_read(f.readline().split()[0])
        for i in range(5):
            tower['TwFAM2Sh'][i] = float_read(f.readline().split()[0])
        f.readline()
        tower['TwSSM1Sh'] = [None] * 5
        tower['TwSSM2Sh'] = [None] * 5
        for i in range(5):
            tower['TwSSM1Sh'][i] = float_read(f.readline().split()[0])
        for i in range(5):
            tower['TwSSM2Sh'][i] = float_read(f.readline().split()[0])

        f.close()
        return tower

    # ------------------------------------------------------------------
    # WRITE
    # ------------------------------------------------------------------

    def write(self, data: dict, file_path: Path, base_dir: Path, outlist: dict = None) -> None:
        file_path = Path(file_path)
        base_dir = Path(base_dir)

        # Accept either {'ElastoDyn': {...}} or a flat ED dict
        if 'ElastoDyn' in data:
            ed = data['ElastoDyn']
            blades_data = data.get('ElastoDynBlade', None)
            tower_data = data.get('ElastoDynTower', None)
        else:
            ed = data
            blades_data = None
            tower_data = None

        f = open(file_path, 'w')

        f.write('------- ELASTODYN INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')

        # Simulation Control
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(ed['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['Method'], 'Method', '- Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['DT'], 'DT', '- Integration time step (s)\n'))

        # Degrees of Freedom
        f.write('---------------------- DEGREES OF FREEDOM --------------------------------------\n')
        for key in ['FlapDOF1', 'FlapDOF2', 'EdgeDOF', 'PitchDOF', 'TeetDOF', 'DrTrDOF',
                     'GenDOF', 'YawDOF', 'TwFADOF1', 'TwFADOF2', 'TwSSDOF1', 'TwSSDOF2',
                     'PtfmSgDOF', 'PtfmSwDOF', 'PtfmHvDOF', 'PtfmRDOF', 'PtfmPDOF', 'PtfmYDOF']:
            f.write('{!s:<22} {:<11} {:}'.format(ed[key], key, f'- {key} (flag)\n'))

        # Initial Conditions
        f.write('---------------------- INITIAL CONDITIONS --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['OoPDefl'], 'OoPDefl', '- Initial out-of-plane blade-tip displacement (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['IPDefl'], 'IPDefl', '- Initial in-plane blade-tip deflection (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['BlPitch1'], 'BlPitch(1)', '- Blade 1 initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['BlPitch2'], 'BlPitch(2)', '- Blade 2 initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['BlPitch3'], 'BlPitch(3)', '- Blade 3 initial pitch (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TeetDefl'], 'TeetDefl', '- Initial or fixed teeter angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['Azimuth'], 'Azimuth', '- Initial azimuth angle for blade 1 (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['RotSpeed'], 'RotSpeed', '- Initial or fixed rotor speed (rpm)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['NacYaw'], 'NacYaw', '- Initial or fixed nacelle-yaw angle (degrees)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TTDspFA'], 'TTDspFA', '- Initial fore-aft tower-top displacement (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TTDspSS'], 'TTDspSS', '- Initial side-to-side tower-top displacement (meters)\n'))
        for key in ['PtfmSurge', 'PtfmSway', 'PtfmHeave', 'PtfmRoll', 'PtfmPitch', 'PtfmYaw']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- Initial or fixed platform {key.replace("Ptfm","")} displacement\n'))

        # Turbine Configuration
        f.write('---------------------- TURBINE CONFIGURATION -----------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['NumBl'], 'NumBl', '- Number of blades (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TipRad'], 'TipRad', '- The distance from the rotor apex to the blade tip (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['HubRad'], 'HubRad', '- The distance from the rotor apex to the blade root (meters)\n'))
        for i in range(1, 4):
            f.write('{:<22} {:<11} {:}'.format(ed[f'PreCone({i})'], f'PreCone({i})', f'- Blade {i} cone angle (degrees)\n'))
        for key in ['HubCM', 'UndSling', 'Delta3', 'AzimB1Up', 'OverHang', 'ShftGagL', 'ShftTilt']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key}\n'))
        for key in ['NacCMxn', 'NacCMyn', 'NacCMzn', 'NcIMUxn', 'NcIMUyn', 'NcIMUzn']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key} (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['Twr2Shft'], 'Twr2Shft', '- Vertical distance from the tower-top to the rotor shaft (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TowerHt'], 'TowerHt', '- Height of tower above ground level (meters)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TowerBsHt'], 'TowerBsHt', '- Height of tower base above ground level (meters)\n'))
        for key in ['PtfmCMxt', 'PtfmCMyt', 'PtfmCMzt', 'PtfmRefxt', 'PtfmRefyt', 'PtfmRefzt']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key} (meters)\n'))

        # Mass and Inertia
        f.write('---------------------- MASS AND INERTIA ----------------------------------------\n')
        for i in range(1, 4):
            f.write('{:<22} {:<11} {:}'.format(ed[f'TipMass({i})'], f'TipMass({i})', f'- Tip-brake mass, blade {i} (kg)\n'))
        for i in range(1, 4):
            f.write('{:<22} {:<11} {:}'.format(ed[f'PBrIner({i})'], f'PBrIner({i})', f'- Pitch bearing inertia, blade {i} (kg m^2)\n'))
        for i in range(1, 4):
            f.write('{:<22} {:<11} {:}'.format(ed[f'BlPIner({i})'], f'BlPIner({i})', f'- Blade pitch inertia, blade {i} (kg m^2)\n'))
        for key in ['HubMass', 'HubIner', 'HubIner_Teeter', 'GenIner', 'NacMass', 'NacYIner',
                     'YawBrMass', 'PtfmMass', 'PtfmRIner', 'PtfmPIner', 'PtfmYIner',
                     'PtfmXYIner', 'PtfmYZIner', 'PtfmXZIner']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key}\n'))

        # Blade
        f.write('---------------------- BLADE ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['BldNodes'], 'BldNodes', '- Number of blade nodes (per blade) used for analysis (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['BldFile1']+'"', 'BldFile(1)', '- Name of file containing properties for blade 1\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['BldFile2']+'"', 'BldFile(2)', '- Name of file containing properties for blade 2\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['BldFile3']+'"', 'BldFile(3)', '- Name of file containing properties for blade 3\n'))

        # Rotor-Teeter
        f.write('---------------------- ROTOR-TEETER --------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['TeetMod'], 'TeetMod', '- Rotor-teeter spring/damper model (switch)\n'))
        for key in ['TeetDmpP', 'TeetDmp', 'TeetCDmp', 'TeetSStP', 'TeetHStP', 'TeetSSSp', 'TeetHSSp']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key}\n'))

        # Yaw Friction
        f.write('---------------------- YAW-FRICTION --------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['YawFrctMod'], 'YawFrctMod', '- Yaw-friction model (switch)\n'))
        for key in ['M_CSmax', 'M_FCSmax', 'M_MCSmax', 'M_CD', 'M_FCD', 'M_MCD', 'sig_v', 'sig_v2', 'OmgCut']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key}\n'))

        # Drivetrain
        f.write('---------------------- DRIVETRAIN ----------------------------------------------\n')
        for key in ['GBoxEff', 'GBRatio', 'DTTorSpr', 'DTTorDmp']:
            f.write('{:<22} {:<11} {:}'.format(ed[key], key, f'- {key}\n'))

        # Furling
        f.write('---------------------- FURLING -------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(ed['Furling'], 'Furling', '- Read in additional model properties for furling turbine (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['FurlFile']+'"', 'FurlFile', '- Name of file containing furling properties\n'))

        # Tower
        f.write('---------------------- TOWER ---------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(ed['TwrNodes'], 'TwrNodes', '- Number of tower nodes used for analysis (-)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['TwrFile']+'"', 'TwrFile', '- Name of file containing tower properties\n'))

        # Output
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(ed['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['OutFile'], 'OutFile', '- Switch to determine where output will be placed (currently unused)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(ed['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+ed['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['TStart'], 'TStart', '- Time to begin tabular output (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['DecFact'], 'DecFact', '- Decimation factor for tabular output (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['NTwGages'], 'NTwGages', '- Number of tower nodes that have strain gages for output (-)\n'))
        if ed['TwrGagNd'] != 0:
            f.write('{:<22} {:<11} {:}'.format(', '.join(['%d'%int(i) for i in ed['TwrGagNd']]), 'TwrGagNd', '- List of tower nodes that have strain gages\n'))
        else:
            f.write('{:<22} {:<11} {:}'.format('', 'TwrGagNd', '- List of tower nodes that have strain gages\n'))
        f.write('{:<22} {:<11} {:}'.format(ed['NBlGages'], 'NBlGages', '- Number of blade nodes that have strain gages for output (-)\n'))
        if ed['BldGagNd'] != 0:
            f.write('{:<22} {:<11} {:}'.format(', '.join(['%d'%int(i) for i in ed['BldGagNd']]), 'BldGagNd', '- List of blade nodes that have strain gages\n'))
        else:
            f.write('{:<22} {:<11} {:}'.format('', 'BldGagNd', '- List of blade nodes that have strain gages\n'))

        # OutList — emit the captured channels (mirrors baseline write)
        f.write('                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n')
        if outlist is not None:
            emit_outlist(f, outlist, 'ElastoDyn')
        f.write('END of OutList section (the word "END" must appear in the first 3 columns of the last OutList line)\n')

        # Optional nodal output
        if 'BldNd_BladesOut' in ed:
            f.write('====== Outputs for all blade stations =========================== [optional section]\n')
            f.write('{:<22d} {:<11} {:}'.format(ed['BldNd_BladesOut'], 'BldNd_BladesOut', '- Number of blades to output all node information at (-)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ed['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output (-)\n'))
            f.write('                   OutList     - The next line(s) contains a list of output parameters.\n')
            f.write('END (the word "END" must appear in the first 3 columns of this last OutList line in the optional nodal output section)\n')

        f.write('---------------------------------------------------------------------------------------\n')
        f.flush()
        os.fsync(f)
        f.close()

        # Write blade files
        if blades_data is not None and blades_data:
            if isinstance(blades_data, dict):
                # Single dict = all blades identical
                for i in range(1, 4):
                    bld_file = ed.get(f'BldFile{i}', '')
                    if bld_file:
                        self._write_blade(blades_data, base_dir / bld_file)
            elif isinstance(blades_data, list):
                for i, blade in enumerate(blades_data):
                    bld_file = ed.get(f'BldFile{i+1}', '')
                    if bld_file and blade:
                        self._write_blade(blade, base_dir / bld_file)

        # Write tower file
        if tower_data is not None and tower_data:
            twr_file = ed.get('TwrFile', '')
            if twr_file:
                self._write_tower(tower_data, base_dir / twr_file)

    @staticmethod
    def _write_blade(blade: dict, blade_path: Path) -> None:
        f = open(blade_path, 'w')

        f.write('------- ELASTODYN INDIVIDUAL BLADE INPUT FILE --------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- BLADE PARAMETERS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(blade['NBlInpSt'], 'NBlInpSt', '- Number of blade input stations (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['BldFlDmp1'], 'BldFlDmp(1)', '- Blade flap mode #1 structural damping in percent of critical (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['BldFlDmp2'], 'BldFlDmp(2)', '- Blade flap mode #2 structural damping in percent of critical (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['BldEdDmp1'], 'BldEdDmp(1)', '- Blade edge mode #1 structural damping in percent of critical (%)\n'))
        f.write('---------------------- BLADE ADJUSTMENT FACTORS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(blade['FlStTunr1'], 'FlStTunr(1)', '- Blade flapwise modal stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['FlStTunr2'], 'FlStTunr(2)', '- Blade flapwise modal stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['AdjBlMs'], 'AdjBlMs', '- Factor to adjust blade mass density (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['AdjFlSt'], 'AdjFlSt', '- Factor to adjust blade flap stiffness (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(blade['AdjEdSt'], 'AdjEdSt', '- Factor to adjust blade edge stiffness (-)\n'))
        f.write('---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------\n')
        f.write('    BlFract      StrcTwst       BMassDen        FlpStff        EdgStff\n')
        f.write('      (-)         (deg)          (kg/m)         (Nm^2)         (Nm^2)\n')
        for i in range(blade['NBlInpSt']):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(
                blade['BlFract'][i], blade['StrcTwst'][i], blade['BMassDen'][i],
                blade['FlpStff'][i], blade['EdgStff'][i]))
        f.write('---------------------- BLADE MODE SHAPES ---------------------------------------\n')
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(blade['BldFl1Sh'][i], f'BldFl1Sh({i+2})', f'- Flap mode 1, coeff of x^{i+2}\n'))
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(blade['BldFl2Sh'][i], f'BldFl2Sh({i+2})', f'- Flap mode 2, coeff of x^{i+2}\n'))
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(blade['BldEdgSh'][i], f'BldEdgSh({i+2})', f'- Edge mode 1, coeff of x^{i+2}\n'))

        f.flush()
        os.fsync(f)
        f.close()

    @staticmethod
    def _write_tower(tower: dict, tower_path: Path) -> None:
        f = open(tower_path, 'w')

        f.write('------- ELASTODYN TOWER INPUT FILE -------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- TOWER PARAMETERS ----------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(tower['NTwInpSt'], 'NTwInpSt', '- Number of input stations to specify tower geometry\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['TwrFADmp1'], 'TwrFADmp(1)', '- Tower 1st fore-aft mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['TwrFADmp2'], 'TwrFADmp(2)', '- Tower 2nd fore-aft mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['TwrSSDmp1'], 'TwrSSDmp(1)', '- Tower 1st side-to-side mode structural damping ratio (%)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['TwrSSDmp2'], 'TwrSSDmp(2)', '- Tower 2nd side-to-side mode structural damping ratio (%)\n'))
        f.write('---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(tower['FAStTunr1'], 'FAStTunr(1)', '- Tower fore-aft modal stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['FAStTunr2'], 'FAStTunr(2)', '- Tower fore-aft modal stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['SSStTunr1'], 'SSStTunr(1)', '- Tower side-to-side stiffness tuner, 1st mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['SSStTunr2'], 'SSStTunr(2)', '- Tower side-to-side stiffness tuner, 2nd mode (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['AdjTwMa'], 'AdjTwMa', '- Factor to adjust tower mass density (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['AdjFASt'], 'AdjFASt', '- Factor to adjust tower fore-aft stiffness (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(tower['AdjSSSt'], 'AdjSSSt', '- Factor to adjust tower side-to-side stiffness (-)\n'))
        f.write('---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------\n')
        f.write('  HtFract       TMassDen         TwFAStif       TwSSStif\n')
        f.write('   (-)           (kg/m)           (Nm^2)         (Nm^2)\n')
        for i in range(tower['NTwInpSt']):
            f.write('{: 2.15e} {: 2.15e} {: 2.15e} {: 2.15e}\n'.format(
                tower['HtFract'][i], tower['TMassDen'][i],
                tower['TwFAStif'][i], tower['TwSSStif'][i]))
        f.write('---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n')
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(tower['TwFAM1Sh'][i], f'TwFAM1Sh({i+2})', f'- Mode 1, coefficient of x^{i+2} term\n'))
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(tower['TwFAM2Sh'][i], f'TwFAM2Sh({i+2})', f'- Mode 2, coefficient of x^{i+2} term\n'))
        f.write('---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n')
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(tower['TwSSM1Sh'][i], f'TwSSM1Sh({i+2})', f'- Mode 1, coefficient of x^{i+2} term\n'))
        for i in range(5):
            f.write('{:<22} {:<11} {:}'.format(tower['TwSSM2Sh'][i], f'TwSSM2Sh({i+2})', f'- Mode 2, coefficient of x^{i+2} term\n'))

        f.flush()
        os.fsync(f)
        f.close()
