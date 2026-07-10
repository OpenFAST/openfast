"""OpenFAST driver — orchestrates module-level IO classes to read/write complete simulation decks.

This module provides:
    - init_fst_vt(): Initialize the fst_vt variable tree structure
    - OpenFASTDriver: Compose module IOs to read/write .fst + all referenced files
"""
import os, copy
from pathlib import Path

from ..io.elastodyn import ElastoDynIO
from ..io.aerodyn import AeroDynIO
from ..io.inflowwind import InflowWindIO
from ..io.beamdyn import BeamDynIO
from ..io.servodyn import ServoDynIO
from ..io.hydrodyn import HydroDynIO
from ..io.seastate import SeaStateIO
from ..io.subdyn import SubDynIO
from ..io.moordyn import MoorDynIO
from ..io.map_io import MAPIO
from ..io.extptfm import ExtPtfmIO
from ..io.simple_elastodyn import SimpleElastoDynIO
from ..io.aerodisk import AeroDiskIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
    fix_path,
)

try:
    from openfast_io.FAST_vars_out import FstOutput
except ImportError:
    FstOutput = {}

from openfast_io.outlist import OutList, capture_outlist, emit_outlist


def init_fst_vt() -> dict:
    """Initialize the fst_vt structure.

    Key names and structure match the legacy InputReader_OpenFAST.fst_vt.
    WEIS and downstream code access these keys directly — do not rename or remove.

    Notes:
        - ElastoDynBlade, AeroDynBlade, BeamDynBlade default to [{},{},{}] (list
          of per-blade dicts). The reader's blade-dedup logic may collapse to a
          single dict when all blades are identical.
        - BStC/NStC/TStC/SStC are lists of StC parameter dicts, empty by default.
        - _path keys are set dynamically during read(), not here.
    """
    return {
        'Fst': {},
        # Mirror the legacy openfast_io reader exactly: start from FstOutput (its
        # built-in defaults) and ADD the deck's channels during read() via
        # capture_outlist. The legacy reader does deepcopy(FstOutput) + set_outlist
        # per module; reproducing that init is required for behavior parity (the prior
        # bug was skipping the per-module capture entirely, leaving ONLY the defaults).
        'outlist': copy.deepcopy(FstOutput) if FstOutput else {},
        'description': '',
        'ElastoDyn': {},
        'SimpleElastoDyn': {},
        'ElastoDynBlade': [{}, {}, {}],
        'ElastoDynTower': {},
        'InflowWind': {},
        'AeroDyn': {},
        'AeroDisk': {},
        'AeroDynBlade': [{}, {}, {}],
        'AeroDynPolar': [],
        'ServoDyn': {},
        'DISCON_in': {},
        'spd_trq': {},
        'BStC': [],
        'NStC': [],
        'TStC': [],
        'SStC': [],
        'HydroDyn': {},
        'SeaState': {},
        'MoorDyn': {},
        'SubDyn': {},
        'ExtPtfm': {},
        'MAP': {},
        'BeamDyn': [{}, {}, {}],
        'BeamDynBlade': [{}, {}, {}],
        'WaterKin': {},
        'SoilDyn': {},
    }


class OpenFASTDriver:
    """Reads and writes a complete OpenFAST simulation deck.

    Composes module-level IO objects. Returns fst_vt with the same
    structure as the legacy InputReader_OpenFAST.fst_vt.
    """

    def __init__(self):
        self._elastodyn = ElastoDynIO()
        self._aerodyn = AeroDynIO()
        self._inflowwind = InflowWindIO()
        self._beamdyn = BeamDynIO()
        self._servodyn = ServoDynIO()
        self._hydrodynamics = HydroDynIO()
        self._seastate = SeaStateIO()
        self._subdyn = SubDynIO()
        self._moordyn = MoorDynIO()
        self._map = MAPIO()
        self._extptfm = ExtPtfmIO()
        self._simple_elastodyn = SimpleElastoDynIO()
        self._aerodisk = AeroDiskIO()

    def read(self, fst_path: Path) -> dict:
        """Read .fst and all referenced module files. Returns complete fst_vt."""
        fst_path = Path(fst_path)
        fst_vt = init_fst_vt()
        base_dir = fst_path.parent

        # Callback that captures a module's OutList section into fst_vt['outlist'],
        # mirroring legacy openfast_io read_outlist/set_outlist. Modules that take a
        # read_outlist_fn use this; freeform modules (SubDyn/SeaState) pass freeform=True.
        def _cap(f, module, freeform=False):
            return capture_outlist(f, fst_vt['outlist'], module, freeform=freeform)

        def _cap_ff(f, module):
            return capture_outlist(f, fst_vt['outlist'], module, freeform=True)

        fst_vt['Fst'] = self._read_main_input(fst_path, base_dir)

        n_rotors = fst_vt['Fst'].get('NRotors', 1)
        if n_rotors > 1:
            raise ValueError(
                'openfast_io does not currently support multi-rotor turbines (NRotors > 1), '
                'this feature will be added in a future release'
            )

        ed_file = os.path.join(str(base_dir), fst_vt['Fst'].get('EDFile', ''))
        fastdir = str(base_dir)

        # ------- ElastoDyn -------
        comp_elast = fst_vt['Fst'].get('CompElast', 1)
        if comp_elast == 3:
            sed_rel = fst_vt['Fst'].get('EDFile', '')
            sed_file = os.path.normpath(os.path.join(fastdir, sed_rel))
            if os.path.isfile(sed_file):
                sed_data = self._simple_elastodyn.read(sed_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['SimpleElastoDyn'] = sed_data.get('SimpleElastoDyn', {})
        elif comp_elast in (1, 2):
            if os.path.isfile(ed_file):
                fst_vt['Fst']['EDFile_path'] = os.path.split(fst_vt['Fst']['EDFile'])[0]
                ed_data = self._elastodyn.read(Path(ed_file), Path(os.path.dirname(ed_file)), outlist=fst_vt['outlist'])

                fst_vt['ElastoDyn'] = ed_data.get('ElastoDyn', {})
                fst_vt['ElastoDynTower'] = ed_data.get('ElastoDynTower', {})

                # Blade deduplication logic — match legacy behavior
                blades = ed_data.get('ElastoDynBlade', [{}, {}, {}])
                bldFile1 = fst_vt['ElastoDyn'].get('BldFile1', '')
                bldFile2 = fst_vt['ElastoDyn'].get('BldFile2', '')
                bldFile3 = fst_vt['ElastoDyn'].get('BldFile3', '')
                num_bl = fst_vt['ElastoDyn'].get('NumBl', 3)

                if bldFile1 == bldFile2 and bldFile1 == bldFile3:
                    fst_vt['ElastoDynBlade'] = blades[0] if blades else {}
                elif num_bl == 2 and bldFile1 == bldFile2:
                    fst_vt['ElastoDynBlade'] = blades[0] if blades else {}
                elif num_bl == 1:
                    fst_vt['ElastoDynBlade'] = blades[0] if blades else {}
                else:
                    fst_vt['ElastoDynBlade'] = blades

        # ------- BeamDyn (per-blade) -------
        # Match the legacy openfast_io FAST_reader: read BeamDyn whenever BDBldFile(1)
        # exists on disk, regardless of CompElast. The legacy reader reads (and captures
        # the OutList of) an inactive BeamDyn when its blade file is present; reproduce that.
        _bd1 = os.path.normpath(os.path.join(fastdir, fst_vt['Fst'].get('BDBldFile(1)', '')))
        if comp_elast == 2 or os.path.isfile(_bd1):
            num_bl = fst_vt['ElastoDyn'].get('NumBl', 3)
            bd_blades = []
            bd_blade_data = []
            for i in range(min(num_bl, 3)):
                bd_file_key = f'BDBldFile({i+1})'
                bd_rel = fst_vt['Fst'].get(bd_file_key, '')
                bd_file = os.path.normpath(os.path.join(fastdir, bd_rel))
                if os.path.isfile(bd_file):
                    bd_data = self._beamdyn.read(bd_file, base_dir=os.path.dirname(bd_file), outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                    bd_blades.append(bd_data.get('BeamDyn', {}))
                    bd_blade_data.append(bd_data.get('BeamDynBlade', {}))
                else:
                    bd_blades.append({})
                    bd_blade_data.append({})
            # Blade dedup — match legacy (FAST_reader): collapse to a single dict when
            # the BeamDyn blade files are identical. Mirrors the ElastoDyn block above.
            bd1 = fst_vt['Fst'].get('BDBldFile(1)', '')
            bd2 = fst_vt['Fst'].get('BDBldFile(2)', '')
            bd3 = fst_vt['Fst'].get('BDBldFile(3)', '')
            if (bd1 == bd2 == bd3) or num_bl == 1 or (num_bl == 2 and bd1 == bd2):
                fst_vt['BeamDyn'] = bd_blades[0] if bd_blades else {}
                fst_vt['BeamDynBlade'] = bd_blade_data[0] if bd_blade_data else {}
            else:
                fst_vt['BeamDyn'] = bd_blades
                fst_vt['BeamDynBlade'] = bd_blade_data

        # ------- InflowWind -------
        comp_inflow = fst_vt['Fst'].get('CompInflow', 0)
        if comp_inflow == 1:
            ifw_rel = fst_vt['Fst'].get('InflowFile', '')
            ifw_file = os.path.normpath(os.path.join(fastdir, ifw_rel))
            if os.path.isfile(ifw_file):
                ifw_data = self._inflowwind.read(ifw_file, base_dir=os.path.dirname(ifw_file), outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['InflowWind'] = ifw_data.get('InflowWind', {})

        # ------- AeroDyn -------
        comp_aero = fst_vt['Fst'].get('CompAero', 0)
        if comp_aero == 2:
            aero_rel = fst_vt['Fst'].get('AeroFile', '')
            aero_file = os.path.normpath(os.path.join(fastdir, aero_rel))
            if os.path.isfile(aero_file):
                num_bl = fst_vt['ElastoDyn'].get('NumBl', 3)
                ad_data = self._aerodyn.read(
                    aero_file,
                    base_dir=os.path.dirname(aero_file),
                    num_blades=num_bl,
                    aero_file_path=fst_vt['Fst'].get('AeroFile_path', ''),
                    outlist=fst_vt['outlist'],
                    read_outlist_fn=_cap,
                )
                fst_vt['AeroDyn'] = ad_data.get('AeroDyn', {})

                # AeroDynBlade dedup (same pattern as ElastoDyn)
                ad_blade = ad_data.get('AeroDynBlade', {})
                if isinstance(ad_blade, dict):
                    fst_vt['AeroDynBlade'] = ad_blade  # already collapsed
                else:
                    fst_vt['AeroDynBlade'] = ad_blade

                fst_vt['AeroDynPolar'] = ad_data.get('AeroDynPolar', [])
        elif comp_aero == 1:  # AeroDisk {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}
            aero_rel = fst_vt['Fst'].get('AeroFile', '')
            aero_file = os.path.normpath(os.path.join(fastdir, aero_rel))
            if os.path.isfile(aero_file):
                adsk_data = self._aerodisk.read(aero_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['AeroDisk'] = adsk_data.get('AeroDisk', {})

        # ------- ServoDyn -------
        comp_servo = fst_vt['Fst'].get('CompServo', 0)
        if comp_servo == 1:
            sd_rel = fst_vt['Fst'].get('ServoFile', '')
            sd_file = os.path.normpath(os.path.join(fastdir, sd_rel))
            if os.path.isfile(sd_file):
                sd_data = self._servodyn.read(
                    sd_file,
                    base_dir=fastdir,
                    servo_file_rel=sd_rel,
                    outlist=fst_vt['outlist'],
                    read_outlist_fn=_cap,
                )
                fst_vt['ServoDyn'] = sd_data.get('ServoDyn', {})
                fst_vt['BStC'] = sd_data.get('BStC', [])
                fst_vt['NStC'] = sd_data.get('NStC', [])
                fst_vt['TStC'] = sd_data.get('TStC', [])
                fst_vt['SStC'] = sd_data.get('SStC', [])
                if 'DISCON_in' in sd_data:
                    fst_vt['DISCON_in'] = sd_data['DISCON_in']
                if 'spd_trq' in sd_data:
                    fst_vt['spd_trq'] = sd_data['spd_trq']

        # ------- HydroDyn -------
        comp_hydro = fst_vt['Fst'].get('CompHydro', 0)
        if comp_hydro == 1:
            hd_rel = fst_vt['Fst'].get('HydroFile', '')
            hd_file = os.path.normpath(os.path.join(fastdir, hd_rel))
            if os.path.isfile(hd_file):
                fst_vt['Fst']['HydroFile_path'] = os.path.split(hd_rel)[0]
                hd_data = self._hydrodynamics.read(hd_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['HydroDyn'] = hd_data.get('HydroDyn', {})

        # ------- SeaState -------
        comp_seast = fst_vt['Fst'].get('CompSeaSt', 0)
        if comp_seast == 1:
            ss_rel = fst_vt['Fst'].get('SeaStFile', '')
            ss_file = os.path.normpath(os.path.join(fastdir, ss_rel))
            if os.path.isfile(ss_file):
                ss_data = self._seastate.read(ss_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap_ff)
                fst_vt['SeaState'] = ss_data.get('SeaState', {})

        # ------- SubDyn / ExtPtfm -------
        comp_sub = fst_vt['Fst'].get('CompSub', 0)
        if comp_sub == 1:
            sub_rel = fst_vt['Fst'].get('SubFile', '')
            sub_file = os.path.normpath(os.path.join(fastdir, sub_rel))
            if os.path.isfile(sub_file):
                fst_vt['Fst']['SubFile_path'] = os.path.split(sub_rel)[0]
                sub_data = self._subdyn.read(sub_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap_ff)
                fst_vt['SubDyn'] = sub_data.get('SubDyn', {})
        elif comp_sub == 2:
            sub_rel = fst_vt['Fst'].get('SubFile', '')
            sub_file = os.path.normpath(os.path.join(fastdir, sub_rel))
            if os.path.isfile(sub_file):
                fst_vt['Fst']['SubFile_path'] = os.path.split(sub_rel)[0]
                ep_data = self._extptfm.read(sub_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['ExtPtfm'] = ep_data.get('ExtPtfm', {})

        # ------- MoorDyn / MAP -------
        comp_mooring = fst_vt['Fst'].get('CompMooring', 0)
        if comp_mooring == 1:  # MAP
            moor_rel = fst_vt['Fst'].get('MooringFile', '')
            moor_file = os.path.normpath(os.path.join(fastdir, moor_rel))
            if os.path.isfile(moor_file):
                fst_vt['Fst']['MooringFile_path'] = os.path.split(moor_rel)[0]
                map_data = self._map.read(moor_file)
                fst_vt['MAP'] = map_data.get('MAP', {})
        elif comp_mooring == 3:  # MoorDyn
            moor_rel = fst_vt['Fst'].get('MooringFile', '')
            moor_file = os.path.normpath(os.path.join(fastdir, moor_rel))
            if os.path.isfile(moor_file):
                fst_vt['Fst']['MooringFile_path'] = os.path.split(moor_rel)[0]
                md_data = self._moordyn.read(moor_file, outlist=fst_vt['outlist'], read_outlist_fn=_cap)
                fst_vt['MoorDyn'] = md_data.get('MoorDyn', {})

        return fst_vt

    def write(self, fst_vt: dict, output_dir: Path, case_name: str) -> list:
        """Write all enabled module files + .fst. Returns list of written paths.

        Parameters
        ----------
        fst_vt : dict
            Complete variable tree (same shape as returned by read()).
        output_dir : Path
            Directory where all files will be written.
        case_name : str
            Base name for the .fst and sub-files (e.g. "5MW_Land").
        """
        import numpy as np

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written = []
        fst = fst_vt['Fst']

        # ------- ElastoDyn / SimpleElastoDyn -------
        comp_elast = fst.get('CompElast', 1)
        if comp_elast == 3:
            sed_name = case_name + '_SimpleElastoDyn.dat'
            fst['EDFile'] = sed_name
            sed_path = str(output_dir / sed_name)
            self._simple_elastodyn.write({'SimpleElastoDyn': fst_vt.get('SimpleElastoDyn', {})}, sed_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(sed_path)
        elif comp_elast in (1, 2):
            # ElastoDyn blade(s)
            edblade = fst_vt.get('ElastoDynBlade', {})
            if isinstance(edblade, list):
                for i, bld in enumerate(edblade):
                    bld_name = '{}_{}_ElastoDynBlade_{}.dat'.format(case_name, 'ElastoDyn', i + 1)
                    fst_vt['ElastoDyn']['BldFile{}'.format(i + 1)] = bld_name
                    bld_path = str(output_dir / bld_name)
                    self._elastodyn._write_blade(bld, bld_path)
                    written.append(bld_path)
            elif isinstance(edblade, dict) and edblade:
                bld_name = case_name + '_ElastoDynBlade.dat'
                fst_vt['ElastoDyn']['BldFile1'] = bld_name
                fst_vt['ElastoDyn']['BldFile2'] = bld_name
                fst_vt['ElastoDyn']['BldFile3'] = bld_name
                bld_path = str(output_dir / bld_name)
                self._elastodyn._write_blade(edblade, bld_path)
                written.append(bld_path)

            # Tower
            twr_name = case_name + '_ElastoDynTower.dat'
            fst_vt['ElastoDyn']['TwrFile'] = twr_name
            twr_path = str(output_dir / twr_name)
            self._elastodyn._write_tower(fst_vt.get('ElastoDynTower', {}), twr_path)
            written.append(twr_path)

            # Main ED
            ed_name = case_name + '_ElastoDyn.dat'
            fst['EDFile'] = ed_name
            ed_path = str(output_dir / ed_name)
            self._elastodyn.write(
                {'ElastoDyn': fst_vt.get('ElastoDyn', {}),
                 'ElastoDynTower': fst_vt.get('ElastoDynTower', {}),
                 'ElastoDynBlade': fst_vt.get('ElastoDynBlade', {})},
                ed_path,
                base_dir=str(output_dir),
                outlist=fst_vt.get('outlist'),
            )
            written.append(ed_path)

        # ------- BeamDyn -------
        # Write whenever fst_vt['BeamDyn'] is populated (symmetric with the read guard,
        # which reads BeamDyn whenever the blade file exists) — not gated on CompElast.
        # Handle both the collapsed-dict shape (identical blades) and the list shape
        # (distinct blades), and assign a UNIQUE per-blade BldFile so the writer does not
        # send every blade to the same path (silent blade-property collision).
        bd_data = fst_vt.get('BeamDyn')
        bd_blade_all = fst_vt.get('BeamDynBlade')
        num_bl_w = fst_vt.get('ElastoDyn', {}).get('NumBl', 3)

        def _write_bd_blade(idx, bd_src, blade, main_name, blade_name):
            bd = dict(bd_src)
            bd['BldFile'] = blade_name
            bd_path = str(output_dir / main_name)
            self._beamdyn.write(
                {'BeamDyn': bd, 'BeamDynBlade': blade},
                bd_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'),
            )
            written.append(bd_path)

        if isinstance(bd_data, dict) and bd_data:
            # Collapsed identical blades: one file, all BDBldFile(i) point at it.
            blade = bd_blade_all if isinstance(bd_blade_all, dict) else \
                (bd_blade_all[0] if isinstance(bd_blade_all, list) and bd_blade_all else {})
            main_name = case_name + '_BeamDyn.dat'
            _write_bd_blade(0, bd_data, blade, main_name, case_name + '_BeamDyn_Blade.dat')
            for k in range(num_bl_w):
                fst['BDBldFile({})'.format(k + 1)] = main_name
        elif isinstance(bd_data, list):
            for i, bd in enumerate(bd_data):
                if bd:
                    main_name = case_name + '_BeamDyn_{}.dat'.format(i + 1)
                    fst['BDBldFile({})'.format(i + 1)] = main_name
                    blade = bd_blade_all[i] if isinstance(bd_blade_all, list) and i < len(bd_blade_all) else {}
                    _write_bd_blade(i, bd, blade, main_name, case_name + '_BeamDyn_Blade_{}.dat'.format(i + 1))

        # ------- InflowWind -------
        if fst.get('CompInflow', 0) == 1:
            ifw_name = case_name + '_InflowWind.dat'
            fst['InflowFile'] = ifw_name
            ifw_path = str(output_dir / ifw_name)
            self._inflowwind.write({'InflowWind': fst_vt.get('InflowWind', {})}, ifw_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(ifw_path)

        # ------- AeroDyn / AeroDisk -------
        comp_aero = fst.get('CompAero', 0)
        if comp_aero == 2:
            ad_name = case_name + '_AeroDyn.dat'
            fst['AeroFile'] = ad_name
            ad_path = str(output_dir / ad_name)
            self._aerodyn.write(
                {'AeroDyn': fst_vt.get('AeroDyn', {}),
                 'AeroDynBlade': fst_vt.get('AeroDynBlade', {}),
                 'AeroDynPolar': fst_vt.get('AeroDynPolar', [])},
                ad_path,
                base_dir=str(output_dir),
                outlist=fst_vt.get('outlist'),
            )
            written.append(ad_path)
        elif comp_aero == 1:  # AeroDisk {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}
            adsk_name = case_name + '_AeroDisk.dat'
            fst['AeroFile'] = adsk_name
            adsk_path = str(output_dir / adsk_name)
            self._aerodisk.write({'AeroDisk': fst_vt.get('AeroDisk', {})}, adsk_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(adsk_path)

        # ------- ServoDyn -------
        if fst.get('CompServo', 0) == 1:
            sd_name = case_name + '_ServoDyn.dat'
            fst['ServoFile'] = sd_name
            sd_path = str(output_dir / sd_name)
            self._servodyn.write(
                {'ServoDyn': fst_vt.get('ServoDyn', {}),
                 'BStC': fst_vt.get('BStC', []),
                 'NStC': fst_vt.get('NStC', []),
                 'TStC': fst_vt.get('TStC', []),
                 'SStC': fst_vt.get('SStC', []),
                 'spd_trq': fst_vt.get('spd_trq')},
                sd_path,
                base_dir=str(output_dir),
                outlist=fst_vt.get('outlist'),
            )
            written.append(sd_path)

        # ------- SeaState -------
        if fst.get('CompSeaSt', 0) == 1:
            ss_name = case_name + '_SeaState.dat'
            fst['SeaStFile'] = ss_name
            ss_path = str(output_dir / ss_name)
            self._seastate.write({'SeaState': fst_vt.get('SeaState', {})}, ss_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(ss_path)

        # ------- HydroDyn -------
        if fst.get('CompHydro', 0) == 1:
            hd_name = case_name + '_HydroDyn.dat'
            fst['HydroFile'] = hd_name
            hd_path = str(output_dir / hd_name)
            self._hydrodynamics.write({'HydroDyn': fst_vt.get('HydroDyn', {})}, hd_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(hd_path)

        # ------- SubDyn / ExtPtfm -------
        comp_sub = fst.get('CompSub', 0)
        if comp_sub == 1:
            sub_name = case_name + '_SubDyn.dat'
            fst['SubFile'] = sub_name
            sub_path = str(output_dir / sub_name)
            self._subdyn.write({'SubDyn': fst_vt.get('SubDyn', {})}, sub_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(sub_path)
        elif comp_sub == 2:
            ep_name = case_name + '_ExtPtfm.dat'
            fst['SubFile'] = ep_name
            ep_path = str(output_dir / ep_name)
            self._extptfm.write({'ExtPtfm': fst_vt.get('ExtPtfm', {})}, ep_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(ep_path)

        # ------- MoorDyn / MAP -------
        comp_mooring = fst.get('CompMooring', 0)
        if comp_mooring == 1:
            map_name = case_name + '_MAP.dat'
            fst['MooringFile'] = map_name
            map_path = str(output_dir / map_name)
            self._map.write({'MAP': fst_vt.get('MAP', {})}, map_path, base_dir=str(output_dir))
            written.append(map_path)
        elif comp_mooring == 3:
            md_name = case_name + '_MoorDyn.dat'
            fst['MooringFile'] = md_name
            md_path = str(output_dir / md_name)
            self._moordyn.write({'MoorDyn': fst_vt.get('MoorDyn', {})}, md_path, base_dir=str(output_dir), outlist=fst_vt.get('outlist'))
            written.append(md_path)

        # ------- Main .fst file -------
        fst_path = str(output_dir / (case_name + '.fst'))
        self._write_main_input(fst_vt, fst_path)
        written.append(fst_path)

        return written

    def _write_main_input(self, fst_vt: dict, fst_path: str) -> None:
        """Write the main .fst input file."""
        import numpy as np
        fst = fst_vt['Fst']

        with open(fst_path, 'w') as f:
            f.write('------- OpenFAST INPUT FILE -------------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')

            # Simulation Control
            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(fst['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + str(fst.get('AbortLevel', 'FATAL')) + '"', 'AbortLevel', '- Error level when simulation should abort\n'))
            f.write('{:<22} {:<11} {:}'.format(fst['TMax'], 'TMax', '- Total run time (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst['DT'], 'DT', '- Recommended module time step (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('ModCoupling', 1), 'ModCoupling', '- Module coupling method (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('InterpOrder', 1), 'InterpOrder', '- Interpolation order\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('NumCrctn', 0), 'NumCrctn', '- Numerical damping parameter\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('RhoInf', 1), 'RhoInf', '- Convergence iteration error tolerance\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('ConvTol', 1e-5), 'ConvTol', '- Convergence tolerance\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('MaxConvIter', 4), 'MaxConvIter', '- Maximum number of convergence iterations\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('DT_UJac', 'default'), 'DT_UJac', '- Time between Jacobian calls\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('UJacSclFact', 1e6), 'UJacSclFact', '- Scaling factor for Jacobians\n'))

            # Feature Switches
            f.write('---------------------- FEATURE SWITCHES AND FLAGS ------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(fst.get('NRotors', 1), 'NRotors', '- Number of rotors\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompElast', 1), 'CompElast', '- Compute structural dynamics\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompInflow', 0), 'CompInflow', '- Compute inflow wind velocities\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompAero', 0), 'CompAero', '- Compute aerodynamic loads\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompServo', 0), 'CompServo', '- Compute control and electrical-drive dynamics\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompSeaSt', 0), 'CompSeaSt', '- Compute sea state information\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompHydro', 0), 'CompHydro', '- Compute hydrodynamic loads\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompSub', 0), 'CompSub', '- Compute sub-structural dynamics\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompMooring', 0), 'CompMooring', '- Compute mooring system\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompIce', 0), 'CompIce', '- Compute ice loads\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('CompSoil', 0), 'CompSoil', '- Compute soil dynamics\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('MHK', 0), 'MHK', '- MHK turbine type\n'))
            mirror = fst.get('MirrorRotor', [False])
            f.write('{:<22} {:<11} {:}'.format(' '.join([str(b)[0] for b in np.array(mirror, dtype=bool)]), 'MirrorRotor', '- List of rotor rotation directions\n'))

            # Environmental Conditions
            f.write('---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(fst.get('Gravity', 9.81), 'Gravity', '- Gravitational acceleration (m/s^2)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('AirDens', 1.225), 'AirDens', '- Air density (kg/m^3)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('WtrDens', 1025.0), 'WtrDens', '- Water density (kg/m^3)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('KinVisc', 1.464e-5), 'KinVisc', '- Kinematic viscosity (m^2/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('SpdSound', 335.0), 'SpdSound', '- Speed of sound (m/s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('Patm', 103500.0), 'Patm', '- Atmospheric pressure (Pa)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('Pvap', 1700.0), 'Pvap', '- Vapour pressure (Pa)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('WtrDpth', 0.0), 'WtrDpth', '- Water depth (m)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('MSL2SWL', 0.0), 'MSL2SWL', '- Offset between still-water level and mean sea level (m)\n'))

            # Input Files
            f.write('---------------------- INPUT FILES ---------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('EDFile', 'unused') + '"', 'EDFile', '- ElastoDyn input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('BDBldFile(1)', 'unused') + '"', 'BDBldFile(1)', '- BeamDyn blade 1 input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('BDBldFile(2)', 'unused') + '"', 'BDBldFile(2)', '- BeamDyn blade 2 input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('BDBldFile(3)', 'unused') + '"', 'BDBldFile(3)', '- BeamDyn blade 3 input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('InflowFile', 'unused') + '"', 'InflowFile', '- InflowWind input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('AeroFile', 'unused') + '"', 'AeroFile', '- AeroDyn input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('ServoFile', 'unused') + '"', 'ServoFile', '- ServoDyn input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('SeaStFile', 'unused') + '"', 'SeaStFile', '- SeaState input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('HydroFile', 'unused') + '"', 'HydroFile', '- HydroDyn input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('SubFile', 'unused') + '"', 'SubFile', '- SubDyn/ExtPtfm input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('MooringFile', 'unused') + '"', 'MooringFile', '- MoorDyn/MAP input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('IceFile', 'unused') + '"', 'IceFile', '- Ice input file\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + fst.get('SoilFile', 'unused') + '"', 'SoilFile', '- Soil input file\n'))

            # Output
            f.write('---------------------- OUTPUT --------------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('SumPrint', False), 'SumPrint', '- Print summary data\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('SttsTime', 10.0), 'SttsTime', '- Screen status message interval (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('ChkptTime', 99999.9), 'ChkptTime', '- Checkpoint interval (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('DT_Out', 'default'), 'DT_Out', '- Time step for tabular output (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('TStart', 0.0), 'TStart', '- Time to begin tabular output (s)\n'))
            f.write('{:<22d} {:<11} {:}'.format(fst.get('OutFileFmt', 2), 'OutFileFmt', '- Output file format\n'))
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('TabDelim', True), 'TabDelim', '- Tab delimited output\n'))
            f.write('{:<22} {:<11} {:}'.format('"' + str(fst.get('OutFmt', 'ES10.3E2')) + '"', 'OutFmt', '- Format for text tabular output\n'))

            # Linearization
            f.write('---------------------- LINEARIZATION -------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('Linearize', False), 'Linearize', '- Linearization analysis\n'))
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('CalcSteady', False), 'CalcSteady', '- Calculate steady-state periodic operating point\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('TrimCase', 3), 'TrimCase', '- Controller parameter to trim\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('TrimTol', 0.001), 'TrimTol', '- Trim tolerance\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('TrimGain', 0.001), 'TrimGain', '- Trim gain\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('Twr_Kdmp', 0), 'Twr_Kdmp', '- Tower damping factor\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('Bld_Kdmp', 0), 'Bld_Kdmp', '- Blade damping factor\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('NLinTimes', 2), 'NLinTimes', '- Number of linearization times\n'))
            lin_times = fst.get('LinTimes', [30.0, 60.0])
            f.write('{:<22} {:<11} {:}'.format(', '.join(['{:f}'.format(t) for t in np.array(lin_times, dtype=float)]), 'LinTimes', '- Linearization times (s)\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('LinInputs', 1), 'LinInputs', '- Inputs in linearization\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('LinOutputs', 1), 'LinOutputs', '- Outputs in linearization\n'))
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('LinOutJac', False), 'LinOutJac', '- Include full Jacobians\n'))
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('LinOutMod', False), 'LinOutMod', '- Write module-level linearization\n'))

            # Visualization
            f.write('---------------------- VISUALIZATION ------------------------------------------\n')
            f.write('{:<22} {:<11} {:}'.format(fst.get('WrVTK', 0), 'WrVTK', '- VTK visualization data output\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('VTK_type', 1), 'VTK_type', '- Type of VTK visualization data\n'))
            f.write('{!s:<22} {:<11} {:}'.format(fst.get('VTK_fields', False), 'VTK_fields', '- Write mesh fields to VTK\n'))
            f.write('{:<22} {:<11} {:}'.format(fst.get('VTK_fps', 1), 'VTK_fps', '- Frame rate for VTK output\n'))

    def _read_main_input(self, fst_path: Path, base_dir: Path) -> dict:
        """Read main .fst file. Logic from FAST_reader.read_MainInput()."""
        fst = {}
        f = open(fst_path)

        # Header
        f.readline()
        # We store the description but it goes in fst_vt['description'], not fst dict
        # For now, skip — the driver's read() can set it separately
        f.readline()

        # Simulation Control
        f.readline()
        fst['Echo'] = bool_read(f.readline().split()[0])
        fst['AbortLevel'] = quoted_read(f.readline().split()[0])
        fst['TMax'] = float_read(f.readline().split()[0])
        fst['DT'] = float_read(f.readline().split()[0])
        fst['ModCoupling'] = int(f.readline().split()[0])
        fst['InterpOrder'] = int(f.readline().split()[0])
        fst['NumCrctn'] = int(f.readline().split()[0])
        fst['RhoInf'] = float_read(f.readline().split()[0])
        fst['ConvTol'] = float_read(f.readline().split()[0])
        fst['MaxConvIter'] = int(f.readline().split()[0])
        fst['DT_UJac'] = float_read(f.readline().split()[0])
        fst['UJacSclFact'] = float_read(f.readline().split()[0])

        # Feature Switches and Flags
        f.readline()
        fst['NRotors'] = int(f.readline().split()[0])
        fst['CompElast'] = int(f.readline().split()[0])
        fst['CompInflow'] = int(f.readline().split()[0])
        fst['CompAero'] = int(f.readline().split()[0])
        fst['CompServo'] = int(f.readline().split()[0])
        fst['CompSeaSt'] = int(f.readline().split()[0])
        fst['CompHydro'] = int(f.readline().split()[0])
        fst['CompSub'] = int(f.readline().split()[0])
        fst['CompMooring'] = int(f.readline().split()[0])
        fst['CompIce'] = int(f.readline().split()[0])
        fst['CompSoil'] = int(f.readline().split()[0])
        fst['MHK'] = int(f.readline().split()[0])
        fst['MirrorRotor'] = read_array(f, fst['NRotors'], array_type=bool)

        # Environmental conditions
        f.readline()
        fst['Gravity'] = float_read(f.readline().split()[0])
        fst['AirDens'] = float_read(f.readline().split()[0])
        fst['WtrDens'] = float_read(f.readline().split()[0])
        fst['KinVisc'] = float_read(f.readline().split()[0])
        fst['SpdSound'] = float_read(f.readline().split()[0])
        fst['Patm'] = float_read(f.readline().split()[0])
        fst['Pvap'] = float_read(f.readline().split()[0])
        fst['WtrDpth'] = float_read(f.readline().split()[0])
        fst['MSL2SWL'] = float_read(f.readline().split()[0])

        # Input Files
        f.readline()
        fst['EDFile'] = quoted_read(f.readline().split()[0])
        fst['BDBldFile(1)'] = quoted_read(f.readline().split()[0])
        fst['BDBldFile(2)'] = quoted_read(f.readline().split()[0])
        fst['BDBldFile(3)'] = quoted_read(f.readline().split()[0])
        fst['InflowFile'] = quoted_read(f.readline().split()[0])
        fst['AeroFile'] = quoted_read(f.readline().split()[0])
        fst['ServoFile'] = quoted_read(f.readline().split()[0])
        fst['SeaStFile'] = quoted_read(f.readline().split()[0])
        fst['HydroFile'] = quoted_read(f.readline().split()[0])
        fst['SubFile'] = quoted_read(f.readline().split()[0])
        fst['MooringFile'] = quoted_read(f.readline().split()[0])
        fst['IceFile'] = quoted_read(f.readline().split()[0])
        fst['SoilFile'] = quoted_read(f.readline().split()[0])

        # Output Parameters
        f.readline()
        fst['SumPrint'] = bool_read(f.readline().split()[0])
        fst['SttsTime'] = float_read(f.readline().split()[0])
        fst['ChkptTime'] = float_read(f.readline().split()[0])
        fst['DT_Out'] = float_read(f.readline().split()[0])
        fst['TStart'] = float_read(f.readline().split()[0])
        fst['OutFileFmt'] = int(f.readline().split()[0])
        fst['TabDelim'] = bool_read(f.readline().split()[0])
        fst['OutFmt'] = quoted_read(f.readline().split()[0])

        # Linearization
        f.readline()
        fst['Linearize'] = f.readline().split()[0]
        fst['CalcSteady'] = f.readline().split()[0]
        fst['TrimCase'] = f.readline().split()[0]
        fst['TrimTol'] = f.readline().split()[0]
        fst['TrimGain'] = f.readline().split()[0]
        fst['Twr_Kdmp'] = f.readline().split()[0]
        fst['Bld_Kdmp'] = f.readline().split()[0]
        fst['NLinTimes'] = int(f.readline().split()[0])
        fst['LinTimes'] = read_array(f, fst['NLinTimes'], array_type=float)
        fst['LinInputs'] = f.readline().split()[0]
        fst['LinOutputs'] = f.readline().split()[0]
        fst['LinOutJac'] = f.readline().split()[0]
        fst['LinOutMod'] = f.readline().split()[0]

        # Visualization
        f.readline()
        fst['WrVTK'] = int(f.readline().split()[0])
        fst['VTK_type'] = int(f.readline().split()[0])
        fst['VTK_fields'] = bool_read(f.readline().split()[0])
        fst['VTK_fps'] = float_read(f.readline().split()[0])

        f.close()

        # File paths
        fst['EDFile_path'] = os.path.split(fst['EDFile'])[0]
        fst['BDBldFile(1_path)'] = os.path.split(fst['BDBldFile(1)'])[0]
        fst['BDBldFile(2_path)'] = os.path.split(fst['BDBldFile(2)'])[0]
        fst['BDBldFile(3_path)'] = os.path.split(fst['BDBldFile(3)'])[0]
        fst['InflowFile_path'] = os.path.split(fst['InflowFile'])[0]
        fst['AeroFile_path'] = os.path.split(fst['AeroFile'])[0]
        fst['ServoFile_path'] = os.path.split(fst['ServoFile'])[0]
        fst['HydroFile_path'] = os.path.split(fst['HydroFile'])[0]
        fst['SubFile_path'] = os.path.split(fst['SubFile'])[0]
        fst['MooringFile_path'] = os.path.split(fst['MooringFile'])[0]
        fst['IceFile_path'] = os.path.split(fst['IceFile'])[0]

        return fst

    @staticmethod
    def _resolve(base_dir: Path, rel_path: str) -> Path:
        return base_dir / rel_path.strip('"').strip("'")
