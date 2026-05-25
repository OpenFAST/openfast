"""
Standalone driver tests: read, roundtrip (read → write → re-read → compare),
and key-field smoke checks.

Test case names mirror the CTestList.cmake entries for each module driver so
this file is the Python counterpart to the CMake regression suite.  Each
parametrize list intentionally matches the *same* identifiers used in:

    reg_tests/CTestList.cmake

Pattern:
  1. Read the driver input file from the r-test directory.
  2. Write to a temporary directory.
  3. Re-read from the temporary directory.
  4. Compare original and re-read dicts with compare_fst_vt.

Only the driver dict keys are compared (not simulation outputs).
External data files that are only *referenced* (e.g. csv timeseries, wind
files, binary outb files) are not copied to the temp dir because they are not
written by the driver — the comparison uses removeFileRef=True so file-name
strings are excluded.

Smoke tests (``test_*_read_smoke``) verify a handful of known scalar values
from specific r-test cases to confirm parsing correctness beyond structural
roundtrip equality.
"""
from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Locate r-test root  (same logic as conftest.py: parents[2] = openfast repo)
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parents[2]   # tests/ → openfast_io (pkg) → openfast_io (proj) → openfast
_MODULES_DIR = _REPO_ROOT / 'reg_tests' / 'r-test' / 'modules'

# ---------------------------------------------------------------------------
# Import driver classes and comparison helper
# ---------------------------------------------------------------------------
from openfast_io.drivers.aerodyn_driver import AeroDynStandaloneDriver
from openfast_io.drivers.beamdyn_driver import BeamDynStandaloneDriver
from openfast_io.drivers.hydrodyn_driver import HydroDynStandaloneDriver
from openfast_io.drivers.subdyn_driver import SubDynStandaloneDriver
from openfast_io.drivers.inflowwind_driver import InflowWindStandaloneDriver
from openfast_io.drivers.seastate_driver import SeaStateStandaloneDriver
from openfast_io.drivers.moordyn_driver import MoorDynStandaloneDriver
from openfast_io.drivers.aerodisk_driver import AeroDiskStandaloneDriver
from openfast_io.drivers.simple_elastodyn_driver import SimpleElastoDynStandaloneDriver
from openfast_io.drivers.unsteadyaero_driver import UnsteadyAeroStandaloneDriver
from openfast_io.FileTools import compare_fst_vt


# ---------------------------------------------------------------------------
# Generic roundtrip helper
# ---------------------------------------------------------------------------

def _roundtrip(driver, dvr_path: Path, case_dir: Path) -> dict:
    """Read → write to tmpdir → re-read.  Returns diff dict (empty = pass)."""

    # --- READ ---
    data_orig = driver.read(dvr_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        out_dvr = tmp_path / dvr_path.name

        # Copy external data files that the driver references but does not
        # write (timeseries csv, wind files, polar tables, motion files ...).
        # We copy *all* non-module-input files from the case dir so that the
        # re-read can resolve every reference without modification.
        for src in case_dir.iterdir():
            if src.is_file():
                shutil.copy2(src, tmp_path / src.name)
        # Copy any sub-directories (e.g. Airfoils/)
        for src in case_dir.iterdir():
            if src.is_dir():
                shutil.copytree(src, tmp_path / src.name, dirs_exist_ok=True)

        # --- WRITE ---
        driver.write(data_orig, out_dvr)

        # --- RE-READ ---
        data_reread = driver.read(out_dvr)

    # --- COMPARE ---
    diff = compare_fst_vt(
        data_orig, data_reread,
        ignoreVars=['OutRootName', 'TMax', 'TStart', 'OutFileFmt'],
        removeFileRef=True,
        removeArrayProps=True,
        print_diff=False,
    )
    return diff


# ---------------------------------------------------------------------------
# Helper: skip if r-test case dir doesn't exist
# ---------------------------------------------------------------------------

def _case_dir(module: str, case_name: str) -> Path:
    return _MODULES_DIR / module / case_name


def _dvr_file(case_dir: Path, filename: str) -> Path:
    return case_dir / filename


# ===========================================================================
# AeroDyn  (mirrors ad_regression calls in CTestList.cmake)
# ===========================================================================

# Each entry: (case_name, driver_filename)
_AD_CASES = [
    ('ad_timeseries_shutdown',         'ad_driver.dvr'),
    ('ad_BAR_SineMotion',              'ad_driver.dvr'),
    ('ad_BAR_SineMotion_UA4_DBEMT3',   'ad_driver.dvr'),
    ('ad_BAR_RNAMotion',               'ad_driver.dvr'),
    ('ad_MHK_RM1_Fixed',               'ad_driver.dvr'),
    ('ad_MHK_RM1_Fixed_IfW',           'ad_driver.dvr'),
    ('ad_MHK_RM1_Floating',            'ad_driver.dvr'),
    ('ad_MultipleHAWT',                'ad_driver.dvr'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _AD_CASES, ids=[c[0] for c in _AD_CASES])
def test_aerodyn_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('aerodyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = AeroDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'AeroDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# BeamDyn  (mirrors bd_regression calls in CTestList.cmake)
# ===========================================================================

_BD_CASES = [
    ('bd_5MW_dynamic',               'bd_driver.inp'),
    ('bd_5MW_dynamic_gravity_Az00',  'bd_driver.inp'),
    ('bd_5MW_dynamic_gravity_Az90',  'bd_driver.inp'),
    ('bd_5MW_dynamic_modal_damping', 'bd_driver.inp'),
    ('bd_curved_beam',               'bd_driver.inp'),
    ('bd_isotropic_rollup',          'bd_driver.inp'),
    ('bd_static_cantilever_beam',    'bd_driver.inp'),
    ('bd_static_twisted_with_k1',    'bd_driver.inp'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _BD_CASES, ids=[c[0] for c in _BD_CASES])
def test_beamdyn_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('beamdyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = BeamDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'BeamDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# HydroDyn  (mirrors hd_regression calls in CTestList.cmake)
# ===========================================================================

_HD_CASES = [
    ('hd_5MW_ITIBarge_DLL_WTurb_WavesIrr',       'hd_driver.inp'),
    ('hd_5MW_OC3Spar_DLL_WTurb_WavesIrr',        'hd_driver.inp'),
    ('hd_5MW_OC4Semi_WSt_WavesWN',               'hd_driver.inp'),
    ('hd_5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti', 'hd_driver.inp'),
    ('hd_NBodyMod1',                             'hd_driver.inp'),
    ('hd_NBodyMod2',                             'hd_driver.inp'),
    ('hd_NBodyMod3',                             'hd_driver.inp'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _HD_CASES)
def test_hydrodyn_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('hydrodyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = HydroDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'HydroDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# SubDyn  (mirrors sd_regression calls in CTestList.cmake)
# ===========================================================================

def _sd_dvr_name(case_name: str) -> str:
    """SubDyn driver files are named <case_name>.dvr."""
    return case_name + '.dvr'


_SD_CASES = [
    'SD_Cable_5Joints',
    'SD_PendulumDamp',
    'SD_Rigid',
    'SD_SparHanging',
    'SD_AnsysComp2_Cable',
    'SD_Spring_Case1',
    'SD_Spring_Case2',
    'SD_Spring_Case3',
    'SD_Revolute_Joint',
    'SD_2Beam_Spring',
    'SD_2Beam_Cantilever',
    'SD_CantileverBeam_Rectangular',
]


@pytest.mark.parametrize('case_name', _SD_CASES)
def test_subdyn_driver_roundtrip(case_name):
    cdir = _case_dir('subdyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, _sd_dvr_name(case_name))
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = SubDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'SubDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# InflowWind  (mirrors ifw_regression calls in CTestList.cmake)
# ===========================================================================

_IFW_CASES = [
    ('ifw_turbsimff',    'ifw_driver.inp'),
    ('ifw_uniform',      'ifw_driver.inp'),
    ('ifw_BoxExceed',    'ifw_driver.inp'),
    ('ifw_BoxExceedTwr', 'ifw_driver.inp'),
    ('ifw_HAWC',         'ifw_driver.inp'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _IFW_CASES, ids=[c[0] for c in _IFW_CASES])
def test_inflowwind_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('inflowwind', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = InflowWindStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'InflowWind roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# SeaState  (mirrors seast_regression calls in CTestList.cmake)
# ===========================================================================

_SS_CASES = [
    ('seastate_1',               'seastate_driver.inp'),
    ('seastate_wr_kin1',         'seastate_driver.inp'),
    ('seastate_CNW1',            'seastate_driver.inp'),
    ('seastate_CNW2',            'seastate_driver.inp'),
    ('seastate_WaveMod7_WaveStMod1', 'seastate_driver.inp'),
    ('seastate_WaveMod7_WaveStMod2', 'seastate_driver.inp'),
    ('seastate_WaveMod7_WaveStMod3', 'seastate_driver.inp'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _SS_CASES, ids=[c[0] for c in _SS_CASES])
def test_seastate_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('seastate', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = SeaStateStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'SeaState roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# MoorDyn  (mirrors md_regression calls in CTestList.cmake)
# ===========================================================================

_MD_CASES = [
    ('md_5MW_OC4Semi',  'md_driver.inp'),
    ('md_BodiesAndRods','md_driver.inp'),
    ('md_bodyDrag',     'md_driver.inp'),
    ('md_cable',        'md_driver.inp'),
    ('md_case2',        'md_driver.inp'),
    ('md_case5',        'md_driver.inp'),
    ('md_float',        'md_driver.inp'),
    ('md_horizontal',   'md_driver.inp'),
    ('md_no_line',      'md_driver.inp'),
    ('md_vertical',     'md_driver.inp'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _MD_CASES, ids=[c[0] for c in _MD_CASES])
def test_moordyn_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('moordyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = MoorDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'MoorDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# AeroDisk  (mirrors adsk_regression calls in CTestList.cmake)
# ===========================================================================

_ADSK_CASES = [
    ('adsk_timeseries_shutdown', 'adsk_driver.dvr'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _ADSK_CASES, ids=[c[0] for c in _ADSK_CASES])
def test_aerodisk_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('aerodisk', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = AeroDiskStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'AeroDisk roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# SimpleElastoDyn  (mirrors sed_regression calls in CTestList.cmake)
# ===========================================================================

_SED_CASES = [
    ('sed_test_HSSbrk',    'sed_driver.dvr'),
    ('sed_test_freewheel', 'sed_driver.dvr'),
]


@pytest.mark.parametrize('case_name,dvr_filename', _SED_CASES, ids=[c[0] for c in _SED_CASES])
def test_simple_elastodyn_driver_roundtrip(case_name, dvr_filename):
    cdir = _case_dir('simple-elastodyn', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')
    dvr_path = _dvr_file(cdir, dvr_filename)
    if not dvr_path.is_file():
        pytest.skip(f'driver file not found: {dvr_path}')

    driver = SimpleElastoDynStandaloneDriver()
    diff = _roundtrip(driver, dvr_path, cdir)
    assert diff == {}, f'SimpleElastoDyn roundtrip diff for {case_name}:\n{diff}'


# ===========================================================================
# UnsteadyAero  (mirrors ua_regression calls in CTestList.cmake)
# ===========================================================================

def _ua_dvr_files(case_dir: Path):
    """UnsteadyAero cases can have multiple .dvr files (UA2.dvr, UA3.dvr, ...)."""
    return sorted(case_dir.glob('UA*.dvr'))


_UA_CASES = ['ua_redfreq', 'ua_elast']


@pytest.mark.parametrize('case_name', _UA_CASES)
def test_unsteadyaero_driver_roundtrip(case_name):
    cdir = _case_dir('unsteadyaero', case_name)
    if not cdir.is_dir():
        pytest.skip(f'r-test dir not found: {cdir}')

    dvr_files = _ua_dvr_files(cdir)
    if not dvr_files:
        pytest.skip(f'No UA*.dvr files found in {cdir}')

    driver = UnsteadyAeroStandaloneDriver()
    for dvr_path in dvr_files:
        diff = _roundtrip(driver, dvr_path, cdir)
        assert diff == {}, f'UnsteadyAero roundtrip diff for {case_name}/{dvr_path.name}:\n{diff}'


# ===========================================================================
# Smoke tests — verify known scalar values from specific r-test cases.
# These confirm parsing correctness beyond structural roundtrip equality.
# ===========================================================================

class TestDriverReadSmoke:
    """Spot-check key fields from known r-test driver files."""

    def test_aerodyn_read_smoke(self):
        cdir = _case_dir('aerodyn', 'ad_BAR_SineMotion')
        if not cdir.is_dir():
            pytest.skip('ad_BAR_SineMotion r-test not available')
        dvr = AeroDynStandaloneDriver().read(cdir / 'ad_driver.dvr')['AeroDynDriver']
        assert dvr['AnalysisType'] == 1
        assert dvr['NumTurbines'] == 1
        assert abs(dvr['TMax'] - 7.0) < 0.1
        turb = dvr['Turbines'][0]
        assert turb['NumBlades'] == 3
        assert turb['BasicHAWTFormat'] is False

    def test_beamdyn_read_smoke(self):
        cdir = _case_dir('beamdyn', 'bd_5MW_dynamic')
        if not cdir.is_dir():
            pytest.skip('bd_5MW_dynamic r-test not available')
        dvr = BeamDynStandaloneDriver().read(cdir / 'bd_driver.inp')['BeamDynDriver']
        assert dvr['DynamicSolve'] is True
        assert dvr['t_final'] == 30.0
        assert abs(dvr['dt'] - 0.002) < 1e-6
        assert abs(dvr['Gy'] - (-9.8)) < 0.01
        assert len(dvr['GlbDCM']) == 3

    def test_hydrodyn_read_smoke(self):
        cdir = _case_dir('hydrodyn', 'hd_5MW_OC4Semi_WSt_WavesWN')
        if not cdir.is_dir():
            pytest.skip('hd_5MW_OC4Semi r-test not available')
        dvr = HydroDynStandaloneDriver().read(cdir / 'hd_driver.inp')['HydroDynDriver']
        assert abs(dvr['Gravity'] - 9.80665) < 0.001
        assert dvr['WtrDens'] == 1025
        assert dvr['WtrDpth'] == 200
        assert dvr['PRPInputsMod'] == 2

    def test_subdyn_read_smoke(self):
        cdir = _case_dir('subdyn', 'SD_Rigid')
        if not cdir.is_dir():
            pytest.skip('SD_Rigid r-test not available')
        dvr = SubDynStandaloneDriver().read(cdir / 'SD_Rigid.dvr')['SubDynDriver']
        assert abs(dvr['Gravity'] - 9.81) < 0.01
        assert dvr['NSteps'] == 2000
        assert dvr['TP_RefPoint_Z'] == 30.0

    def test_inflowwind_read_smoke(self):
        cdir = _case_dir('inflowwind', 'ifw_uniform')
        if not cdir.is_dir():
            pytest.skip('ifw_uniform r-test not available')
        dvr = InflowWindStandaloneDriver().read(cdir / 'ifw_driver.inp')['InflowWindDriver']
        assert dvr['IfWFileName'] == 'ifw_primary.inp'
        assert dvr['NumTSteps'] == 8
        assert dvr['DT'] == 0.1

    def test_seastate_read_smoke(self):
        cdir = _case_dir('seastate', 'seastate_1')
        if not cdir.is_dir():
            pytest.skip('seastate_1 r-test not available')
        dvr = SeaStateStandaloneDriver().read(cdir / 'seastate_driver.inp')['SeaStateDriver']
        assert abs(dvr['Gravity'] - 9.80665) < 0.001
        assert dvr['WtrDens'] == 1025
        assert dvr['NSteps'] == 801

    def test_moordyn_read_smoke(self):
        cdir = _case_dir('moordyn', 'md_5MW_OC4Semi')
        if not cdir.is_dir():
            pytest.skip('md_5MW_OC4Semi r-test not available')
        dvr = MoorDynStandaloneDriver().read(cdir / 'md_driver.inp')['MoorDynDriver']
        assert abs(dvr['Gravity'] - 9.80665) < 0.001
        assert dvr['rhoW'] == 1025.0
        assert dvr['TMax'] == 60

    def test_aerodisk_read_smoke(self):
        cdir = _case_dir('aerodisk', 'adsk_timeseries_shutdown')
        if not cdir.is_dir():
            pytest.skip('adsk_timeseries_shutdown r-test not available')
        dvr = AeroDiskStandaloneDriver().read(cdir / 'adsk_driver.dvr')['AeroDiskDriver']
        assert dvr['AirDens'] == 1.225
        assert dvr['RotorRad'] == 63.0

    def test_simple_elastodyn_read_smoke(self):
        cdir = _case_dir('simple-elastodyn', 'sed_test_freewheel')
        if not cdir.is_dir():
            pytest.skip('sed_test_freewheel r-test not available')
        dvr = SimpleElastoDynStandaloneDriver().read(cdir / 'sed_driver.dvr')['SimpleElastoDynDriver']
        assert dvr['SEDiptFile'] == 'sed_primary.inp'
        assert dvr['TimeseriesFile'] == 'Free.csv'

    def test_unsteadyaero_read_smoke(self):
        cdir = _case_dir('unsteadyaero', 'ua_redfreq')
        if not cdir.is_dir():
            pytest.skip('ua_redfreq r-test not available')
        dvr_files = sorted(cdir.glob('UA*.dvr'))
        if not dvr_files:
            pytest.skip('No UA*.dvr files found')
        dvr = UnsteadyAeroStandaloneDriver().read(dvr_files[0])['UnsteadyAeroDriver']
        assert dvr['FldDens'] == 1.225
        assert dvr['UAMod'] == 2
        assert dvr['Chord'] == 3.5
        assert len(dvr['MassMatrix']) == 3

