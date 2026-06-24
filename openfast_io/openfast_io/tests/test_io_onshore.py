"""Edge-case IO tests for onshore modules.

Consolidated from the former per-module test files. Read/write/OutList
round-trip integration is covered by ``test_driver_openfast.py`` and the
external differential harness; only unique edge cases survive here.

Modules: ElastoDyn, SimpleElastoDyn, BeamDyn, AeroDyn, AeroDisk,
InflowWind, ServoDyn.
"""
import os
import tempfile

import pytest
from pathlib import Path

from openfast_io.io.aerodyn import AeroDynIO
from openfast_io.io.elastodyn import ElastoDynIO
from openfast_io.io.simple_elastodyn import SimpleElastoDynIO
from openfast_io.io.beamdyn import BeamDynIO
from openfast_io.io.inflowwind import InflowWindIO
from openfast_io.io.aerodisk import AeroDiskIO
from openfast_io.io.servodyn import ServoDynIO


# ── Fixtures / paths ─────────────────────────────────────────────────────────

R_TEST_5MW = os.path.join(
    os.path.dirname(__file__), '..', '..', '..',
    'reg_tests', 'r-test', 'glue-codes', 'openfast', '5MW_Land_DLL_WTurb'
)


@pytest.fixture
def r_test_5mw_dir():
    d = os.path.normpath(R_TEST_5MW)
    if not os.path.isdir(d):
        pytest.skip('r-test 5MW data not found')
    return d


@pytest.fixture
def ad_io():
    return AeroDynIO()


R_TEST_BASE = os.path.normpath(os.path.join(
    os.path.dirname(__file__), '..', '..', '..',
    'reg_tests', 'r-test', 'glue-codes', 'openfast'
))


@pytest.fixture
def baseline_dir():
    d = os.path.join(R_TEST_BASE, '5MW_Baseline')
    if not os.path.isdir(d):
        pytest.skip('r-test 5MW_Baseline data not found')
    return d


# AeroDisk / SimpleElastoDyn r-test paths
_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                      'glue-codes', 'openfast')
_ADSK_SED_DIR = os.path.join(_RTEST, '5MW_Land_DLL_WTurb_ADsk_SED')
_SED_FILE = os.path.join(_ADSK_SED_DIR, 'NREL_5MW_Simplified-ElastoDyn.dat')
_ADSK_FILE = os.path.join(_ADSK_SED_DIR, 'NRELOffshrBsline5MW_Onshore_AeroDisk.dat')

_has_sed = os.path.isfile(_SED_FILE)
_has_adsk = os.path.isfile(_ADSK_FILE)
skipsed = pytest.mark.skipif(not _has_sed, reason='SimpleElastoDyn r-test data missing')
skipadsk = pytest.mark.skipif(not _has_adsk, reason='AeroDisk r-test data missing')

# ServoDyn r-test paths
_5MW_DIR = os.path.join(_RTEST, '5MW_Land_DLL_WTurb')
_5MW_SD = os.path.join(_5MW_DIR, 'NRELOffshrBsline5MW_Onshore_ServoDyn.dat')
_STC_DIR = os.path.join(_RTEST, 'StC_test_OC4Semi')
_STC_SD = os.path.join(_STC_DIR, 'ServoDyn_with_StC.dat')
_HAVE_5MW = os.path.isfile(_5MW_SD)
_HAVE_STC = os.path.isfile(_STC_SD)


# ── AeroDyn edge cases ───────────────────────────────────────────────────────

def test_aerodyn_polars(ad_io, r_test_5mw_dir):
    """Verify airfoil polar data is read."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')

    ad = result['AeroDyn']
    assert 'af_data' in ad
    assert len(ad['af_data']) == ad['NumAFfiles']

    # Each af_data entry should be a list of tab dicts
    for afi, tabs in enumerate(ad['af_data']):
        assert isinstance(tabs, list)
        assert len(tabs) > 0
        polar = tabs[0]
        assert 'Alpha' in polar
        assert 'Cl' in polar
        assert 'Cd' in polar
        assert polar['NumAlf'] > 0


def test_aerodyn_coords(ad_io, r_test_5mw_dir):
    """Verify airfoil coordinate data is read."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')

    ad = result['AeroDyn']
    assert 'af_coord' in ad
    assert 'ac' in ad
    assert len(ad['af_coord']) == ad['NumAFfiles']


def test_aerodyn_bem_options(ad_io, r_test_5mw_dir):
    """BEM-specific parameters should be present."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    ad = result['AeroDyn']

    # BEM params
    assert 'Skew_Mod' in ad
    assert 'TipLoss' in ad
    assert 'HubLoss' in ad
    assert 'TanInd' in ad
    assert 'MaxIter' in ad
    assert 'SectAvg' in ad


def test_aerodyn_nacelle_and_hub(ad_io, r_test_5mw_dir):
    """Hub and nacelle parameters should be present."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    ad = result['AeroDyn']

    assert 'VolHub' in ad
    assert 'HubCenBx' in ad
    assert 'VolNac' in ad
    assert 'NacCenB' in ad
    assert 'TFinAero' in ad


# ── ElastoDyn edge cases ─────────────────────────────────────────────────────

def test_elastodyn_blade_data(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    blades = result['ElastoDynBlade']
    assert isinstance(blades, list)
    assert len(blades) == 3
    # All 3 blades reference the same file, so all should have data
    assert blades[0]['NBlInpSt'] == 6
    assert len(blades[0]['BlFract']) == 6


def test_elastodyn_tower_data(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    tower = result['ElastoDynTower']
    assert tower['NTwInpSt'] == 3
    assert len(tower['HtFract']) == 3


# ── BeamDyn edge cases ───────────────────────────────────────────────────────

def test_beamdyn_blade(baseline_dir):
    io = BeamDynIO()
    bd_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_BeamDyn.dat')
    result = io.read(bd_file, baseline_dir)

    bld = result['BeamDynBlade']
    assert bld['station_total'] > 0
    assert len(bld['radial_stations']) == bld['station_total']
    assert bld['beam_stiff'].shape == (bld['station_total'], 6, 6)
    assert bld['beam_inertia'].shape == (bld['station_total'], 6, 6)


def test_beamdyn_member_geometry(baseline_dir):
    io = BeamDynIO()
    bd_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_BeamDyn.dat')
    result = io.read(bd_file, baseline_dir)

    bd = result['BeamDyn']
    for mem in bd['members']:
        assert 'kp_xr' in mem
        assert 'kp_yr' in mem
        assert 'kp_zr' in mem
        assert 'initial_twist' in mem
        assert len(mem['kp_xr']) == len(mem['kp_yr'])


# ── InflowWind edge cases ────────────────────────────────────────────────────

def test_inflowwind_steady(baseline_dir):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_Steady8mps.dat')
    result = io.read(ifw_file, baseline_dir)
    ifw = result['InflowWind']
    assert ifw['WindType'] == 1
    assert ifw['HWindSpeed'] == 8.0


def test_inflowwind_hawc_params(baseline_dir):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_12mps.dat')
    result = io.read(ifw_file, baseline_dir)
    ifw = result['InflowWind']

    # HAWC params should be present even if not used
    assert 'ScaleMethod' in ifw
    assert 'URef' in ifw
    assert 'WindProfile' in ifw


# ── AeroDisk edge cases ──────────────────────────────────────────────────────

@skipadsk
def test_aerodisk_disk_table():
    io = AeroDiskIO()
    ad = io.read(_ADSK_FILE)['AeroDisk']
    assert 'actuatorDiskTable' in ad
    tbl = ad['actuatorDiskTable']
    assert 'attr' in tbl
    assert 'data' in tbl
    assert len(tbl['data']) > 0


# ── SimpleElastoDyn edge cases ───────────────────────────────────────────────

@skipsed
def test_sed_outlist():
    io = SimpleElastoDynIO()
    sed = io.read(_SED_FILE)['SimpleElastoDyn']
    assert '_outlist' in sed
    assert len(sed['_outlist']) > 0
    assert 'BlPitch1' in sed['_outlist']


# ── ServoDyn edge cases ──────────────────────────────────────────────────────

@pytest.mark.skipif(not _HAVE_5MW, reason='r-test data not available')
class TestServoDyn5MW:

    def test_servodyn_generator_torque(self):
        io = ServoDynIO()
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'VSContrl' in sd
        assert 'GenEff' in sd

    def test_servodyn_bladed_interface(self):
        io = ServoDynIO()
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'DLL_FileName' in sd
        assert 'DLL_NumTrq' in sd


@pytest.mark.skipif(not _HAVE_STC, reason='r-test StC data not available')
class TestServoDynStC:

    def test_servodyn_stc_counts(self):
        io = ServoDynIO()
        result = io.read(_STC_SD, base_dir=_STC_DIR,
                         servo_file_rel='ServoDyn_with_StC.dat')
        sd = result['ServoDyn']
        # This test case should have StC files defined
        # Verify counts match the loaded lists
        assert len(result['BStC']) == sd['NumBStC']
        assert len(result['NStC']) == sd['NumNStC']
        assert len(result['TStC']) == sd['NumTStC']
        assert len(result['SStC']) == sd['NumSStC']

    def test_servodyn_stc_fields(self):
        io = ServoDynIO()
        result = io.read(_STC_SD, base_dir=_STC_DIR,
                         servo_file_rel='ServoDyn_with_StC.dat')
        # Check at least one StC list has entries and fields are present
        all_stc = result['BStC'] + result['NStC'] + result['TStC'] + result['SStC']
        assert len(all_stc) > 0, 'Expected at least one StC file in test case'
        stc = all_stc[0]
        assert 'StC_DOF_MODE' in stc
        assert 'StC_X_M' in stc
        assert 'SpringForceTable' in stc


class TestServoDynWrite:

    @staticmethod
    def _make_minimal_sd() -> dict:
        """Construct a minimal ServoDyn data dict for writing."""
        sd = {
            'Echo': False, 'DT': 0.005,
            'PCMode': 0, 'TPCOn': 0.0,
        }
        for idx in range(1, 4):
            sd[f'PitNeut({idx})'] = 0.0
            sd[f'PitSpr({idx})'] = 0.0
            sd[f'PitDamp({idx})'] = 0.0
            sd[f'TPitManS({idx})'] = 9999.9
            sd[f'PitManRat({idx})'] = 2.0
            sd[f'BlPitchF({idx})'] = 0.0
        sd.update({
            'VSContrl': 5, 'GenModel': 1, 'GenEff': 94.4,
            'GenTiStr': True, 'GenTiStp': True,
            'SpdGenOn': 0.0, 'TimGenOn': 0.0, 'TimGenOf': 9999.9,
            'VS_RtGnSp': 0.0, 'VS_RtTq': 0.0, 'VS_Rgn2K': 0.0, 'VS_SlPc': 0.0,
            'SIG_SlPc': 0.0, 'SIG_SySp': 0.0, 'SIG_RtTq': 0.0, 'SIG_PORt': 0.0,
            'TEC_Freq': 0.0, 'TEC_NPol': 0, 'TEC_SRes': 0.0, 'TEC_RRes': 0.0,
            'TEC_VLL': 0.0, 'TEC_SLR': 0.0, 'TEC_RLR': 0.0, 'TEC_MR': 0.0,
            'HSSBrMode': 0, 'THSSBrDp': 0.0, 'HSSBrDT': 0.0, 'HSSBrTqF': 0.0,
            'YCMode': 0, 'TYCOn': 0.0, 'YawNeut': 0.0, 'YawSpr': 0.0,
            'YawDamp': 0.0, 'TYawManS': 0.0, 'YawManRat': 0.0, 'NacYawF': 0.0,
            'AfCmode': 0, 'AfC_Mean': 0.0, 'AfC_Amp': 0.0, 'AfC_Phase': 0.0,
            'NumBStC': 0, 'BStCfiles': [],
            'NumNStC': 0, 'NStCfiles': [],
            'NumTStC': 0, 'TStCfiles': [],
            'NumSStC': 0, 'SStCfiles': [],
            'CCmode': 0,
            'DLL_FileName': 'libdiscon.so', 'DLL_InFile': 'DISCON.IN',
            'DLL_ProcName': 'DISCON', 'DLL_DT': 'default',
            'DLL_Ramp': False, 'BPCutoff': 0.0, 'NacYaw_North': 0.0,
            'Ptch_Cntrl': 1, 'Ptch_SetPnt': 0.0, 'Ptch_Min': 0.0,
            'Ptch_Max': 90.0, 'PtchRate_Min': -8.0, 'PtchRate_Max': 8.0,
            'Gain_OM': 0.0, 'GenSpd_MinOM': 0.0, 'GenSpd_MaxOM': 0.0,
            'GenSpd_Dem': 0.0, 'GenTrq_Dem': 0.0, 'GenPwr_Dem': 0.0,
            'DLL_NumTrq': 0, 'GenSpd_TLU': [], 'GenTrq_TLU': [],
            'SumPrint': False, 'OutFile': 1, 'TabDelim': True,
            'OutFmt': 'ES10.3E2', 'TStart': 0.0,
        })
        return {'ServoDyn': sd, 'BStC': [], 'NStC': [], 'TStC': [], 'SStC': []}

    def test_servodyn_write_reread_synthetic(self):
        io = ServoDynIO()
        data = self._make_minimal_sd()
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, 'ServoDyn.dat')
            io.write(data, out)
            result = io.read(out, base_dir=tmp)
        sd_orig = data['ServoDyn']
        sd_read = result['ServoDyn']
        assert sd_read['PCMode'] == sd_orig['PCMode']
        assert sd_read['VSContrl'] == sd_orig['VSContrl']
        assert sd_read['GenEff'] == pytest.approx(sd_orig['GenEff'])
