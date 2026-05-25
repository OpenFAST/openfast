"""Tests for AeroDynIO read/write extraction."""
import os
import pytest
from pathlib import Path

from openfast_io.io.aerodyn import AeroDynIO


# ── Fixtures ───────────────────────────────────────────────────────────────

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


# ── Unit tests ─────────────────────────────────────────────────────────────

def test_aerodyn_io_is_module_io():
    from openfast_io.io.base import ModuleIO
    assert issubclass(AeroDynIO, ModuleIO)


def test_read_5mw_aerodyn(ad_io, r_test_5mw_dir):
    """Read the 5MW AeroDyn file and verify key parameters."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    
    assert 'AeroDyn' in result
    assert 'AeroDynBlade' in result
    
    ad = result['AeroDyn']
    assert ad['Echo'] == False
    assert isinstance(ad['Wake_Mod'], int)
    assert ad['NumAFfiles'] > 0
    assert len(ad['AFNames']) == ad['NumAFfiles']
    
    # Tower data
    assert ad['NumTwrNds'] > 0
    assert len(ad['TwrElev']) == ad['NumTwrNds']
    assert len(ad['TwrDiam']) == ad['NumTwrNds']


def test_read_5mw_aerodyn_blade(ad_io, r_test_5mw_dir):
    """Read and verify blade data."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    
    ad_blade = result['AeroDynBlade']
    # 5MW has identical blades → should be collapsed to single dict
    assert isinstance(ad_blade, dict), "Identical blades should be collapsed to dict"
    assert ad_blade['NumBlNds'] > 0
    assert len(ad_blade['BlSpn']) == ad_blade['NumBlNds']
    assert len(ad_blade['BlChord']) == ad_blade['NumBlNds']


def test_read_5mw_aerodyn_polars(ad_io, r_test_5mw_dir):
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


def test_read_5mw_aerodyn_coords(ad_io, r_test_5mw_dir):
    """Verify airfoil coordinate data is read."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    
    ad = result['AeroDyn']
    assert 'af_coord' in ad
    assert 'ac' in ad
    assert len(ad['af_coord']) == ad['NumAFfiles']


def test_write_and_reread_5mw(ad_io, r_test_5mw_dir, tmp_path):
    """Write 5MW AeroDyn data to tmp, re-read, verify key values match."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    original = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    
    # Write it out
    out_file = tmp_path / 'test_AeroDyn.dat'
    ad_io.write(original, out_file, tmp_path, naming_out='test')
    
    assert out_file.exists()
    
    # Re-read
    reread = ad_io.read(str(out_file), str(tmp_path), num_blades=3, aero_file_path='')
    
    # Compare key scalar values
    for key in ('Wake_Mod', 'TwrPotent', 'TwrShadow', 'NumAFfiles', 'NumTwrNds',
                'BEM_Mod', 'DBEMT_Mod', 'AFTabMod', 'UA_Mod', 'NBlOuts', 'NTwOuts'):
        assert original['AeroDyn'][key] == reread['AeroDyn'][key], f"Mismatch on {key}"
    
    # Compare blade data
    orig_bld = original['AeroDynBlade']
    new_bld = reread['AeroDynBlade']
    if isinstance(orig_bld, dict):
        assert isinstance(new_bld, dict)
        assert orig_bld['NumBlNds'] == new_bld['NumBlNds']
    
    # Compare polar count
    assert len(original['AeroDyn']['af_data']) == len(reread['AeroDyn']['af_data'])


def test_read_returns_bem_options(ad_io, r_test_5mw_dir):
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


def test_read_nacelle_and_hub(ad_io, r_test_5mw_dir):
    """Hub and nacelle parameters should be present."""
    ad_file = os.path.join(r_test_5mw_dir, 'NRELOffshrBsline5MW_Onshore_AeroDyn.dat')
    result = ad_io.read(ad_file, r_test_5mw_dir, num_blades=3, aero_file_path='')
    ad = result['AeroDyn']
    
    assert 'VolHub' in ad
    assert 'HubCenBx' in ad
    assert 'VolNac' in ad
    assert 'NacCenB' in ad
    assert 'TFinAero' in ad
