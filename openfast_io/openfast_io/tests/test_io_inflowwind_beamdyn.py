"""Tests for InflowWindIO and BeamDynIO extraction."""
import os
import pytest
from pathlib import Path

from openfast_io.io.inflowwind import InflowWindIO
from openfast_io.io.beamdyn import BeamDynIO


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


# ── InflowWind Tests ──────────────────────────────────────────────────────

def test_inflowwind_is_module_io():
    from openfast_io.io.base import ModuleIO
    assert issubclass(InflowWindIO, ModuleIO)


def test_read_inflowwind(baseline_dir):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_12mps.dat')
    result = io.read(ifw_file, baseline_dir)
    
    assert 'InflowWind' in result
    ifw = result['InflowWind']
    
    assert 'WindType' in ifw
    assert isinstance(ifw['WindType'], int)
    assert 'HWindSpeed' in ifw
    assert 'RefHt' in ifw
    assert 'PLExp' in ifw
    assert 'SensorType' in ifw
    assert 'SumPrint' in ifw


def test_read_inflowwind_steady(baseline_dir):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_Steady8mps.dat')
    result = io.read(ifw_file, baseline_dir)
    ifw = result['InflowWind']
    assert ifw['WindType'] == 1
    assert ifw['HWindSpeed'] == 8.0


def test_write_reread_inflowwind(baseline_dir, tmp_path):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_Steady8mps.dat')
    original = io.read(ifw_file, baseline_dir)
    
    out_file = tmp_path / 'test_InflowWind.dat'
    io.write(original, out_file, tmp_path)
    
    assert out_file.exists()
    reread = io.read(str(out_file), str(tmp_path))
    
    for key in ('WindType', 'HWindSpeed', 'RefHt', 'PLExp', 'SensorType'):
        assert original['InflowWind'][key] == reread['InflowWind'][key], f"Mismatch: {key}"


def test_inflowwind_hawc_params(baseline_dir):
    io = InflowWindIO()
    ifw_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_InflowWind_12mps.dat')
    result = io.read(ifw_file, baseline_dir)
    ifw = result['InflowWind']
    
    # HAWC params should be present even if not used
    assert 'ScaleMethod' in ifw
    assert 'URef' in ifw
    assert 'WindProfile' in ifw


# ── BeamDyn Tests ─────────────────────────────────────────────────────────

def test_beamdyn_is_module_io():
    from openfast_io.io.base import ModuleIO
    assert issubclass(BeamDynIO, ModuleIO)


def test_read_beamdyn(baseline_dir):
    io = BeamDynIO()
    bd_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_BeamDyn.dat')
    result = io.read(bd_file, baseline_dir)
    
    assert 'BeamDyn' in result
    assert 'BeamDynBlade' in result
    
    bd = result['BeamDyn']
    assert 'member_total' in bd
    assert bd['member_total'] > 0
    assert 'members' in bd
    assert len(bd['members']) == bd['member_total']
    assert 'order_elem' in bd
    assert 'BldFile' in bd


def test_read_beamdyn_blade(baseline_dir):
    io = BeamDynIO()
    bd_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_BeamDyn.dat')
    result = io.read(bd_file, baseline_dir)
    
    bld = result['BeamDynBlade']
    assert bld['station_total'] > 0
    assert len(bld['radial_stations']) == bld['station_total']
    assert bld['beam_stiff'].shape == (bld['station_total'], 6, 6)
    assert bld['beam_inertia'].shape == (bld['station_total'], 6, 6)


def test_write_reread_beamdyn(baseline_dir, tmp_path):
    io = BeamDynIO()
    bd_file = os.path.join(baseline_dir, 'NRELOffshrBsline5MW_BeamDyn.dat')
    original = io.read(bd_file, baseline_dir)

    # Set the BldFile name for output
    original['BeamDyn']['BldFile'] = 'test_BD_Blade.dat'
    
    out_file = tmp_path / 'test_BeamDyn.dat'
    io.write(original, out_file, tmp_path)
    
    assert out_file.exists()
    assert (tmp_path / 'test_BD_Blade.dat').exists()
    
    reread = io.read(str(out_file), str(tmp_path))
    
    assert original['BeamDyn']['member_total'] == reread['BeamDyn']['member_total']
    assert original['BeamDyn']['kp_total'] == reread['BeamDyn']['kp_total']
    assert original['BeamDynBlade']['station_total'] == reread['BeamDynBlade']['station_total']


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
