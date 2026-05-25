from pathlib import Path
import pytest
from openfast_io.drivers.openfast import OpenFASTDriver, init_fst_vt


def test_init_fst_vt_has_expected_keys():
    """fst_vt structure contains all required keys."""
    fst_vt = init_fst_vt()
    for key in ['Fst', 'outlist', 'ElastoDyn', 'AeroDyn', 'ServoDyn',
                'HydroDyn', 'SeaState', 'MoorDyn', 'SubDyn', 'BeamDyn',
                'ElastoDynBlade', 'AeroDynBlade', 'BeamDynBlade',
                'BStC', 'NStC', 'TStC', 'SStC', 'DISCON_in', 'spd_trq',
                'AeroDynPolar', 'SimpleElastoDyn', 'description',
                'WaterKin', 'SoilDyn',
                'ExtPtfm', 'MAP', 'AeroDisk']:
        assert key in fst_vt, f"fst_vt missing key: {key}"


def test_init_fst_vt_stc_are_lists():
    fst_vt = init_fst_vt()
    for key in ['BStC', 'NStC', 'TStC', 'SStC']:
        assert isinstance(fst_vt[key], list)


def test_driver_reads_main_input(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['Fst']['TMax'] == 60.0
    assert fst_vt['Fst']['CompElast'] == 1


def test_driver_populates_elastodyn(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['ElastoDyn']
    assert fst_vt['ElastoDyn']['NumBl'] == 3


def test_elastodyn_blade_collapsed_when_identical(r_test_5mw_dir):
    """When all blade files are the same, ElastoDynBlade is collapsed to a single dict."""
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    # Should be a single dict (collapsed), not a list
    assert isinstance(fst_vt['ElastoDynBlade'], dict)


def test_tower_data_populated(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['ElastoDynTower']
    assert 'NTwInpSt' in fst_vt['ElastoDynTower']


# --- InflowWind integration tests ---

def test_driver_populates_inflowwind(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['InflowWind'], "InflowWind should be populated when CompInflow=1"
    assert 'WindType' in fst_vt['InflowWind']


# --- AeroDyn integration tests ---

def test_driver_populates_aerodyn(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['AeroDyn'], "AeroDyn should be populated when CompAero=2"
    assert 'Wake_Mod' in fst_vt['AeroDyn']


def test_driver_populates_aerodyn_blade(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    # AeroDynBlade can be a single dict (collapsed) or list
    assert fst_vt['AeroDynBlade']


# --- ServoDyn integration tests ---

def test_driver_populates_servodyn(r_test_5mw_dir):
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert fst_vt['ServoDyn'], "ServoDyn should be populated when CompServo=1"
    assert 'PCMode' in fst_vt['ServoDyn']
    assert 'VSContrl' in fst_vt['ServoDyn']


def test_driver_stc_lists(r_test_5mw_dir):
    """5MW_Land has no StC files, so the lists should be empty."""
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    assert isinstance(fst_vt['BStC'], list)
    assert isinstance(fst_vt['NStC'], list)


# --- BeamDyn integration tests ---

def test_driver_populates_beamdyn(r_test_dir):
    """5MW_Land_BD_DLL_WTurb has CompElast=2 (BeamDyn enabled)."""
    bd_dir = r_test_dir / "glue-codes" / "openfast" / "5MW_Land_BD_DLL_WTurb"
    if not bd_dir.exists():
        pytest.skip("5MW_Land_BD_DLL_WTurb r-test case not found")
    driver = OpenFASTDriver()
    fst_path = bd_dir / "5MW_Land_BD_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)
    # BeamDyn should be a list of per-blade dicts
    assert isinstance(fst_vt['BeamDyn'], list)
    assert len(fst_vt['BeamDyn']) > 0
    # First blade should have content
    if fst_vt['BeamDyn'][0]:
        assert 'member_total' in fst_vt['BeamDyn'][0]


# ---------------------------------------------------------------------------
# HydroDyn / SeaState via OC4 Semi case
# ---------------------------------------------------------------------------

def test_driver_populates_hydrodynamics(r_test_dir):
    """StC_test_OC4Semi has CompHydro=1 (HydroDyn enabled)."""
    oc4_dir = r_test_dir / "glue-codes" / "openfast" / "StC_test_OC4Semi"
    if not oc4_dir.exists():
        pytest.skip("StC_test_OC4Semi r-test case not found")
    driver = OpenFASTDriver()
    fst_path = oc4_dir / "StC_test_OC4Semi.fst"
    fst_vt = driver.read(fst_path)
    hd = fst_vt['HydroDyn']
    assert isinstance(hd, dict)
    assert hd.get('NBody', 0) >= 1
    assert hd.get('NJoints', 0) > 0
    assert hd.get('NMembers', 0) > 0


def test_driver_populates_seastate(r_test_dir):
    """StC_test_OC4Semi has CompSeaSt=1 (SeaState enabled)."""
    oc4_dir = r_test_dir / "glue-codes" / "openfast" / "StC_test_OC4Semi"
    if not oc4_dir.exists():
        pytest.skip("StC_test_OC4Semi r-test case not found")
    driver = OpenFASTDriver()
    fst_path = oc4_dir / "StC_test_OC4Semi.fst"
    fst_vt = driver.read(fst_path)
    ss = fst_vt['SeaState']
    assert isinstance(ss, dict)
    assert isinstance(ss.get('WaveMod'), int)


# ---------------------------------------------------------------------------
# SubDyn via OC3 Monopile case (CompSub=1)
# ---------------------------------------------------------------------------

def test_driver_populates_subdyn(r_test_dir):
    """5MW_OC3Mnpl_DLL_WTurb_WavesIrr has CompSub=1 (SubDyn enabled)."""
    oc3_dir = r_test_dir / "glue-codes" / "openfast" / "5MW_OC3Mnpl_DLL_WTurb_WavesIrr"
    if not oc3_dir.exists():
        pytest.skip("5MW_OC3Mnpl r-test case not found")
    driver = OpenFASTDriver()
    fst_path = oc3_dir / "5MW_OC3Mnpl_DLL_WTurb_WavesIrr.fst"
    fst_vt = driver.read(fst_path)
    sd = fst_vt['SubDyn']
    assert isinstance(sd, dict)
    assert sd.get('NJoints', 0) > 0
    assert sd.get('NMembers', 0) > 0


# ---------------------------------------------------------------------------
# MoorDyn via OC4 Semi Wave case (CompMooring=3)
# ---------------------------------------------------------------------------

def test_driver_populates_moordyn(r_test_dir):
    """5MW_OC4Semi_WSt_WavesWN has CompMooring=3 (MoorDyn enabled)."""
    oc4w_dir = r_test_dir / "glue-codes" / "openfast" / "5MW_OC4Semi_WSt_WavesWN"
    if not oc4w_dir.exists():
        pytest.skip("5MW_OC4Semi_WSt_WavesWN r-test case not found")
    driver = OpenFASTDriver()
    fst_path = oc4w_dir / "5MW_OC4Semi_WSt_WavesWN.fst"
    fst_vt = driver.read(fst_path)
    md = fst_vt['MoorDyn']
    assert isinstance(md, dict)
    assert 'Name' in md
    assert len(md['Name']) >= 1


# ---------------------------------------------------------------------------
# MAP via OC3 Spar case (CompMooring=1)
# ---------------------------------------------------------------------------

def test_driver_populates_map(r_test_dir):
    """5MW_OC3Spar_DLL_WTurb_WavesIrr has CompMooring=1 (MAP++ enabled)."""
    spar_dir = r_test_dir / "glue-codes" / "openfast" / "5MW_OC3Spar_DLL_WTurb_WavesIrr"
    if not spar_dir.exists():
        pytest.skip("5MW_OC3Spar_DLL_WTurb_WavesIrr r-test case not found")
    driver = OpenFASTDriver()
    fst_path = spar_dir / "5MW_OC3Spar_DLL_WTurb_WavesIrr.fst"
    fst_vt = driver.read(fst_path)
    m = fst_vt['MAP']
    assert isinstance(m, dict)
    assert 'LineType' in m
    assert len(m['LineType']) >= 1
    assert 'Node' in m
    assert len(m['Node']) > 0


# ---------------------------------------------------------------------------
# Driver write roundtrip
# ---------------------------------------------------------------------------

def test_driver_write_roundtrip(r_test_5mw_dir, tmp_path):
    """Read a 5MW case, write it, re-read. Key scalars should survive."""
    driver = OpenFASTDriver()
    fst_path = r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"
    fst_vt = driver.read(fst_path)

    written = driver.write(fst_vt, tmp_path, "roundtrip_test")
    assert len(written) > 0

    # Re-read the written .fst
    new_fst = tmp_path / "roundtrip_test.fst"
    assert new_fst.exists()
    fst_vt2 = driver.read(new_fst)

    # Check key scalars survived
    assert fst_vt2['Fst']['TMax'] == fst_vt['Fst']['TMax']
    assert fst_vt2['Fst']['CompElast'] == fst_vt['Fst']['CompElast']
    assert fst_vt2['Fst']['CompAero'] == fst_vt['Fst']['CompAero']
    assert fst_vt2['ElastoDyn']['NumBl'] == fst_vt['ElastoDyn']['NumBl']
