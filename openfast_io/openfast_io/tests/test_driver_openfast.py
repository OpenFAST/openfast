"""OpenFASTDriver read/write tests.

Strengthened 2026-06-23 after a differential audit found the prior round-trip test
(4 scalars only) passed while OutList read/write was broken on 63/77 decks, ServoDyn/
HydroDyn/SeaState writers dumped the full channel registry, and BeamDyn blades collided.
The round-trip checks here now assert OutList preservation, write->reread idempotency,
and absence of private-key pollution — the properties those bugs violated.
"""
import copy
import tempfile
from pathlib import Path

import pytest

from openfast_io.drivers.openfast import OpenFASTDriver, init_fst_vt


def _case(r_test_dir, name):
    d = r_test_dir / "glue-codes" / "openfast" / name
    if not d.exists():
        pytest.skip(f"r-test case not found: {name}")
    return d / f"{name}.fst"


# ---------------------------------------------------------------------------
# Helpers — the invariants the audit-found bugs violated
# ---------------------------------------------------------------------------

def outlist_sets(fst_vt):
    """{module: frozenset(channels set True)} from fst_vt['outlist']."""
    def true_ch(d):
        out = set()
        for k, v in d.items():
            if isinstance(v, dict):
                out |= true_ch(v)
            elif v is True:
                out.add(k)
        return out
    return {m: frozenset(true_ch(v)) for m, v in fst_vt.get("outlist", {}).items()
            if isinstance(v, dict) and true_ch(v)}


def assert_no_outlist_pollution(fst_vt):
    """No module data dict should carry a private '_outlist' key (channels live only in
    fst_vt['outlist']). Regression guard for the ExtPtfm/AeroDisk/SED pollution bug."""
    for mod, val in fst_vt.items():
        if mod == "outlist":
            continue
        for d in (val if isinstance(val, list) else [val]):
            if isinstance(d, dict):
                assert "_outlist" not in d, f"fst_vt['{mod}'] leaks a private '_outlist' key"


def roundtrip(driver, fst_path, tmp_path):
    """read -> write -> reread, returning (src, reread). Asserts the invariants that the
    OutList / full-registry-emit / spd_trq / pollution bugs all broke."""
    src = driver.read(fst_path)
    driver.write(src, tmp_path, "rt")
    reread = driver.read(tmp_path / "rt.fst")  # must not raise (spd_trq regression)

    # The exact requested OutList channel set must survive per module — catches both
    # silent drops AND the full-registry explosion (e.g. ServoDyn 2 -> 518).
    assert outlist_sets(reread) == outlist_sets(src), (
        f"OutList not preserved on round-trip\n src={outlist_sets(src)}\n re ={outlist_sets(reread)}"
    )
    assert_no_outlist_pollution(src)
    assert_no_outlist_pollution(reread)
    return src, reread


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

def test_init_fst_vt_has_expected_keys():
    fst_vt = init_fst_vt()
    for key in ['Fst', 'outlist', 'ElastoDyn', 'AeroDyn', 'ServoDyn', 'HydroDyn',
                'SeaState', 'MoorDyn', 'SubDyn', 'BeamDyn', 'ElastoDynBlade',
                'AeroDynBlade', 'BeamDynBlade', 'BStC', 'NStC', 'TStC', 'SStC',
                'DISCON_in', 'spd_trq', 'AeroDynPolar', 'SimpleElastoDyn', 'description',
                'WaterKin', 'SoilDyn', 'ExtPtfm', 'MAP', 'AeroDisk']:
        assert key in fst_vt, f"fst_vt missing key: {key}"


def test_init_fst_vt_stc_are_lists():
    fst_vt = init_fst_vt()
    for key in ['BStC', 'NStC', 'TStC', 'SStC']:
        assert isinstance(fst_vt[key], list)


# ---------------------------------------------------------------------------
# Read population (per module)
# ---------------------------------------------------------------------------

def test_driver_reads_main_input(r_test_5mw_dir):
    fst_vt = OpenFASTDriver().read(r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst")
    assert fst_vt['Fst']['TMax'] == 60.0
    assert fst_vt['Fst']['CompElast'] == 1


def test_driver_populates_elastodyn(r_test_5mw_dir):
    fst_vt = OpenFASTDriver().read(r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst")
    assert fst_vt['ElastoDyn']['NumBl'] == 3
    assert isinstance(fst_vt['ElastoDynBlade'], dict), "identical blades collapse to a dict"
    assert fst_vt['ElastoDynTower'] and 'NTwInpSt' in fst_vt['ElastoDynTower']


def test_driver_populates_inflowwind_aerodyn_servodyn(r_test_5mw_dir):
    fst_vt = OpenFASTDriver().read(r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst")
    assert 'WindType' in fst_vt['InflowWind']
    assert 'Wake_Mod' in fst_vt['AeroDyn']
    assert fst_vt['AeroDynBlade']
    assert 'VSContrl' in fst_vt['ServoDyn']
    assert isinstance(fst_vt['BStC'], list) and isinstance(fst_vt['NStC'], list)


def test_driver_populates_hydrodynamics(r_test_dir):
    fst_vt = OpenFASTDriver().read(_case(r_test_dir, "StC_test_OC4Semi"))
    hd = fst_vt['HydroDyn']
    assert hd.get('NBody', 0) >= 1 and hd.get('NJoints', 0) > 0 and hd.get('NMembers', 0) > 0
    assert isinstance(fst_vt['SeaState'].get('WaveMod'), int)


def test_driver_populates_subdyn(r_test_dir):
    sd = OpenFASTDriver().read(_case(r_test_dir, "5MW_OC3Mnpl_DLL_WTurb_WavesIrr"))['SubDyn']
    assert sd.get('NJoints', 0) > 0 and sd.get('NMembers', 0) > 0


def test_driver_populates_moordyn(r_test_dir):
    md = OpenFASTDriver().read(_case(r_test_dir, "5MW_OC4Semi_WSt_WavesWN"))['MoorDyn']
    assert 'Name' in md and len(md['Name']) >= 1


def test_driver_populates_map(r_test_dir):
    m = OpenFASTDriver().read(_case(r_test_dir, "5MW_OC3Spar_DLL_WTurb_WavesIrr"))['MAP']
    assert 'LineType' in m and len(m['LineType']) >= 1
    assert 'Node' in m and len(m['Node']) > 0


# ---------------------------------------------------------------------------
# OutList capture — regression guard for the "ElastoDyn-only" read bug
# ---------------------------------------------------------------------------

def test_outlist_captured_for_every_active_module(r_test_5mw_dir):
    """The read bug dropped every module's OutList except ElastoDyn. Assert the deck's
    actually-requested channels are captured for InflowWind and ServoDyn too."""
    sets = outlist_sets(OpenFASTDriver().read(r_test_5mw_dir / "5MW_Land_DLL_WTurb.fst"))
    assert len(sets.get('ElastoDyn', ())) > 10
    assert len(sets.get('InflowWind', ())) >= 1, "InflowWind OutList dropped on read"
    assert len(sets.get('ServoDyn', ())) >= 1, "ServoDyn OutList dropped on read"


# ---------------------------------------------------------------------------
# BeamDyn — collapse contract + the HIGH-severity distinct-blade collision
# ---------------------------------------------------------------------------

def test_beamdyn_identical_blades_collapse_to_dict(r_test_dir):
    """5MW_Land_BD has 3 identical blade files -> fst_vt['BeamDyn'] must be a single dict
    (the fst_vt list->dict contract WEIS/WISDEM rely on), mirroring ElastoDynBlade."""
    fst_vt = OpenFASTDriver().read(_case(r_test_dir, "5MW_Land_BD_DLL_WTurb"))
    assert isinstance(fst_vt['BeamDyn'], dict), "identical BeamDyn blades must collapse to a dict"
    assert 'member_total' in fst_vt['BeamDyn']
    assert isinstance(fst_vt['BeamDynBlade'], dict)


def test_beamdyn_distinct_blades_do_not_collide(r_test_dir):
    """HIGH-severity regression: the writer must give each distinct blade its own file.
    Previously every blade was written to the same path, so blades 1&2 silently inherited
    blade 3's properties. The corpus has no distinct-blade deck, so synthesize one."""
    driver = OpenFASTDriver()
    fst_vt = driver.read(_case(r_test_dir, "5MW_Land_BD_DLL_WTurb"))
    bd0 = fst_vt['BeamDyn'] if isinstance(fst_vt['BeamDyn'], dict) else fst_vt['BeamDyn'][0]
    bl0 = fst_vt['BeamDynBlade'] if isinstance(fst_vt['BeamDynBlade'], dict) else fst_vt['BeamDynBlade'][0]

    def mark(val):
        b = copy.deepcopy(bl0)
        b['beam_stiff'][0][0][0] = val  # (nstation,6,6) ndarray — a round-trippable field
        return b

    fst_vt['BeamDyn'] = [copy.deepcopy(bd0) for _ in range(3)]
    fst_vt['BeamDynBlade'] = [mark(1000.0), mark(2000.0), mark(3000.0)]
    for i in range(3):
        fst_vt['Fst'][f'BDBldFile({i + 1})'] = f'blade_{i + 1}.dat'

    out = Path(tempfile.mkdtemp(prefix="bd_distinct_"))
    driver.write(fst_vt, out, "case")
    blades = driver.read(out / "case.fst")['BeamDynBlade']
    assert not isinstance(blades, dict), "3 distinct blades must not collapse to one dict"
    markers = sorted(b['beam_stiff'][0][0][0] for b in blades)
    assert markers == [1000.0, 2000.0, 3000.0], f"blade collision — got {markers}"


# ---------------------------------------------------------------------------
# Round-trip — OutList preservation + idempotency across deck types.
# Each deck exercises a different bug the weak 4-scalar test missed.
# ---------------------------------------------------------------------------

ROUNDTRIP_CASES = [
    "5MW_Land_DLL_WTurb",            # land: ElastoDyn/AeroDyn/ServoDyn/InflowWind
    "5MW_Land_BD_DLL_WTurb",         # active BeamDyn
    "5MW_Land_DLL_WTurb_ADsk",       # AeroDisk private-_outlist path
    "5MW_Land_DLL_WTurb_SED",        # SimpleElastoDyn private-_outlist path
    "5MW_OC4Semi_WSt_WavesWN",       # HydroDyn/SeaState/MoorDyn — full-registry-emit bug
    "5MW_OC4Jckt_ExtPtfm",           # ExtPtfm pollution path
    "SWRT_YFree_VS_WTurb",           # VSContrl=3 -> spd_trq.dat round-trip (R2)
]


@pytest.mark.parametrize("name", ROUNDTRIP_CASES)
def test_driver_roundtrip_preserves_outlist_and_scalars(r_test_dir, name, tmp_path):
    src, reread = roundtrip(OpenFASTDriver(), _case(r_test_dir, name), tmp_path)
    for k in ('TMax', 'CompElast', 'CompAero'):
        assert reread['Fst'][k] == src['Fst'][k]


def test_servodyn_outlist_does_not_explode_on_write(r_test_dir, tmp_path):
    """Direct guard for the full-registry-emit bug: ServoDyn round-tripped to 518 channels
    (the whole registry) instead of the ~2 the deck requests."""
    src, reread = roundtrip(OpenFASTDriver(), _case(r_test_dir, "5MW_OC4Semi_WSt_WavesWN"), tmp_path)
    for mod in ('ServoDyn', 'HydroDyn', 'SeaState'):
        n_src, n_re = len(outlist_sets(src).get(mod, ())), len(outlist_sets(reread).get(mod, ()))
        assert n_re == n_src, f"{mod} OutList exploded on write: {n_src} -> {n_re}"
