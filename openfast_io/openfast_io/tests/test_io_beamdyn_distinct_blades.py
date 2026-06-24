"""Regression test for the BeamDyn distinct-blade write collision.

The r-test corpus has no deck with 3 *distinct* BeamDyn blades, so the historical
collision — the driver never assigned a per-blade ``BldFile``, sending every blade to
the same path so blades 1&2 silently inherited blade 3's properties — was invisible to
the integration suite. This test synthesizes the distinct-blade case from a real BeamDyn
deck and asserts the writer emits 3 distinct blade files that round-trip to distinct data.

Pure I/O (no OpenFAST binary). Skips if the r-test submodule is not checked out.
"""
import copy
import tempfile
from pathlib import Path

import pytest

from openfast_io.drivers.openfast import OpenFASTDriver

# repo root = .../openfast_io/openfast_io/tests/ -> parents[3]
REPO = Path(__file__).resolve().parents[3]
DECK = REPO / "reg_tests/r-test/glue-codes/openfast/5MW_Land_BD_DLL_WTurb/5MW_Land_BD_DLL_WTurb.fst"

pytestmark = pytest.mark.skipif(not DECK.is_file(), reason="reg_tests/r-test submodule not checked out")


def _mark(blade, val):
    b = copy.deepcopy(blade)
    b["beam_stiff"][0][0][0] = val  # (nstation, 6, 6) ndarray — a real, round-trippable field
    return b


def test_distinct_beamdyn_blades_do_not_collide():
    driver = OpenFASTDriver()
    fst_vt = driver.read(DECK)

    bd = fst_vt["BeamDyn"]
    bb = fst_vt["BeamDynBlade"]
    bd0 = bd if isinstance(bd, dict) else bd[0]
    bl0 = bb if isinstance(bb, dict) else bb[0]

    # expand the (identical) blades into 3 DISTINCT blades
    fst_vt["BeamDyn"] = [copy.deepcopy(bd0) for _ in range(3)]
    fst_vt["BeamDynBlade"] = [_mark(bl0, 1000.0), _mark(bl0, 2000.0), _mark(bl0, 3000.0)]
    for i in range(3):
        fst_vt["Fst"][f"BDBldFile({i + 1})"] = f"blade_{i + 1}.dat"

    out = Path(tempfile.mkdtemp(prefix="bd_distinct_"))
    driver.write(fst_vt, out, "case")

    fst2 = driver.read(out / "case.fst")
    blades = fst2["BeamDynBlade"]
    assert not isinstance(blades, dict), "3 distinct blades must NOT collapse to one dict"
    markers = [b["beam_stiff"][0][0][0] for b in blades]
    assert sorted(markers) == [1000.0, 2000.0, 3000.0], (
        f"BeamDyn blade collision — blades did not round-trip to distinct values: {markers}"
    )
