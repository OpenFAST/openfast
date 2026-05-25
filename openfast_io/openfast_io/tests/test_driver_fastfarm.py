"""Tests for FASTFarmDriver."""

import pytest
from pathlib import Path

from openfast_io.drivers import FASTFarmDriver


@pytest.fixture
def ff_tsinflow_dir(r_test_fastfarm_dir):
    """Path to TSinflow test case."""
    candidate = r_test_fastfarm_dir / "TSinflow"
    if not candidate.exists():
        pytest.skip("TSinflow FAST.Farm case not found")
    return candidate


@pytest.fixture
def ff_data(ff_tsinflow_dir):
    """Read the TSinflow FAST.Farm file and return the result dict."""
    drv = FASTFarmDriver()
    return drv.read(ff_tsinflow_dir / "FAST.Farm.fstf")


class TestFASTFarmDriverInit:
    def test_instantiates(self):
        drv = FASTFarmDriver()
        assert drv is not None

    def test_has_read_method(self):
        drv = FASTFarmDriver()
        assert callable(getattr(drv, "read", None))


class TestFASTFarmRead:
    def test_read_returns_dict_with_keys(self, ff_data):
        assert "FASTFarm" in ff_data
        assert "Turbines" in ff_data

    def test_farm_scalars(self, ff_data):
        ff = ff_data["FASTFarm"]
        assert ff["TMax"] == pytest.approx(90.0)
        assert ff["NumTurbines"] == 2
        assert ff["Mod_AmbWind"] in (1, 2, 3)

    def test_turbine_count_matches(self, ff_data):
        assert len(ff_data["Turbines"]) == ff_data["FASTFarm"]["NumTurbines"]

    def test_turbine_rows_positions(self, ff_data):
        rows = ff_data["FASTFarm"]["TurbineRows"]
        assert len(rows) == 2
        # Each row should have WT_X, WT_Y, WT_Z, WT_FASTInFile
        for row in rows:
            assert "WT_X" in row
            assert "WT_Y" in row
            assert "WT_Z" in row
            assert "WT_FASTInFile" in row

    def test_each_turbine_has_fst_vt_structure(self, ff_data):
        """Each turbine entry should contain OpenFAST driver output keys."""
        for turb in ff_data["Turbines"]:
            # Should have at least ElastoDyn (CompElast=1) and ServoDyn
            assert "ElastoDyn" in turb, "Turbine missing ElastoDyn"
            assert "ServoDyn" in turb, "Turbine missing ServoDyn"

    def test_turbine_farm_positions(self, ff_data):
        """Each turbine has _farm_position metadata."""
        for turb in ff_data["Turbines"]:
            pos = turb.get("_farm_position")
            assert pos is not None
            assert "WT_X" in pos
            assert "WT_Y" in pos

    def test_shared_mooring_fields(self, ff_data):
        ff = ff_data["FASTFarm"]
        assert "SharedMoorFile" in ff
        assert "DT_Mooring" in ff

    def test_ambient_wind_inflowwind_fields(self, ff_data):
        ff = ff_data["FASTFarm"]
        assert "DT_Low" in ff
        assert "DT_High" in ff
        assert "NX_Low" in ff
        assert "InflowFile" in ff

    def test_remaining_text_captured(self, ff_data):
        ff = ff_data["FASTFarm"]
        # _remaining should capture wake dynamics, curled-wake, output sections
        assert "_remaining" in ff
        assert len(ff["_remaining"]) > 0
