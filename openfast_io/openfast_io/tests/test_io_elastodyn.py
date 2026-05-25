from pathlib import Path
import pytest
from openfast_io.io.elastodyn import ElastoDynIO


def test_elastodyn_io_is_module_io():
    from openfast_io.io.base import ModuleIO
    assert issubclass(ElastoDynIO, ModuleIO)


def test_read_returns_dict(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    assert isinstance(result, dict)
    assert 'ElastoDyn' in result
    assert 'ElastoDynBlade' in result
    assert 'ElastoDynTower' in result


def test_read_flapdof1_is_bool(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    ed = result['ElastoDyn']
    assert isinstance(ed['FlapDOF1'], bool)
    assert ed['FlapDOF1'] is True


def test_read_rotspeed_is_float(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    ed = result['ElastoDyn']
    assert isinstance(ed['RotSpeed'], float)
    assert abs(ed['RotSpeed'] - 12.1) < 1e-6


def test_read_numbl(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    assert result['ElastoDyn']['NumBl'] == 3


def test_read_towerht(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    assert abs(result['ElastoDyn']['TowerHt'] - 87.6) < 0.1


def test_read_blade_data(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    blades = result['ElastoDynBlade']
    assert isinstance(blades, list)
    assert len(blades) == 3
    # All 3 blades reference the same file, so all should have data
    assert blades[0]['NBlInpSt'] == 6
    assert len(blades[0]['BlFract']) == 6


def test_read_tower_data(sample_ed_file, tmp_path):
    io = ElastoDynIO()
    result = io.read(sample_ed_file, tmp_path)
    tower = result['ElastoDynTower']
    assert tower['NTwInpSt'] == 3
    assert len(tower['HtFract']) == 3


def test_roundtrip(sample_ed_file, tmp_path):
    """Read → write → read produces consistent dict."""
    io = ElastoDynIO()
    original = io.read(sample_ed_file, tmp_path)
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    # Copy blade/tower files to output dir so write can find them
    import shutil
    shutil.copy(tmp_path / "test_blade.dat", out_dir / "test_blade.dat")
    shutil.copy(tmp_path / "test_tower.dat", out_dir / "test_tower.dat")

    out_path = out_dir / "ElastoDyn_out.dat"
    io.write(original, out_path, out_dir)
    reread = io.read(out_path, out_dir)

    ed_orig = original['ElastoDyn']
    ed_re = reread['ElastoDyn']
    assert ed_orig['FlapDOF1'] == ed_re['FlapDOF1']
    assert abs(ed_orig['RotSpeed'] - ed_re['RotSpeed']) < 1e-6
    assert ed_orig['NumBl'] == ed_re['NumBl']
    assert abs(ed_orig['TowerHt'] - ed_re['TowerHt']) < 0.01


def test_read_real_5mw_file(r_test_5mw_dir):
    """Integration test against the NREL 5MW r-test case."""
    ed_path = r_test_5mw_dir / "NRELOffshrBsline5MW_Onshore_ElastoDyn.dat"
    if not ed_path.exists():
        pytest.skip("ElastoDyn file not found in r-test")
    io = ElastoDynIO()
    result = io.read(ed_path, r_test_5mw_dir)
    ed = result['ElastoDyn']
    assert ed['NumBl'] == 3
    assert abs(ed['TowerHt'] - 87.6) < 0.1
