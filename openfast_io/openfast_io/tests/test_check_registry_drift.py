from pathlib import Path
import pytest
from openfast_io.tools.check_registry_drift import parse_registry, scan_reader_params, RegistryParam


def test_parse_registry_returns_list_of_params(tmp_path):
    """parse_registry correctly extracts ED_InputFile entries."""
    reg_content = """\
# ElastoDyn Registry test
typedef\t^\tED_InputFile\tLOGICAL\tFlapDOF1\t-\t-\t-\t"First flapwise blade mode DOF"\t-
typedef\t^\tED_InputFile\tReKi\tRotSpeed\t-\t-\t-\t"Initial rotor speed"\trad/s
typedef\t^\tContinuousStateType\tR8Ki\tQT\t{:}\t-\t-\t"Displacement DOF vector"\t-
"""
    reg_path = tmp_path / "ElastoDyn_Registry.txt"
    reg_path.write_text(reg_content)

    params = parse_registry(reg_path, 'ED_InputFile')
    names = [p.name for p in params]

    assert 'FlapDOF1' in names
    assert 'RotSpeed' in names
    assert 'QT' not in names  # ContinuousStateType, not InputFile


def test_parse_registry_captures_units(tmp_path):
    reg_content = "typedef\t^\tED_InputFile\tReKi\tRotSpeed\t-\t-\t-\t\"Initial rotor speed\"\trad/s\n"
    reg_path = tmp_path / "ED_Registry.txt"
    reg_path.write_text(reg_content)
    params = parse_registry(reg_path, 'ED_InputFile')
    assert params[0].units == 'rad/s'


def test_scan_reader_params_finds_fst_vt_assignments(tmp_path):
    reader_content = """\
def read_ElastoDyn(self):
    fst_vt = self.fst_vt
    fst_vt['ElastoDyn']['FlapDOF1'] = True
    fst_vt['ElastoDyn']['RotSpeed'] = 12.1
    fst_vt['AeroDyn']['TwrAero'] = False
"""
    reader_path = tmp_path / "FAST_reader.py"
    reader_path.write_text(reader_content)

    params = scan_reader_params(reader_path, 'ElastoDyn')
    assert 'FlapDOF1' in params
    assert 'RotSpeed' in params
    assert 'TwrAero' not in params  # AeroDyn key, not ElastoDyn


def test_parse_registry_handles_empty_file(tmp_path):
    reg_path = tmp_path / "empty.txt"
    reg_path.write_text("")
    params = parse_registry(reg_path, 'ED_InputFile')
    assert params == []


def test_parse_registry_skips_comments(tmp_path):
    reg_content = """\
# This is a comment
! This is also a comment
typedef\t^\tED_InputFile\tReKi\tDT\t-\t-\t-\t"Time step"\ts
"""
    reg_path = tmp_path / "reg.txt"
    reg_path.write_text(reg_content)
    params = parse_registry(reg_path, 'ED_InputFile')
    assert len(params) == 1
    assert params[0].name == 'DT'
