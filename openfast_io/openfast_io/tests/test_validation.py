import pytest
from openfast_io.validation import validate_fst_vt, ValidationIssue


def test_valid_deck_returns_no_errors():
    fst_vt = {
        'Fst': {'CompElast': 1, 'CompAero': 2, 'CompServo': 0,
                'CompHydro': 0, 'NRotors': 1},
        'ElastoDyn': {'NumBl': 3},
        'AeroDynBlade': [{}, {}, {}],
        'outlist': {},
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=False)
    errors = [i for i in issues if i.severity == 'ERROR']
    assert len(errors) == 0


def test_blade_count_mismatch_is_error():
    fst_vt = {
        'Fst': {'CompElast': 1, 'CompAero': 2, 'CompServo': 0,
                'CompHydro': 0, 'NRotors': 1},
        'ElastoDyn': {'NumBl': 3},
        'AeroDynBlade': [{'data': True}, {'data': True}],  # only 2 for a 3-blade turbine
        'outlist': {},
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=False)
    errors = [i for i in issues if i.severity == 'ERROR']
    assert any('blade' in i.message.lower() for i in errors)


def test_removed_param_produces_warning():
    fst_vt = {
        'Fst': {'CompElast': 1, 'CompAero': 2, 'NRotors': 1},
        'AeroDyn': {'Buoyancy': True},
        'ElastoDyn': {'NumBl': 3},
        'AeroDynBlade': [{}, {}, {}],
        'outlist': {},
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=False)
    warnings = [i for i in issues if i.severity == 'WARNING']
    assert any('Buoyancy' in i.message for i in warnings)


def test_hydrodyn_info_when_comp_zero():
    fst_vt = {
        'Fst': {'CompHydro': 0, 'CompServo': 0},
        'HydroDyn': {'WaveMod': 1},
        'ElastoDyn': {},
        'AeroDynBlade': [],
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=False)
    infos = [i for i in issues if i.severity == 'INFO']
    assert any('HydroDyn' in i.message for i in infos)


def test_validation_issue_dataclass():
    issue = ValidationIssue(severity='ERROR', modules=['Fst'], parameter='TMax',
                            message='TMax must be positive')
    assert issue.severity == 'ERROR'
    assert issue.modules == ['Fst']


def test_check_files_with_base_dir_existing_ref_no_issue(tmp_path):
    """A file reference stored relative to the case dir (base_dir) should
    resolve cleanly when it actually exists there, even though it does not
    exist relative to the current working directory."""
    blade_file = tmp_path / "blade1.dat"
    blade_file.write_text("dummy blade file\n")

    fst_vt = {
        'Fst': {'CompElast': 1, 'CompAero': 2, 'CompServo': 0,
                'CompHydro': 0, 'NRotors': 1},
        'ElastoDyn': {'NumBl': 3, 'BldFile1': 'blade1.dat'},
        'AeroDynBlade': [{}, {}, {}],
        'outlist': {},
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=True,
                              base_dir=tmp_path)
    errors = [i for i in issues if i.severity == 'ERROR']
    assert not any('BldFile1' in (i.parameter or '') for i in errors)


def test_check_files_with_base_dir_missing_ref_is_error(tmp_path):
    """A file reference that does not exist under base_dir should still be
    flagged as an ERROR."""
    fst_vt = {
        'Fst': {'CompElast': 1, 'CompAero': 2, 'CompServo': 0,
                'CompHydro': 0, 'NRotors': 1},
        'ElastoDyn': {'NumBl': 3, 'BldFile1': 'does_not_exist.dat'},
        'AeroDynBlade': [{}, {}, {}],
        'outlist': {},
    }
    issues = validate_fst_vt(fst_vt, version='5.0.0', check_files=True,
                              base_dir=tmp_path)
    errors = [i for i in issues if i.severity == 'ERROR']
    assert any('BldFile1' in (i.parameter or '') for i in errors)
    assert any('does_not_exist.dat' in i.message for i in errors)
