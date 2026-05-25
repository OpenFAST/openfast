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
