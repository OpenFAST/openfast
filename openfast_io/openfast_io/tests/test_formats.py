import json
import pytest
from openfast_io.formats import fst_vt_to_json, fst_vt_from_json, fst_vt_to_yaml, fst_vt_from_yaml


def test_json_roundtrip_preserves_values():
    fst_vt = {'ElastoDyn': {'NumBl': 3, 'RotSpeed': 12.1, 'FlapDOF1': True}}
    s = fst_vt_to_json(fst_vt)
    recovered = fst_vt_from_json(s)
    assert recovered['ElastoDyn']['NumBl'] == 3
    assert abs(recovered['ElastoDyn']['RotSpeed'] - 12.1) < 1e-9
    assert recovered['ElastoDyn']['FlapDOF1'] is True


def test_json_is_valid_json():
    fst_vt = {'Fst': {'TMax': 60.0}, 'ElastoDyn': {'NumBl': 3}}
    s = fst_vt_to_json(fst_vt)
    parsed = json.loads(s)
    assert parsed['Fst']['TMax'] == 60.0


def test_yaml_roundtrip():
    fst_vt = {'ElastoDyn': {'NumBl': 3, 'RotSpeed': 12.1}}
    s = fst_vt_to_yaml(fst_vt)
    recovered = fst_vt_from_yaml(s)
    assert recovered['ElastoDyn']['NumBl'] == 3


def test_json_handles_empty():
    fst_vt = {'Fst': {}, 'ElastoDyn': {}}
    s = fst_vt_to_json(fst_vt)
    recovered = fst_vt_from_json(s)
    assert recovered == fst_vt


def test_yaml_handles_lists():
    fst_vt = {'AeroDynPolar': [{'Re': 1e6}, {'Re': 2e6}]}
    s = fst_vt_to_yaml(fst_vt)
    recovered = fst_vt_from_yaml(s)
    assert len(recovered['AeroDynPolar']) == 2
