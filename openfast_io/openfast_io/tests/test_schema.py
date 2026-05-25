from openfast_io.schema import get_schema, get_param_info, FILE_REF_PARAMS


def test_get_schema_returns_dict():
    s = get_schema('ElastoDyn', '5.0.0')
    assert isinstance(s, dict)
    assert len(s) > 10


def test_known_param_has_required_fields():
    s = get_schema('ElastoDyn', '5.0.0')
    p = s['FlapDOF1']
    assert 'type' in p
    assert 'desc' in p
    assert p['type'] == bool


def test_get_param_info_returns_description():
    info = get_param_info('ElastoDyn', 'FlapDOF1', '5.0.0')
    assert 'desc' in info
    assert 'flapwise' in info['desc'].lower()


def test_file_ref_params_marks_bldfile():
    s = get_schema('ElastoDyn', '5.0.0')
    assert s.get('BldFile', {}).get('is_file_ref') is True or \
           s.get('TwrFile', {}).get('is_file_ref') is True


def test_removed_param_not_in_v5_schema():
    s = get_schema('AeroDyn', '5.0.0')
    assert 'Buoyancy' not in s


def test_removed_param_in_v4_schema():
    s = get_schema('AeroDyn', '4.0.0')
    assert 'Buoyancy' in s


def test_file_ref_params_dict_has_fst_key():
    assert 'Fst' in FILE_REF_PARAMS
    assert 'EDFile' in FILE_REF_PARAMS['Fst']


def test_fst_schema_has_tmax():
    s = get_schema('Fst', '5.0.0')
    assert 'TMax' in s
    assert s['TMax']['type'] == float
    assert s['TMax']['units'] == 's'


def test_get_param_info_unknown_param():
    info = get_param_info('ElastoDyn', 'NonExistentParam', '5.0.0')
    assert info == {}
