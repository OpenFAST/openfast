import io
import pytest
from openfast_io.outlist import OutList


def test_enable_and_query():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1'])
    ol.enable('AeroDyn', ['RtAeroCp'])
    assert ol.enabled('ElastoDyn') == {'RotSpeed', 'BldPitch1'}
    result = ol.all_enabled()
    assert 'RtAeroCp' in result
    assert 'BldPitch1' in result
    assert 'RotSpeed' in result


def test_disable():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1', 'GenPwr'])
    ol.disable('ElastoDyn', ['GenPwr'])
    assert 'GenPwr' not in ol.enabled('ElastoDyn')
    assert ol.is_enabled('RotSpeed')
    assert not ol.is_enabled('GenPwr')


def test_enabled_by_module():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1'])
    ol.enable('AeroDyn', ['RtAeroCp'])
    by_mod = ol.enabled_by_module()
    assert 'ElastoDyn' in by_mod
    assert 'AeroDyn' in by_mod
    assert by_mod['ElastoDyn'] == ['BldPitch1', 'RotSpeed']


def test_roundtrip_fst_output():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1'])
    legacy = ol.to_fst_output()
    # Round-trip back
    ol2 = OutList.from_fst_output(legacy)
    assert ol2.enabled('ElastoDyn') == ol.enabled('ElastoDyn')


def test_validate_catches_unknown_channels():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'TotallyBogusChannel'])
    unknown = ol.validate()
    assert any('TotallyBogusChannel' in u for u in unknown)
    assert not any('RotSpeed' in u for u in unknown)


def test_read_from_file():
    content = '"RotSpeed"\n"BldPitch1", "GenPwr"\nEND of OutList\n'
    f = io.StringIO(content)
    channels = OutList.read_from_file(f, 'ElastoDyn')
    assert 'RotSpeed' in channels
    assert 'BldPitch1' in channels
    assert 'GenPwr' in channels


def test_write_to_file():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed', 'BldPitch1'])
    buf = io.StringIO()
    ol.write_to_file(buf, 'ElastoDyn')
    text = buf.getvalue()
    assert '"BldPitch1"' in text
    assert '"RotSpeed"' in text
    assert 'END' in text


def test_clear():
    ol = OutList()
    ol.enable('ElastoDyn', ['RotSpeed'])
    ol.enable('AeroDyn', ['RtAeroCp'])
    ol.clear('ElastoDyn')
    assert ol.enabled('ElastoDyn') == set()
    assert ol.enabled('AeroDyn') == {'RtAeroCp'}
    ol.clear()
    assert ol.all_enabled() == []


def test_is_enabled_checks_all_modules():
    ol = OutList()
    ol.enable('AeroDyn', ['RtAeroCp'])
    assert ol.is_enabled('RtAeroCp')
    assert not ol.is_enabled('RotSpeed')
