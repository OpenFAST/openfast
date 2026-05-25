"""Tests for SimpleElastoDynIO and AeroDiskIO."""
import os
import tempfile

import pytest

from openfast_io.io.simple_elastodyn import SimpleElastoDynIO
from openfast_io.io.aerodisk import AeroDiskIO

_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                       'glue-codes', 'openfast')

_ADSK_SED_DIR = os.path.join(_RTEST, '5MW_Land_DLL_WTurb_ADsk_SED')
_SED_FILE = os.path.join(_ADSK_SED_DIR, 'NREL_5MW_Simplified-ElastoDyn.dat')
_ADSK_FILE = os.path.join(_ADSK_SED_DIR, 'NRELOffshrBsline5MW_Onshore_AeroDisk.dat')

_has_sed = os.path.isfile(_SED_FILE)
_has_adsk = os.path.isfile(_ADSK_FILE)

skipsed = pytest.mark.skipif(not _has_sed, reason='SimpleElastoDyn r-test data missing')
skipadsk = pytest.mark.skipif(not _has_adsk, reason='AeroDisk r-test data missing')


# =================================================================
# SimpleElastoDynIO
# =================================================================
class TestSimpleElastoDynIO:
    @skipsed
    def test_read_returns_dict(self):
        io = SimpleElastoDynIO()
        result = io.read(_SED_FILE)
        assert 'SimpleElastoDyn' in result
        assert isinstance(result['SimpleElastoDyn'], dict)

    @skipsed
    def test_read_scalars(self):
        io = SimpleElastoDynIO()
        sed = io.read(_SED_FILE)['SimpleElastoDyn']
        assert sed['Echo'] == True
        assert sed['IntMethod'] == 3
        assert sed['NumBl'] == 3
        assert sed['TipRad'] == 63.0
        assert sed['GBoxRatio'] == 97.0

    @skipsed
    def test_read_outlist(self):
        io = SimpleElastoDynIO()
        sed = io.read(_SED_FILE)['SimpleElastoDyn']
        assert '_outlist' in sed
        assert len(sed['_outlist']) > 0
        assert 'BlPitch1' in sed['_outlist']

    @skipsed
    def test_write_then_reread(self):
        io = SimpleElastoDynIO()
        data = io.read(_SED_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            s1 = data['SimpleElastoDyn']
            s2 = data2['SimpleElastoDyn']
            assert s2['NumBl'] == s1['NumBl']
            assert s2['TipRad'] == s1['TipRad']
            assert s2['GBoxRatio'] == s1['GBoxRatio']
            assert len(s2['_outlist']) == len(s1['_outlist'])
        finally:
            os.unlink(tmp_path)

    @skipsed
    def test_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(SimpleElastoDynIO, ModuleIO)


# =================================================================
# AeroDiskIO
# =================================================================
class TestAeroDiskIO:
    @skipadsk
    def test_read_returns_dict(self):
        io = AeroDiskIO()
        result = io.read(_ADSK_FILE)
        assert 'AeroDisk' in result
        assert isinstance(result['AeroDisk'], dict)

    @skipadsk
    def test_read_scalars(self):
        io = AeroDiskIO()
        ad = io.read(_ADSK_FILE)['AeroDisk']
        assert ad['AirDens'] == 1.225
        assert ad['RotorRad'] == 63.0
        assert 'TSR' in ad['InColNames']

    @skipadsk
    def test_read_disk_table(self):
        io = AeroDiskIO()
        ad = io.read(_ADSK_FILE)['AeroDisk']
        assert 'actuatorDiskTable' in ad
        tbl = ad['actuatorDiskTable']
        assert 'attr' in tbl
        assert 'data' in tbl
        assert len(tbl['data']) > 0

    @skipadsk
    def test_read_outlist(self):
        io = AeroDiskIO()
        ad = io.read(_ADSK_FILE)['AeroDisk']
        assert '_outlist' in ad
        assert len(ad['_outlist']) > 0

    @skipadsk
    def test_write_then_reread(self):
        io = AeroDiskIO()
        data = io.read(_ADSK_FILE)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, 'AeroDisk.dat')
            io.write(data, out_path)
            data2 = io.read(out_path)
            ad1 = data['AeroDisk']
            ad2 = data2['AeroDisk']
            assert ad2['AirDens'] == ad1['AirDens']
            assert ad2['RotorRad'] == ad1['RotorRad']
            assert ad2['InColNames'] == ad1['InColNames']

    @skipadsk
    def test_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(AeroDiskIO, ModuleIO)
