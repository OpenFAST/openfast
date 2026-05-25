"""Tests for facade classes (backward-compat wrappers).

Tests both import paths:
  - from openfast_io.FAST_reader import InputReader_OpenFAST  (canonical)
  - from openfast_io.facade    import InputReader_Facade       (alias)
"""
import os
import tempfile

import pytest

from openfast_io.FAST_reader import InputReader_OpenFAST
from openfast_io.FAST_writer import InputWriter_OpenFAST
from openfast_io.facade import InputReader_Facade, InputWriter_Facade

_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                       'glue-codes', 'openfast')

_5MW_DIR = os.path.join(_RTEST, '5MW_Land_DLL_WTurb')
_FST_FILE = '5MW_Land_DLL_WTurb.fst'

_has_5mw = os.path.isfile(os.path.join(_5MW_DIR, _FST_FILE))
skip5mw = pytest.mark.skipif(not _has_5mw, reason='5MW r-test data missing')


class TestReaderFacade:
    @skip5mw
    def test_reader_execute(self):
        reader = InputReader_Facade()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()
        assert reader.fst_vt['Fst']['TMax'] == 60.0
        assert reader.fst_vt['Fst']['CompElast'] == 1

    @skip5mw
    def test_reader_populates_elastodyn(self):
        reader = InputReader_Facade()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()
        assert 'NumBl' in reader.fst_vt['ElastoDyn']

    @skip5mw
    def test_reader_has_fst_vt_attribute(self):
        reader = InputReader_Facade()
        assert hasattr(reader, 'fst_vt')
        assert 'Fst' in reader.fst_vt


class TestWriterFacade:
    @skip5mw
    def test_writer_execute(self):
        reader = InputReader_Facade()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()

        with tempfile.TemporaryDirectory() as tmpdir:
            writer = InputWriter_Facade()
            writer.fst_vt = reader.fst_vt
            writer.FAST_namingOut = 'test_write'
            writer.FAST_runDirectory = tmpdir
            writer.execute()

            fst_out = os.path.join(tmpdir, 'test_write.fst')
            assert os.path.isfile(fst_out)

    @skip5mw
    def test_writer_roundtrip(self):
        """Write then re-read preserves key scalars."""
        reader = InputReader_Facade()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()
        original = reader.fst_vt

        with tempfile.TemporaryDirectory() as tmpdir:
            writer = InputWriter_Facade()
            writer.fst_vt = original
            writer.FAST_namingOut = 'rt'
            writer.FAST_runDirectory = tmpdir
            writer.execute()

            reader2 = InputReader_Facade()
            reader2.FAST_InputFile = 'rt.fst'
            reader2.FAST_directory = tmpdir
            reader2.execute()

            assert reader2.fst_vt['Fst']['TMax'] == original['Fst']['TMax']
            assert reader2.fst_vt['Fst']['CompElast'] == original['Fst']['CompElast']

    def test_writer_update(self):
        writer = InputWriter_Facade()
        writer.fst_vt = {'Fst': {'TMax': 60.0, 'DT': 0.005}}
        writer.update({'Fst': {'TMax': 120.0}})
        assert writer.fst_vt['Fst']['TMax'] == 120.0
        assert writer.fst_vt['Fst']['DT'] == 0.005


class TestLegacyImportPaths:
    """Verify that the canonical FAST_reader / FAST_writer import paths work."""

    def test_facade_aliases_are_same_class(self):
        assert InputReader_Facade is InputReader_OpenFAST
        assert InputWriter_Facade is InputWriter_OpenFAST

    @skip5mw
    def test_legacy_reader_execute(self):
        reader = InputReader_OpenFAST()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()
        assert reader.fst_vt['Fst']['TMax'] == 60.0
        assert 'NumBl' in reader.fst_vt['ElastoDyn']

    @skip5mw
    def test_legacy_writer_roundtrip(self):
        reader = InputReader_OpenFAST()
        reader.FAST_InputFile = _FST_FILE
        reader.FAST_directory = _5MW_DIR
        reader.execute()

        with tempfile.TemporaryDirectory() as tmpdir:
            writer = InputWriter_OpenFAST()
            writer.fst_vt = reader.fst_vt
            writer.FAST_namingOut = 'legacy_rt'
            writer.FAST_runDirectory = tmpdir
            writer.execute()

            reader2 = InputReader_OpenFAST()
            reader2.FAST_InputFile = 'legacy_rt.fst'
            reader2.FAST_directory = tmpdir
            reader2.execute()

            assert reader2.fst_vt['Fst']['TMax'] == reader.fst_vt['Fst']['TMax']

    def test_parsing_re_exports(self):
        """External code importing parsing helpers from FAST_reader still works."""
        from openfast_io.FAST_reader import bool_read, float_read, int_read, quoted_read
        assert bool_read('True') is True
        assert float_read('3.14') == pytest.approx(3.14)
        assert int_read('42') == 42

    def test_writer_helper_functions(self):
        """External code importing helper functions from FAST_writer still works."""
        from openfast_io.FAST_writer import auto_format, float_default_out, int_default_out
        assert callable(auto_format)
        assert '3.1400' in float_default_out(3.14, trim=True)
