"""Tests for ExtPtfmIO."""
import os
import tempfile

import numpy as np
import pytest

from openfast_io.io.extptfm import ExtPtfmIO

_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                       'glue-codes', 'openfast')

_EP_DIR = os.path.join(_RTEST, '5MW_OC4Jckt_ExtPtfm')
_EP_FILE = os.path.join(_EP_DIR, 'ExtPtfm.dat')

_has_ep = os.path.isfile(_EP_FILE)
skipep = pytest.mark.skipif(not _has_ep, reason='ExtPtfm r-test data missing')


class TestExtPtfmIO:
    @skipep
    def test_read_returns_dict(self):
        io = ExtPtfmIO()
        result = io.read(_EP_FILE)
        assert 'ExtPtfm' in result
        assert isinstance(result['ExtPtfm'], dict)

    @skipep
    def test_read_scalars(self):
        io = ExtPtfmIO()
        ep = io.read(_EP_FILE)['ExtPtfm']
        assert ep['Echo'] == False
        assert ep['IntMethod'] == 3
        assert ep['RBMod'] == 0
        assert ep['NActiveDOFList'] == -1
        assert ep['HasConnections'] == False
        assert ep['HasUserForcing'] == True

    @skipep
    def test_read_superelement(self):
        io = ExtPtfmIO()
        ep = io.read(_EP_FILE)['ExtPtfm']
        flex = ep['FlexASCII']
        assert flex['nDOF'] == 31
        assert flex['MassMatrix'].shape == (31, 31)
        assert flex['StiffnessMatrix'].shape == (31, 31)
        assert flex['DampingMatrix'].shape == (31, 31)
        assert flex['WeightConstant'].shape == (1, 31)
        assert flex['WeightStiffness'].shape == (31, 31)

    @skipep
    def test_read_user_forcing(self):
        io = ExtPtfmIO()
        ep = io.read(_EP_FILE)['ExtPtfm']
        assert 'UserForcing' in ep
        uf = ep['UserForcing']
        assert uf['nSteps'] == 501
        # time + 31 DOF columns
        assert uf['ForceTimeSeries'].shape == (501, 32)

    @skipep
    def test_read_outlist(self):
        io = ExtPtfmIO()
        ep = io.read(_EP_FILE)['ExtPtfm']
        assert '_outlist' in ep
        assert len(ep['_outlist']) > 0
        assert 'IntrfFx' in ep['_outlist']

    @skipep
    def test_write_then_reread(self):
        io = ExtPtfmIO()
        data = io.read(_EP_FILE)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, 'ExtPtfm.dat')
            io.write(data, out_path)
            data2 = io.read(out_path)
            ep1 = data['ExtPtfm']
            ep2 = data2['ExtPtfm']
            assert ep2['IntMethod'] == ep1['IntMethod']
            assert ep2['RBMod'] == ep1['RBMod']
            assert ep2['FlexASCII']['nDOF'] == ep1['FlexASCII']['nDOF']
            np.testing.assert_allclose(
                ep2['FlexASCII']['MassMatrix'],
                ep1['FlexASCII']['MassMatrix'],
                rtol=1e-6,
            )

    @skipep
    def test_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(ExtPtfmIO, ModuleIO)
