"""Tests for HydroDynIO and SeaStateIO."""
import os
import tempfile

import numpy as np
import pytest

from openfast_io.io.hydrodyn import HydroDynIO
from openfast_io.io.seastate import SeaStateIO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                       'glue-codes', 'openfast')
_OC4_DIR = os.path.join(_RTEST, 'StC_test_OC4Semi')
_HD_FILE = os.path.join(_OC4_DIR,
                        'NRELOffshrBsline5MW_OC4DeepCwindSemi_HydroDyn.dat')
_SS_FILE = os.path.join(_OC4_DIR, 'SeaState.dat')

_has_oc4 = os.path.isfile(_HD_FILE) and os.path.isfile(_SS_FILE)
skipoc4 = pytest.mark.skipif(not _has_oc4, reason='r-test data missing')


# ===================================================================
# HydroDynIO tests
# ===================================================================
class TestHydroDynIO:
    """Test HydroDynIO read and write against OC4 r-test data."""

    @skipoc4
    def test_read_returns_dict(self):
        io = HydroDynIO()
        result = io.read(_HD_FILE)
        assert 'HydroDyn' in result
        hd = result['HydroDyn']
        assert isinstance(hd, dict)

    @skipoc4
    def test_read_scalars(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        assert isinstance(hd['PotMod'], int)
        assert isinstance(hd['RdtnTMax'], float)
        assert hd['NBody'] >= 1

    @skipoc4
    def test_read_potfile_paths(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        assert isinstance(hd['PotFile'], list)
        assert len(hd['PotFile']) == hd['NBody']

    @skipoc4
    def test_read_arrays(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        assert len(hd['WAMITULEN']) == hd['NBody']
        assert len(hd['PtfmRefxt']) == hd['NBody']

    @skipoc4
    def test_read_matrices(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        assert hd['AddCLin'].shape[0] == 6
        assert hd['AddBLin'].shape[0] == 6
        assert hd['AddBQuad'].shape[0] == 6

    @skipoc4
    def test_read_joints(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        n = hd['NJoints']
        assert len(hd['JointID']) == n
        assert len(hd['Jointxi']) == n

    @skipoc4
    def test_read_members(self):
        io = HydroDynIO()
        hd = io.read(_HD_FILE)['HydroDyn']
        n = hd['NMembers']
        assert len(hd['MemberID']) == n
        assert len(hd['PropPot']) == n

    @skipoc4
    def test_write_then_reread(self):
        io = HydroDynIO()
        data = io.read(_HD_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            hd1 = data['HydroDyn']
            hd2 = data2['HydroDyn']
            # Scalar fields
            assert hd2['PotMod'] == hd1['PotMod']
            assert hd2['NBody'] == hd1['NBody']
            assert hd2['NJoints'] == hd1['NJoints']
            assert hd2['NMembers'] == hd1['NMembers']
            # Matrices
            np.testing.assert_allclose(hd2['AddCLin'], hd1['AddCLin'], atol=1e-10)
            np.testing.assert_allclose(hd2['AddBLin'], hd1['AddBLin'], atol=1e-10)
        finally:
            os.unlink(tmp_path)


# ===================================================================
# SeaStateIO tests
# ===================================================================
class TestSeaStateIO:
    """Test SeaStateIO read and write against OC4 r-test data."""

    @skipoc4
    def test_read_returns_dict(self):
        io = SeaStateIO()
        result = io.read(_SS_FILE)
        assert 'SeaState' in result
        ss = result['SeaState']
        assert isinstance(ss, dict)

    @skipoc4
    def test_read_scalars(self):
        io = SeaStateIO()
        ss = io.read(_SS_FILE)['SeaState']
        # WtrDens/WtrDpth/MSL2SWL may be "default" strings
        assert isinstance(ss['WtrDens'], (float, str))
        assert isinstance(ss['WtrDpth'], (float, str))
        assert isinstance(ss['WaveMod'], int)

    @skipoc4
    def test_read_wave_params(self):
        io = SeaStateIO()
        ss = io.read(_SS_FILE)['SeaState']
        assert isinstance(ss['WaveHs'], float)
        assert isinstance(ss['WaveTp'], float)
        assert isinstance(ss['WaveSeed1'], int)

    @skipoc4
    def test_read_current(self):
        io = SeaStateIO()
        ss = io.read(_SS_FILE)['SeaState']
        assert isinstance(ss['CurrMod'], int)

    @skipoc4
    def test_write_then_reread(self):
        io = SeaStateIO()
        data = io.read(_SS_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            ss1 = data['SeaState']
            ss2 = data2['SeaState']
            assert ss2['WtrDens'] == ss1['WtrDens']
            assert ss2['WtrDpth'] == ss1['WtrDpth']
            assert ss2['WaveMod'] == ss1['WaveMod']
            assert ss2['CurrMod'] == ss1['CurrMod']
        finally:
            os.unlink(tmp_path)


# ===================================================================
# Synthetic data test
# ===================================================================
class TestSynthetic:
    def test_hydrodynamics_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(HydroDynIO, ModuleIO)
        assert issubclass(SeaStateIO, ModuleIO)
