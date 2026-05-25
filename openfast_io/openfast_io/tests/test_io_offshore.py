"""Tests for SubDynIO, MoorDynIO, and MAPIO."""
import os
import tempfile

import numpy as np
import pytest

from openfast_io.io.subdyn import SubDynIO
from openfast_io.io.moordyn import MoorDynIO
from openfast_io.io.map_io import MAPIO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                       'glue-codes', 'openfast')

_OC3_DIR = os.path.join(_RTEST, '5MW_OC3Mnpl_DLL_WTurb_WavesIrr')
_SD_FILE = os.path.join(_OC3_DIR, 'NRELOffshrBsline5MW_OC3Monopile_SubDyn.dat')

_OC4_DIR = os.path.join(_RTEST, 'StC_test_OC4Semi')
_MD_FILE = os.path.join(_OC4_DIR, 'NRELOffshrBsline5MW_OC4DeepCwindSemi_MoorDyn.dat')

_MAP_DIR = os.path.join(_RTEST, '5MW_OC4Semi_WSt_WavesWN')
_MAP_FILE = os.path.join(_MAP_DIR, 'NRELOffshrBsline5MW_OC4DeepCwindSemi_MAP.dat')

_has_sd = os.path.isfile(_SD_FILE)
_has_md = os.path.isfile(_MD_FILE)
_has_map = os.path.isfile(_MAP_FILE)

skipsd = pytest.mark.skipif(not _has_sd, reason='SubDyn r-test data missing')
skipmd = pytest.mark.skipif(not _has_md, reason='MoorDyn r-test data missing')
skipmap = pytest.mark.skipif(not _has_map, reason='MAP r-test data missing')


# =================================================================
# SubDynIO
# =================================================================
class TestSubDynIO:
    @skipsd
    def test_read_returns_dict(self):
        io = SubDynIO()
        result = io.read(_SD_FILE)
        assert 'SubDyn' in result
        assert isinstance(result['SubDyn'], dict)

    @skipsd
    def test_read_scalars(self):
        io = SubDynIO()
        sd = io.read(_SD_FILE)['SubDyn']
        assert isinstance(sd['FEMMod'], int)
        assert isinstance(sd['NJoints'], int)
        assert sd['NJoints'] > 0

    @skipsd
    def test_read_joints(self):
        io = SubDynIO()
        sd = io.read(_SD_FILE)['SubDyn']
        assert len(sd['JointID']) == sd['NJoints']
        assert len(sd['JointXss']) == sd['NJoints']

    @skipsd
    def test_read_members(self):
        io = SubDynIO()
        sd = io.read(_SD_FILE)['SubDyn']
        assert len(sd['MemberID']) == sd['NMembers']

    @skipsd
    def test_read_prop_sets(self):
        io = SubDynIO()
        sd = io.read(_SD_FILE)['SubDyn']
        assert len(sd['PropSetID1']) == sd['NPropSetsCyl']

    @skipsd
    def test_write_then_reread(self):
        io = SubDynIO()
        data = io.read(_SD_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            sd1 = data['SubDyn']
            sd2 = data2['SubDyn']
            assert sd2['NJoints'] == sd1['NJoints']
            assert sd2['NMembers'] == sd1['NMembers']
            assert sd2['NPropSetsCyl'] == sd1['NPropSetsCyl']
        finally:
            os.unlink(tmp_path)


# =================================================================
# MoorDynIO
# =================================================================
class TestMoorDynIO:
    @skipmd
    def test_read_returns_dict(self):
        io = MoorDynIO()
        result = io.read(_MD_FILE)
        assert 'MoorDyn' in result
        assert isinstance(result['MoorDyn'], dict)

    @skipmd
    def test_read_line_types(self):
        io = MoorDynIO()
        md = io.read(_MD_FILE)['MoorDyn']
        assert 'Name' in md
        assert len(md['Name']) >= 1

    @skipmd
    def test_read_points(self):
        io = MoorDynIO()
        md = io.read(_MD_FILE)['MoorDyn']
        assert 'Point_ID' in md
        assert len(md['Point_ID']) > 0

    @skipmd
    def test_read_lines(self):
        io = MoorDynIO()
        md = io.read(_MD_FILE)['MoorDyn']
        assert 'Line_ID' in md
        assert len(md['Line_ID']) > 0

    @skipmd
    def test_write_then_reread(self):
        io = MoorDynIO()
        data = io.read(_MD_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            md1 = data['MoorDyn']
            md2 = data2['MoorDyn']
            assert len(md2['Name']) == len(md1['Name'])
            assert len(md2['Point_ID']) == len(md1['Point_ID'])
            assert len(md2['Line_ID']) == len(md1['Line_ID'])
        finally:
            os.unlink(tmp_path)


# =================================================================
# MAPIO
# =================================================================
class TestMAPIO:
    @skipmap
    def test_read_returns_dict(self):
        io = MAPIO()
        result = io.read(_MAP_FILE)
        assert 'MAP' in result
        assert isinstance(result['MAP'], dict)

    @skipmap
    def test_read_line_types(self):
        io = MAPIO()
        m = io.read(_MAP_FILE)['MAP']
        assert 'LineType' in m
        assert len(m['LineType']) >= 1

    @skipmap
    def test_read_nodes(self):
        io = MAPIO()
        m = io.read(_MAP_FILE)['MAP']
        assert 'Node' in m
        assert len(m['Node']) > 0

    @skipmap
    def test_write_then_reread(self):
        io = MAPIO()
        data = io.read(_MAP_FILE)
        with tempfile.NamedTemporaryFile(suffix='.dat', mode='w', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            io.write(data, tmp_path)
            data2 = io.read(tmp_path)
            m1 = data['MAP']
            m2 = data2['MAP']
            assert len(m2['LineType']) == len(m1['LineType'])
            assert len(m2['Node']) == len(m1['Node'])
        finally:
            os.unlink(tmp_path)


# =================================================================
# Implements ModuleIO
# =================================================================
class TestModuleIOSubclass:
    def test_subdyn_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(SubDynIO, ModuleIO)

    def test_moordyn_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(MoorDynIO, ModuleIO)

    def test_map_implements_module_io(self):
        from openfast_io.io.base import ModuleIO
        assert issubclass(MAPIO, ModuleIO)
