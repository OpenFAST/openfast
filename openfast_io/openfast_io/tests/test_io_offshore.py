"""Edge-case IO tests for offshore modules.

Consolidated from the former per-module test files. Read/write/OutList
round-trip integration is covered by ``test_driver_openfast.py`` and the
external differential harness; only unique edge cases survive here.

Modules: HydroDyn, SeaState, SubDyn, MoorDyn, MAP, ExtPtfm.
"""
import os
import tempfile

import numpy as np
import pytest

from openfast_io.io.hydrodyn import HydroDynIO
from openfast_io.io.seastate import SeaStateIO
from openfast_io.io.subdyn import SubDynIO
from openfast_io.io.moordyn import MoorDynIO
from openfast_io.io.map_io import MAPIO
from openfast_io.io.extptfm import ExtPtfmIO

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

_OC3_DIR = os.path.join(_RTEST, '5MW_OC3Mnpl_DLL_WTurb_WavesIrr')
_SD_FILE = os.path.join(_OC3_DIR, 'NRELOffshrBsline5MW_OC3Monopile_SubDyn.dat')

_MD_FILE = os.path.join(_OC4_DIR, 'NRELOffshrBsline5MW_OC4DeepCwindSemi_MoorDyn.dat')

_MAP_DIR = os.path.join(_RTEST, '5MW_OC4Semi_WSt_WavesWN')
_MAP_FILE = os.path.join(_MAP_DIR, 'NRELOffshrBsline5MW_OC4DeepCwindSemi_MAP.dat')

_EP_DIR = os.path.join(_RTEST, '5MW_OC4Jckt_ExtPtfm')
_EP_FILE = os.path.join(_EP_DIR, 'ExtPtfm.dat')

_has_oc4 = os.path.isfile(_HD_FILE) and os.path.isfile(_SS_FILE)
_has_sd = os.path.isfile(_SD_FILE)
_has_md = os.path.isfile(_MD_FILE)
_has_map = os.path.isfile(_MAP_FILE)
_has_ep = os.path.isfile(_EP_FILE)

skipoc4 = pytest.mark.skipif(not _has_oc4, reason='r-test data missing')
skipsd = pytest.mark.skipif(not _has_sd, reason='SubDyn r-test data missing')
skipmd = pytest.mark.skipif(not _has_md, reason='MoorDyn r-test data missing')
skipmap = pytest.mark.skipif(not _has_map, reason='MAP r-test data missing')
skipep = pytest.mark.skipif(not _has_ep, reason='ExtPtfm r-test data missing')


# ===================================================================
# HydroDyn edge cases
# ===================================================================
@skipoc4
def test_hydrodyn_potfile_paths():
    io = HydroDynIO()
    hd = io.read(_HD_FILE)['HydroDyn']
    assert isinstance(hd['PotFile'], list)
    assert len(hd['PotFile']) == hd['NBody']


@skipoc4
def test_hydrodyn_arrays():
    io = HydroDynIO()
    hd = io.read(_HD_FILE)['HydroDyn']
    assert len(hd['WAMITULEN']) == hd['NBody']
    assert len(hd['PtfmRefxt']) == hd['NBody']


@skipoc4
def test_hydrodyn_matrices():
    io = HydroDynIO()
    hd = io.read(_HD_FILE)['HydroDyn']
    assert hd['AddCLin'].shape[0] == 6
    assert hd['AddBLin'].shape[0] == 6
    assert hd['AddBQuad'].shape[0] == 6


@skipoc4
def test_hydrodyn_joints():
    io = HydroDynIO()
    hd = io.read(_HD_FILE)['HydroDyn']
    n = hd['NJoints']
    assert len(hd['JointID']) == n
    assert len(hd['Jointxi']) == n


@skipoc4
def test_hydrodyn_members():
    io = HydroDynIO()
    hd = io.read(_HD_FILE)['HydroDyn']
    n = hd['NMembers']
    assert len(hd['MemberID']) == n
    assert len(hd['PropPot']) == n


# ===================================================================
# SeaState edge cases
# ===================================================================
@skipoc4
def test_seastate_wave_params():
    io = SeaStateIO()
    ss = io.read(_SS_FILE)['SeaState']
    assert isinstance(ss['WaveHs'], float)
    assert isinstance(ss['WaveTp'], float)
    assert isinstance(ss['WaveSeed1'], int)


@skipoc4
def test_seastate_current():
    io = SeaStateIO()
    ss = io.read(_SS_FILE)['SeaState']
    assert isinstance(ss['CurrMod'], int)


# ===================================================================
# SubDyn edge cases
# ===================================================================
@skipsd
def test_subdyn_joints():
    io = SubDynIO()
    sd = io.read(_SD_FILE)['SubDyn']
    assert len(sd['JointID']) == sd['NJoints']
    assert len(sd['JointXss']) == sd['NJoints']


@skipsd
def test_subdyn_members():
    io = SubDynIO()
    sd = io.read(_SD_FILE)['SubDyn']
    assert len(sd['MemberID']) == sd['NMembers']


@skipsd
def test_subdyn_prop_sets():
    io = SubDynIO()
    sd = io.read(_SD_FILE)['SubDyn']
    assert len(sd['PropSetID1']) == sd['NPropSetsCyl']


# ===================================================================
# MoorDyn edge cases
# ===================================================================
@skipmd
def test_moordyn_line_types():
    io = MoorDynIO()
    md = io.read(_MD_FILE)['MoorDyn']
    assert 'Name' in md
    assert len(md['Name']) >= 1


@skipmd
def test_moordyn_points():
    io = MoorDynIO()
    md = io.read(_MD_FILE)['MoorDyn']
    assert 'Point_ID' in md
    assert len(md['Point_ID']) > 0


@skipmd
def test_moordyn_lines():
    io = MoorDynIO()
    md = io.read(_MD_FILE)['MoorDyn']
    assert 'Line_ID' in md
    assert len(md['Line_ID']) > 0


# ===================================================================
# MAP edge cases
# ===================================================================
@skipmap
def test_map_line_types():
    io = MAPIO()
    m = io.read(_MAP_FILE)['MAP']
    assert 'LineType' in m
    assert len(m['LineType']) >= 1


@skipmap
def test_map_nodes():
    io = MAPIO()
    m = io.read(_MAP_FILE)['MAP']
    assert 'Node' in m
    assert len(m['Node']) > 0


# ===================================================================
# ExtPtfm edge cases
# ===================================================================
@skipep
def test_extptfm_superelement():
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
def test_extptfm_user_forcing():
    io = ExtPtfmIO()
    ep = io.read(_EP_FILE)['ExtPtfm']
    assert 'UserForcing' in ep
    uf = ep['UserForcing']
    assert uf['nSteps'] == 501
    # time + 31 DOF columns
    assert uf['ForceTimeSeries'].shape == (501, 32)
