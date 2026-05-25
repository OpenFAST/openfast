"""Tests for ServoDynIO (read / write ServoDyn + StC sub-files)."""
import os
import tempfile
import shutil

import pytest

from openfast_io.io.servodyn import ServoDynIO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(__file__)
_RTEST = os.path.join(_HERE, '..', '..', '..', 'reg_tests', 'r-test',
                      'glue-codes', 'openfast')

# 5MW_Land_DLL_WTurb – has ServoDyn but NO StC files
_5MW_DIR = os.path.join(_RTEST, '5MW_Land_DLL_WTurb')
_5MW_SD  = os.path.join(_5MW_DIR, 'NRELOffshrBsline5MW_Onshore_ServoDyn.dat')

# StC_test_OC4Semi – has ServoDyn WITH StC files
_STC_DIR = os.path.join(_RTEST, 'StC_test_OC4Semi')
_STC_SD  = os.path.join(_STC_DIR, 'ServoDyn_with_StC.dat')

_HAVE_5MW = os.path.isfile(_5MW_SD)
_HAVE_STC = os.path.isfile(_STC_SD)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
io = ServoDynIO()


# ---------------------------------------------------------------------------
# Tests – 5MW basic ServoDyn (no StC)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not _HAVE_5MW, reason='r-test data not available')
class Test5MWServoDyn:

    def test_read_returns_dict(self):
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        assert 'ServoDyn' in result
        sd = result['ServoDyn']
        assert isinstance(sd, dict)

    def test_read_pitch_control(self):
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'PCMode' in sd
        assert 'PitNeut(1)' in sd

    def test_read_generator_torque(self):
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'VSContrl' in sd
        assert 'GenEff' in sd

    def test_read_bladed_interface(self):
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'DLL_FileName' in sd
        assert 'DLL_NumTrq' in sd

    def test_read_output_section(self):
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        sd = result['ServoDyn']
        assert 'SumPrint' in sd
        assert 'OutFmt' in sd

    def test_stc_lists_empty(self):
        """5MW_Land has no StC files – lists should be empty."""
        result = io.read(_5MW_SD, base_dir=_5MW_DIR)
        assert result['BStC'] == []
        assert result['NStC'] == []
        assert result['TStC'] == []
        assert result['SStC'] == []


# ---------------------------------------------------------------------------
# Tests – StC_test_OC4Semi (ServoDyn with StC sub-files)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not _HAVE_STC, reason='r-test StC data not available')
class TestStCServoDyn:

    def test_read_stc_counts(self):
        result = io.read(_STC_SD, base_dir=_STC_DIR,
                         servo_file_rel='ServoDyn_with_StC.dat')
        sd = result['ServoDyn']
        # This test case should have StC files defined 
        # Verify counts match the loaded lists
        assert len(result['BStC']) == sd['NumBStC']
        assert len(result['NStC']) == sd['NumNStC']
        assert len(result['TStC']) == sd['NumTStC']
        assert len(result['SStC']) == sd['NumSStC']

    def test_read_stc_fields(self):
        result = io.read(_STC_SD, base_dir=_STC_DIR,
                         servo_file_rel='ServoDyn_with_StC.dat')
        # Check at least one StC list has entries and fields are present
        all_stc = result['BStC'] + result['NStC'] + result['TStC'] + result['SStC']
        assert len(all_stc) > 0, 'Expected at least one StC file in test case'
        stc = all_stc[0]
        assert 'StC_DOF_MODE' in stc
        assert 'StC_X_M' in stc
        assert 'SpringForceTable' in stc

    def test_write_and_reread(self):
        """Round-trip: read → write → reread → compare key fields."""
        result = io.read(_STC_SD, base_dir=_STC_DIR,
                         servo_file_rel='ServoDyn_with_StC.dat')

        with tempfile.TemporaryDirectory() as tmp:
            out_path = os.path.join(tmp, 'ServoDyn.dat')
            io.write(result, out_path, base_dir=tmp, run_dir=tmp)

            # Reread
            result2 = io.read(out_path, base_dir=tmp,
                              servo_file_rel='ServoDyn.dat')

        sd1 = result['ServoDyn']
        sd2 = result2['ServoDyn']

        # Spot-check a few fields
        for key in ['PCMode', 'VSContrl', 'GenEff', 'HSSBrMode',
                    'YCMode', 'AfCmode', 'CCmode', 'DLL_NumTrq',
                    'NumBStC', 'NumNStC', 'NumTStC', 'NumSStC']:
            assert sd1[key] == sd2[key], f'Mismatch on {key}: {sd1[key]} vs {sd2[key]}'

        # StC counts should match
        assert len(result['BStC']) == len(result2['BStC'])
        assert len(result['NStC']) == len(result2['NStC'])
        assert len(result['TStC']) == len(result2['TStC'])
        assert len(result['SStC']) == len(result2['SStC'])

        # If StC entries exist, spot-check one field
        all_stc1 = result['BStC'] + result['NStC'] + result['TStC'] + result['SStC']
        all_stc2 = result2['BStC'] + result2['NStC'] + result2['TStC'] + result2['SStC']
        for s1, s2 in zip(all_stc1, all_stc2):
            assert s1['StC_DOF_MODE'] == s2['StC_DOF_MODE']
            assert s1['StC_X_M'] == pytest.approx(s2['StC_X_M'])


# ---------------------------------------------------------------------------
# Tests – write only (unit test for writer with synthetic data)
# ---------------------------------------------------------------------------
class TestServoDynWrite:

    @staticmethod
    def _make_minimal_sd() -> dict:
        """Construct a minimal ServoDyn data dict for writing."""
        sd = {
            'Echo': False, 'DT': 0.005,
            'PCMode': 0, 'TPCOn': 0.0,
        }
        for idx in range(1, 4):
            sd[f'PitNeut({idx})'] = 0.0
            sd[f'PitSpr({idx})'] = 0.0
            sd[f'PitDamp({idx})'] = 0.0
            sd[f'TPitManS({idx})'] = 9999.9
            sd[f'PitManRat({idx})'] = 2.0
            sd[f'BlPitchF({idx})'] = 0.0
        sd.update({
            'VSContrl': 5, 'GenModel': 1, 'GenEff': 94.4,
            'GenTiStr': True, 'GenTiStp': True,
            'SpdGenOn': 0.0, 'TimGenOn': 0.0, 'TimGenOf': 9999.9,
            'VS_RtGnSp': 0.0, 'VS_RtTq': 0.0, 'VS_Rgn2K': 0.0, 'VS_SlPc': 0.0,
            'SIG_SlPc': 0.0, 'SIG_SySp': 0.0, 'SIG_RtTq': 0.0, 'SIG_PORt': 0.0,
            'TEC_Freq': 0.0, 'TEC_NPol': 0, 'TEC_SRes': 0.0, 'TEC_RRes': 0.0,
            'TEC_VLL': 0.0, 'TEC_SLR': 0.0, 'TEC_RLR': 0.0, 'TEC_MR': 0.0,
            'HSSBrMode': 0, 'THSSBrDp': 0.0, 'HSSBrDT': 0.0, 'HSSBrTqF': 0.0,
            'YCMode': 0, 'TYCOn': 0.0, 'YawNeut': 0.0, 'YawSpr': 0.0,
            'YawDamp': 0.0, 'TYawManS': 0.0, 'YawManRat': 0.0, 'NacYawF': 0.0,
            'AfCmode': 0, 'AfC_Mean': 0.0, 'AfC_Amp': 0.0, 'AfC_Phase': 0.0,
            'NumBStC': 0, 'BStCfiles': [],
            'NumNStC': 0, 'NStCfiles': [],
            'NumTStC': 0, 'TStCfiles': [],
            'NumSStC': 0, 'SStCfiles': [],
            'CCmode': 0,
            'DLL_FileName': 'libdiscon.so', 'DLL_InFile': 'DISCON.IN',
            'DLL_ProcName': 'DISCON', 'DLL_DT': 'default',
            'DLL_Ramp': False, 'BPCutoff': 0.0, 'NacYaw_North': 0.0,
            'Ptch_Cntrl': 1, 'Ptch_SetPnt': 0.0, 'Ptch_Min': 0.0,
            'Ptch_Max': 90.0, 'PtchRate_Min': -8.0, 'PtchRate_Max': 8.0,
            'Gain_OM': 0.0, 'GenSpd_MinOM': 0.0, 'GenSpd_MaxOM': 0.0,
            'GenSpd_Dem': 0.0, 'GenTrq_Dem': 0.0, 'GenPwr_Dem': 0.0,
            'DLL_NumTrq': 0, 'GenSpd_TLU': [], 'GenTrq_TLU': [],
            'SumPrint': False, 'OutFile': 1, 'TabDelim': True,
            'OutFmt': 'ES10.3E2', 'TStart': 0.0,
        })
        return {'ServoDyn': sd, 'BStC': [], 'NStC': [], 'TStC': [], 'SStC': []}

    def test_write_creates_file(self):
        data = self._make_minimal_sd()
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, 'ServoDyn.dat')
            io.write(data, out)
            assert os.path.isfile(out)

    def test_write_reread_synthetic(self):
        data = self._make_minimal_sd()
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, 'ServoDyn.dat')
            io.write(data, out)
            result = io.read(out, base_dir=tmp)
        sd_orig = data['ServoDyn']
        sd_read = result['ServoDyn']
        assert sd_read['PCMode'] == sd_orig['PCMode']
        assert sd_read['VSContrl'] == sd_orig['VSContrl']
        assert sd_read['GenEff'] == pytest.approx(sd_orig['GenEff'])
