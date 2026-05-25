"""
HydroDynIO -- read / write HydroDyn input files.

Produces ``{'HydroDyn': hd}``
"""
from __future__ import annotations

import os
from typing import Any, Dict, Optional, Callable

import numpy as np

from .base import ModuleIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
)


def _get_outlist(outlist_dict: dict, keys: list) -> list:
    out = []
    for key in keys:
        if key in outlist_dict:
            out.append(outlist_dict[key])
    return out


# ======================================================================
# HydroDynIO
# ======================================================================
class HydroDynIO(ModuleIO):
    """Read / write HydroDyn input files."""

    # ------------------------------------------------------------------
    def read(
        self,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        read_outlist_fn: Optional[Callable] = None,
        **kwargs,
    ) -> dict:
        hd: Dict[str, Any] = {}
        hd_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        f = open(hd_file)

        f.readline()
        f.readline()

        hd['Echo'] = bool_read(f.readline().split()[0])

        # FLOATING PLATFORM
        f.readline()
        hd['PotMod']      = int_read(f.readline().split()[0])
        hd['ExctnMod']    = int_read(f.readline().split()[0])
        hd['ExctnDisp']   = int_read(f.readline().split()[0])
        hd['ExctnCutOff'] = int_read(f.readline().split()[0])
        hd['PtfmYMod']    = int_read(f.readline().split()[0])
        hd['PtfmRefY']    = float_read(f.readline().split()[0])
        hd['PtfmYCutOff'] = float_read(f.readline().split()[0])
        hd['NExctnHdg']   = int_read(f.readline().split()[0])
        hd['RdtnMod']     = int_read(f.readline().split()[0])
        hd['RdtnTMax']    = float_read(f.readline().split()[0])
        hd['RdtnDT']      = float_read(f.readline().split()[0])
        hd['NBody']       = int_read(f.readline().split()[0])
        hd['NBodyMod']    = int_read(f.readline().split()[0])

        pot_strings = read_array(f, hd['NBody'], str)
        pot_strings = [os.path.normpath(os.path.join(os.path.split(hd_file)[0], ps)) for ps in pot_strings]
        hd['PotFile']      = pot_strings
        hd['WAMITULEN']    = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmRefxt']    = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmRefyt']    = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmRefzt']    = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmRefztRot'] = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmVol0']     = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmCOBxt']    = read_array(f, hd['NBody'], array_type=float)
        hd['PtfmCOByt']    = read_array(f, hd['NBody'], array_type=float)
        hd['NAddDOF']      = read_array(f, hd['NBody'], array_type=int)

        # 2ND-ORDER FLOATING PLATFORM FORCES
        f.readline()
        hd['MnDrift']   = int_read(f.readline().split()[0])
        hd['NewmanApp'] = int_read(f.readline().split()[0])
        hd['DiffQTF']   = int_read(f.readline().split()[0])
        hd['SumQTF']    = int_read(f.readline().split()[0])

        # PLATFORM ADDITIONAL STIFFNESS AND DAMPING
        f.readline()
        NBody = hd['NBody']
        if hd['NBodyMod'] == 1:
            hd['AddF0'] = [float(f.readline().strip().split()[0]) for _ in range(6 * NBody)]
        elif hd['NBodyMod'] > 1:
            hd['AddF0'] = [[float(idx) for idx in f.readline().strip().split()[:NBody]] for _ in range(6)]
        else:
            raise Exception("Invalid value for NBodyMod")

        _mat_rows = 6 * NBody if hd['NBodyMod'] == 1 else 6
        hd['AddCLin']  = np.array([[float(idx) for idx in f.readline().strip().split()[:6 * NBody]] for _ in range(_mat_rows)])
        hd['AddBLin']  = np.array([[float(idx) for idx in f.readline().strip().split()[:6 * NBody]] for _ in range(_mat_rows)])
        hd['AddBQuad'] = np.array([[float(idx) for idx in f.readline().strip().split()[:6 * NBody]] for _ in range(_mat_rows)])

        # STRIP THEORY OPTIONS
        f.readline()
        hd['WaveDisp'] = int_read(f.readline().split()[0])
        hd['AMMod']    = int_read(f.readline().split()[0])
        hd['HstMod']   = int_read(f.readline().split()[0])

        # AXIAL COEFFICIENTS
        f.readline()
        hd['NAxCoef'] = int_read(f.readline().split()[0])
        n = hd['NAxCoef']
        hd['AxCoefID']  = [None] * n
        hd['AxCd']      = [None] * n
        hd['AxCa']      = [None] * n
        hd['AxCp']      = [None] * n
        hd['AxFDMod']   = [None] * n
        hd['AxVnCOff']  = [None] * n
        hd['AxFDLoFSc'] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['AxCoefID'][i]  = int(ln[0])
            hd['AxCd'][i]      = float_read(ln[1])
            hd['AxCa'][i]      = float_read(ln[2])
            hd['AxCp'][i]      = float_read(ln[3])
            hd['AxFDMod'][i]   = float_read(ln[4])
            hd['AxVnCOff'][i]  = float_read(ln[5])
            hd['AxFDLoFSc'][i] = float_read(ln[6])

        # MEMBER JOINTS
        f.readline()
        hd['NJoints'] = int_read(f.readline().split()[0])
        n = hd['NJoints']
        hd['JointID']    = [None] * n
        hd['Jointxi']    = [None] * n
        hd['Jointyi']    = [None] * n
        hd['Jointzi']    = [None] * n
        hd['JointAxID']  = [None] * n
        hd['JointOvrlp'] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['JointID'][i]    = int(ln[0])
            hd['Jointxi'][i]    = float(ln[1])
            hd['Jointyi'][i]    = float(ln[2])
            hd['Jointzi'][i]    = float(ln[3])
            hd['JointAxID'][i]  = int(ln[4])
            hd['JointOvrlp'][i] = int(ln[5])

        # CIRCULAR MEMBER CROSS-SECTION PROPERTIES
        f.readline()
        hd['NPropSetsCyl'] = int_read(f.readline().split()[0])
        n = hd['NPropSetsCyl']
        hd['CylPropSetID'] = [None] * n
        hd['CylPropD']     = [None] * n
        hd['CylPropThck']  = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['CylPropSetID'][i] = int(ln[0])
            hd['CylPropD'][i]     = float(ln[1])
            hd['CylPropThck'][i]  = float(ln[2])

        # RECTANGULAR MEMBER CROSS-SECTION PROPERTIES
        f.readline()
        hd['NPropSetsRec'] = int_read(f.readline().split()[0])
        n = hd['NPropSetsRec']
        hd['RecPropSetID'] = [None] * n
        hd['RecPropA']     = [None] * n
        hd['RecPropB']     = [None] * n
        hd['RecPropThck']  = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['RecPropSetID'][i] = int(ln[0])
            hd['RecPropA'][i]     = float(ln[1])
            hd['RecPropB'][i]     = float(ln[2])
            hd['RecPropThck'][i]  = float(ln[3])

        # SIMPLE CIRCULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline(); f.readline(); f.readline()
        ln = f.readline().split()
        for j, key in enumerate(['CylSimplCd', 'CylSimplCdMG', 'CylSimplCa', 'CylSimplCaMG',
                                  'CylSimplCp', 'CylSimplCpMG', 'CylSimplAxCd', 'CylSimplAxCdMG',
                                  'CylSimplAxCa', 'CylSimplAxCaMG', 'CylSimplAxCp', 'CylSimplAxCpMG',
                                  'CylSimplCb', 'CylSimplCbMG']):
            hd[key] = float_read(ln[j])

        # SIMPLE RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline(); f.readline(); f.readline()
        ln = f.readline().split()
        for j, key in enumerate(['RecSimplCdA', 'RecSimplCdAMG', 'RecSimplCdB', 'RecSimplCdBMG',
                                  'RecSimplCaA', 'RecSimplCaAMG', 'RecSimplCaB', 'RecSimplCaBMG',
                                  'RecSimplCp', 'RecSimplCpMG', 'RecSimplAxCd', 'RecSimplAxCdMG',
                                  'RecSimplAxCa', 'RecSimplAxCaMG', 'RecSimplAxCp', 'RecSimplAxCpMG',
                                  'RecSimplCb', 'RecSimplCbMG']):
            hd[key] = float_read(ln[j])

        # DEPTH-BASED CIRCULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline()
        hd['NCoefDpthCyl'] = int_read(f.readline().split()[0])
        n = hd['NCoefDpthCyl']
        cyl_dpth_keys = ['CylDpth', 'CylDpthCd', 'CylDpthCdMG', 'CylDpthCa', 'CylDpthCaMG',
                         'CylDpthCp', 'CylDpthCpMG', 'CylDpthAxCd', 'CylDpthAxCdMG',
                         'CylDpthAxCa', 'CylDpthAxCaMG', 'CylDpthAxCp', 'CylDpthAxCpMG',
                         'CylDpthCb', 'CylDpthCbMG']
        for k in cyl_dpth_keys:
            hd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            for j, k in enumerate(cyl_dpth_keys):
                hd[k][i] = float_read(ln[j])

        # DEPTH-BASED RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline()
        hd['NCoefDpthRec'] = int_read(f.readline().split()[0])
        n = hd['NCoefDpthRec']
        rec_dpth_keys = ['RecDpth', 'RecDpthCdA', 'RecDpthCdAMG', 'RecDpthCdB', 'RecDpthCdBMG',
                         'RecDpthCaA', 'RecDpthCaAMG', 'RecDpthCaB', 'RecDpthCaBMG',
                         'RecDpthCp', 'RecDpthCpMG', 'RecDpthAxCd', 'RecDpthAxCdMG',
                         'RecDpthAxCa', 'RecDpthAxCaMG', 'RecDpthAxCp', 'RecDpthAxCpMG',
                         'RecDpthCb', 'RecDpthCbMG']
        for k in rec_dpth_keys:
            hd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            for j, k in enumerate(rec_dpth_keys):
                hd[k][i] = float_read(ln[j])

        # MEMBER-BASED CIRCULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline()
        hd['NCoefMembersCyl'] = int_read(f.readline().split()[0])
        n = hd['NCoefMembersCyl']
        cyl_mem_keys = ['MemberID_HydCCyl',
                        'CylMemberCd1', 'CylMemberCd2', 'CylMemberCdMG1', 'CylMemberCdMG2',
                        'CylMemberCa1', 'CylMemberCa2', 'CylMemberCaMG1', 'CylMemberCaMG2',
                        'CylMemberCp1', 'CylMemberCp2', 'CylMemberCpMG1', 'CylMemberCpMG2',
                        'CylMemberAxCd1', 'CylMemberAxCd2', 'CylMemberAxCdMG1', 'CylMemberAxCdMG2',
                        'CylMemberAxCa1', 'CylMemberAxCa2', 'CylMemberAxCaMG1', 'CylMemberAxCaMG2',
                        'CylMemberAxCp1', 'CylMemberAxCp2', 'CylMemberAxCpMG1', 'CylMemberAxCpMG2',
                        'CylMemberCb1', 'CylMemberCb2', 'CylMemberCbMG1', 'CylMemberCbMG2']
        for k in cyl_mem_keys:
            hd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['MemberID_HydCCyl'][i] = int(ln[0])
            for j, k in enumerate(cyl_mem_keys[1:], 1):
                hd[k][i] = float_read(ln[j])

        # MEMBER-BASED RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS
        f.readline()
        hd['NCoefMembersRec'] = int_read(f.readline().split()[0])
        n = hd['NCoefMembersRec']
        rec_mem_keys = ['MemberID_HydCRec',
                        'RecMemberCdA1', 'RecMemberCdA2', 'RecMemberCdAMG1', 'RecMemberCdAMG2',
                        'RecMemberCdB1', 'RecMemberCdB2', 'RecMemberCdBMG1', 'RecMemberCdBMG2',
                        'RecMemberCaA1', 'RecMemberCaA2', 'RecMemberCaAMG1', 'RecMemberCaAMG2',
                        'RecMemberCaB1', 'RecMemberCaB2', 'RecMemberCaBMG1', 'RecMemberCaBMG2',
                        'RecMemberCp1', 'RecMemberCp2', 'RecMemberCpMG1', 'RecMemberCpMG2',
                        'RecMemberAxCd1', 'RecMemberAxCd2', 'RecMemberAxCdMG1', 'RecMemberAxCdMG2',
                        'RecMemberAxCa1', 'RecMemberAxCa2', 'RecMemberAxCaMG1', 'RecMemberAxCaMG2',
                        'RecMemberAxCp1', 'RecMemberAxCp2', 'RecMemberAxCpMG1', 'RecMemberAxCpMG2',
                        'RecMemberCb1', 'RecMemberCb2', 'RecMemberCbMG1', 'RecMemberCbMG2']
        for k in rec_mem_keys:
            hd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['MemberID_HydCRec'][i] = int(ln[0])
            for j, k in enumerate(rec_mem_keys[1:], 1):
                hd[k][i] = float_read(ln[j])

        # MEMBERS
        f.readline()
        hd['NMembers'] = int_read(f.readline().split()[0])
        n = hd['NMembers']
        mem_keys = ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2',
                    'MSecGeom', 'MSpinOrient', 'MDivSize', 'MCoefMod', 'MHstLMod', 'PropPot']
        for k in mem_keys:
            hd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['MemberID'][i]    = int(ln[0])
            hd['MJointID1'][i]   = int(ln[1])
            hd['MJointID2'][i]   = int(ln[2])
            hd['MPropSetID1'][i] = int(ln[3])
            hd['MPropSetID2'][i] = int(ln[4])
            hd['MSecGeom'][i]    = int(ln[5])
            hd['MSpinOrient'][i] = float(ln[6])
            hd['MDivSize'][i]    = float(ln[7])
            hd['MCoefMod'][i]    = int(ln[8])
            hd['MHstLMod'][i]    = int(ln[9])
            hd['PropPot'][i]     = bool_read(ln[10])

        # FILLED MEMBERS
        f.readline()
        hd['NFillGroups'] = int_read(f.readline().split()[0])
        n = hd['NFillGroups']
        hd['FillNumM']  = [None] * n
        hd['FillMList'] = [None] * n
        hd['FillFSLoc'] = [None] * n
        hd['FillDens']  = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            n_fill = int(ln[0])
            hd['FillNumM'][i]  = n_fill
            hd['FillMList'][i] = [int(j) for j in ln[1:1 + n_fill]]
            hd['FillFSLoc'][i] = float(ln[n_fill + 1])
            if ln[n_fill + 2] == 'DEFAULT':
                hd['FillDens'][i] = 'DEFAULT'
            else:
                hd['FillDens'][i] = float(ln[n_fill + 2])

        # MARINE GROWTH
        f.readline()
        hd['NMGDepths'] = int_read(f.readline().split()[0])
        n = hd['NMGDepths']
        hd['MGDpth'] = [None] * n
        hd['MGThck'] = [None] * n
        hd['MGDens'] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['MGDpth'][i] = float(ln[0])
            hd['MGThck'][i] = float(ln[1])
            hd['MGDens'][i] = float(ln[2])

        # MEMBER OUTPUT LIST
        f.readline()
        hd['NMOutputs'] = int_read(f.readline().split()[0])
        n = hd['NMOutputs']
        hd['MemberID_out'] = [None] * n
        hd['NOutLoc']      = [None] * n
        hd['NodeLocs']     = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            hd['MemberID_out'][i] = int(ln[0])
            hd['NOutLoc'][i]      = int(ln[1])
            hd['NodeLocs'][i]     = float(ln[2])

        # JOINT OUTPUT LIST
        f.readline()
        hd['NJOutputs'] = int_read(f.readline().split()[0])
        if int(hd['NJOutputs']) > 0:
            hd['JOutLst'] = [int(idx.strip()) for idx in f.readline().split('JOutLst')[0].split(',')]
        else:
            f.readline()
            hd['JOutLst'] = [0]

        # OUTPUT
        f.readline()
        hd['HDSum']    = bool_read(f.readline().split()[0])
        hd['OutAll']   = bool_read(f.readline().split()[0])
        hd['OutSwtch'] = int_read(f.readline().split()[0])
        hd['OutFmt']   = quoted_read(f.readline().split()[0])
        hd['OutSFmt']  = quoted_read(f.readline().split()[0])

        # Outlist
        f.readline()
        if read_outlist_fn is not None:
            read_outlist_fn(f, 'HydroDyn')

        f.close()
        return {'HydroDyn': hd}

    # ------------------------------------------------------------------
    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        **kwargs,
    ) -> None:
        hd = data['HydroDyn']

        with open(file_path, 'w') as f:
            f.write('------- HydroDyn Input File --------------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('{!s:<22} {:<11} {:}'.format(hd['Echo'], 'Echo', '- Echo the input file data (flag)\n'))
            f.write('---------------------- FLOATING PLATFORM --------------------------------------- [unused with WaveMod=6]\n')
            f.write('{:<22d} {:<11} {:}'.format(hd['PotMod'], 'PotMod', '- Potential-flow model {0: none=no potential flow, 1: frequency-to-time-domain transforms based on WAMIT output, 2: fluid-impulse theory (FIT)} (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['ExctnMod'], 'ExctnMod', '- Wave-excitation model (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['ExctnDisp'], 'ExctnDisp', '- Method of computing Wave Excitation (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['ExctnCutOff'], 'ExctnCutOff', '- Cutoff frequency (Hz)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['PtfmYMod'], 'PtfmYMod', '- Model for large platform yaw offset (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['PtfmRefY'], 'PtfmRefY', '- Platform reference yaw offset (deg)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['PtfmYCutOff'], 'PtfmYCutOff', '- Cutoff frequency for PRP yaw filtering (Hz)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['NExctnHdg'], 'NExctnHdg', '- Number of platform yaw/heading angles (-)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['RdtnMod'], 'RdtnMod', '- Radiation memory-effect model (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['RdtnTMax'], 'RdtnTMax', '- Analysis time for wave radiation kernel calculations (sec)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['RdtnDT'], 'RdtnDT', '- Time step for wave radiation kernel calculations (sec)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['NBody'], 'NBody', '- Number of WAMIT bodies (-)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['NBodyMod'], 'NBodyMod', '- Body coupling model (switch)\n'))
            pot = hd['PotFile']
            f.write('{:<22} {:<11} {:}'.format('"{}"'.format('", "'.join(pot) if isinstance(pot, list) else pot), 'PotFile', '- Root name of potential-flow model data\n'))

            for arr_key, label, desc in [
                ('WAMITULEN', 'WAMITULEN', '- Characteristic body length scale (m)'),
                ('PtfmRefxt', 'PtfmRefxt', '- xt offset of body reference point (m)'),
                ('PtfmRefyt', 'PtfmRefyt', '- yt offset of body reference point (m)'),
                ('PtfmRefzt', 'PtfmRefzt', '- zt offset of body reference point (m)'),
                ('PtfmRefztRot', 'PtfmRefztRot', '- Rotation about zt (deg)'),
                ('PtfmVol0', 'PtfmVol0', '- Displaced volume (m^3)'),
                ('PtfmCOBxt', 'PtfmCOBxt', '- xt offset of COB (m)'),
                ('PtfmCOByt', 'PtfmCOByt', '- yt offset of COB (m)'),
                ('NAddDOF', 'NAddDOF', '- Number of additional DOF (-)'),
            ]:
                f.write('{:<22} {:<11} {:}'.format(', '.join([f'{val}' for val in hd[arr_key]]), label, desc + '\n'))

            f.write('---------------------- 2ND-ORDER FLOATING PLATFORM FORCES ----------------------\n')
            f.write('{:<22} {:<11} {:}'.format(hd['MnDrift'], 'MnDrift', '- Mean-drift 2nd-order forces\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['NewmanApp'], 'NewmanApp', '- Newman approximation\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['DiffQTF'], 'DiffQTF', '- Full difference-frequency QTF\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['SumQTF'], 'SumQTF', '- Full summation-frequency QTF\n'))

            f.write('---------------------- PLATFORM ADDITIONAL STIFFNESS AND DAMPING  --------------\n')
            _n_f0 = 6 * hd['NBody'] if hd['NBodyMod'] == 1 else 6
            for j in range(_n_f0):
                val = hd['AddF0'][j]
                if isinstance(val, float):
                    ln = '{:14}   '.format(val)
                elif isinstance(val, (list, np.ndarray)):
                    ln = '{:14}   '.format(' '.join([f'{v}' for v in val]))
                else:
                    ln = '{:14}   '.format(val)
                if j == 0:
                    ln += 'AddF0    - Additional preload (N, N-m)\n'
                else:
                    ln += '\n'
                f.write(ln)
            for mat_key, mat_label in [('AddCLin', 'AddCLin'), ('AddBLin', 'AddBLin'), ('AddBQuad', 'AddBQuad')]:
                n_mat_rows = len(hd[mat_key])
                for j in range(n_mat_rows):
                    try:
                        ln = " ".join(['{:14}'.format(i) for i in hd[mat_key][j, :]])
                    except Exception:
                        ln = " ".join(['{:14}'.format(i) for i in hd[mat_key][j]])
                    if j == 0:
                        ln += f"   {mat_label}\n"
                    else:
                        ln += "\n"
                    f.write(ln)

            f.write('---------------------- STRIP THEORY OPTIONS --------------------------------------\n')
            f.write('{:<22d} {:<11} {:}'.format(hd['WaveDisp'], 'WaveDisp', '- Wave kinematics method (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['AMMod'], 'AMMod', '- Added-mass force method (switch)\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['HstMod'], 'HstMod', '- Hydrostatic loads method (switch)\n'))

            # AXIAL COEFFICIENTS
            f.write('---------------------- AXIAL COEFFICIENTS --------------------------------------\n')
            f.write('{:<22d} {:<11} {:}'.format(hd['NAxCoef'], 'NAxCoef', '- Number of axial coefficients (-)\n'))
            ax_hdrs = ['AxCoefID', 'AxCd', 'AxCa', 'AxCp', 'AxFDMod', 'AxVnCOff', 'AxFDLoFSc']
            f.write(" ".join(['{:^11s}'.format(h) for h in ax_hdrs]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in ax_hdrs]) + '\n')
            for i in range(hd['NAxCoef']):
                ln = ['{:^11d}'.format(hd['AxCoefID'][i])]
                for k in ax_hdrs[1:]:
                    ln.append('{:^11}'.format(hd[k][i]))
                f.write(" ".join(ln) + '\n')

            # MEMBER JOINTS
            f.write('---------------------- MEMBER JOINTS -------------------------------------------\n')
            f.write('{:<22d} {:<11} {:}'.format(hd['NJoints'], 'NJoints', '- Number of joints (-)\n'))
            jt_hdrs = ['JointID', 'Jointxi', 'Jointyi', 'Jointzi', 'JointAxID', 'JointOvrlp']
            f.write(" ".join(['{:^11s}'.format(h) for h in jt_hdrs]) + '\n')
            f.write(" ".join(['{:^11s}'.format(u) for u in ['(-)', '(m)', '(m)', '(m)', '(-)', '(switch)']]) + '\n')
            for i in range(hd['NJoints']):
                ln = ['{:^11d}'.format(hd['JointID'][i]),
                      '{:^11}'.format(hd['Jointxi'][i]), '{:^11}'.format(hd['Jointyi'][i]),
                      '{:^11}'.format(hd['Jointzi'][i]),
                      '{:^11d}'.format(hd['JointAxID'][i]), '{:^11d}'.format(hd['JointOvrlp'][i])]
                f.write(" ".join(ln) + '\n')

            # CYLINDRICAL MEMBER CROSS-SECTION PROPERTIES
            f.write('---------------------- CYLINDRICAL MEMBER CROSS-SECTION PROPERTIES -------------------------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NPropSetsCyl'], 'NPropSetsCyl', '- Number of cylindrical member property sets (-)\n'))
            f.write(" ".join(['{:^11s}'.format(h) for h in ['PropSetID', 'PropD', 'PropThck']]) + '\n')
            f.write(" ".join(['{:^11s}'.format(h) for h in ['(-)', '(m)', '(m)']]) + '\n')
            for i in range(hd['NPropSetsCyl']):
                f.write(" ".join(['{:^11d}'.format(hd['CylPropSetID'][i]),
                                  '{:^11}'.format(hd['CylPropD'][i]),
                                  '{:^11}'.format(hd['CylPropThck'][i])]) + '\n')

            # RECTANGULAR MEMBER CROSS-SECTION PROPERTIES
            f.write('---------------------- RECTANGULAR MEMBER CROSS-SECTION PROPERTIES -------------------------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NPropSetsRec'], 'NPropSetsRec', '- Number of rectangular member property sets (-)\n'))
            f.write(" ".join(['{:^11s}'.format(h) for h in ['PropSetID', 'PropA', 'PropB', 'PropThck']]) + '\n')
            f.write(" ".join(['{:^11s}'.format(h) for h in ['(-)', '(m)', '(m)', '(m)']]) + '\n')
            for i in range(hd['NPropSetsRec']):
                f.write(" ".join(['{:^11d}'.format(hd['RecPropSetID'][i]),
                                  '{:^11}'.format(hd['RecPropA'][i]),
                                  '{:^11}'.format(hd['RecPropB'][i]),
                                  '{:^11}'.format(hd['RecPropThck'][i])]) + '\n')

            # SIMPLE CYLINDRICAL HYDRO COEFFICIENTS
            f.write('---------------------- SIMPLE CYLINDRICAL-MEMBER HYDRODYNAMIC COEFFICIENTS (model 1) --------------\n')
            cyl_simpl_keys = ['CylSimplCd', 'CylSimplCdMG', 'CylSimplCa', 'CylSimplCaMG', 'CylSimplCp', 'CylSimplCpMG',
                              'CylSimplAxCd', 'CylSimplAxCdMG', 'CylSimplAxCa', 'CylSimplAxCaMG',
                              'CylSimplAxCp', 'CylSimplAxCpMG', 'CylSimplCb', 'CylSimplCbMG']
            hdr_names = ['SimplCd', 'SimplCdMG', 'SimplCa', 'SimplCaMG', 'SimplCp', 'SimplCpMG',
                         'SimplAxCd', 'SimplAxCdMG', 'SimplAxCa', 'SimplAxCaMG', 'SimplAxCp', 'SimplAxCpMG', 'SimplCb', 'SimplCbMG']
            f.write(" ".join(['{:^11s}'.format(h) for h in hdr_names]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in hdr_names]) + '\n')
            f.write(" ".join(['{:^11}'.format(hd[k]) for k in cyl_simpl_keys]) + '\n')

            # SIMPLE RECTANGULAR HYDRO COEFFICIENTS
            f.write('---------------------- SIMPLE RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS (model 1) --------------\n')
            rec_simpl_keys = ['RecSimplCdA', 'RecSimplCdAMG', 'RecSimplCdB', 'RecSimplCdBMG',
                              'RecSimplCaA', 'RecSimplCaAMG', 'RecSimplCaB', 'RecSimplCaBMG',
                              'RecSimplCp', 'RecSimplCpMG', 'RecSimplAxCd', 'RecSimplAxCdMG',
                              'RecSimplAxCa', 'RecSimplAxCaMG', 'RecSimplAxCp', 'RecSimplAxCpMG',
                              'RecSimplCb', 'RecSimplCbMG']
            rec_hdr = ['SimplCdA', 'SimplCdAMG', 'SimplCdB', 'SimplCdBMG', 'SimplCaA', 'SimplCaAMG',
                       'SimplCaB', 'SimplCaBMG', 'SimplCp', 'SimplCpMG', 'SimplAxCd', 'SimplAxCdMG',
                       'SimplAxCa', 'SimplAxCaMG', 'SimplAxCp', 'SimplAxCpMG', 'SimplCb', 'SimplCbMG']
            f.write(" ".join(['{:^11s}'.format(h) for h in rec_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in rec_hdr]) + '\n')
            f.write(" ".join(['{:^11}'.format(hd[k]) for k in rec_simpl_keys]) + '\n')

            # DEPTH-BASED CYL COEFFICIENTS
            f.write('---------------------- DEPTH-BASED CYLINDRICAL-MEMBER HYDRODYNAMIC COEFFICIENTS (model 2) ---------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NCoefDpthCyl'], 'NCoefDpthCyl', '- Number of depth-dependent cylindrical-member coefficients (-)\n'))
            cyl_dpth_hdr = ['Dpth', 'DpthCd', 'DpthCdMG', 'DpthCa', 'DpthCaMG', 'DpthCp', 'DpthCpMG',
                            'DpthAxCd', 'DpthAxCdMG', 'DpthAxCa', 'DpthAxCaMG', 'DpthAxCp', 'DpthAxCpMG', 'DpthCb', 'DpthCbMG']
            f.write(" ".join(['{:^11s}'.format(h) for h in cyl_dpth_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in cyl_dpth_hdr]) + '\n')
            cyl_dpth_keys_list = ['CylDpth', 'CylDpthCd', 'CylDpthCdMG', 'CylDpthCa', 'CylDpthCaMG',
                                  'CylDpthCp', 'CylDpthCpMG', 'CylDpthAxCd', 'CylDpthAxCdMG',
                                  'CylDpthAxCa', 'CylDpthAxCaMG', 'CylDpthAxCp', 'CylDpthAxCpMG',
                                  'CylDpthCb', 'CylDpthCbMG']
            for i in range(hd['NCoefDpthCyl']):
                f.write(" ".join(['{:^11}'.format(hd[k][i]) for k in cyl_dpth_keys_list]) + '\n')

            # DEPTH-BASED REC COEFFICIENTS
            f.write('---------------------- DEPTH-BASED RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS (model 2) ---------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NCoefDpthRec'], 'NCoefDpthRec', '- Number of depth-dependent rectangular-member coefficients (-)\n'))
            rec_dpth_hdr = ['Dpth', 'DpthCdA', 'DpthCdAMG', 'DpthCdB', 'DpthCdBMG',
                            'DpthCaA', 'DpthCaAMG', 'DpthCaB', 'DpthCaBMG', 'DpthCp', 'DpthCpMG',
                            'DpthAxCd', 'DpthAxCdMG', 'DpthAxCa', 'DpthAxCaMG', 'DpthAxCp', 'DpthAxCpMG', 'DpthCb', 'DpthCbMG']
            f.write(" ".join(['{:^11s}'.format(h) for h in rec_dpth_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in rec_dpth_hdr]) + '\n')
            rec_dpth_keys_list = ['RecDpth', 'RecDpthCdA', 'RecDpthCdAMG', 'RecDpthCdB', 'RecDpthCdBMG',
                                  'RecDpthCaA', 'RecDpthCaAMG', 'RecDpthCaB', 'RecDpthCaBMG',
                                  'RecDpthCp', 'RecDpthCpMG', 'RecDpthAxCd', 'RecDpthAxCdMG',
                                  'RecDpthAxCa', 'RecDpthAxCaMG', 'RecDpthAxCp', 'RecDpthAxCpMG',
                                  'RecDpthCb', 'RecDpthCbMG']
            for i in range(hd['NCoefDpthRec']):
                f.write(" ".join(['{:^11}'.format(hd[k][i]) for k in rec_dpth_keys_list]) + '\n')

            # MEMBER-BASED CYL COEFFICIENTS
            f.write('---------------------- MEMBER-BASED CYLINDRICAL-MEMBER HYDRODYNAMIC COEFFICIENTS (model 3) --------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NCoefMembersCyl'], 'NCoefMembersCyl', '- Number of member-based cylindrical-member coefficients (-)\n'))
            cyl_mem_hdr = ['MemberID', 'MemberCd1', 'MemberCd2', 'MemberCdMG1', 'MemberCdMG2',
                           'MemberCa1', 'MemberCa2', 'MemberCaMG1', 'MemberCaMG2',
                           'MemberCp1', 'MemberCp2', 'MemberCpMG1', 'MemberCpMG2',
                           'MemberAxCd1', 'MemberAxCd2', 'MemberAxCdMG1', 'MemberAxCdMG2',
                           'MemberAxCa1', 'MemberAxCa2', 'MemberAxCaMG1', 'MemberAxCaMG2',
                           'MemberAxCp1', 'MemberAxCp2', 'MemberAxCpMG1', 'MemberAxCpMG2',
                           'MemberCb1', 'MemberCb2', 'MemberCbMG1', 'MemberCbMG2']
            f.write(" ".join(['{:^11s}'.format(h) for h in cyl_mem_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in cyl_mem_hdr]) + '\n')
            cyl_mem_w_keys = ['MemberID_HydCCyl',
                        'CylMemberCd1', 'CylMemberCd2', 'CylMemberCdMG1', 'CylMemberCdMG2',
                        'CylMemberCa1', 'CylMemberCa2', 'CylMemberCaMG1', 'CylMemberCaMG2',
                        'CylMemberCp1', 'CylMemberCp2', 'CylMemberCpMG1', 'CylMemberCpMG2',
                        'CylMemberAxCd1', 'CylMemberAxCd2', 'CylMemberAxCdMG1', 'CylMemberAxCdMG2',
                        'CylMemberAxCa1', 'CylMemberAxCa2', 'CylMemberAxCaMG1', 'CylMemberAxCaMG2',
                        'CylMemberAxCp1', 'CylMemberAxCp2', 'CylMemberAxCpMG1', 'CylMemberAxCpMG2',
                        'CylMemberCb1', 'CylMemberCb2', 'CylMemberCbMG1', 'CylMemberCbMG2']
            for i in range(hd['NCoefMembersCyl']):
                ln = ['{:^11d}'.format(hd['MemberID_HydCCyl'][i])]
                for k in cyl_mem_w_keys[1:]:
                    ln.append('{:^11}'.format(hd[k][i]))
                f.write(" ".join(ln) + '\n')

            # MEMBER-BASED REC COEFFICIENTS
            f.write('---------------------- MEMBER-BASED RECTANGULAR-MEMBER HYDRODYNAMIC COEFFICIENTS (model 3) --------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NCoefMembersRec'], 'NCoefMembersRec', '- Number of member-based rectangular-member coefficients (-)\n'))
            rec_mem_hdr = ['MemberID', 'MemberCdA1', 'MemberCdA2', 'MemberCdAMG1', 'MemberCdAMG2',
                           'MemberCdB1', 'MemberCdB2', 'MemberCdBMG1', 'MemberCdBMG2',
                           'MemberCaA1', 'MemberCaA2', 'MemberCaAMG1', 'MemberCaAMG2',
                           'MemberCaB1', 'MemberCaB2', 'MemberCaBMG1', 'MemberCaBMG2',
                           'MemberCp1', 'MemberCp2', 'MemberCpMG1', 'MemberCpMG2',
                           'MemberAxCd1', 'MemberAxCd2', 'MemberAxCdMG1', 'MemberAxCdMG2',
                           'MemberAxCa1', 'MemberAxCa2', 'MemberAxCaMG1', 'MemberAxCaMG2',
                           'MemberAxCp1', 'MemberAxCp2', 'MemberAxCpMG1', 'MemberAxCpMG2',
                           'MemberCb1', 'MemberCb2', 'MemberCbMG1', 'MemberCbMG2']
            f.write(" ".join(['{:^11s}'.format(h) for h in rec_mem_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in rec_mem_hdr]) + '\n')
            rec_mem_w_keys = ['MemberID_HydCRec',
                        'RecMemberCdA1', 'RecMemberCdA2', 'RecMemberCdAMG1', 'RecMemberCdAMG2',
                        'RecMemberCdB1', 'RecMemberCdB2', 'RecMemberCdBMG1', 'RecMemberCdBMG2',
                        'RecMemberCaA1', 'RecMemberCaA2', 'RecMemberCaAMG1', 'RecMemberCaAMG2',
                        'RecMemberCaB1', 'RecMemberCaB2', 'RecMemberCaBMG1', 'RecMemberCaBMG2',
                        'RecMemberCp1', 'RecMemberCp2', 'RecMemberCpMG1', 'RecMemberCpMG2',
                        'RecMemberAxCd1', 'RecMemberAxCd2', 'RecMemberAxCdMG1', 'RecMemberAxCdMG2',
                        'RecMemberAxCa1', 'RecMemberAxCa2', 'RecMemberAxCaMG1', 'RecMemberAxCaMG2',
                        'RecMemberAxCp1', 'RecMemberAxCp2', 'RecMemberAxCpMG1', 'RecMemberAxCpMG2',
                        'RecMemberCb1', 'RecMemberCb2', 'RecMemberCbMG1', 'RecMemberCbMG2']
            for i in range(hd['NCoefMembersRec']):
                ln = ['{:^11d}'.format(hd['MemberID_HydCRec'][i])]
                for k in rec_mem_w_keys[1:]:
                    ln.append('{:^11}'.format(hd[k][i]))
                f.write(" ".join(ln) + '\n')

            # MEMBERS
            f.write('-------------------- MEMBERS -------------------------------------------------\n')
            f.write('{:<11d} {:<11} {:}'.format(hd['NMembers'], 'NMembers', '- Number of members (-)\n'))
            mem_hdr = ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2',
                       'MSecGeom', 'MSpinOrient', 'MDivSize', 'MCoefMod', 'MHstLMod', 'PropPot']
            f.write(" ".join(['{:^11s}'.format(h) for h in mem_hdr]) + '\n')
            f.write(" ".join(['{:^11s}'.format(u) for u in ['(-)', '(-)', '(-)', '(-)', '(-)', '(switch)', '(deg)', '(m)', '(switch)', '(switch)', '(flag)']]) + '\n')
            for i in range(hd['NMembers']):
                ln = ['{:^11d}'.format(hd['MemberID'][i]),
                      '{:^11d}'.format(hd['MJointID1'][i]), '{:^11d}'.format(hd['MJointID2'][i]),
                      '{:^11d}'.format(hd['MPropSetID1'][i]), '{:^11d}'.format(hd['MPropSetID2'][i]),
                      '{:^11d}'.format(hd['MSecGeom'][i]), '{:^11}'.format(hd['MSpinOrient'][i]),
                      '{:^11}'.format(hd['MDivSize'][i]),
                      '{:^11d}'.format(hd['MCoefMod'][i]), '{:^11d}'.format(hd['MHstLMod'][i]),
                      '{!s:^11}'.format(hd['PropPot'][i])]
                f.write(" ".join(ln) + '\n')

            # FILLED MEMBERS
            f.write("---------------------- FILLED MEMBERS ------------------------------------------\n")
            f.write('{:<11d} {:<11} {:}'.format(hd['NFillGroups'], 'NFillGroups', '- Number of filled member groups (-)\n'))
            f.write(" ".join(['{:^11s}'.format(h) for h in ['FillNumM', 'FillMList', 'FillFSLoc', 'FillDens']]) + '\n')
            f.write(" ".join(['{:^11s}'.format(h) for h in ['(-)', '(-)', '(m)', '(kg/m^3)']]) + '\n')
            for i in range(hd['NFillGroups']):
                ln = ['{:^11d}'.format(hd['FillNumM'][i]),
                      " ".join(['%d' % j for j in hd['FillMList'][i]]),
                      '{:^11}'.format(hd['FillFSLoc'][i]),
                      '{:^11}'.format(hd['FillDens'][i])]
                f.write(" ".join(ln) + '\n')

            # MARINE GROWTH
            f.write("---------------------- MARINE GROWTH -------------------------------------------\n")
            f.write('{:<11d} {:<11} {:}'.format(hd['NMGDepths'], 'NMGDepths', '- Number of marine-growth depths specified (-)\n'))
            f.write(" ".join(['{:^11s}'.format(h) for h in ['MGDpth', 'MGThck', 'MGDens']]) + '\n')
            f.write(" ".join(['{:^11s}'.format(h) for h in ['(m)', '(m)', '(kg/m^3)']]) + '\n')
            for i in range(hd['NMGDepths']):
                f.write(" ".join(['{:^11}'.format(hd['MGDpth'][i]),
                                  '{:^11}'.format(hd['MGThck'][i]),
                                  '{:^11}'.format(hd['MGDens'][i])]) + '\n')

            # MEMBER OUTPUT LIST
            f.write("---------------------- MEMBER OUTPUT LIST --------------------------------------\n")
            f.write('{:<11d} {:<11} {:}'.format(hd['NMOutputs'], 'NMOutputs', '- Number of member outputs (-)\n'))
            f.write(" ".join(['{:^11s}'.format(h) for h in ['MemberID_out', 'NOutLoc', 'NodeLocs']]) + '\n')
            f.write(" ".join(['{:^11s}'.format('(-)') for _ in range(3)]) + '\n')
            for i in range(hd['NMOutputs']):
                f.write(" ".join(['{:^11d}'.format(hd['MemberID_out'][i]),
                                  '{:^11d}'.format(hd['NOutLoc'][i]),
                                  '{:^11}'.format(hd['NodeLocs'][i])]) + '\n')

            # JOINT OUTPUT LIST
            f.write("---------------------- JOINT OUTPUT LIST ---------------------------------------\n")
            f.write('{:<22d} {:<11} {:}'.format(hd['NJOutputs'], 'NJOutputs', '- Number of joint outputs\n'))
            f.write('{:<22} {:<11} {:}'.format(" ".join(["%d" % i for i in hd['JOutLst']]), 'JOutLst', '- List of JointIDs\n'))

            # OUTPUT
            f.write("---------------------- OUTPUT --------------------------------------------------\n")
            f.write('{!s:<22} {:<11} {:}'.format(hd['HDSum'], 'HDSum', '- Output a summary file [flag]\n'))
            f.write('{!s:<22} {:<11} {:}'.format(hd['OutAll'], 'OutAll', '- Output all member and joint loads [flag]\n'))
            f.write('{:<22d} {:<11} {:}'.format(hd['OutSwtch'], 'OutSwtch', '- Output channels to file (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['OutFmt'], 'OutFmt', '- Output format\n'))
            f.write('{:<22} {:<11} {:}'.format(hd['OutSFmt'], 'OutSFmt', '- Output format for header strings\n'))

            f.write('---------------------- OUTPUT CHANNELS -----------------------------------------\n')
            if outlist is not None:
                ol = _get_outlist(outlist, ['HydroDyn'])
                for channel_list in ol:
                    for ch in channel_list:
                        f.write('"' + ch + '"\n')
            f.write('END of output channels and end of file.\n')
