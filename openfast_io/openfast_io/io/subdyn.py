"""
SubDynIO – read / write SubDyn input files.

Produces ``{'SubDyn': sd}``
"""
from __future__ import annotations

import os
from typing import Any, Dict, Optional, Callable

import numpy as np

from .base import ModuleIO
from ..outlist import emit_outlist, capture_outlist
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
    read_array,
    readline_filterComments,
)


class SubDynIO(ModuleIO):
    """Read / write SubDyn input files."""

    def read(
        self,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        read_outlist_fn: Optional[Callable] = None,
        **kwargs,
    ) -> dict:
        sd: Dict[str, Any] = {}
        sd_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        f = open(sd_file)
        f.readline()
        f.readline()
        f.readline()

        # SIMULATION CONTROL
        sd['Echo']      = bool_read(f.readline().split()[0])
        sd['SDdeltaT']  = float_read(f.readline().split()[0])
        sd['IntMethod'] = int_read(f.readline().split()[0])
        sd['SttcSolve'] = bool_read(f.readline().split()[0])
        f.readline()

        # FEA and CRAIG-BAMPTON PARAMETERS
        sd['FEMMod']       = int_read(f.readline().split()[0])
        sd['NDiv']         = int_read(f.readline().split()[0])
        sd['Nmodes']       = int_read(f.readline().split()[0])
        sd['JDampings']    = float_read(f.readline().split()[0])
        sd['GuyanDampMod'] = int_read(f.readline().split()[0])
        sd['RayleighDamp'] = read_array(f, 2, array_type=float)
        sd['GuyanDampSize'] = int_read(f.readline().split()[0])
        sd['GuyanDamp'] = np.array([
            [float_read(idx.strip(',')) for idx in f.readline().strip().split()[:sd['GuyanDampSize']]]
            for _ in range(sd['GuyanDampSize'])
        ])

        f.readline(); f.readline(); f.readline()

        # INITIAL RIGID-BODY POSITION
        ln = f.readline().split()
        sd['RBSurge'] = float(ln[0])
        sd['RBSway']  = float(ln[1])
        sd['RBHeave'] = float(ln[2])
        sd['RBRoll']  = float(ln[3])
        sd['RBPitch'] = float(ln[4])
        sd['RBYaw']   = float(ln[5])

        f.readline()

        # STRUCTURE JOINTS
        sd['NJoints'] = int_read(f.readline().split()[0])
        n = sd['NJoints']
        for k in ['JointID', 'JointXss', 'JointYss', 'JointZss', 'JointType',
                   'JointDirX', 'JointDirY', 'JointDirZ', 'JointStiff']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['JointID'][i]    = int(ln[0])
            sd['JointXss'][i]   = float(ln[1])
            sd['JointYss'][i]   = float(ln[2])
            sd['JointZss'][i]   = float(ln[3])
            sd['JointType'][i]  = int(ln[4])
            sd['JointDirX'][i]  = float(ln[5])
            sd['JointDirY'][i]  = float(ln[6])
            sd['JointDirZ'][i]  = float(ln[7])
            sd['JointStiff'][i] = float(ln[8])

        f.readline()

        # BASE REACTION JOINTS
        sd['NReact'] = int_read(f.readline().split()[0])
        n = sd['NReact']
        for k in ['RJointID', 'RctTDXss', 'RctTDYss', 'RctTDZss',
                   'RctRDXss', 'RctRDYss', 'RctRDZss', 'Rct_SoilFile']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['RJointID'][i] = int(ln[0])
            sd['RctTDXss'][i] = int(ln[1])
            sd['RctTDYss'][i] = int(ln[2])
            sd['RctTDZss'][i] = int(ln[3])
            sd['RctRDXss'][i] = int(ln[4])
            sd['RctRDYss'][i] = int(ln[5])
            sd['RctRDZss'][i] = int(ln[6])
            sd['Rct_SoilFile'][i] = ln[7] if len(ln) == 8 else 'None'

        f.readline()

        # INTERFACE JOINTS
        sd['NInterf'] = int_read(f.readline().split()[0])
        n = sd['NInterf']
        for k in ['IJointID', 'TPID', 'ItfTDXss', 'ItfTDYss', 'ItfTDZss',
                   'ItfRDXss', 'ItfRDYss', 'ItfRDZss']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['IJointID'][i] = int(ln[0])
            sd['TPID'][i]     = int(ln[1])
            sd['ItfTDXss'][i] = int(ln[2])
            sd['ItfTDYss'][i] = int(ln[3])
            sd['ItfTDZss'][i] = int(ln[4])
            sd['ItfRDXss'][i] = int(ln[5])
            sd['ItfRDYss'][i] = int(ln[6])
            sd['ItfRDZss'][i] = int(ln[7])

        f.readline()

        # MEMBERS
        sd['NMembers'] = int_read(f.readline().split()[0])
        n = sd['NMembers']
        for k in ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2',
                   'MType', 'M_Spin', 'M_COSMID']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['MemberID'][i]    = int(ln[0])
            sd['MJointID1'][i]   = int(ln[1])
            sd['MJointID2'][i]   = int(ln[2])
            sd['MPropSetID1'][i] = int(ln[3])
            sd['MPropSetID2'][i] = int(ln[4])
            if ln[5].lower() == '1c':
                sd['MType'][i] = 1
            elif ln[5].lower() == '1r':
                sd['MType'][i] = -1
            else:
                sd['MType'][i] = int(ln[5])
            if sd['MType'][i] == 5:
                sd['M_Spin'][i]   = 0.
                sd['M_COSMID'][i] = int(ln[6])
            else:
                sd['M_Spin'][i]   = float(ln[6])
                sd['M_COSMID'][i] = -1

        f.readline()

        # MEMBER X-SECTION PROPERTY DATA 1/3 — Circular
        sd['NPropSetsCyl'] = int_read(f.readline().split()[0])
        n = sd['NPropSetsCyl']
        for k in ['PropSetID1', 'YoungE1', 'ShearG1', 'MatDens1', 'XsecD', 'XsecT']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['PropSetID1'][i] = int(ln[0])
            sd['YoungE1'][i]    = float(ln[1])
            sd['ShearG1'][i]    = float(ln[2])
            sd['MatDens1'][i]   = float(ln[3])
            sd['XsecD'][i]      = float(ln[4])
            sd['XsecT'][i]      = float(ln[5])

        f.readline()

        # MEMBER X-SECTION PROPERTY DATA 2/3 — Rectangular
        sd['NPropSetsRec'] = int_read(f.readline().split()[0])
        n = sd['NPropSetsRec']
        for k in ['PropSetID2', 'YoungE2', 'ShearG2', 'MatDens2', 'XsecSa', 'XsecSb', 'XsecT2']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['PropSetID2'][i] = int(ln[0])
            sd['YoungE2'][i]    = float(ln[1])
            sd['ShearG2'][i]    = float(ln[2])
            sd['MatDens2'][i]   = float(ln[3])
            sd['XsecSa'][i]     = float(ln[4])
            sd['XsecSb'][i]     = float(ln[5])
            sd['XsecT2'][i]     = float(ln[6])

        f.readline()

        # MEMBER X-SECTION PROPERTY DATA 3/3 — Arbitrary
        sd['NXPropSets'] = int_read(f.readline().split()[0])
        n = sd['NXPropSets']
        arb_keys = ['PropSetID3', 'YoungE3', 'ShearG3', 'MatDens3', 'XsecA',
                    'XsecAsx', 'XsecAsy', 'XsecJxx', 'XsecJyy', 'XsecJ0', 'XsecJt']
        for k in arb_keys:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['PropSetID3'][i] = int(ln[0])
            for j, k in enumerate(arb_keys[1:], 1):
                sd[k][i] = float(ln[j])

        # CABLE PROPERTIES
        f.readline()
        sd['NCablePropSets'] = int_read(f.readline().split()[0])
        n = sd['NCablePropSets']
        for k in ['CablePropSetID', 'CableEA', 'CableMatDens', 'CableT0']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['CablePropSetID'][i] = int(ln[0])
            sd['CableEA'][i]        = float(ln[1])
            sd['CableMatDens'][i]   = float(ln[2])
            sd['CableT0'][i]        = float(ln[3])

        # RIGID LINK PROPERTIES
        f.readline()
        sd['NRigidPropSets'] = int_read(f.readline().split()[0])
        n = sd['NRigidPropSets']
        for k in ['RigidPropSetID', 'RigidMatDens']:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['RigidPropSetID'][i] = int(ln[0])
            sd['RigidMatDens'][i]   = float(ln[1])

        # SPRING ELEMENT PROPERTIES
        f.readline()
        sd['NSpringPropSets'] = int_read(f.readline().split()[0])
        n = sd['NSpringPropSets']
        spring_list = ['k11','k12','k13','k14','k15','k16',
                       'k22','k23','k24','k25','k26',
                       'k33','k34','k35','k36',
                       'k44','k45','k46',
                       'k55','k56',
                       'k66']
        sd['SpringPropSetID'] = [None] * n
        for sl in spring_list:
            sd[sl] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['SpringPropSetID'][i] = int(ln[0])
            for j, sl in enumerate(spring_list):
                sd[sl][i] = ln[j + 1]

        # MEMBER COSINE MATRICES
        f.readline()
        sd['NCOSMs'] = int_read(f.readline().split()[0])
        n = sd['NCOSMs']
        cosm_keys = ['COSMID', 'COSM11', 'COSM12', 'COSM13',
                     'COSM21', 'COSM22', 'COSM23',
                     'COSM31', 'COSM32', 'COSM33']
        for k in cosm_keys:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['COSMID'][i] = int(ln[0])
            for j in range(1, 10):
                sd[cosm_keys[j]][i] = float(ln[j])

        f.readline()

        # JOINT ADDITIONAL CONCENTRATED MASSES
        sd['NCmass'] = int_read(f.readline().split()[0])
        n = sd['NCmass']
        mass_keys = ['CMJointID', 'JMass', 'JMXX', 'JMYY', 'JMZZ',
                     'JMXY', 'JMXZ', 'JMYZ', 'MCGX', 'MCGY', 'MCGZ']
        for k in mass_keys:
            sd[k] = [None] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split()
            sd['CMJointID'][i] = int(ln[0])
            for j, k in enumerate(mass_keys[1:], 1):
                sd[k][i] = float(ln[j])

        f.readline()

        # OUTPUT
        sd['SumPrint'] = bool_read(f.readline().split()[0])

        # Optional OutCBModes / OutFEMModes
        file_pos = f.tell()
        line = f.readline()
        if 'OutCBModes' in line:
            sd['OutCBModes'] = int_read(line.split()[0])
        else:
            f.seek(file_pos)
        file_pos = f.tell()
        line = f.readline()
        if 'OutFEMModes' in line:
            sd['OutFEMModes'] = int_read(line.split()[0])
        else:
            f.seek(file_pos)

        sd['OutCOSM']  = bool_read(f.readline().split()[0])
        sd['OutAll']   = bool_read(f.readline().split()[0])
        sd['OutSwtch'] = int_read(f.readline().split()[0])
        sd['TabDelim'] = bool_read(f.readline().split()[0])
        sd['OutDec']   = int_read(f.readline().split()[0])
        sd['OutFmt']   = quoted_read(f.readline().split()[0])
        sd['OutSFmt']  = quoted_read(f.readline().split()[0])

        f.readline()

        # MEMBER OUTPUT LIST
        sd['NMOutputs'] = int_read(f.readline().split()[0])
        n = sd['NMOutputs']
        sd['MemberID_out'] = [None] * n
        sd['NOutCnt']      = [None] * n
        sd['NodeCnt']      = [[None]] * n
        f.readline(); f.readline()
        for i in range(n):
            ln = f.readline().split('!')[0].split()
            sd['MemberID_out'][i] = int(ln[0])
            sd['NOutCnt'][i]      = int(ln[1])
            sd['NodeCnt'][i]      = [int(node) for node in ln[2:]]

        f.readline()

        # SSOutList
        if read_outlist_fn is not None and outlist is not None:
            read_outlist_fn(f, 'SubDyn')
        else:
            # Standalone use with no shared registry — consume the section
            # safely and stash the found channels so write() can still emit
            # them (mirrors the aerodisk.py fallback pattern).
            channels = capture_outlist(f, None, 'SubDyn', freeform=True)
            sd['_outlist'] = {ch: True for ch in channels}

        f.close()
        return {'SubDyn': sd}

    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        **kwargs,
    ) -> None:
        sd = data['SubDyn']
        with open(file_path, 'w') as f:
            f.write('----------- SubDyn MultiMember Support Structure Input File ------------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('-------------------------- SIMULATION CONTROL ---------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sd['Echo'], 'Echo', '- Echo input data to "<rootname>.SD.ech" (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['SDdeltaT'], 'SDdeltaT', '- Local Integration Step. Use "default" to base it on the glue-code time step.\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['IntMethod'], 'IntMethod', '- Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2].\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['SttcSolve'], 'SttcSolve', '- Solve dynamics about static equilibrium point\n'))
            f.write('--- FEA and CRAIG-BAMPTON PARAMETERS ---\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['FEMMod'], 'FEMMod', '- FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B); 2=Timoshenko; 3= E-B w/ shear; 4=MITC3+]\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['NDiv'], 'NDiv', '- Number of sub-elements per member\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['Nmodes'], 'Nmodes', '- Number of internal modes to retain.\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['JDampings'], 'JDampings', '- Damping Ratios for each retained mode.\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['GuyanDampMod'], 'GuyanDampMod', '- Guyan damping (switch)\n'))
            f.write('{:<22} {:<11} {:}'.format(', '.join([str(v) for v in sd['RayleighDamp']]), 'RayleighDamp', '- Mass and stiffness proportional damping coefficients\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['GuyanDampSize'], 'GuyanDampSize', '- Guyan damping matrix size (square)\n'))
            for i in range(sd['GuyanDampSize']):
                try:
                    f.write(' '.join(['{:14}'.format(v) for v in sd['GuyanDamp'][i, :]]) + '\n')
                except Exception:
                    f.write(' '.join(['{:14}'.format(v) for v in sd['GuyanDamp'][i]]) + '\n')
            f.write('------- INITIAL RIGID-BODY POSITION [used only for floating structure with more than one transition pieces] -------\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['RBSurge', 'RBSway', 'RBHeave', 'RBRoll', 'RBPitch', 'RBYaw']]) + '\n')
            f.write(" ".join(['{:^11s}'.format(i) for i in ['(m)', '(m)', '(m)', '(deg)', '(deg)', '(deg)']]) + '\n')
            f.write(" ".join(['{:^11}'.format(sd[k]) for k in ['RBSurge', 'RBSway', 'RBHeave', 'RBRoll', 'RBPitch', 'RBYaw']]) + '\n')
            f.write('---- STRUCTURE JOINTS ----\n')

            # Helper for writing a table section
            def _write_table(keys, header_keys, units, count_key, count_label, desc):
                f.write('{:<22d} {:<11} {:}'.format(sd[count_key], count_label, desc + '\n'))
                f.write(' '.join(['{:^11s}'.format(h) for h in header_keys]) + '\n')
                f.write(' '.join(['{:^11s}'.format(u) for u in units]) + '\n')

            # Joints
            _write_table([], ['JointID', 'JointXss', 'JointYss', 'JointZss', 'JointType', 'JointDirX', 'JointDirY', 'JointDirZ', 'JointStiff'],
                         ['(-)', '(m)', '(m)', '(m)', '(-)', '(-)', '(-)', '(-)', '(Nm/rad)'],
                         'NJoints', 'NJoints', '- Number of joints (-)')
            for i in range(sd['NJoints']):
                f.write(' '.join(['{:^11d}'.format(sd['JointID'][i]),
                                  '{:^11}'.format(sd['JointXss'][i]), '{:^11}'.format(sd['JointYss'][i]), '{:^11}'.format(sd['JointZss'][i]),
                                  '{:^11d}'.format(sd['JointType'][i]),
                                  '{:^11}'.format(sd['JointDirX'][i]), '{:^11}'.format(sd['JointDirY'][i]), '{:^11}'.format(sd['JointDirZ'][i]),
                                  '{:^11}'.format(sd['JointStiff'][i])]) + '\n')

            # Base Reactions
            f.write('---- BASE REACTION JOINTS ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NReact'], 'NReact', '- Number of joints with reaction forces\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['RJointID', 'RctTDXss', 'RctTDYss', 'RctTDZss', 'RctRDXss', 'RctRDYss', 'RctRDZss', 'Rct_SoilFile']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(8)]) + '\n')
            for i in range(sd['NReact']):
                f.write(' '.join(['{:^11d}'.format(sd['RJointID'][i]),
                                  '{:^11d}'.format(sd['RctTDXss'][i]), '{:^11d}'.format(sd['RctTDYss'][i]), '{:^11d}'.format(sd['RctTDZss'][i]),
                                  '{:^11d}'.format(sd['RctRDXss'][i]), '{:^11d}'.format(sd['RctRDYss'][i]), '{:^11d}'.format(sd['RctRDZss'][i]),
                                  '{:^11}'.format(sd['Rct_SoilFile'][i])]) + '\n')

            # Interface Joints
            f.write('---- INTERFACE JOINTS ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NInterf'], 'NInterf', '- Number of interface joints\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['IJointID', 'TPID', 'ItfTDXss', 'ItfTDYss', 'ItfTDZss', 'ItfRDXss', 'ItfRDYss', 'ItfRDZss']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(8)]) + '\n')
            for i in range(sd['NInterf']):
                f.write(' '.join(['{:^11d}'.format(sd['IJointID'][i]), '{:^11d}'.format(sd['TPID'][i]),
                                  '{:^11d}'.format(sd['ItfTDXss'][i]), '{:^11d}'.format(sd['ItfTDYss'][i]), '{:^11d}'.format(sd['ItfTDZss'][i]),
                                  '{:^11d}'.format(sd['ItfRDXss'][i]), '{:^11d}'.format(sd['ItfRDYss'][i]), '{:^11d}'.format(sd['ItfRDZss'][i])]) + '\n')

            # Members
            f.write('---- MEMBERS ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NMembers'], 'NMembers', '- Number of frame members\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['MemberID', 'MJointID1', 'MJointID2', 'MPropSetID1', 'MPropSetID2', 'MType', 'COSMID/Spin']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(7)]) + '\n')
            for i in range(sd['NMembers']):
                if sd['MType'][i] == 5:
                    spin_val = '{:^11d}'.format(sd['M_COSMID'][i])
                else:
                    spin_val = '{:^11}'.format(sd['M_Spin'][i])
                f.write(' '.join(['{:^11d}'.format(sd['MemberID'][i]),
                                  '{:^11d}'.format(sd['MJointID1'][i]), '{:^11d}'.format(sd['MJointID2'][i]),
                                  '{:^11d}'.format(sd['MPropSetID1'][i]), '{:^11d}'.format(sd['MPropSetID2'][i]),
                                  '{:^11d}'.format(sd['MType'][i]), spin_val]) + '\n')

            # Cylindrical Section Properties
            f.write('---- MEMBER X-SECTION PROPERTY data 1/3 [isotropic material for circular cross-sections] ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NPropSetsCyl'], 'NPropSetsCyl', '- Number of circular member property sets\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecD', 'XsecT']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(6)]) + '\n')
            for i in range(sd['NPropSetsCyl']):
                f.write(' '.join(['{:^11d}'.format(sd['PropSetID1'][i]),
                                  '{:^11}'.format(sd['YoungE1'][i]), '{:^11}'.format(sd['ShearG1'][i]),
                                  '{:^11}'.format(sd['MatDens1'][i]),
                                  '{:^11}'.format(sd['XsecD'][i]), '{:^11}'.format(sd['XsecT'][i])]) + '\n')

            # Rectangular 
            f.write('---- MEMBER X-SECTION PROPERTY data 2/3 [isotropic material for rectangular cross-sections] ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NPropSetsRec'], 'NPropSetsRec', '- Number of rectangular member property sets\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecSa', 'XsecSb', 'XsecT']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(7)]) + '\n')
            for i in range(sd['NPropSetsRec']):
                f.write(' '.join(['{:^11d}'.format(sd['PropSetID2'][i]),
                                  '{:^11}'.format(sd['YoungE2'][i]), '{:^11}'.format(sd['ShearG2'][i]),
                                  '{:^11}'.format(sd['MatDens2'][i]),
                                  '{:^11}'.format(sd['XsecSa'][i]), '{:^11}'.format(sd['XsecSb'][i]),
                                  '{:^11}'.format(sd['XsecT2'][i])]) + '\n')

            # Arbitrary
            f.write('---- MEMBER X-SECTION PROPERTY data 3/3 [arbitrary cross-sections] ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NXPropSets'], 'NXPropSets', '- Number of arbitrary member property sets\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecA', 'XsecAsx', 'XsecAsy', 'XsecJxx', 'XsecJyy', 'XsecJ0', 'XsecJt']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(11)]) + '\n')
            arb_keys = ['PropSetID3', 'YoungE3', 'ShearG3', 'MatDens3', 'XsecA',
                        'XsecAsx', 'XsecAsy', 'XsecJxx', 'XsecJyy', 'XsecJ0', 'XsecJt']
            for i in range(sd['NXPropSets']):
                ln = ['{:^11d}'.format(int(sd['PropSetID3'][i]))]
                for k in arb_keys[1:]:
                    ln.append('{:^11}'.format(sd[k][i]))
                f.write(' '.join(ln) + '\n')

            # Cable
            f.write('---- CABLE PROPERTIES ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NCablePropSets'], 'NCablePropSets', '- Number of cable property sets\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID', 'EA', 'MatDens', 'T0']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(4)]) + '\n')
            for i in range(sd['NCablePropSets']):
                f.write(' '.join(['{:^11d}'.format(sd['CablePropSetID'][i]),
                                  '{:^11}'.format(sd['CableEA'][i]),
                                  '{:^11}'.format(sd['CableMatDens'][i]),
                                  '{:^11}'.format(sd['CableT0'][i])]) + '\n')

            # Rigid Link
            f.write('---- RIGID LINK PROPERTIES ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NRigidPropSets'], 'NRigidPropSets', '- Number of rigid link property sets\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID', 'MatDens']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(2)]) + '\n')
            for i in range(sd['NRigidPropSets']):
                f.write(' '.join(['{:^11d}'.format(sd['RigidPropSetID'][i]),
                                  '{:^11}'.format(sd['RigidMatDens'][i])]) + '\n')

            # Spring
            f.write('---- SPRING ELEMENT PROPERTIES ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NSpringPropSets'], 'NSpringPropSets', '- Number of spring element property sets\n'))
            spring_list = ['k11','k12','k13','k14','k15','k16',
                           'k22','k23','k24','k25','k26',
                           'k33','k34','k35','k36',
                           'k44','k45','k46',
                           'k55','k56','k66']
            f.write(' '.join(['{:^11s}'.format(h) for h in ['PropSetID'] + spring_list]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(len(spring_list) + 1)]) + '\n')
            for i in range(sd['NSpringPropSets']):
                ln = ['{:^11d}'.format(sd['SpringPropSetID'][i])]
                for sl in spring_list:
                    ln.append('{:^11}'.format(sd[sl][i]))
                f.write(' '.join(ln) + '\n')

            # COSM
            f.write('---- MEMBER COSINE MATRICES ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NCOSMs'], 'NCOSMs', '- Number of unique cosine matrices\n'))
            cosm_keys = ['COSMID', 'COSM11', 'COSM12', 'COSM13', 'COSM21', 'COSM22', 'COSM23', 'COSM31', 'COSM32', 'COSM33']
            f.write(' '.join(['{:^11s}'.format(h) for h in cosm_keys]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in cosm_keys]) + '\n')
            for i in range(sd['NCOSMs']):
                ln = ['{:^11d}'.format(sd['COSMID'][i])]
                for k in cosm_keys[1:]:
                    ln.append('{:^11}'.format(sd[k][i]))
                f.write(' '.join(ln) + '\n')

            # Concentrated Masses
            f.write('---- JOINT ADDITIONAL CONCENTRATED MASSES ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NCmass'], 'NCmass', '- Number of joints with concentrated masses\n'))
            mass_keys = ['CMJointID', 'JMass', 'JMXX', 'JMYY', 'JMZZ', 'JMXY', 'JMXZ', 'JMYZ', 'MCGX', 'MCGY', 'MCGZ']
            f.write(' '.join(['{:^11s}'.format(h) for h in mass_keys]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in mass_keys]) + '\n')
            for i in range(sd['NCmass']):
                ln = ['{:^11d}'.format(sd['CMJointID'][i])]
                for k in mass_keys[1:]:
                    ln.append('{:^11}'.format(sd[k][i]))
                f.write(' '.join(ln) + '\n')

            # Output
            f.write('---- OUTPUT ------------------------------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(sd['SumPrint'], 'SumPrint', '- Output summary file\n'))
            if 'OutCBModes' in sd:
                f.write('{:<22d} {:<11} {:}'.format(sd['OutCBModes'], 'OutCBModes', '- Output CB modes\n'))
            if 'OutFEMModes' in sd:
                f.write('{:<22d} {:<11} {:}'.format(sd['OutFEMModes'], 'OutFEMModes', '- Output FEM modes\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['OutCOSM'], 'OutCOSM', '- Output cosine matrices (flag)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['OutAll'], 'OutAll', '- Output all member/joint loads (flag)\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['OutSwtch'], 'OutSwtch', '- Output channels (switch)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(sd['TabDelim'], 'TabDelim', '- Tab-delimited output (flag)\n'))
            f.write('{:<22d} {:<11} {:}'.format(sd['OutDec'], 'OutDec', '- Decimation of output\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['OutFmt'], 'OutFmt', '- Output format\n'))
            f.write('{:<22} {:<11} {:}'.format(sd['OutSFmt'], 'OutSFmt', '- Output format for header strings\n'))

            # Member Output List
            f.write('---- MEMBER OUTPUT LIST ----\n')
            f.write('{:<22d} {:<11} {:}'.format(sd['NMOutputs'], 'NMOutputs', '- Number of member outputs\n'))
            f.write(' '.join(['{:^11s}'.format(h) for h in ['MemberID', 'NOutCnt', 'NodeCnt']]) + '\n')
            f.write(' '.join(['{:^11s}'.format('(-)') for _ in range(3)]) + '\n')
            for i in range(sd['NMOutputs']):
                f.write(' '.join(['{:^11d}'.format(sd['MemberID_out'][i]),
                                  '{:^11d}'.format(sd['NOutCnt'][i]),
                                  ' '.join([str(nc) for nc in sd['NodeCnt'][i]])]) + '\n')

            f.write('---- SSOutList ----\n')
            if outlist is not None:
                emit_outlist(f, outlist, 'SubDyn')
            else:
                for ch in sd.get('_outlist', {}):  # standalone fallback (no shared registry)
                    f.write('"' + ch + '"\n')
            f.write('END of output channels and end of file.\n')
