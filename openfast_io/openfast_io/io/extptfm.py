"""
ExtPtfmIO – read / write ExtPtfm (External Platform) input files.

Produces ``{'ExtPtfm': ep}``

ExtPtfm reads the main input file plus up to four sub-files:
  - Superelement (Guyan/Craig-Bampton reduction matrices)
  - Connections (optional)
  - UserForcing (optional)
  - ConnForcing (optional)
"""
from __future__ import annotations

import os
import re
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


class ExtPtfmIO(ModuleIO):
    """Read / write ExtPtfm input files."""

    # ------------------------------------------------------------------
    # READ
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
        ep: Dict[str, Any] = {}
        ep_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path
        ep_dir = os.path.dirname(ep_file)

        f = open(ep_file)
        f.readline()
        f.readline()
        f.readline()

        # Simulation Control
        ep['Echo'] = bool_read(f.readline().split()[0])
        ep['DT'] = float_read(f.readline().split()[0])
        ep['IntMethod'] = int_read(f.readline().split()[0])
        f.readline()

        # Reduction inputs
        ep['RBMod'] = int_read(f.readline().split()[0])
        ep['Red_FileName'] = os.path.join(ep_dir, quoted_read(f.readline().split()[0]))
        ep['NActiveDOFList'] = int_read(f.readline().split()[0])
        ep['ActiveDOFList'] = read_array(f, None, split_val='ActiveDOFList', array_type=int)
        ep['NInitPosList'] = int_read(f.readline().split()[0])
        ep['InitPosList'] = read_array(f, None, split_val='InitPosList', array_type=float)
        ep['NInitVelList'] = int_read(f.readline().split()[0])
        ep['InitVelList'] = read_array(f, None, split_val='InitVelList', array_type=float)
        f.readline()

        # Connection inputs
        ep['HasConnections'] = bool_read(f.readline().split()[0])
        ep['Conn_FileName'] = os.path.join(ep_dir, quoted_read(f.readline().split()[0]))
        f.readline()

        # User forcing inputs
        ep['HasUserForcing'] = bool_read(f.readline().split()[0])
        ep['Force_FileName'] = os.path.join(ep_dir, quoted_read(f.readline().split()[0]))
        ep['HasConnForcing'] = bool_read(f.readline().split()[0])
        ep['FConn_FileName'] = os.path.join(ep_dir, quoted_read(f.readline().split()[0]))
        f.readline()

        # Output
        ep['SumPrint'] = bool_read(f.readline().split()[0])
        ep['OutFile'] = int_read(f.readline().split()[0])
        ep['TabDelim'] = bool_read(f.readline().split()[0])
        ep['OutFmt'] = quoted_read(f.readline().split()[0])
        ep['TStart'] = float_read(f.readline().split()[0])

        # Output channels
        f.readline()
        data = f.readline()
        while data.split().__len__() == 0:
            data = f.readline()

        ep['_outlist'] = {}
        while data.split()[0] != 'END':
            if data.find('"') >= 0:
                channels = data.split('"')
                channel_list = channels[1].split(',')
            else:
                row_string = data.split(',')
                if len(row_string) == 1:
                    channel_list = row_string[0].split('\n')[0]
                else:
                    channel_list = row_string
            if isinstance(channel_list, list):
                for ch in channel_list:
                    ch = ch.strip()
                    if ch:
                        ep['_outlist'][ch] = True
            else:
                ch = channel_list.strip()
                if ch:
                    ep['_outlist'][ch] = True
            data = f.readline()

        # Populate the shared registry (freeform — ExtPtfm CBD*/CBF* channels are
        # not in the standard channel registry). data was pre-read, so we copy the
        # parsed set rather than re-reading the section.
        if outlist is not None:
            outlist.setdefault('ExtPtfm', {})
            for ch in ep['_outlist']:
                outlist['ExtPtfm'][ch] = True
            # Channels now live in the shared registry; don't pollute fst_vt['ExtPtfm']
            # with a private '_outlist' key baseline never creates.
            ep.pop('_outlist', None)

        f.close()

        # Read sub-files
        ep['FlexASCII'] = {}
        self._read_superelement(ep['Red_FileName'], ep['FlexASCII'])

        if ep['HasConnections']:
            ep['Connections'] = {}
            self._read_connections(ep['Conn_FileName'], ep['FlexASCII']['nDOF'], ep['Connections'])

        if ep['HasUserForcing']:
            ep['UserForcing'] = {}
            self._read_user_forcing(ep['Force_FileName'], ep['FlexASCII']['nDOF'], ep['UserForcing'])

        if ep['HasConnForcing']:
            ep['ConnForcing'] = {}
            nConn = ep['Connections']['nConn']
            self._read_conn_forcing(ep['FConn_FileName'], nConn, ep['ConnForcing'])

        return {'ExtPtfm': ep}

    # ------------------------------------------------------------------
    # Sub-file readers
    # ------------------------------------------------------------------

    @staticmethod
    def _readmat(n, m, lines, iStart):
        M = np.zeros((n, m))
        for j in range(n):
            M[j, :] = np.array(lines[iStart + j].split()).astype(float)
        return M

    def _read_superelement(self, se_file, flex):
        """Read superelement file into *flex* dict."""
        with open(se_file) as fh:
            lines = fh.read().splitlines()

        nDOF = -1
        i = 0
        while i < len(lines):
            lo = lines[i].lower()
            if lo.find('!mass') == 0:
                flex['MassMatrix'] = self._readmat(nDOF, nDOF, lines, i + 1)
                i += nDOF
            elif lo.find('!stiffness') == 0:
                flex['StiffnessMatrix'] = self._readmat(nDOF, nDOF, lines, i + 1)
                i += nDOF
            elif lo.find('!damping') == 0:
                flex['DampingMatrix'] = self._readmat(nDOF, nDOF, lines, i + 1)
                i += nDOF
            elif lo.find('!weight constant') == 0:
                flex['WeightConstant'] = self._readmat(1, nDOF, lines, i + 1)
                i += 1
            elif lo.find('!weight stiffness') == 0:
                flex['WeightStiffness'] = self._readmat(nDOF, nDOF, lines, i + 1)
                i += nDOF
            elif len(lo) > 0 and lo[0] == '!' and lo.find('!dimension') == 0:
                flex['nDOF'] = int(lo.split(':')[1])
                nDOF = flex['nDOF']
            i += 1

    def _read_connections(self, conn_file, nDOF, conn):
        """Read connections file into *conn* dict."""
        with open(conn_file) as fh:
            lines = fh.read().splitlines()

        nConn = -1
        i = 0
        while i < len(lines):
            lo = lines[i].lower()
            if lo.find('!connection') == 0:
                conn['Position'] = self._readmat(nConn, 3, lines, i + 1)
                i += nConn
            elif lo.find('!displacement') == 0:
                conn['Displacement'] = self._readmat(3 * nConn, nDOF, lines, i + 1)
                i += 3 * nConn
            elif len(lo) > 0 and lo[0] == '!' and lo.find('!nconn') == 0:
                conn['nConn'] = int(lo.split(':')[1])
                nConn = conn['nConn']
            i += 1

    def _read_user_forcing(self, force_file, nDOF, uf):
        """Read user forcing file into *uf* dict."""
        with open(force_file) as fh:
            lines = fh.read().splitlines()

        nSteps = -1
        i = 0
        while i < len(lines):
            lo = lines[i].lower()
            if lo.find('!forcing') == 0:
                uf['ForceTimeSeries'] = self._readmat(nSteps, 1 + nDOF, lines, i + 1)
                i += nSteps
            elif len(lo) > 0 and lo[0] == '!' and lo.find('!nsteps') == 0:
                uf['nSteps'] = int(lo.split(':')[1])
                nSteps = uf['nSteps']
            i += 1

    def _read_conn_forcing(self, fconn_file, nConn, cf):
        """Read connection forcing file into *cf* dict."""
        with open(fconn_file) as fh:
            lines = fh.read().splitlines()

        nSteps = -1
        i = 0
        while i < len(lines):
            lo = lines[i].lower()
            if lo.find('!forcing') == 0:
                cf['ForceTimeSeries'] = self._readmat(nSteps, 1 + 3 * nConn, lines, i + 1)
                i += nSteps
            elif len(lo) > 0 and lo[0] == '!' and lo.find('!nsteps') == 0:
                cf['nSteps'] = int(lo.split(':')[1])
                nSteps = cf['nSteps']
            i += 1

    # ------------------------------------------------------------------
    # WRITE
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
        ep = data['ExtPtfm']
        out_dir = os.path.dirname(file_path) or '.'

        # Write sub-files first
        se_name = os.path.basename(ep.get('Red_FileName', 'ExtPtfm_SE.dat'))
        se_path = os.path.join(out_dir, se_name)
        self._write_superelement(ep['FlexASCII'], se_path)

        if ep.get('HasConnections') and 'Connections' in ep:
            conn_name = os.path.basename(ep.get('Conn_FileName', 'ExtPtfm_Conn.dat'))
            conn_path = os.path.join(out_dir, conn_name)
            self._write_connections(ep['Connections'], conn_path)

        if ep.get('HasUserForcing') and 'UserForcing' in ep:
            frc_name = os.path.basename(ep.get('Force_FileName', 'ExtPtfm_UserFrc.dat'))
            frc_path = os.path.join(out_dir, frc_name)
            self._write_user_forcing(ep['UserForcing'], frc_path)

        if ep.get('HasConnForcing') and 'ConnForcing' in ep:
            fcn_name = os.path.basename(ep.get('FConn_FileName', 'ExtPtfm_ConnFrc.dat'))
            fcn_path = os.path.join(out_dir, fcn_name)
            self._write_conn_forcing(ep['ConnForcing'], fcn_path)

        with open(file_path, 'w') as f:
            f.write('---------------------- EXTPTFM INPUT FILE --------------------------------------\n')
            f.write('Comment describing the model\n')
            f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ep['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ep['DT'], 'DT', '- Communication interval for controllers (s) (or "default")\n'))
            f.write('{:<22d} {:<11} {:}'.format(ep['IntMethod'], 'IntMethod', '- Integration Method {1:RK4; 2:AB4, 3:ABM4} (switch)\n'))
            f.write('---------------------- REDUCTION INPUTS ----------------------------------------\n')
            f.write('{:<22d} {:<11} {:}'.format(ep['RBMod'], 'RBMod', '- Method for handling rigid-body motion\n'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + se_name + '"', 'Red_FileName', '- Path of superelement file\n'))
            f.write('{:<22d} {:<11} {:}'.format(ep['NActiveDOFList'], 'NActiveDOFList', '- Number of active CB modes\n'))
            f.write('{:<22} {:<11} {:}'.format(', '.join([str(v) for v in ep['ActiveDOFList']]), 'ActiveDOFList', '- List of active CB mode indices\n'))
            f.write('{:<22d} {:<11} {:}'.format(ep['NInitPosList'], 'NInitPosList', '- Number of initial positions\n'))
            f.write('{:<22} {:<11} {:}'.format(', '.join([str(v) for v in ep['InitPosList']]), 'InitPosList', '- List of initial positions\n'))
            f.write('{:<22d} {:<11} {:}'.format(ep['NInitVelList'], 'NInitVelList', '- Number of initial velocities\n'))
            f.write('{:<22} {:<11} {:}'.format(', '.join([str(v) for v in ep['InitVelList']]), 'InitVelList', '- List of initial velocities\n'))
            f.write('---------------------- CONNECTION INPUTS ---------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ep['HasConnections'], 'Connections', '- Flag for connection points\n'))
            conn_name_w = os.path.basename(ep.get('Conn_FileName', 'none'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + conn_name_w + '"', 'Conn_FileName', '- Path of connection file\n'))
            f.write('---------------------- USER FORCING INPUTS -------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ep['HasUserForcing'], 'UserForcing', '- Flag for user forcing\n'))
            frc_name_w = os.path.basename(ep.get('Force_FileName', 'none'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + frc_name_w + '"', 'Force_FileName', '- Path of user forcing file\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ep['HasConnForcing'], 'ConnForcing', '- Flag for connection forcing\n'))
            fcn_name_w = os.path.basename(ep.get('FConn_FileName', 'none'))
            f.write('{!s:<22} {:<11} {:}'.format('"' + fcn_name_w + '"', 'FConn_FileName', '- Path of connection forcing file\n'))
            f.write('---------------------- OUTPUT --------------------------------------------------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ep['SumPrint'], 'SumPrint', '- Print summary data\n'))
            f.write('{:<22d} {:<11} {:}'.format(ep['OutFile'], 'OutFile', '- Output switch\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ep['TabDelim'], 'TabDelim', '- Tab delimited output\n'))
            f.write('{!s:<22} {:<11} {:}'.format(ep['OutFmt'], 'OutFmt', '- Output format\n'))
            f.write('{:<22f} {:<11} {:}'.format(ep['TStart'], 'TStart', '- Time to begin output\n'))
            f.write('                    OutList\n')
            out_channels = ep.get('_outlist', {})
            if outlist and 'ExtPtfm' in outlist:
                out_channels = outlist['ExtPtfm']
            for ch in out_channels:
                f.write('"' + ch + '"\n')
            f.write('END of input file\n')

    # ------------------------------------------------------------------
    # Sub-file writers
    # ------------------------------------------------------------------

    @staticmethod
    def _mat_to_string(M):
        return '\n'.join(''.join('{:16.8e}'.format(x) for x in row) for row in M)

    def _write_superelement(self, flex, path):
        with open(path, 'w') as f:
            f.write('!Dimension: {}\n'.format(flex['nDOF']))
            f.write('\n!Mass Matrix\n')
            f.write(self._mat_to_string(flex['MassMatrix']))
            f.write('\n\n!Stiffness Matrix\n')
            f.write(self._mat_to_string(flex['StiffnessMatrix']))
            f.write('\n\n!Damping Matrix\n')
            f.write(self._mat_to_string(flex['DampingMatrix']))
            f.write('\n\n!Weight constant\n')
            f.write(self._mat_to_string(flex['WeightConstant']))
            f.write('\n\n!Weight stiffness matrix\n')
            f.write(self._mat_to_string(flex['WeightStiffness']))
            f.write('\n')

    def _write_connections(self, conn, path):
        with open(path, 'w') as f:
            f.write('!nConn: {}\n'.format(conn['nConn']))
            f.write('\n!Connections\n')
            f.write(self._mat_to_string(conn['Position']))
            f.write('\n\n!Displacement\n')
            f.write(self._mat_to_string(conn['Displacement']))
            f.write('\n')

    def _write_user_forcing(self, uf, path):
        with open(path, 'w') as f:
            f.write('!nSteps: {}\n'.format(uf['nSteps']))
            f.write('\n!Forcing:\n')
            f.write(self._mat_to_string(uf['ForceTimeSeries']))
            f.write('\n')

    def _write_conn_forcing(self, cf, path):
        with open(path, 'w') as f:
            f.write('!nSteps: {}\n'.format(cf['nSteps']))
            f.write('\n!Forcing:\n')
            f.write(self._mat_to_string(cf['ForceTimeSeries']))
            f.write('\n')
