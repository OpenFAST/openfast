"""BeamDyn module IO — reads and writes BeamDyn input files.

Extracted from FAST_reader.py and FAST_writer.py.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from .base import ModuleIO
from ..parsing import (
    bool_read,
    float_read,
    int_read,
    quoted_read,
)


def _float_default_out(val, trim=False):
    if isinstance(val, float):
        return '{:.4f}'.format(val) if trim else '{: 22f}'.format(val)
    else:
        return '{:}'.format(val) if trim else '{:<22}'.format(val)


def _int_default_out(val, trim=False):
    if isinstance(val, int):
        return '{:22d}'.format(val) if not trim else '{:d}'.format(val)
    else:
        return '{:}'.format(val) if trim else '{:<22}'.format(val)


def _get_outlist(outlist_dict, channel_list):
    def loop_dict(vartree, outlist_i):
        for var in vartree.keys():
            if isinstance(vartree[var], dict):
                loop_dict(vartree[var], outlist_i)
            else:
                if vartree[var]:
                    outlist_i.append(var)
        return outlist_i

    if not channel_list:
        channel_list = outlist_dict.keys()
    outlist = []
    for var in channel_list:
        var = var.replace(' ', '')
        outlist_i = loop_dict(outlist_dict[var], [])
        if outlist_i:
            outlist.append(sorted(outlist_i))
    return outlist


class BeamDynIO(ModuleIO):
    """Reads and writes BeamDyn input files (per blade).

    read() returns::

        {
            'BeamDyn': { ... main BD params ... },
            'BeamDynBlade': { ... blade material props ... },
        }

    Note: Unlike the monolithic reader which stores per-blade arrays indexed by
    blade number, this IO handles one blade at a time.
    """

    def read(self, file_path: Path, base_dir: Path, *,
             outlist: dict | None = None,
             read_outlist_fn=None) -> dict:
        bd = {}
        file_path = str(file_path)

        f = open(file_path)
        f.readline(); f.readline(); f.readline()

        # Simulation control
        bd['Echo']             = bool_read(f.readline().split()[0])
        bd['QuasiStaticInit']  = bool_read(f.readline().split()[0])
        bd['rhoinf']           = float_read(f.readline().split()[0])
        bd['quadrature']       = int_read(f.readline().split()[0])
        bd['refine']           = int_read(f.readline().split()[0])
        bd['n_fact']           = int_read(f.readline().split()[0])
        bd['DTBeam']           = float_read(f.readline().split()[0])
        bd['load_retries']     = int_read(f.readline().split()[0])
        bd['NRMax']            = int_read(f.readline().split()[0])
        bd['stop_tol']         = float_read(f.readline().split()[0])
        bd['tngt_stf_fd']      = bool_read(f.readline().split()[0])
        bd['tngt_stf_comp']    = bool_read(f.readline().split()[0])
        bd['tngt_stf_pert']    = float_read(f.readline().split()[0])
        bd['tngt_stf_difftol'] = float_read(f.readline().split()[0])
        bd['RotStates']        = bool_read(f.readline().split()[0])
        f.readline()

        # Geometry
        bd['member_total'] = int_read(f.readline().split()[0])
        bd['kp_total']     = int_read(f.readline().split()[0])
        bd['members']      = []
        for i in range(bd['member_total']):
            ln = f.readline().split()
            n_pts_i = int(ln[1])
            member_i = {}
            member_i['kp_xr']         = [None] * n_pts_i
            member_i['kp_yr']         = [None] * n_pts_i
            member_i['kp_zr']         = [None] * n_pts_i
            member_i['initial_twist'] = [None] * n_pts_i
            f.readline(); f.readline()
            for j in range(n_pts_i):
                ln = f.readline().split()
                member_i['kp_xr'][j]         = float(ln[0])
                member_i['kp_yr'][j]         = float(ln[1])
                member_i['kp_zr'][j]         = float(ln[2])
                member_i['initial_twist'][j] = float(ln[3])
            bd['members'].append(member_i)

        # Mesh
        f.readline()
        bd['order_elem'] = int_read(f.readline().split()[0])

        # Material
        f.readline()
        bd['BldFile'] = f.readline().split()[0].replace('"', '').replace("'", '')

        # Outputs
        f.readline()
        bd['SumPrint']  = bool_read(f.readline().split()[0])
        bd['OutFmt']    = quoted_read(f.readline().split()[0])
        bd['NNodeOuts'] = int_read(f.readline().split()[0])
        bd['OutNd']     = [idx.strip() for idx in f.readline().split('OutNd')[0].split(',')]

        # OutList
        f.readline()
        if read_outlist_fn is not None and outlist is not None:
            read_outlist_fn(f, 'BeamDyn')
        else:
            line = f.readline()
            while line and 'END' not in line.split('!')[0].upper()[:3]:
                line = f.readline()

        # Optional nodal output
        try:
            f.readline()
            bd['BldNd_BlOutNd'] = f.readline().split()[0]
            f.readline()
            if read_outlist_fn is not None and outlist is not None:
                read_outlist_fn(f, 'BeamDyn_Nodes')
            else:
                line = f.readline()
                while line and 'END' not in line.split('!')[0].upper()[:3]:
                    line = f.readline()
        except Exception:
            pass

        f.close()

        # Read blade material file
        blade_file = os.path.join(os.path.dirname(file_path), bd['BldFile'])
        bd_blade = self._read_blade(blade_file)

        return {'BeamDyn': bd, 'BeamDynBlade': bd_blade}

    def write(self, data: dict, file_path: Path, base_dir: Path, *,
              naming_out: str = 'openfast',
              outlist: dict | None = None,
              fst_bd_filename: str | None = None) -> None:
        bd = data['BeamDyn']
        bd_blade = data['BeamDynBlade']
        run_dir = str(base_dir)

        # Write blade material file first
        blade_file = os.path.abspath(os.path.join(run_dir, bd['BldFile']))
        self._write_blade(bd_blade, blade_file)

        f = open(str(file_path), 'w')

        f.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(bd['Echo'], 'Echo', '- Echo input data to "<RootName>.ech" (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd['QuasiStaticInit'], 'QuasiStaticInit', '- Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd['rhoinf'], 'rhoinf', '- Numerical damping parameter for generalized-alpha integrator\n'))
        f.write('{:<22d} {:<11} {:}'.format(bd['quadrature'], 'quadrature', '- Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n'))
        f.write('{:<22} {:<11} {:}'.format(bd['refine'], 'refine', '- Refinement factor for trapezoidal quadrature (-) [DEFAULT = 1; used only when quadrature=2]\n'))
        f.write('{:<22} {:<11} {:}'.format(bd['n_fact'], 'n_fact', '- Factorization frequency for the Jacobian in N-R iteration(-) [DEFAULT = 5]\n'))
        f.write(_float_default_out(bd['DTBeam']) + '   {:<11} {:}'.format('DTBeam', '- Time step size (s).\n'))
        f.write(_int_default_out(bd['load_retries']) + '   {:<11} {:}'.format('load_retries', '- Number of factored load retries before quitting the aimulation [DEFAULT = 20]\n'))
        f.write(_int_default_out(bd['NRMax']) + '   {:<11} {:}'.format('NRMax', '- Max number of iterations in Newton-Raphson algorithm (-). [DEFAULT = 10]\n'))
        f.write(_float_default_out(bd['stop_tol']) + '   {:<11} {:}'.format('stop_tol', '- Tolerance for stopping criterion (-) [DEFAULT = 1E-5]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd['tngt_stf_fd'], 'tngt_stf_fd', '- Use finite differenced tangent stiffness matrix? (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd['tngt_stf_comp'], 'tngt_stf_comp', '- Compare analytical finite differenced tangent stiffness matrix? (flag)\n'))
        f.write(_float_default_out(bd['tngt_stf_pert']) + '   {:<11} {:}'.format('tngt_stf_pert', '- Perturbation size for finite differencing (-) [DEFAULT = 1E-6]\n'))
        f.write(_float_default_out(bd['tngt_stf_difftol']) + '   {:<11} {:}'.format('tngt_stf_difftol', '- Maximum allowable relative difference between analytical and fd tangent stiffness (-); [DEFAULT = 0.1]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(bd['RotStates'], 'RotStates', '- Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n'))
        f.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(bd['member_total'], 'member_total', '- Total number of members (-)\n'))
        f.write('{:<22d} {:<11} {:}'.format(bd['kp_total'], 'kp_total', '- Total number of key points (-) [must be at least 3]\n'))
        for i in range(bd['member_total']):
            mem = bd['members'][i]
            f.write('{:<22} {:<11} {:}'.format(' '.join(['%d' % (i + 1), '%d' % len(mem['kp_xr'])]), '', '- Member number; Number of key points in this member\n'))
            f.write(" ".join(['{:^21s}'.format(h) for h in ['kp_xr', 'kp_yr', 'kp_zr', 'initial_twist']]) + '\n')
            f.write(" ".join(['{:^21s}'.format(h) for h in ['(m)', '(m)', '(m)', '(deg)']]) + '\n')
            for j in range(len(mem['kp_xr'])):
                ln = ['{: 2.14e}'.format(mem['kp_xr'][j]),
                      '{: 2.14e}'.format(mem['kp_yr'][j]),
                      '{: 2.14e}'.format(mem['kp_zr'][j]),
                      '{: 2.14e}'.format(mem['initial_twist'][j])]
                f.write(" ".join(ln) + '\n')
        f.write('---------------------- MESH PARAMETER ------------------------------------------\n')
        f.write('{:<22d} {:<11} {:}'.format(bd['order_elem'], 'order_elem', '- Order of interpolation (basis) function (-)\n'))
        f.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format('"' + bd['BldFile'] + '"', 'BldFile', '- Name of file containing properties for blade (quoted string)\n'))
        f.write('---------------------- OUTPUTS -------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(bd['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"' + bd['OutFmt'] + '"', 'OutFmt', '- Format used for text tabular output, excluding the time channel.\n'))
        f.write('{:<22} {:<11} {:}'.format(bd['NNodeOuts'], 'NNodeOuts', '- Number of nodes to output to file [0 - 9] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(bd['OutNd']), 'OutNd', '- Nodes whose values will be output  (-)\n'))
        f.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')

        if outlist is not None:
            ol = _get_outlist(outlist, ['BeamDyn'])
            for channel_list in ol:
                for ch in channel_list:
                    f.write('"' + ch + '"\n')
        f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')

        if 'BldNd_BlOutNd' in bd:
            f.write('====== Outputs for all blade stations (same ending as above for B1N1.... =========================== [optional section]\n')
            f.write('{!s:<22} {:<11} {:}'.format(bd['BldNd_BlOutNd'], 'BldNd_BlOutNd', '- Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)\n'))
            f.write('                   OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, BeamDyn_Nodes tab for a listing of available output channels, (-)\n')
            if outlist is not None:
                opt_ol = _get_outlist(outlist, ['BeamDyn_Nodes'])
                for opt_channel_list in opt_ol:
                    for ch in opt_channel_list:
                        f.write('"' + ch + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')

        f.write('---------------------------------------------------------------------------------------')
        f.flush()
        os.fsync(f)
        f.close()

    # ------------------------------------------------------------------
    # PRIVATE: Blade material read/write
    # ------------------------------------------------------------------

    @staticmethod
    def _read_blade(blade_file):
        bld = {}
        f = open(blade_file)
        f.readline(); f.readline(); f.readline()

        bld['station_total'] = int_read(f.readline().split()[0])
        bld['damp_type']     = int_read(f.readline().split()[0])
        f.readline(); f.readline(); f.readline()

        # Stiffness-proportional damping
        ln = f.readline().split()
        bld['mu1'] = float(ln[0])
        bld['mu2'] = float(ln[1])
        bld['mu3'] = float(ln[2])
        bld['mu4'] = float(ln[3])
        bld['mu5'] = float(ln[4])
        bld['mu6'] = float(ln[5])
        f.readline()

        # Modal damping
        n_modes = int(f.readline().split()[0])
        bld['n_modes'] = n_modes
        bld['zeta'] = np.array(f.readline().strip().replace(',', ' ').split()[:n_modes], dtype=float).tolist()
        f.readline()

        # Distributed properties
        bld['radial_stations'] = np.zeros(bld['station_total'])
        bld['beam_stiff']      = np.zeros((bld['station_total'], 6, 6))
        bld['beam_inertia']    = np.zeros((bld['station_total'], 6, 6))
        for i in range(bld['station_total']):
            bld['radial_stations'][i] = float_read(f.readline().split()[0])
            for j in range(6):
                bld['beam_stiff'][i, j, :] = np.array([float(val) for val in f.readline().strip().split()])
            f.readline()
            for j in range(6):
                bld['beam_inertia'][i, j, :] = np.array([float(val) for val in f.readline().strip().split()])
            f.readline()

        f.close()
        return bld

    @staticmethod
    def _write_blade(bld, blade_file):
        from pathlib import Path as _Path
        _Path(blade_file).parent.mkdir(parents=True, exist_ok=True)
        f = open(blade_file, 'w')

        f.write('------- BEAMDYN INDIVIDUAL BLADE INPUT FILE --------------------------\n')
        f.write('Generated with OpenFAST_IO\n')
        f.write('------ Blade Parameters --------------------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(bld['station_total'], 'station_total', '- Number of blade input stations (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(bld['damp_type'], 'damp_type', '- Damping type (switch) {0: none, 1: stiffness-proportional, 2: modal}\n'))
        f.write('------ Stiffness-Proportional Damping [used only if damp_type=1] ---------------\n')
        f.write(" ".join(['{:^11s}'.format(h) for h in ['mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6']]) + '\n')
        f.write(" ".join(['{:^11s}'.format(h) for h in ['(-)', '(-)', '(-)', '(-)', '(-)', '(-)']]) + '\n')
        mu = [bld['mu1'], bld['mu2'], bld['mu3'], bld['mu4'], bld['mu5'], bld['mu6']]
        f.write(" ".join(['{:^11f}'.format(m) for m in mu]) + '\n')
        f.write('------ Modal Damping [used only if damp_type=2] --------------------------------\n')
        f.write('{:<22} {:<11} {:}\n'.format(bld['n_modes'], 'n_modes', '- Number of modal damping coefficients (-)'))
        f.write('{:<22} {:<11} {:}\n'.format(" ".join([repr(v) for v in bld['zeta']]), 'zeta', ' - Damping coefficients for mode 1 through n_modes'))
        f.write('------ Distributed Properties --------------------------------------------------\n')
        for i in range(len(bld['radial_stations'])):
            f.write('{: 2.15e}\n'.format(bld['radial_stations'][i]))
            for j in range(6):
                f.write(" ".join(['{: 2.15e}'.format(v) for v in bld['beam_stiff'][i, j, :]]) + '\n')
            f.write('\n')
            for j in range(6):
                f.write(" ".join(['{: 2.15e}'.format(v) for v in bld['beam_inertia'][i, j, :]]) + '\n')
            f.write('\n')

        f.write('\n')
        f.flush()
        os.fsync(f)
        f.close()
