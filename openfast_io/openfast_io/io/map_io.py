"""
MAPIO – read / write MAP++ mooring input files.

Produces ``{'MAP': m}``
"""
from __future__ import annotations

import os
from typing import Any, Dict, Optional, Callable

from .base import ModuleIO
from ..parsing import float_read


class MAPIO(ModuleIO):
    """Read / write MAP++ input files."""

    def read(
        self,
        file_path: str,
        base_dir: str = '',
        **kwargs,
    ) -> dict:
        m: Dict[str, Any] = {}
        map_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        f = open(map_file)

        # --- LINE DICTIONARY ---
        f.readline()  # dashed header line
        f.readline()  # column names
        f.readline()  # units
        for k in ['LineType', 'Diam', 'MassDenInAir', 'EA', 'CB', 'CIntDamp', 'Ca', 'Cdn', 'Cdt']:
            m[k] = []
        data_line = f.readline().strip().split()
        while data_line and data_line[0][:3] != '---':
            m['LineType'].append(str(data_line[0]))
            m['Diam'].append(float_read(data_line[1]))
            m['MassDenInAir'].append(float_read(data_line[2]))
            m['EA'].append(float_read(data_line[3]))
            m['CB'].append(float_read(data_line[4]))
            m['CIntDamp'].append(float_read(data_line[5]))
            m['Ca'].append(float_read(data_line[6]))
            m['Cdn'].append(float_read(data_line[7]))
            m['Cdt'].append(float_read(data_line[8]))
            data_line = f.readline().strip().split()

        # --- NODE PROPERTIES ---
        f.readline()  # column names
        f.readline()  # units
        for k in ['Node', 'Type', 'X', 'Y', 'Z', 'M', 'B', 'FX', 'FY', 'FZ']:
            m[k] = []
        data_node = f.readline().strip().split()
        while data_node and data_node[0][:3] != '---':
            m['Node'].append(int(data_node[0]))
            m['Type'].append(str(data_node[1]))
            m['X'].append(float_read(data_node[2]))
            m['Y'].append(float_read(data_node[3]))
            m['Z'].append(float_read(data_node[4]))
            m['M'].append(float_read(data_node[5]))
            m['B'].append(float_read(data_node[6]))
            m['FX'].append(float_read(data_node[7]))
            m['FY'].append(float_read(data_node[8]))
            m['FZ'].append(float_read(data_node[9]))
            data_node = f.readline().strip().split()

        # --- LINE PROPERTIES ---
        f.readline()  # column names
        f.readline()  # units
        for k in ['Line', 'LineType_prop', 'UnstrLen', 'NodeAnch', 'NodeFair', 'Flags']:
            m[k] = []
        data_lp = f.readline().strip().split()
        while data_lp and data_lp[0][:3] != '---':
            m['Line'].append(int(data_lp[0]))
            m['LineType_prop'].append(str(data_lp[1]))
            m['UnstrLen'].append(float_read(data_lp[2]))
            m['NodeAnch'].append(int(data_lp[3]))
            m['NodeFair'].append(int(data_lp[4]))
            m['Flags'].append([str(val) for val in data_lp[5:]])
            data_lp = f.readline().strip().split()

        # --- SOLVER OPTIONS ---
        f.readline()  # column names (Option)
        f.readline()  # units (-)
        m['Option'] = []
        data_solver = f.readline().strip().split()
        while len(data_solver) > 0:
            m['Option'].append([str(val) for val in data_solver])
            data_solver = f.readline().strip().split()

        f.close()
        return {'MAP': m}

    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        **kwargs,
    ) -> None:
        m = data['MAP']
        with open(file_path, 'w') as f:
            f.write('---------------------- LINE DICTIONARY ---------------------------------------\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['LineType', 'Diam', 'MassDenInAir', 'EA', 'CB', 'CIntDamp', 'Ca', 'Cdn', 'Cdt']]) + '\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(m)', '(kg/m)', '(N)', '(-)', '(Pa-s)', '(-)', '(-)', '(-)']]) + '\n')
            for i in range(len(m.get('Diam', []))):
                ln = []
                ln.append('{:^11}'.format(m['LineType'][i]))
                ln.append('{:^11}'.format(m['Diam'][i]))
                ln.append('{:^11}'.format(m['MassDenInAir'][i]))
                ln.append('{:^11}'.format(m['EA'][i]))
                ln.append('{:<11}'.format(m['CB'][i]))
                ln.append('{:<11}'.format(m['CIntDamp'][i]))
                ln.append('{:<11}'.format(m['Ca'][i]))
                ln.append('{:<11}'.format(m['Cdn'][i]))
                ln.append('{:<11}'.format(m['Cdt'][i]))
                f.write(" ".join(ln) + '\n')

            f.write('---------------------- NODE PROPERTIES ---------------------------------------\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['Node', 'Type', 'X', 'Y', 'Z', 'M', 'B', 'FX', 'FY', 'FZ']]) + '\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(m)', '(m)', '(kg)', '(m^3)', '(N)', '(N)', '(N)']]) + '\n')
            for i in range(len(m.get('Node', []))):
                ln = []
                ln.append('{:<11}'.format(m['Node'][i]))
                ln.append('{:<11}'.format(m['Type'][i]))
                ln.append('{:<11}'.format(m['X'][i]))
                ln.append('{:<11}'.format(m['Y'][i]))
                ln.append('{:<11}'.format(m['Z'][i]))
                ln.append('{:<11}'.format(m['M'][i]))
                ln.append('{:<11}'.format(m['B'][i]))
                ln.append('{:<11}'.format(m['FX'][i]))
                ln.append('{:<11}'.format(m['FY'][i]))
                ln.append('{:<11}'.format(m['FZ'][i]))
                f.write(" ".join(ln) + '\n')

            f.write('---------------------- LINE PROPERTIES ---------------------------------------\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['Line', 'LineType', 'UnstrLen', 'NodeAnch', 'NodeFair', 'Flags']]) + '\n')
            f.write(" ".join(['{:<11s}'.format(i) for i in ['(-)', '(-)', '(m)', '(-)', '(-)', '(-)']]) + '\n')
            for i in range(len(m.get('Line', []))):
                ln = []
                ln.append('{:^11d}'.format(m['Line'][i]))
                ln.append('{:^11}'.format(m['LineType_prop'][i]))
                ln.append('{:^11}'.format(m['UnstrLen'][i]))
                ln.append('{:^11d}'.format(m['NodeAnch'][i]))
                ln.append('{:^11d}'.format(m['NodeFair'][i]))
                ln.append('{:<11}'.format(" ".join(m['Flags'][i])))
                f.write(" ".join(ln) + '\n')

            f.write('---------------------- SOLVER OPTIONS-----------------------------------------\n')
            f.write('{:<11s}'.format('Option') + '\n')
            f.write('{:<11s}'.format('(-)') + '\n')
            for opt in m.get('Option', []):
                f.write(" ".join(opt) + '\n')
            f.write('\n')
