"""
AeroDiskIO – read / write AeroDisk input files.

Produces ``{'AeroDisk': ad}``

AeroDisk references an external CSV file containing CpCtCq tables.
"""
from __future__ import annotations

import os
from typing import Any, Callable, Dict, Optional

from .base import ModuleIO
from ..parsing import bool_read, float_read, quoted_read

try:
    from openfast_io.FAST_output_reader import load_ascii_output
except ImportError:
    load_ascii_output = None


class AeroDiskIO(ModuleIO):
    """Read / write AeroDisk input files."""

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
        ad: Dict[str, Any] = {}
        ad_dir = os.path.dirname(file_path) or '.'

        f = open(file_path)
        f.readline()
        f.readline()
        f.readline()

        # Simulation Control
        ad['Echo'] = bool_read(f.readline().split()[0])
        ad['DT'] = float_read(f.readline().split()[0])

        # Environmental Conditions
        f.readline()
        ad['AirDens'] = float_read(f.readline().split()[0])

        # Actuator Disk Properties
        f.readline()
        ad['RotorRad'] = float_read(f.readline().split()[0])

        # InColNames
        ad['InColNames'] = [x.strip() for x in quoted_read(f.readline().split('InColNames')[0]).split(',')]

        # InColDims
        ad['InColDims'] = [int(x) for x in f.readline().split('InColDims')[0].split(',')]

        # CSV file reference (starts with '@')
        line = f.readline()
        if line.strip().startswith('@'):
            csv_rel = line.strip()[1:].strip()
            ad['actuatorDiskFile'] = os.path.join(ad_dir, csv_rel)

            if load_ascii_output is not None:
                data, info = load_ascii_output(
                    ad['actuatorDiskFile'],
                    headerLines=3,
                    descriptionLine=0,
                    attributeLine=1,
                    unitLine=2,
                    delimiter=',',
                )
                ad['actuatorDiskTable'] = {
                    'dsc': info['description'],
                    'attr': info['attribute_names'],
                    'units': info['attribute_units'],
                    'data': data,
                }
            else:
                ad['actuatorDiskTable'] = {}
        else:
            raise Exception('Expecting a file reference to the actuator disk CSV file')

        # Output
        f.readline()
        f.readline()

        # Read output list
        ad['_outlist'] = {}
        data_line = f.readline()
        while data_line.split().__len__() == 0:
            data_line = f.readline()
        while data_line.split()[0] != 'END':
            if data_line.find('"') >= 0:
                channels = data_line.split('"')
                channel_list = channels[1].split(',')
            else:
                row_string = data_line.split(',')
                if len(row_string) == 1:
                    channel_list = row_string[0].split('\n')[0]
                else:
                    channel_list = row_string
            if isinstance(channel_list, list):
                for ch in channel_list:
                    ch = ch.strip()
                    if ch:
                        ad['_outlist'][ch] = True
            else:
                ch = channel_list.strip()
                if ch:
                    ad['_outlist'][ch] = True
            data_line = f.readline()

        f.close()
        return {'AeroDisk': ad}

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
        ad = data['AeroDisk']
        out_dir = os.path.dirname(file_path) or '.'

        # Write CSV table first
        csv_name = os.path.basename(ad.get('actuatorDiskFile', 'AeroDiskProp.csv'))
        csv_path = os.path.join(out_dir, csv_name)
        if 'actuatorDiskTable' in ad and ad['actuatorDiskTable']:
            tbl = ad['actuatorDiskTable']
            with open(csv_path, 'w') as cf:
                cf.write('{}\n'.format(tbl.get('dsc', '')))
                cf.write('{}\n'.format(', '.join(str(a) for a in tbl['attr'])))
                cf.write('{}\n'.format(', '.join(str(u) for u in tbl['units'])))
                for row in tbl['data']:
                    cf.write('{}\n'.format(', '.join('{:.6f}'.format(v) for v in row)))

        with open(file_path, 'w') as f:
            f.write('--- AERO DISK INPUT FILE -------\n')
            f.write('Generated with OpenFAST_IO\n')
            f.write('--- SIMULATION CONTROL ---------\n')
            f.write('{!s:<22} {:<11} {:}'.format(ad['Echo'], 'Echo', '- Echo input data to "<RootName>.ADsk.ech" (flag)\n'))
            f.write('{:<22} {:<11} {:}'.format(ad['DT'], 'DT', '- Integration time step (s)\n'))
            f.write('--- ENVIRONMENTAL CONDITIONS ---\n')
            f.write('{:<22f} {:<11} {:}'.format(ad['AirDens'], 'AirDens', '- Air density (kg/m^3) (or "default")\n'))
            f.write('--- ACTUATOR DISK PROPERTIES ---\n')
            f.write('{:<22f} {:<11} {:}'.format(ad['RotorRad'], 'RotorRad', '- Rotor radius (m) (or "default")\n'))
            f.write('"{}" {:<11} {:}'.format(
                ', '.join(ad['InColNames']),
                'InColNames',
                '- Input column headers\n',
            ))
            f.write('{:<22} {:<11} {:}'.format(
                ', '.join(str(d) for d in ad['InColDims']),
                'InColDims',
                '- Number of unique values in each column\n',
            ))
            f.write('@{}\n'.format(csv_name))
            f.write('--- OUTPUTS --------------------\n')
            f.write('{:<22} {:<11} {:}'.format('OutList', 'OutList', '- The next line(s) contains a list of output parameters.\n'))
            out_channels = ad.get('_outlist', {})
            if outlist and 'AeroDisk' in outlist:
                out_channels = outlist['AeroDisk']
            for ch in out_channels:
                f.write('"' + ch + '"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of the last OutList line)\n')
            f.write('---------------------------------------------------------------------------------------\n')
