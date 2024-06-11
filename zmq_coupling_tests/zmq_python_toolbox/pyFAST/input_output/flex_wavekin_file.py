import numpy as np
import pandas as pd
import os
import re
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass
from .csv_file import CSVFile


class FLEXWaveKinFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.wko'] #'.001 etc..'

    @staticmethod
    def formatName():
        return 'FLEX WaveKin file'

    def _read(self):
        numeric_const_pattern = r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
        rx = re.compile(numeric_const_pattern, re.VERBOSE)
        def extract_floats(s):
            v=np.array(rx.findall(s))
            return v


        try:
            csv = CSVFile(self.filename, sep=' ', commentLines=list(np.arange(11)),detectColumnNames=False)
        except:
            raise WrongFormatError('Unable to parse Flex WaveKin file as CSV with 11 header lines')

        header = csv.header 
        self['header'] = csv.header[0:2]
        self['data'] = csv.data
        try:
            self['MaxCrestHeight'] = float(extract_floats(header[2])[0])
            self['MaxLongiVel']    = float(extract_floats(header[3])[0])
            self['MaxLongiAcc']    = float(extract_floats(header[4])[0])
            dat = extract_floats(header[5]).astype(float)
            self['WaterDepth'] = dat[0]
            self['Hs']         = dat[1]
            self['Tp']         = dat[2]
            self['SpecType']   = dat[3]
        except:
            raise BrokenFormatError('Unable to parse floats from header lines 3-6')

        try:
            nDisp = int(extract_floats(header[6])[0])
            nRelD = int(extract_floats(header[8])[0])
        except:
            raise BrokenFormatError('Unable to parse int from header lines 7 and 9')

        try:
            displ = extract_floats(header[7]).astype(float)
            depth = extract_floats(header[9]).astype(float)
        except:
            raise BrokenFormatError('Unable to parse displacements or depths from header lines 8 and 10')
        if len(displ)!=nDisp:
            print(displ)
            raise BrokenFormatError('Number of displacements ({}) does not match number provided ({})'.format(nDisp, len(displ)))
        if len(depth)!=nRelD:
            print(depth)
            raise BrokenFormatError('Number of rel depth ({}) does not match number provided ({})'.format(nRelD, len(depth)))

        self['RelDepth']      = depth
        self['Displacements'] = displ

        cols=['Time_[s]', 'WaveElev_[m]']
        for j,x in enumerate(displ):
            for i,z in enumerate(depth):
                cols+=['u_z={:.1f}_x={:.1f}_[m/s]'.format(z*self['WaterDepth']*-1,x)]
            for i,z in enumerate(depth):
                cols+=['a_z={:.1f}_x={:.1f}_[m/s^2]'.format(z*self['WaterDepth'],x)]

        if len(cols)!=len(self['data'].columns):
            raise BrokenFormatError('Number of columns not valid')
        self['data'].columns = cols

#     def _write(self):
#         with open(self.filename,'w') as f:
#             f.write(self.toString)

    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)

        for k in ['MaxCrestHeight','MaxLongiVel','MaxLongiAcc','WaterDepth','Hs','Tp','SpecType','RelDepth','Displacements']:
            s += '{:15s}: {}\n'.format(k,self[k])
        if len(self['header'])>0:
            s += 'header        : '+ '  ,'.join(self['header'])+'\n'
        if len(self['data'])>0:
            s += 'data size     : {}x{}'.format(self['data'].shape[0],self['data'].shape[1])
        return s

    def _toDataFrame(self):
        return self['data']

