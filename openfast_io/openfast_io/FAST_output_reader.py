#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
Created on 03/09/2015
@author: MMPE
Copied from https://github.com/WISDEM/AeroelasticSE/tree/openmdao1/src/AeroelasticSE/old_files on 15 Aug 2016 by Ganesh Vijayakumar
'''
import os
import numpy as np
import struct

class FASTOutputFile():
    """ 
    Read an OpenFAST ouput file (.out, .outb, .elev).

    Main methods
    ------------
    - read, write, toDataFrame
    """

    def __init__(self, filename=None, method='numpy'):
        """
        Load a FAST binary or ascii output file

        Parameters
        ----------
        filename : str
            filename

        Returns
        -------
        data: ndarray
            data values
        info: dict
            info containing:
                - name: filename
                - description: description of dataset
                - attribute_names: list of attribute names
                - attribute_units: list of attribute units
        """

        assert os.path.isfile(filename), "File, %s, does not exists" % filename

        ext = os.path.splitext(filename)[1]
        if ext in ['.out']:
            with open(filename, 'r') as f:
                try:
                    f.readline()
                except UnicodeDecodeError:
                    self.data, self.info = load_binary_output(filename)

            self.data, self.info = load_ascii_output(filename)

        elif ext == '.outb':
            self.data, self.info, self.pack = load_binary_output(filename)

        else:
            raise ValueError(f'File extension, {ext}, not supported')
        
        if method == 'pandas':
            self.toDataFrame()
            


    def toDataFrame(self):
        import pandas as pd
        """
        Returns object into one DataFrame, or a dictionary of DataFrames
        Borrowed from openfast_toolbox
        """
        # --- Example (returning one DataFrame):
        #  return pd.DataFrame(data=np.zeros((10,2)),columns=['Col1','Col2'])
        if self.info['attribute_units'] is not None:
            if len(self.info['attribute_names'])!=len(self.info['attribute_units']):
                cols=self.info['attribute_names']
                print('[WARN] not all columns have units! Skipping units')
            else:
                cols=[n+'_['+u.replace('sec','s')+']' for n,u in zip(self.info['attribute_names'],self.info['attribute_units'])]
        else:
            cols=self.info['attribute_names']
        if isinstance(self.data, pd.DataFrame):
            df= self.data
            df.columns=cols
        else:
            if len(cols)!=self.data.shape[1]:
                raise Warning('Inconstistent number of columns between headers ({}) and data ({}) for file {}'.format(len(cols), self.data.shape[1], self.filename))
            df = pd.DataFrame(data=self.data,columns=cols)

        return df


def load_ascii_output(filename, headerLines = 8, descriptionLine = 4, attributeLine = 6, unitLine = 7, delimiter = None):
    with open(filename) as f:
        info = {}
        info['name'] = os.path.splitext(os.path.basename(filename))[0]
        header = [f.readline() for _ in range(headerLines)]
        info['description'] = header[descriptionLine].strip()
        info['attribute_names'] = header[attributeLine].strip().split(delimiter)
        info['attribute_units'] = [unit[1:-1] for unit in header[unitLine].strip().split(delimiter)]  #removing "()"
        data = np.array([line.split(delimiter) for line in f.readlines()], dtype=float)
        return data, info

def load_binary_output(filename):
    """
    Ported from ReadFASTbinary.m by Mads M Pedersen, DTU Wind
    Info about ReadFASTbinary.m:
    Author: Bonnie Jonkman, National Renewable Energy Laboratory
    (c) 2012, National Renewable Energy Laboratory
    Edited for FAST v7.02.00b-bjj  22-Oct-2012
    """

    def fread(fid, n, type):
        fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
        return struct.unpack(fmt * n, fid.read(nbytes * n))
    
    FileFmtID_WithTime = 1    # File identifiers used in FAST
    FileFmtID_WithoutTime = 2
    FileFmtID_NoCompressWithoutTime = 3
    FileFmtID_ChanLen_In = 4
    
    with open(filename, 'rb') as fid:
        FileID = fread(fid, 1, 'int16')[0]       # FAST output file format, INT(2)
        
        if FileID == FileFmtID_ChanLen_In:
            LenName  = fread(fid, 1, 'int16')[0] # Number of characters in channel names and units
        else:
            LenName = 10                         # default number of characters per channel name

        
        NumOutChans = fread(fid, 1, 'int32')[0]  # The number of output channels, INT(4)
        NT = fread(fid, 1, 'int32')[0]           # The number of time steps, INT(4)
        
        if FileID == FileFmtID_WithTime:
            TimeScl = fread(fid, 1, 'float64')   # The time slopes for scaling, REAL(8)
            TimeOff = fread(fid, 1, 'float64')   # The time offsets for scaling, REAL(8)
        else:
            TimeOut1 = fread(fid, 1, 'float64')  # The first time in the time series, REAL(8)
            TimeIncr = fread(fid, 1, 'float64')  # The time increment, REAL(8)

        if FileID != FileFmtID_NoCompressWithoutTime:
            ColScl = fread(fid, NumOutChans, 'float32')  # The channel slopes for scaling, REAL(4)
            ColOff = fread(fid, NumOutChans, 'float32')  # The channel offsets for scaling, REAL(4)

        LenDesc = fread(fid, 1, 'int32')[0]          # The number of characters in the description string, INT(4)
        DescStrASCII = fread(fid, LenDesc, 'uint8')  # DescStr converted to ASCII
        DescStr = "".join(map(chr, DescStrASCII)).strip()

        ChanName = []                                     # initialize the ChanName cell array
        for iChan in range(NumOutChans + 1):
            ChanNameASCII = fread(fid, LenName, 'uint8')  # ChanName converted to numeric ASCII
            ChanName.append("".join(map(chr, ChanNameASCII)).strip())

        ChanUnit = []                                     # initialize the ChanUnit cell array
        for iChan in range(NumOutChans + 1):
            ChanUnitASCII = fread(fid, LenName, 'uint8')  # ChanUnit converted to numeric ASCII
            ChanUnit.append("".join(map(chr, ChanUnitASCII)).strip()[1:-1])

        # get the channel time series
        nPts = NT * NumOutChans                   # number of data points in the file
        if FileID == FileFmtID_WithTime:
            PackedTime = fread(fid, NT, 'int32')  # read the time data
            cnt = len(PackedTime)
            if cnt < NT:
                raise Exception('Could not read entire %s file: read %d of %d time values' % (filename, cnt, NT))
        
        if FileID == FileFmtID_NoCompressWithoutTime:
            PackedData = fread(fid, nPts, 'float64')    # read the channel data
        else:
            PackedData = fread(fid, nPts, 'int16')    # read the channel data
            
        cnt = len(PackedData)
        if cnt < nPts:
            raise Exception('Could not read entire %s file: read %d of %d values' % (filename, cnt, nPts))

    if FileID == FileFmtID_NoCompressWithoutTime:
        pack = np.array(PackedData).reshape(NT, NumOutChans)
        data = pack
    else:
        # Scale the packed binary to real data
        pack = np.array(PackedData).reshape(NT, NumOutChans)
        data = (pack - ColOff) / ColScl

    if FileID == FileFmtID_WithTime:
        time = (np.array(PackedTime) - TimeOff) / TimeScl;
    else:
        time = TimeOut1 + TimeIncr * np.arange(NT)

    data = np.concatenate([time.reshape(NT, 1), data], 1)
    pack = np.concatenate([time.reshape(NT, 1), pack], 1)

    info = {'name': os.path.splitext(os.path.basename(filename))[0],
            'description': DescStr,
            'attribute_names': ChanName,
            'attribute_units': ChanUnit}
    return data, info, pack



if __name__=="__main__":

    from openfast_io.FileTools import check_rtest_cloned

    parent_dir = os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) + os.sep

    of_outputfile = os.path.join(parent_dir, 'reg_tests', 'r-test', 'glue-codes', 
                                 'openfast', '5MW_Land_BD_DLL_WTurb', '5MW_Land_BD_DLL_WTurb.outb')

    check_rtest_cloned(of_outputfile)

    d,i,p = load_binary_output(of_outputfile)

    print(tuple(i['attribute_names']))
    print(type(d))
    print(i)
    print(len(i['attribute_names']))
    print(np.shape(d))
