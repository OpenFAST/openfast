""" 
This module contains:
  - PolarFile: class to read different kind of polar formats

  - 
"""

import os
import numpy as np
import pandas as pd

# --- Welib polar readers
try:
    from pyFAST.input_output.csv_file import CSVFile
except:
    CSVFile=None
try:
    from pyFAST.input_output.fast_input_file import ADPolarFile
except:
    ADPolarFile=None
try:
    from pyFAST.input_output.hawc2_ae_file import HAWC2AEFile
except:
    HAWC2AEFile=None

class WrongPolarFormatError(Exception): pass
class BrokenPolarFormatError(Exception): pass


# List of columns used to "unify" the dataframes coming out of "PolarFile"
DEFAULT_COLUMNS={'alpha':'Alpha', 'cl':'Cl', 'cd':'Cd', 'cm':'Cm'}
DEFAULT_COLUMNS_EXT={
        'clinv':'Cl_inv', 'clfs':'Cl_fs', 'fs':'fs',
        'cn':'Cn', 'cnpot':'Cn_pot', 'cnz12':'Cn_012', 'cnf':'Cn_f', 'cncd0off':'Cn_Cd0off'
        }


# --------------------------------------------------------------------------------
# --- Small Helper functions
# --------------------------------------------------------------------------------
def _load_txt(filename, commentChars, skiprows=0, **kwargs):
    """ 
    Similar to np.loadtxt but also works if comments are present anywhere in the file (e.g. end of file)
    """
    with open(filename) as f:
        lines = (line for iline, line in enumerate(f) if not line.startswith(commentChars) and iline>=skiprows)
        Lines = list(lines) 
    if len(Lines)==0:
        raise Exception('Zero lines')
    else:
        return np.loadtxt(Lines, **kwargs)


# --------------------------------------------------------------------------------}
# --- Simple classes 
# --------------------------------------------------------------------------------{
class BasePolarFile(dict):
    def __init__(self, filename=None):
        super().__init__()
        self.COMMENT_CHARS=('#','!','%')
        self['header'] = ''
        self['columns'] = []
        self['data']    = np.array([[]])
        self['nPolars']  = 0
        if filename is not None:
            self.read(filename)

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='- header: {}\n'.format(self['header'])
        s+='- columns: {}\n'.format(self['columns'])
        s+='- nPolars:{}\n'.format(self['nPolars'])
        s+='- data: shape {}\n'.format(self['data'].shape)
        s+='        first: {}\n'.format(self['data'][0,:])
        s+='        last:  {}\n'.format(self['data'][-1,:])
        return s

    @staticmethod
    def formatName(): raise NotImplementedError()

    def toDataFrame(self):
        if self['nPolars']==1:
            return pd.DataFrame(data=self['data'], columns=self['columns'])
        else:
            raise NotImplementedError()

class PolarFile_OneLineHeader(BasePolarFile):
    """ Polar file with exatcly one line of header. 
    Column names in header can be separated by spaces or commas.
    Header may start with the following comment characters: ['#','!','%']
    Data may be space or column separated
    """
    @staticmethod
    def formatName(): return 'Polar file one header line'

    def read(self, filename):
        super().__init__()
        with open(filename) as f:
            header = f.readline().strip()
            second = f.readline()
        self['header'] = header
        for c in self.COMMENT_CHARS:
            header = header.lstrip(c)
        sep=',' 
        try:
            self['data'] = np.loadtxt(filename, delimiter=sep, skiprows=1)
        except:
            sep=None
            self['data'] = np.loadtxt(filename, delimiter=sep, skiprows=1)
        self['nPolars']=1

        # --- Detect columns
        nCols = self['data'].shape[1]
        # First, if all values are numeric, abort
        onestring = header.replace(',',' ')
        try:
            vals = np.array(onestring.split()).astype(float) # This should fail
        except:
            pass # Great, it actually failed, the first line is not made of floats
        else:
            raise WrongPolarFormatError('The first line is all numeric, it should contain column names')
        # Then, try to split by commas or space
        colsComma = header.split(',')
        colsSpace = header.split()
        if len(colsComma)==nCols:
            cols = [c.strip() for c in colsComma]
        elif len(colsSpace)==nCols:
            cols = colsSpace
        else:
            raise BrokenPolarFormatError('The number of header columns ({}) does not match the number of columns in the data ({})'.format(len(cols),nCols))
        self['columns'] = cols

class PolarFile_NoHeader(BasePolarFile):
    """ 
    Polar file with no header, or some "meaningless" comments that starts with ['#','!','%']
    Data may be space or column separated
    """
    @staticmethod
    def formatName(): return 'Polar file no header'

    def read(self, filename):
        self['data']    = _load_txt(filename, self.COMMENT_CHARS)
        self['nPolars'] = 1
        # --- Detect columns
        nCols = self['data'].shape[1]
        d = [DEFAULT_COLUMNS['alpha'], DEFAULT_COLUMNS['cl'], DEFAULT_COLUMNS['cd'], DEFAULT_COLUMNS['cm']]
        n2col = {2:d[0:2], 3:d[0:3], 4:d[0:4] }
        if nCols in n2col.keys():
            self['columns'] = n2col[nCols]
        else:
            raise BrokenPolarFormatError('The number of columns in the data ({}) is not amongst the supported ones ({}).'.format(nCols, n2col.keys()))

class PolarFile_AD_Basic(BasePolarFile):
    """ 
    Reads a basic AeroDyn file
    """
    @staticmethod
    def formatName(): return 'Polar AeroDyn file basic'

    def read(self, filename):
        self['data']    = _load_txt(filename, self.COMMENT_CHARS, skiprows = 53)
        self['nPolars'] = 1
        # import pandas as pd
        # df=pd.read_csv(filename, skiprows = 53, header=None, delim_whitespace=True, names=['Alpha','Cl','Cd','Cm']).values
        # --- Detect columns
        nCols = self['data'].shape[1]
        n2col = {2:['Alpha','Cl'], 3:['Alpha','Cl', 'Cm'], 4:['Alpha','Cl', 'Cm', 'Cd'] }
        if nCols in n2col.keys():
            self['columns'] = n2col[nCols]
        else:
            raise BrokenPolarFormatError('The number of columns in the data ({}) is not amongst the supported ones ({}).'.format(nCols, n2col.keys()))


class PolarFile(BasePolarFile):
    """ """
    @staticmethod
    def formatName(): return 'Polar file'


def loadPolarFile(filename, fformat='auto', to_radians=False, standardizeCols=True, verbose=False):
    """ 
    Loads a PolarFile, return a dataFrame
    """
    if not os.path.exists(filename):
        raise Exception('File not found:',filename)
        print('[WARN] Not all file formats supported ')

    allReaders   = [ADPolarFile, PolarFile_OneLineHeader, PolarFile_NoHeader, PolarFile_AD_Basic, CSVFile]
    delimReaders = [PolarFile_OneLineHeader, PolarFile_AD_Basic, CSVFile]

    def tryReading(f, reader):
        if f is not None:
            return f
        if reader is None:
            return None
        try:
            if verbose:
                print('')
                print('PolarFile: trying to read with format: {}'.format(reader.formatName()))
            return reader(filename)
        except:
            if verbose:
                print('>>> PolarFile: Failed to read with format: {}'.format(reader.formatName()))
            pass
    f = None
    Re = np.nan # TODO

    if fformat==None:
        fformat = 'auto'

    if fformat=='ADPolar':
        f = ADPolarFile(filename)

    elif fformat=='delimited':

        for reader in delimReaders:
            f = tryReading(f, reader)
            if f is not None:
                break

    elif fformat=='auto':

        for reader in allReaders:
            f = tryReading(f, reader)
            if f is not None:
                break

    if f is None:
        raise Exception('Unable to read the polar {} using the fileformat {}. Use a supported fileformat'.format(filename, fformat))

    # --- Store in DataFrame
    df = f.toDataFrame()
    if verbose:
        print('PolarFile: Columns before: ',df.columns.values)

    # --- Rename columns - Standardize column names
    if standardizeCols:
        COLS_TODO= {**DEFAULT_COLUMNS,**DEFAULT_COLUMNS_EXT}
        for ic, col in enumerate(df.columns):
            c = col.strip().lower().replace('_','')
            c = c.replace('aoa','alpha')
            c = c.replace('fst','fs')
            c = c.replace('012','z12')
            c = c.replace('cllin','clinv')
            c = c.replace('clpot','clinv')
            known_keys = reversed(sorted(list(COLS_TODO.keys())))
            found=False
            for kk in known_keys:
                if c.startswith(kk):
                    cnew = COLS_TODO.pop(kk)
                    #print('changing {} to {}'.format(c, cnew))
                    df.columns.values[ic] = cnew # rename column
                    found=True
                    break
            if not found:
                print('[WARN] PolarFile: The following column was not understood: {}'.format(col))

        # --- Standardize data
        for k,v in DEFAULT_COLUMNS.items():
            if v not in df.columns:
                 df[v] = np.nan
    if verbose:
        print('PolarFile: Columns after: ',df.columns.values)

    if standardizeCols:
        cAlpha = DEFAULT_COLUMNS['alpha']
        if cAlpha not in df.columns:
            raise Exception('Angle of attack was not detected as part of the columns')
    else:
        cAlpha = df.columns.values[0]

    if to_radians:
        # First, check the data, if the max alpha is above pi, most likely we are in degrees
        _radians = np.mean(np.abs(df[cAlpha])) <= np.pi / 2
        if _radians:
            raise Exception('PolarFile: Asked to convert input to radian, but the data is likely already in radians.')
        df[cAlpha]*=np.pi/180

    Re  = np.nan
    return df, Re 

if __name__ == "__main__":
    from welib.tools.clean_exceptions import *
#     PolarFile_OneLineHeader('data/63-235.csv')
    #f = PolarFile_NoHeader('data/63-235.csv')
#     f = loadPolarFile('data/63-235.csv')
    f = loadPolarFile('data/FFA-W3-241-Re12M.dat', verbose=True)
    #f = loadPolarFile('data/Cylinder.dat')
    #f = loadPolarFile('../../data/NREL5MW/5MW_Baseline/Airfoils/DU21_A17.dat')
    print(f)
    pass
