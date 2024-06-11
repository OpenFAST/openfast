""" 
Tools to read/write OpenFAST output files

Main content:

- class FASTOutputFile()
- data, info = def load_output(filename)
- data, info = def load_ascii_output(filename)
- data, info = def load_binary_output(filename, use_buffer=True)
- def writeDataFrame(df, filename, binary=True)
- def writeBinary(fileName, channels, chanNames, chanUnits, fileID=2, descStr='')

NOTE: 
  - load_binary and writeBinary are not "fully reversible" for now.
      Some small numerical errors are introduced in the conversion.
      Some of the error is likely due to the fact that Python converts to "int" and "float" (double).
      Maybe all the operations should be done in single. I tried but failed. 
      I simply wonder if the operation is perfectly reversible. 
             

"""
from itertools import takewhile
import numpy as np
import pandas as pd
import struct
import ctypes
import os
import re
try:
    from .file import File, WrongFormatError, BrokenReaderError, EmptyFileError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class WrongReaderError(Exception): pass
    class BrokenFormatError(Exception): pass
    class EmptyFileError(Exception): pass
try:
    from .csv_file import CSVFile
except:
    print('CSVFile not available')



FileFmtID_WithTime              = 1 # File identifiers used in FAST
FileFmtID_WithoutTime           = 2
FileFmtID_NoCompressWithoutTime = 3
FileFmtID_ChanLen_In            = 4 # Channel length included in file


# --------------------------------------------------------------------------------}
# --- OUT FILE 
# --------------------------------------------------------------------------------{
class FASTOutputFile(File):
    """ 
    Read an OpenFAST ouput file (.out, .outb, .elev).

    Main methods
    ------------
    - read, write, toDataFrame

    Examples
    --------

        # read an output file, convert it to pandas dataframe, modify it, write it back
        f = FASTOutputFile('5MW.outb')
        df=f.toDataFrame()
        time  = df['Time_[s]']
        Omega = df['RotSpeed_[rpm]'] 
        df['Time_[s]'] -=100
        f.writeDataFrame(df, '5MW_TimeShifted.outb')

    """

    @staticmethod
    def defaultExtensions():
        return ['.out','.outb','.elm','.elev','.dbg','.dbg2']

    @staticmethod
    def formatName():
        return 'FAST output file'

    def __init__(self, filename=None, **kwargs):
        """ Class constructor. If a `filename` is given, the file is read. """
        self.filename = filename
        if filename:
            self.read(**kwargs)

    def read(self, filename=None, **kwargs):
        """ Reads the file self.filename, or `filename` if provided """
        
        # --- Standard tests and exceptions (generic code)
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)
        # --- Calling (children) function to read
        self._read(**kwargs)

    def _read(self):
        def readline(iLine):
            with open(self.filename) as f:
                for i, line in enumerate(f):
                    if i==iLine-1:
                        return line.strip()
                    elif i>=iLine:
                        break

        ext = os.path.splitext(self.filename.lower())[1]
        self.info={}
        self['binary']=False
        try:
            if ext in ['.out','.elev','.dbg','.dbg2']:
                self.data, self.info = load_ascii_output(self.filename)
            elif ext=='.outb':
                self.data, self.info = load_binary_output(self.filename)
                self['binary']=True
            elif ext=='.elm':
                F=CSVFile(filename=self.filename, sep=' ', commentLines=[0,2],colNamesLine=1)
                self.data = F.data
                del F
                self.info['attribute_units']=readline(3).replace('sec','s').split()
                self.info['attribute_names']=self.data.columns.values
            else:
                if isBinary(self.filename):
                    self.data, self.info = load_binary_output(self.filename)
                    self['binary']=True
                else:
                    self.data, self.info = load_ascii_output(self.filename)
                    self['binary']=False
        except MemoryError as e:    
            raise BrokenReaderError('FAST Out File {}: Memory error encountered\n{}'.format(self.filename,e))
        except Exception as e:    
            raise WrongFormatError('FAST Out File {}: {}'.format(self.filename,e.args))
        if self.data.shape[0]==0:
            raise EmptyFileError('This FAST output file contains no data: {}'.format(self.filename))

        if self.info['attribute_units'] is not None:
            self.info['attribute_units'] = [re.sub(r'[()\[\]]','',u) for u in self.info['attribute_units']]


    def _write(self, binary=None, fileID=4): 
        if binary is None:
            binary = self['binary']

        if binary:
            # NOTE: user provide a filename, we allow overwrite
            self.toOUTB(filename=self.filename, fileID=fileID, noOverWrite=False)
        else:
            # ascii output
            with open(self.filename,'w') as f:
                f.write('\t'.join(['{:>10s}'.format(c)         for c in self.info['attribute_names']])+'\n')
                f.write('\t'.join(['{:>10s}'.format('('+u+')') for u in self.info['attribute_units']])+'\n')
                # TODO better..
                f.write('\n'.join(['\t'.join(['{:10.4f}'.format(y[0])]+['{:10.3e}'.format(x) for x in y[1:]]) for y in self.data]))

    def toDataFrame(self):
        """ Returns object into one DataFrame, or a dictionary of DataFrames"""
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
                raise BrokenFormatError('Inconstistent number of columns between headers ({}) and data ({}) for file {}'.format(len(cols), self.data.shape[1], self.filename))
            df = pd.DataFrame(data=self.data,columns=cols)

        return df

    def writeDataFrame(self, df, filename, binary=True):
        writeDataFrame(df, filename, binary=binary)

    def __repr__(self):
        s='<{} object> with attributes:\n'.format(type(self).__name__)
        s+=' - info ({})\n'.format(type(self.info))
        s+=' - data ({})\n'.format(type(self.data))
        s+='and keys: {}\n'.format(self.keys())
        return s

    # --------------------------------------------------------------------------------
    # --- Converters 
    # --------------------------------------------------------------------------------
    def toOUTB(self, filename=None, extension='.outb', fileID=4, noOverWrite=True, **kwargs):
        #NOTE: we override the File class here
        if filename is None:
            base, _ = os.path.splitext(self.filename)
            filename = base + extension
        else:
            base, ext = os.path.splitext(filename)
            if len(ext)!=0:
                extension = ext
        if (filename==self.filename) and noOverWrite:
            raise Exception('Not overwritting {}. Specify a filename or an extension.'.format(filename))
        
        # NOTE: fileID=2 will chop the channels name of long channels use fileID4 instead
        channels = self.data
        chanNames = self.info['attribute_names']
        chanUnits = self.info['attribute_units']
        descStr   = self.info['description']
        if isinstance(descStr, list):
            descStr=(''.join(descStr[:2])).replace('\n','')
        writeBinary(filename, channels, chanNames, chanUnits, fileID=fileID, descStr=descStr)


# --------------------------------------------------------------------------------
# --- Helper low level functions 
# --------------------------------------------------------------------------------
def isBinary(filename):
    with open(filename, 'r') as f:
        try:
            # first try to read as string
            l = f.readline()
            # then look for weird characters
            for c in l:
                code = ord(c)
                if code<10 or (code>14 and code<31):
                    return True
            return False
        except UnicodeDecodeError:
            return True

        



def load_ascii_output(filename, method='numpy', encoding='ascii'):


    if method in ['forLoop','pandas']:
        from .file import numberOfLines
        nLines = numberOfLines(filename, method=2)

    with open(filename, encoding=encoding, errors='ignore') as f:
        info = {}
        info['name'] = os.path.splitext(os.path.basename(filename))[0]
        # Header is whatever is before the keyword `time`
        header = []
        maxHeaderLines=35
        headerRead = False
        for i in range(maxHeaderLines):
            l = f.readline()
            if not l:
                raise Exception('Error finding the end of FAST out file header. Keyword Time missing.')
            # Check for utf-16
            if l[:3] == '\x00 \x00':
                f.close()
                encoding=''
                print('[WARN] Attempt to re-read the file with encoding utf-16')
                return load_ascii_output(filename=filename, method=method, encoding='utf-16')
            first_word = (l+' dummy').lower().split()[0]
            in_header=  (first_word != 'time') and  (first_word != 'alpha')
            if in_header:
                header.append(l)
            else:
                info['description'] = header
                info['attribute_names'] = l.split()
                info['attribute_units'] = [unit[1:-1] for unit in f.readline().split()]
                headerRead=True
                break
        if not headerRead:
            raise WrongFormatError('Could not find the keyword "Time" or "Alpha" in the first {} lines of the file'.format(maxHeaderLines))

        nHeader = len(header)+1
        nCols = len(info['attribute_names'])

        if method=='numpy':
            # The most efficient, and will remove empty lines and the lines that starts with "This"
            #  ("This" is found at the end of some Hydro Out files..)
            data = np.loadtxt(f, comments=('This'))

        elif method =='pandas':
            # Could probably be made more efficient, but 
            f.close()
            nRows = nLines-nHeader
            sep=r'\s+'
            cols= ['C{}'.format(i) for i in range(nCols)]
            df = pd.read_csv(filename, sep=sep, header=0, skiprows=nHeader, names=cols, dtype=float, na_filter=False, nrows=nRows, engine='pyarrow'); print(df)
            data=df.values

        elif method == 'forLoop':
            # The most inefficient
            nRows = nLines-nHeader
            sep=r'\s+'
            data = np.zeros((nRows, nCols))
            for i in range(nRows):
                l = f.readline().strip()
                sp = np.array(l.split()).astype(float)
                data[i,:] = sp[:nCols]

        elif method == 'listCompr':
            # --- Method 4 - List comprehension
            # Data, up to end of file or empty line (potential comment line at the end)
            data = np.array([l.strip().split() for l in takewhile(lambda x: len(x.strip())>0, f.readlines())]).astype(float)
        else:
            raise NotImplementedError()

    return data, info


def load_binary_output(filename, use_buffer=True):
    """
    03/09/15: Ported from ReadFASTbinary.m by Mads M Pedersen, DTU Wind
    24/10/18: Low memory/buffered version by E. Branlard, NREL
    18/01/19: New file format for exctended channels, by E. Branlard, NREL

    Info about ReadFASTbinary.m:
    % Author: Bonnie Jonkman, National Renewable Energy Laboratory
    % (c) 2012, National Renewable Energy Laboratory
    %
    %  Edited for FAST v7.02.00b-bjj  22-Oct-2012
    """
    StructDict = {
            'uint8': ('B', 1, np.uint8), 
            'int16':('h', 2, np.int16), 
            'int32':('i', 4, np.int32), 
            'float32':('f', 4, np.float32),
            'float64':('d', 8, np.float64)}
    def fread(fid, n, dtype):
        fmt, nbytes, npdtype = StructDict[dtype]
        #return np.array(struct.unpack(fmt * n, fid.read(nbytes * n)), dtype=npdtype)
        return struct.unpack(fmt * n, fid.read(nbytes * n))

    def freadRowOrderTableBuffered(fid, n, type_in, nCols, nOff=0, type_out='float64'):
        """ 
        Reads of row-ordered table from a binary file.

        Read `n` data of type `type_in`, assumed to be a row ordered table of `nCols` columns.
        Memory usage is optimized by allocating the data only once.
        Buffered reading is done for improved performances (in particular for 32bit python)

        `nOff` allows for additional column space at the begining of the storage table.
        Typically, `nOff=1`, provides a column at the beginning to store the time vector.

        @author E.Branlard, NREL

        """
        fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type_in]
        nLines          = int(n/nCols)
        GoodBufferSize  = 4096*40
        nLinesPerBuffer = int(GoodBufferSize/nCols)
        BufferSize      = nCols * nLinesPerBuffer
        nBuffer         = int(n/BufferSize)
        # Allocation of data
        data = np.zeros((nLines,nCols+nOff), dtype = type_out)
        # Reading
        try:
            nIntRead   = 0
            nLinesRead = 0
            while nIntRead<n:
                nIntToRead = min(n-nIntRead, BufferSize)
                nLinesToRead = int(nIntToRead/nCols)
                Buffer = np.array(struct.unpack(fmt * nIntToRead, fid.read(nbytes * nIntToRead)))
                Buffer = Buffer.reshape(-1,nCols)
                data[ nLinesRead:(nLinesRead+nLinesToRead),  nOff:(nOff+nCols)  ] = Buffer
                nLinesRead = nLinesRead + nLinesToRead
                nIntRead   = nIntRead   + nIntToRead
        except:
            raise Exception('Read only %d of %d values in file:' % (nIntRead, n, filename))
        return data

    with open(filename, 'rb') as fid:
        #----------------------------        
        # get the header information
        #----------------------------

        FileID = fread(fid, 1, 'int16')[0]  #;             % FAST output file format, INT(2)

        if FileID not in [FileFmtID_WithTime, FileFmtID_WithoutTime, FileFmtID_NoCompressWithoutTime, FileFmtID_ChanLen_In]:
            raise Exception('FileID not supported {}. Is it a FAST binary file?'.format(FileID))

        if FileID == FileFmtID_ChanLen_In: 
            LenName = fread(fid, 1, 'int16')[0] # Number of characters in channel names and units
        else:
            LenName = 10                    # Default number of characters per channel name

        NumOutChans = fread(fid, 1, 'int32')[0]  #;             % The number of output channels, INT(4)
        NT = fread(fid, 1, 'int32')[0]  #;             % The number of time steps, INT(4)

        if FileID == FileFmtID_WithTime:
            TimeScl = fread(fid, 1, 'float64')[0]  # The time slopes for scaling, REAL(8)
            TimeOff = fread(fid, 1, 'float64')[0]  # The time offsets for scaling, REAL(8)
        else:
            TimeOut1 = fread(fid, 1, 'float64')[0]  # The first time in the time series, REAL(8)
            TimeIncr = fread(fid, 1, 'float64')[0]  # The time increment, REAL(8)

        if FileID == FileFmtID_NoCompressWithoutTime:
            ColScl = np.ones ((NumOutChans, 1)) # The channel slopes for scaling, REAL(4)
            ColOff = np.zeros((NumOutChans, 1)) # The channel offsets for scaling, REAL(4)
        else:
            ColScl = fread(fid, NumOutChans, 'float32')  # The channel slopes for scaling, REAL(4)
            ColOff = fread(fid, NumOutChans, 'float32')  # The channel offsets for scaling, REAL(4)

        LenDesc      = fread(fid, 1, 'int32')[0]  #;  % The number of characters in the description string, INT(4)
        DescStrASCII = fread(fid, LenDesc, 'uint8')  #;  % DescStr converted to ASCII
        DescStr      = "".join(map(chr, DescStrASCII)).strip()

        ChanName = []  # initialize the ChanName cell array
        for iChan in range(NumOutChans + 1):
            ChanNameASCII = fread(fid, LenName, 'uint8')  #; % ChanName converted to numeric ASCII
            ChanName.append("".join(map(chr, ChanNameASCII)).strip())

        ChanUnit = []  # initialize the ChanUnit cell array
        for iChan in range(NumOutChans + 1):
            ChanUnitASCII = fread(fid, LenName, 'uint8')  #; % ChanUnit converted to numeric ASCII
            ChanUnit.append("".join(map(chr, ChanUnitASCII)).strip()[1:-1])

        # -------------------------
        #  get the channel time series
        # -------------------------

        nPts = NT * NumOutChans  #;           % number of data points in the file

        if FileID == FileFmtID_WithTime:
            PackedTime = fread(fid, NT, 'int32')  #; % read the time data
            cnt = len(PackedTime)
            if cnt < NT:
                raise Exception('Could not read entire %s file: read %d of %d time values' % (filename, cnt, NT))

        if use_buffer:
            # Reading data using buffers, and allowing an offset for time column (nOff=1)
            if FileID == FileFmtID_NoCompressWithoutTime:
                data = freadRowOrderTableBuffered(fid, nPts, 'float64', NumOutChans, nOff=1, type_out='float64')
            else:
                data = freadRowOrderTableBuffered(fid, nPts, 'int16', NumOutChans, nOff=1, type_out='float64')
        else:
            # NOTE: unpacking huge data not possible on 32bit machines
            if FileID == FileFmtID_NoCompressWithoutTime:
                PackedData = fread(fid, nPts, 'float64')  #; % read the channel data
            else:
                PackedData = fread(fid, nPts, 'int16')  #; % read the channel data

            cnt = len(PackedData)
            if cnt < nPts:
                raise Exception('Could not read entire %s file: read %d of %d values' % (filename, cnt, nPts))
            data = np.array(PackedData).reshape(NT, NumOutChans)
            del PackedData

    if FileID == FileFmtID_WithTime:
        time = (np.array(PackedTime) - TimeOff) / TimeScl;
    else:
        time = TimeOut1 + TimeIncr * np.arange(NT)

    # -------------------------
    #  Scale the packed binary to real data
    # -------------------------
    if use_buffer:
        # Scaling Data
        for iCol in range(NumOutChans):
            if np.isnan(ColScl[iCol]) and np.isnan(ColOff[iCol]):
                data[:,iCol+1] = 0 # probably due to a division by zero in Fortran
            else:
                data[:,iCol+1] = (data[:,iCol+1] - ColOff[iCol]) / ColScl[iCol]
        # Adding time column
        data[:,0] = time
    else:
        # NOTE: memory expensive due to time conversion, and concatenation
        data = (data - ColOff) / ColScl
        data = np.concatenate([time.reshape(NT, 1), data], 1)

    info = {'name': os.path.splitext(os.path.basename(filename))[0],
            'description': DescStr,
            'fileID': FileID,
            'attribute_names': ChanName,
            'attribute_units': ChanUnit}
    return data, info


def writeBinary(fileName, channels, chanNames, chanUnits, fileID=4, descStr=''):
    """
    Write an OpenFAST binary file.

    Based on contributions from
        Hugo Castro, David Schlipf, Hochschule Flensburg
    
    Input:
     FileName      - string: contains file name to open
     Channels      - 2-D array: dimension 1 is time, dimension 2 is channel 
     ChanName      - cell array containing names of output channels
     ChanUnit      - cell array containing unit names of output channels, preferably surrounded by parenthesis
     FileID        - constant that determines if the time is stored in the
                     output, indicating possible non-constant time step
     DescStr       - String describing the file
    """
    # Data sanitization
    chanNames = list(chanNames)
    channels  = np.asarray(channels)
    if chanUnits[0][0]!='(':
        chanUnits = ['('+u+')' for u in chanUnits] # units surrounded by parenthesis to match OpenFAST convention

    nT, nChannelsWithTime = np.shape(channels)
    nChannels             = nChannelsWithTime - 1

    # For FileID =2,  time needs to be present and at the first column
    try:
        iTime = chanNames.index('Time')
    except ValueError:
        raise Exception('`Time` needs to be present in channel names' )
    if iTime!=0:
        raise Exception('`Time` needs to be the first column of `chanName`' )

    time = channels[:,iTime]
    timeStart = time[0]
    timeIncr = (time[-1]-time[0])/(nT-1)
    dataWithoutTime = channels[:,1:]
        
    # Compute data range, scaling and offsets to convert to int16
    #   To use the int16 range to its fullest, the max float is matched to 2^15-1 and the
    #   the min float is matched to -2^15. Thus, we have the to equations we need
    #   to solve to get scaling and offset, see line 120 of ReadFASTbinary:
    #   Int16Max = FloatMax * Scaling + Offset 
    #   Int16Min = FloatMin * Scaling + Offset
    int16Max   = np.single( 32767.0)         # Largest integer represented in 2 bytes,  2**15 - 1
    int16Min   = np.single(-32768.0)         # Smallest integer represented in 2 bytes -2**15
    int16Rng   = np.single(int16Max - int16Min)  # Max Range of 2 byte integer
    mins   = np.min(dataWithoutTime, axis=0)
    ranges = np.single(np.max(dataWithoutTime, axis=0) - mins)
    ranges[ranges==0]=1  # range set to 1 for constant channel. In OpenFAST: /sqrt(epsilon(1.0_SiKi))
    ColScl  = np.single(int16Rng/ranges)
    ColOff  = np.single(int16Min - np.single(mins)*ColScl)
    
    #Just available for fileID 
    if fileID not in [2,4]:
        print("current version just works with fileID = 2 or 4")

    else:
        with open(fileName,'wb') as fid:
            # Notes on struct:
            # @ is used for packing in native byte order
            #  B - unsigned integer 8 bits
            #  h - integer 16 bits
            #  i - integer 32 bits
            #  f - float 32 bits
            #  d - float 64 bits

            # Write header informations
            fid.write(struct.pack('@h',fileID))
            if fileID == FileFmtID_ChanLen_In: 
                maxChanLen = np.max([len(s) for s in chanNames])
                maxUnitLen = np.max([len(s) for s in chanUnits])
                nChar = max(maxChanLen, maxUnitLen)
                fid.write(struct.pack('@h',nChar))
            else:
                nChar = 10

            fid.write(struct.pack('@i',nChannels))
            fid.write(struct.pack('@i',nT))
            fid.write(struct.pack('@d',timeStart))
            fid.write(struct.pack('@d',timeIncr))
            fid.write(struct.pack('@{}f'.format(nChannels), *ColScl))
            fid.write(struct.pack('@{}f'.format(nChannels), *ColOff))
            descStrASCII = [ord(char) for char in descStr]
            fid.write(struct.pack('@i',len(descStrASCII)))
            fid.write(struct.pack('@{}B'.format(len((descStrASCII))), *descStrASCII))

            # Write channel names
            for chan in chanNames:
                chan = chan[:nChar]
                ordchan = [ord(char) for char in chan] + [32]*(nChar-len(chan))
                fid.write(struct.pack('@'+str(nChar)+'B', *ordchan))

            # Write channel units
            for unit in chanUnits:
                unit = unit[:nChar]
                ordunit = [ord(char) for char in unit] + [32]*(nChar-len(unit))
                fid.write(struct.pack('@'+str(nChar)+'B', *ordunit))

            # --- Pack and write data
            # Method 1
            #packedData=np.zeros((nT, nChannels), dtype=np.int16)
            #for iChan in range(nChannels):
            #    packedData[:,iChan] = np.clip( ColScl[iChan]*dataWithoutTime[:,iChan]+ColOff[iChan], int16Min, int16Max)
            #packedData = packedData.ravel()
            ## NOTE: the *packedData converts to a tuple before passing to struct.pack
            ##       which is inefficient
            #fid.write(struct.pack('@{}h'.format(packedData.size), *packedData))

            # --- Method 2
            #packedData=np.zeros((nT, nChannels), dtype=np.int16)
            #for iChan in range(nChannels):
            #    packedData[:,iChan] = np.clip( ColScl[iChan]*dataWithoutTime[:,iChan]+ColOff[iChan], int16Min, int16Max)
            #packedData = packedData.ravel()
            ## Here we use slice assignment
            #buf = (ctypes.c_int16 * len(packedData))()
            #buf[:] = packedData
            #fid.write(buf)

            # --- Method 3 use packedData as slice directly
            packedData = (ctypes.c_int16 * (nT*nChannels))()
            for iChan in range(nChannels):
                packedData[iChan::nChannels] = np.clip( ColScl[iChan]*dataWithoutTime[:,iChan]+ColOff[iChan], int16Min, int16Max).astype(np.int16)
            fid.write(packedData)


            fid.close()

def writeDataFrame(df, filename, binary=True):
    """ write a DataFrame to OpenFAST output format"""
    channels  = df.values
    # attempt to extract units from channel names
    chanNames=[]
    chanUnits=[]
    for c in df.columns:
        c     = c.strip()
        name  = c
        unit = ''
        if c[-1]==']':
            chars=['[',']']
        elif c[-1]==')':
            chars=['(',')']
        else:
            chars=[]
        if len(chars)>0:
            op,cl = chars
            iu=c.rfind(op)
            if iu>1:
                name = c[:iu]
                unit = c[iu+1:].replace(cl,'')
                if name[-1]=='_':
                    name=name[:-1]
                
        chanNames.append(name)
        chanUnits.append(unit)

    if binary:
        writeBinary(filename, channels, chanNames, chanUnits, fileID=FileFmtID_ChanLen_In)
    else:
        NotImplementedError()



if __name__ == "__main__":
    B=FASTOutputFile('tests/example_files/FASTOutBin.outb')
    df=B.toDataFrame()
    B.writeDataFrame(df, 'tests/example_files/FASTOutBin_OUT.outb')
    B.toOUTB(extension='.dat.outb')
    B.toParquet()
    B.toCSV()


