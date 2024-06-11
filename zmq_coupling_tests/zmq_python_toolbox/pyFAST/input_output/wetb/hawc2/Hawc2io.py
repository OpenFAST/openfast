# -*- coding: utf-8 -*-
"""
Author:
    Bjarne S. Kallesoee


Description:
    Reads all HAWC2 output data formats, HAWC2 ascii, HAWC2 binary and FLEX

call ex.:
    # creat data file object, call without extension, but with parth
    file = ReadHawc2("HAWC2ex/tests")
    # if called with ReadOnly = 1 as
    file = ReadHawc2("HAWC2ex/tests",ReadOnly=1)
    # no channels a stored in memory, otherwise read channels are stored for reuse

    # channels are called by a list
    file([0,2,1,1])  => channels 1,3,2,2
    # if empty all channels are returned
    file()  => all channels as 1,2,3,...
    file.t => time vector

1. version: 19/4-2011
2. version: 5/11-2015 fixed columns to get description right, fixed time vector (mmpe@dtu.dk)

Need to be done:
    * add error handling for allmost every thing

"""
import numpy as np
import os

#from wetb import gtsdf

# FIXME: numpy doesn't like io.open binary fid in PY27, why is that? As a hack
# workaround, use opent for PY23 compatibility when handling text files,
# and default open for binary

################################################################################
################################################################################
################################################################################
# Read HAWC2 class
################################################################################
class ReadHawc2(object):
    """
    """
################################################################################
# read *.sel file
    def _ReadSelFile(self):
        """
        Some title
        ==========

        Using docstrings formatted according to the reStructuredText specs
        can be used for automated documentation generation with for instance
        Sphinx: http://sphinx.pocoo.org/.

        Parameters
        ----------
        signal : ndarray
            some description

        Returns
        -------
        output : int
            describe variable
        """

        # read *.sel hawc2 output file for result info
        if self.FileName.lower().endswith('.sel'):
            self.FileName = self.FileName[:-4]
        fid = open(self.FileName + '.sel', 'r')
        Lines = fid.readlines()
        fid.close()
        if Lines[0].lower().find('bhawc')>=0:
            # --- Find line with scan info
            iLine=0
            for i in np.arange(5,10):
                if Lines[i].lower().find('scans')>=0:
                    iLine=i+1
            if iLine==0:
                raise Exception('Cannot find the keyword "scans"')
            temp = Lines[iLine].split()
            self.NrSc = int(temp[0])
            self.NrCh = int(temp[1])
            self.Time = float(temp[2])
            self.Freq = self.NrSc / self.Time
            self.t = np.linspace(0, self.Time, self.NrSc + 1)[1:]
            # --- Find line with channel info
            iLine=0
            for i in np.arange(5,13):
                if Lines[i].lower().find('channel')>=0:
                    iLine=i+1
            if iLine==0:
                raise Exception('Cannot find the keyword "Channel"')

            # reads channel info (name, unit and description)
            Name = []; Unit = []; Description = [];
            for i in range(0, self.NrCh+1):
                if (i+iLine)>=len(Lines):
                    break
                line = Lines[i + iLine].strip()
                if len(line)==0:
                    continue
                # --- removing number and unit
                sp=[sp.strip() for sp in line.split() if len(sp.strip())>0]
                num   = sp[0]
                iNum  = line.find(num)
                line  = line[iNum+len(num)+1:]
                unit  = sp[-1]
                iUnit = line.find(unit)
                line  = line[:iUnit]
                # --- Splitting to find label and description
                sp=[sp.strip() for sp in line.split('\t') if len(sp.strip())>0]
                if len(sp)!=2:
                    for nSpaces in np.arange(2,15):
                        sp=[sp.strip() for sp in line.split(' '*nSpaces) if len(sp.strip())>0]
                        if len(sp)==2:
                            break
                if len(sp)!=2:
                    raise Exception('Dont know how to split the input of the sel file into 4 columns')

                Unit.append(unit)
                Description.append(sp[0]) 
                Name.append(sp[1]) 

            self.ChInfo = [Name, Unit, Description]
            self.FileFormat = 'BHAWC_ASCII'
        else:

            # findes general result info (number of scans, number of channels,
            # simulation time and file format)
            temp = Lines[8].split()
            self.NrSc = int(temp[0])
            self.NrCh = int(temp[1])
            self.Time = float(temp[2])
            self.Freq = self.NrSc / self.Time
            self.t = np.linspace(0, self.Time, self.NrSc + 1)[1:]
            Format = temp[3]
            # reads channel info (name, unit and description)
            Name = []; Unit = []; Description = [];
            for i in range(0, self.NrCh):
                temp = str(Lines[i + 12][12:43]); Name.append(temp.strip())
                temp = str(Lines[i + 12][43:54]); Unit.append(temp.strip())
                temp = str(Lines[i + 12][54:-1]); Description.append(temp.strip())
            self.ChInfo = [Name, Unit, Description]
            # if binary file format, scaling factors are read
            if Format.lower() == 'binary':
                self.ScaleFactor = np.zeros(self.NrCh)
                self.FileFormat = 'HAWC2_BINARY'
                for i in range(0, self.NrCh):
                    self.ScaleFactor[i] = float(Lines[i + 12 + self.NrCh + 2])
            else:
                self.FileFormat = 'HAWC2_ASCII'
################################################################################
# read sensor file for FLEX format
    def _ReadSensorFile(self):
        # read sensor file used if results are saved in FLEX format
        DirName = os.path.dirname(self.FileName)
        try:
            fid = opent(DirName + r"\sensor ", 'r')
        except IOError:
            print("can't finde sensor file for FLEX format")
            return
        Lines = fid.readlines()
        fid.close()
        # reads channel info (name, unit and description)
        self.NrCh = 0
        Name = []
        Unit = []
        Description = []
        for i in range(2, len(Lines)):
            temp = Lines[i]
            if not temp.strip():
                break
            self.NrCh += 1
            temp = str(Lines[i][38:45])
            Unit.append(temp.strip())
            temp = str(Lines[i][45:53])
            Name.append(temp.strip())
            temp = str(Lines[i][53:])
            Description.append(temp.strip())
        self.ChInfo = [Name, Unit, Description]
        # read general info from *.int file
        fid = open(self.FileName, 'rb')
        fid.seek(4 * 19)
        if not np.fromfile(fid, 'int32', 1) == self.NrCh:
            print("number of sensors in sensor file and data file are not consisten")
        fid.seek(4 * (self.NrCh) + 4, 1)
        self.Version = np.fromfile(fid, 'int32',1)[0]
        time_start, time_step = np.fromfile(fid, 'f', 2)
        self.Freq = 1 / time_step
        self.ScaleFactor = np.fromfile(fid, 'f', self.NrCh)
        fid.seek(2 * 4 * self.NrCh + 48 * 2)
        self.NrSc = int(len(np.fromfile(fid, 'int16')) / self.NrCh)
        self.Time = self.NrSc * time_step
        self.t = np.arange(0, self.Time, time_step) + time_start
        fid.close()
################################################################################
# init function, load channel and other general result file info
    def __init__(self, FileName, ReadOnly=0):
        self.FileName = FileName
        self.ReadOnly = ReadOnly
        self.Iknown = []  # to keep track of what has been read all ready
        self.Data = np.zeros(0)
        if FileName.lower().endswith('.sel') or os.path.isfile(FileName + ".sel"):
             self._ReadSelFile()
        elif FileName.lower().endswith('.dat') and os.path.isfile(os.path.splitext(FileName)[0] + ".sel"):
             self.FileName = os.path.splitext(FileName)[0]
             self._ReadSelFile()
        elif FileName.lower().endswith('.int') or FileName.lower().endswith('.res'):
             self.FileFormat = 'FLEX'
             self._ReadSensorFile()
        elif os.path.isfile(self.FileName + ".int"):
             self.FileName = self.FileName + ".int"
             self.FileFormat = 'FLEX'
             self._ReadSensorFile()
        elif os.path.isfile(self.FileName + ".res"):
             self.FileName = self.FileName + ".res"
             self.FileFormat = 'FLEX'
             self._ReadSensorFile()
        elif FileName.lower().endswith('.hdf5') or os.path.isfile(self.FileName + ".hdf5"):
            self.FileFormat = 'GTSDF'
            self.ReadGtsdf()
        else:
            raise Exception("unknown file: " + FileName)
################################################################################
# Read results in binary format
    def ReadBinary(self, ChVec=None):
        ChVec = [] if ChVec is None else ChVec
        if not ChVec:
            ChVec = range(0, self.NrCh)
        with open(self.FileName + '.dat', 'rb') as fid:
            data = np.zeros((self.NrSc, len(ChVec)))
            j = 0
            for i in ChVec:
                fid.seek(i * self.NrSc * 2, 0)
                data[:, j] = np.fromfile(fid, 'int16', self.NrSc) * self.ScaleFactor[i]
                j += 1
        return data
################################################################################
# Read results in ASCII format
    def ReadAscii(self, ChVec=None):
        ChVec = [] if ChVec is None else ChVec
        if not ChVec:
            ChVec = range(0, self.NrCh)
        temp = np.loadtxt(self.FileName + '.dat', usecols=ChVec)
        return temp.reshape((self.NrSc, len(ChVec)))
################################################################################
# Read results in FLEX format
    def ReadFLEX(self, ChVec=None):
        ChVec = [] if ChVec is None else ChVec
        if not ChVec:
            ChVec = range(1, self.NrCh)
        fid = open(self.FileName, 'rb')
        fid.seek(2 * 4 * self.NrCh + 48 * 2)
        temp = np.fromfile(fid, 'int16')
        if self.Version==3:
            temp = temp.reshape(self.NrCh, self.NrSc).transpose()
        else:
            temp = temp.reshape(self.NrSc, self.NrCh)
        fid.close()
        return np.dot(temp[:, ChVec], np.diag(self.ScaleFactor[ChVec]))
################################################################################
# Read results in GTSD format
    def ReadGtsdf(self):
        raise NotImplementedError
        #self.t, data, info = gtsdf.load(self.FileName + '.hdf5')
        #self.Time = self.t
        #self.ChInfo = [['Time'] + info['attribute_names'],
        #               ['s'] + info['attribute_units'],
        #               ['Time'] + info['attribute_descriptions']]
        #self.NrCh = data.shape[1] + 1
        #self.NrSc = data.shape[0]
        #self.Freq = self.NrSc / self.Time
        #self.FileFormat = 'GTSDF'
        #self.gtsdf_description = info['description']
        #data = np.hstack([self.Time[:,np.newaxis], data])
        #return data
################################################################################
# One stop call for reading all data formats
    def ReadAll(self, ChVec=None):
        ChVec = [] if ChVec is None else ChVec
        if not ChVec and not self.FileFormat == 'GTSDF':
            ChVec = range(0, self.NrCh)
        if self.FileFormat == 'HAWC2_BINARY':
            return self.ReadBinary(ChVec)
        elif self.FileFormat == 'HAWC2_ASCII' or self.FileFormat == 'BHAWC_ASCII':
            return self.ReadAscii(ChVec)
        elif self.FileFormat == 'GTSDF':
            return self.ReadGtsdf()
        elif self.FileFormat == 'FLEX':
            return self.ReadFLEX(ChVec)
        else:
            raise Exception('Unknown file format {} for hawc2 out file'.format(self.FileFormat))

################################################################################
# Main read data call, read, save and sort data
    def __call__(self, ChVec=None):
        ChVec = [] if ChVec is None else ChVec
        if not ChVec:
            ChVec = range(0, self.NrCh)
        elif max(ChVec) >= self.NrCh:
            print("to high channel number")
            return
        # if ReadOnly, read data but no storeing in memory
        if self.ReadOnly:
            return self.ReadAll(ChVec)
        # if not ReadOnly, sort in known and new channels, read new channels
        # and return all requested channels
        else:
            # sort into known channels and channels to be read
            I1 = []
            I2 = []  # I1=Channel mapping, I2=Channels to be read
            for i in ChVec:
                try:
                    I1.append(self.Iknown.index(i))
                except:
                    self.Iknown.append(i)
                    I2.append(i)
                    I1.append(len(I1))
            # read new channels
            if I2:
                temp = self.ReadAll(I2)
                # add new channels to Data
                if self.Data.any():
                    self.Data = np.append(self.Data, temp, axis=1)
                # if first call, so Daata is empty
                else:
                    self.Data = temp
            return self.Data[:, tuple(I1)]


################################################################################
################################################################################
################################################################################
# write HAWC2 class, to be implemented
################################################################################

if __name__ == '__main__':
    res_file = ReadHawc2('structure_wind')
    results = res_file.ReadAscii()
    channelinfo = res_file.ChInfo
