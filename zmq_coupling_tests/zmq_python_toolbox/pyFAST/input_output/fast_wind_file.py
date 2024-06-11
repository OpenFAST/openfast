import numpy as np
import pandas as pd
from .csv_file import CSVFile
from .file import isBinary, WrongFormatError

class FASTWndFile(CSVFile):

    @staticmethod
    def defaultExtensions():
        return ['.wnd']

    @staticmethod
    def formatName():
        return 'FAST determ. wind file'

    def __init__(self, *args, **kwargs):
        self.colNames=['Time','WindSpeed','WindDir','VertSpeed','HorizShear','VertShear','LinVShear','GustSpeed']
        self.units=['[s]','[m/s]','[deg]','[m/s]','[-]','[-]','[-]','[m/s]']
        Cols=['{}_{}'.format(c,u) for c,u in zip(self.colNames,self.units)]

        header=[]
        header+=['!Wind file.']
        header+=['!Time  Wind     Wind	Vert.       Horiz.      Vert.       LinV        Gust']
        header+=['!      Speed    Dir    Speed       Shear		Shear       Shear       Speed']

        super(FASTWndFile, self).__init__(sep=' ',commentChar='!',colNames=Cols, header=header, *args, **kwargs)

    def _read(self, *args, **kwargs):
        if isBinary(self.filename):
            raise WrongFormatError('This is a binary file (turbulence file?) not a FAST ascii determinisctic wind file')
        super(FASTWndFile, self)._read(*args, **kwargs)

    def _write(self, *args, **kwargs):
        super(FASTWndFile, self)._write(*args, **kwargs)


    def _toDataFrame(self):
        return self.data


# --------------------------------------------------------------------------------}
# --- Functions specific to file type  
# --------------------------------------------------------------------------------{
    def stepWind(self,WSstep=1,WSmin=3,WSmax=25,tstep=100,dt=0.5,tmin=0,tmax=999):
        """ Set the wind file to a step wind 
        tstep: can be an array of size 2 [tstepmax tstepmin]

        
        """
            
        Steps= np.arange(WSmin,WSmax+WSstep,WSstep)
        if hasattr(tstep,'__len__'):
            tstep = np.around(np.linspace(tstep[0], tstep[1], len(Steps)),0)
        else:
            tstep = len(Steps)*[tstep]
        nCol = len(self.colNames)
        nRow = len(Steps)*2
        M = np.zeros((nRow,nCol));
        M[0,0] = tmin
        M[0,1] = WSmin
        for i,s in enumerate(Steps[:-1]):
            M[2*i+1,0] = tmin + tstep[i]-dt 
            M[2*i+2,0] = tmin + tstep[i]
            tmin +=tstep[i]
            M[2*i+1,1] = Steps[i]
            if i<len(Steps)-1:
                M[2*i+2,1] = Steps[i+1]
            else:
                M[2*i+2,1] = Steps[-1]
        M[-1,0]= max(tmax, tmin+tstep[-1])
        M[-1,1]= WSmax
        self.data=pd.DataFrame(data=M,columns=self.colNames)
