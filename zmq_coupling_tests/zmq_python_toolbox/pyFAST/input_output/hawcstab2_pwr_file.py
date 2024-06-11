import numpy as np
import pandas as pd
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass


class HAWCStab2PwrFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.pwr', '.txt']

    @staticmethod
    def formatName():
        return 'HAWCStab2 power file'

    def _read(self):
        # Reading header line
        with open(self.filename,'r',encoding=self.encoding) as f:
            header = f.readline().strip()
        if len(header)<=0 or header[0]!='#':
            raise WrongFormatError('Pwr File {}: header line does not start with `#`'.format(self.filename)+e.args[0])
        # Extracting column names
        header       = '0 '+header[1:].strip()
        num_and_cols = [s.strip()+']' for s in header.split(']')[:-1]]
        cols         = [(' '.join(col.split(' ')[1:])).strip().replace(' ','_')  for col in num_and_cols]
        # Determining type based on number of columns (NOTE: could use col names as well maybe)
        if len(cols)!=15:
            raise WrongFormatError('Pwr File {}: '.format(self.filename))
        self.colNames=cols
        # Reading numerical data
        try:
            self.data = np.loadtxt(self.filename, skiprows=1)
        except Exception as e:    
            raise BrokenFormatError('Pwr File {}: '.format(self.filename)+e.args[0])

        if self.data.shape[1]!=len(cols):
            raise BrokenFormatError('Pwr File {}: inconsistent number of header columns and data columns.'.format(self.filename)+e.args[0])

    #def _write(self):
        #self.data.to_csv(self.filename,sep=self.false,index=False)

    def _toDataFrame(self):
        return pd.DataFrame(data=self.data, columns=self.colNames)

