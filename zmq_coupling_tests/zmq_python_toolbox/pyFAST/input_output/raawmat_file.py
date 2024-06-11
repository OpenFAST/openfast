# imports
import pandas as pd
import numpy as np
import os


def matfile(path_file,output='pandas',path_key=None):
    ''' Function to parse RAAW .mat historical files with specific data structure

    Parameter
    ---------
    path_file : str
        File path to .mat file to read
    output : str
        Desired output format of 'pandas' or 'xarray'
    path_key : str (optional)
        File path to .csv channel_key for extra checks and consistent formatting
    
    Return
    ------
    ds : dict
        Dict with keys 'data50hz' and 'data1hz' with specified output format
    '''
    from scipy.io import loadmat

    # load in the mat file
    mfile = loadmat(path_file)

    # flatten the mat file
    data = [[row.flat[0] for row in line] for line in mfile['DataCell']]

    # create time column for easy comparison to simulations (0-600 seconds)
    time50hz = np.arange(600,step=0.02)
    time1hz = np.arange(600,step=1)
    names_50hz = ['time']
    names_1hz  = ['time']
    dat_50hz   = [time50hz]
    dat_1hz    = [time1hz]
    # start data reorganization into lists to be converted into pandas
    allnames = []
    for signal in range(mfile['DataCell'].size):
        name = data[signal][0][1][0] # get signal name
        allnames.append(name)
        if data[signal][0][3][0][0] == 0.02:
            names_50hz.append(name)
            dat_50hz.append(data[signal][0][6].squeeze())
        elif data[signal][0][3][0][0] == 1:
            names_1hz.append(name)
            dat_1hz.append(data[signal][0][6].squeeze())
        else:
            print('ERROR: '+name+' ignored!')

    # create time index for each dataframe
    tstamp = np.array2string(data[signal][0][5][0], separator=',')
    tstamp = tstamp.replace(' ','').replace('[','').replace(']','')
    t0 = pd.to_datetime(tstamp,format='%Y,%m,%d,%H,%M,%S')
    tstamp50hz = pd.date_range(start=t0,periods=30000,freq='20ms')
    tstamp1hz = pd.date_range(start=t0,periods=600,freq='1s')
    
    if not path_key==None:
        # load in list of channels
        channel_key = pd.read_csv(path_key)
        list_50hz = channel_key['50hz'].to_list()
        list_1hz = channel_key['1hz'].dropna().to_list()
        # check if file size matches channel list size
        if mfile['DataCell'].size > (len(list_50hz)+len(list_1hz)):
            print('WARNING: channel size in file is greater than size of channel names lists!')
        # find any missing channels from file
        missing50hz = list(set(list_50hz)-set(names_50hz))
        missing1hz = list(set(list_1hz)-set(names_1hz))

    if output == 'xarray':
        import xarray as xr
        # add missing columns filled with NaNs
        nan50 = np.empty(30000)
        nan50[:] = np.nan
        nan1 = np.empty(600)
        nan1[:] = np.nan
        try:
            if missing50hz:
                for x in missing50hz: dat_50hz.append(nan50)
                names_50hz.append(missing50hz)
            if missing1hz:
                for y in missing1hz: dat_1hz.append(nan1)
                names_1hz.append(missing1hz)
        except:
            pass
        # convert lists into xarray
        df50hz = xr.DataArray(np.transpose(dat_50hz), dims=('timestamp','channel'),coords={'timestamp':tstamp50hz,'channel':names_50hz})
        df1hz = xr.DataArray(np.transpose(dat_1hz), dims=('timestamp','channel'),coords={'timestamp':tstamp1hz,'channel':names_1hz})
        # store in dict
        ds = {'data50hz': df50hz, 'data1hz': df1hz}
        return ds

    if output == 'pandas':
        # convert lists into pandas
        df_50hz = pd.DataFrame(dat_50hz).T
        df_50hz.columns = names_50hz
        df_1hz = pd.DataFrame(dat_1hz).T
        df_1hz.columns = names_1hz
        df_50hz.index = tstamp50hz
        df_1hz.index = tstamp1hz
        try:
            # add columns with missing channels filled with NaNs
            for x in missing50hz: df_50hz[x] = np.nan
            for y in missing1hz: df_1hz[y] = np.nan
        except:
            pass
        # store in dict
        ds = {'data50hz': df_50hz, 'data1hz': df_1hz}
        return ds

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    EmptyFileError    = type('EmptyFileError', (Exception,),{})
    WrongFormatError  = type('WrongFormatError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})
    File=dict

class RAAWMatFile(File):
    """ 
    Read a RAAW .mat file. The object behaves as a dictionary.
    
    Main methods
    ------------
    - read, toDataFrame, keys
    
    Examples
    --------
        f = RAAWMatFile('file.mat')
        print(f.keys())
        print(f.toDataFrame().columns)  
    
    """

    @staticmethod
    def defaultExtensions():
        """ List of file extensions expected for this fileformat"""
        return ['.mat']

    @staticmethod
    def formatName():
        """ Short string (~100 char) identifying the file format"""
        return 'RAAW .mat file'

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

    def write(self, filename=None):
        """ Rewrite object to file, or write object to `filename` if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        # Calling (children) function to write
        self._write()

    def _read(self):
        """ Reads self.filename and stores data into self. Self is (or behaves like) a dictionary"""
        self['data'] = matfile(self.filename,output='pandas')

    def _write(self):
        """ Writes to self.filename"""
        # --- Example:
        #with open(self.filename,'w') as f:
        #    f.write(self.toString)
        raise NotImplementedError()

    def toDataFrame(self):
        """ Returns object into one DataFrame, or a dictionary of DataFrames"""
        return self['data']

    # --- Optional functions
    def __repr__(self):
        """ String that is written to screen when the user calls `print()` on the object. 
        Provide short and relevant information to save time for the user. 
        """
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|Main attributes:\n'
        s+='| - filename: {}\n'.format(self.filename)
        # --- Example printing some relevant information for user
        #s+='|Main keys:\n'
        #s+='| - ID: {}\n'.format(self['ID'])
        #s+='| - data : shape {}\n'.format(self['data'].shape)
        s+='|Main methods:\n'
        s+='| - read, toDataFrame, keys\n'
        s+='|Info:\n'
        d1=self['data']['data1hz']
        d5=self['data']['data50hz']
        s+='| - data1hz : {} to {} (n:{}, T:{}, dt:{})\n'.format(d1.index[0], d1.index[-1], len(d1), (d1.index[-1]-d1.index[0]).total_seconds(), (d1.index[1]-d1.index[0]).total_seconds())
        s+='| - data50hz: {} to {} (n:{}, T:{}, dt:{})\n'.format(d5.index[0], d5.index[-1], len(d5), (d5.index[-1]-d5.index[0]).total_seconds(), (d5.index[1]-d5.index[0]).total_seconds())
        return s
