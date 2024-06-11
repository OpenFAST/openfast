import numpy as np
import pandas as pd
import os

try:
    from .file import File, WrongFormatError, BrokenFormatError, OptionalImportError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass
    class OptionalImportError(Exception): pass

class TDMSFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.tdms']

    @staticmethod
    def formatName():
        return 'TDMS file'

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
        try:
            from nptdms import TdmsFile
        except:
            raise OptionalImportError('Install the library nptdms to read this file')

        fh = TdmsFile(self.filename, read_metadata_only=False)
        # --- OLD, using some kind of old version of tdms and probably specific to one file
        #   channels_address = list(fh.objects.keys())
        #   channels_address = [ s.replace("'",'') for s in channels_address]
        #   channel_keys= [ s.split('/')[1:]  for s in channels_address if len(s.split('/'))==3]
        #   # --- Setting up list of signals and times
        #   signals=[]
        #   times=[]
        #   for i,ck in enumerate(channel_keys):
        #       channel = fh.object(ck[0],ck[1])
        #       signals.append(channel.data)
        #       times.append  (channel.time_track())

        #   lenTimes = [len(time) for time in times]
        #   minTimes = [np.min(time) for time in times]
        #   maxTimes = [np.max(time) for time in times]
        #   if len(np.unique(lenTimes))>1:
        #       print(lenTimes)
        #       raise NotImplementedError('Different time length') 
        #       # NOTE: could use fh.as_dataframe
        #   if len(np.unique(minTimes))>1:
        #       print(minTimes)
        #       raise NotImplementedError('Different time span') 
        #   if len(np.unique(maxTimes))>1:
        #       print(maxTimes)
        #       raise NotImplementedError('Different time span') 
        #   # --- Gathering into a data frame with time
        #   time =times[0]
        #   signals = [time]+signals
        #   M = np.column_stack(signals)
        #   colnames = ['Time_[s]'] + [ck[1] for ck in channel_keys]
        #   self['data'] =  pd.DataFrame(data = M, columns=colnames)
        # --- NEW
        self['data'] = fh

        #for group in fh.groups():
        #   for channel in group.channels():
        #       #channel = group['channel name']
        #       print('Group:',group.name , 'Chan:',channel.name)
        #       channel_data = channel[:]
        #       if len(channel_data)>0:
        #           print(' ', type(channel_data))
        #           #print(' ', len(channel_data))
        #           print(' ', channel_data)
        #           print(' ', channel_data[0])
        #           try:
        #               print(channel.time_track())
        #           except KeyError:
        #               print('>>> No time track')

    def write(self, filename=None, df=None):
        """" 
        Write to TDMS file.
        NOTE: for now only using a conversion from dataframe...
        """
        if filename is None:
            filename = self.filename
        if df is None:
            df = self.toDataFrame(split=False)
        writeTDMSFromDataFrame(filename, df)


    def groups(self):
        return self['data'].groups()

    @property
    def groupNames(self):
        return [group.name for group in self['data'].groups()]

    def __repr__(self):
        s ='Class TDMS (key: data)\n'
        s +=' - data: TdmsFile\n'
        s +=' * groupNames: {}\n'.format(self.groupNames)
        #for group in fh.groups():
        #   for channel in group.channels():
        #       print(group.name)
        #       print(channel.name)
        return s

    def toDataFrame(self, split=True):
        """ Export to one (split=False) or several dataframes (split=True)
        Splitting on the group
        """

        def cleanColumns(df):
            # Cleanup columns
            colnames = df.columns
            colnames=[c.replace('\'','') for c in colnames]
            colnames=[c[1:] if c.startswith('/') else c for c in colnames]
            # If there is only one group, we remove the group key
            groupNames = self.groupNames
            if len(groupNames)==1:
                nChar = len(groupNames[0])
                colnames=[c[nChar+1:] for c in colnames] # +1 for the "/"
            df.columns = colnames

        fh = self['data']
        if split:
            # --- One dataframe per group. We skip group that have empty data
            dfs={}
            for group in fh.groups():
                try:
                    df = group.as_dataframe(time_index=True)
                    df.insert(0,'Time_[s]', df.index.values)
                    df.index=np.arange(0,len(df))
                except KeyError:
                    df = group.as_dataframe(time_index=False)
                if len(df)>0:
                    dfs[group.name] = df
            if len(dfs)==1:
                dfs=dfs[group.name]
            return dfs
        else:
            # --- One dataframe with all data
            try:
                df = fh.as_dataframe(time_index=True)
                cleanColumns(df)
                df.insert(0,'Time_[s]', df.index.values)
                df.index=np.arange(0,len(df))
            except KeyError:
                df = fh.as_dataframe(time_index=False)
            return df

def writeTDMSFromDataFrame(filename, df, defaultGroupName='default'):
    """ 
    Write a TDMS file from a pandas dataframe

    Example:
        # --- Create a TDMS file - One group two channels with time track
        time = np.linspace(0,1,20)
        colA = np.random.normal(0,1,20)
        colB = np.random.normal(0,1,20)
        df = pd.DataFrame(data={'Time_[s]':time ,'ColA':colA,'ColB':colB})
        writeTDMSFromDataFrame('out12.tdms', df, defaultGroupName = 'myGroup')

        #--- Create a TDMS file - Two groups, two channels without time track but with timestamp
        TS = np.arange('2010-02', '2010-02-21', dtype='datetime64[D]')
        df = pd.DataFrame(data={'GroupA/ColTime':time,'GroupA/ColA':colA,'GroupB/ColTimestamp': TS,'GroupB/ColA':colB)})
        writeTDMSFromDataFrame('out22.tdms', df)

    """
    from nptdms import TdmsWriter, ChannelObject
    
    defaultGroupName = 'default'

    columns =df.columns

    # Check if first column is time
    if columns[0].lower().find('time')==0:
        t = df.iloc[:,0].values
        n = len(t)
        dt1 = (np.max(t)-np.min(t))/(n-1)
        if n>1:
            dt2 = t[1]-t[0]
        timeProperties = {'wf_increment':dt1, 'wf_start_offset':t[0]}
        columns = columns[1:] # We remove the time column
    else:
        timeProperties = {}

    with TdmsWriter(filename) as tdms_writer:

        channels=[]
        for iCol, col in enumerate(columns):
            sp = col.split('/')
            if len(sp)==2:
                groupName   = sp[0]
                channelName = sp[1]
            else:
                groupName   = defaultGroupName
                channelName = col
            data_array = df[col].values
            channels.append(ChannelObject(groupName, channelName, data_array, timeProperties))
        tdms_writer.write_segment(channels)

if __name__ == '__main__':
    pass
#     f = TDMSFile('TDMS_.tdms')
#     dfs = f.toDataFrame(split=True)
#     print(f)

