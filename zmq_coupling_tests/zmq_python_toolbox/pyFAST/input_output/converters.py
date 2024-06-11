import os


# --------------------------------------------------------------------------------
# --- Writing pandas DataFrame to different formats
# --------------------------------------------------------------------------------
def writeDataFrameToFormat(df, filename, fformat):
    """  
    Write a dataframe to disk based on user-specified fileformat
    - df: pandas dataframe
    - filename: filename 
    - fformat: fileformat in: ['csv', 'outb', 'parquet']
    """

    if fformat=='outb':
        dataFrameToOUTB(df, filename)
    elif fformat=='parquet':
        dataFrameToParquet(df, filename)
    elif fformat=='csv':
        dataFrameToCSV(df, filename, sep=',', index=False)
    else:
        raise Exception('File format not supported for dataframe export `{}`'.format(fformat))

def writeDataFrameAutoFormat(df, filename, fformat=None):
    """ 
    Write a dataframe to disk based on extension
    - df: pandas dataframe
    - filename: filename 
    """
    if fformat is not None:
        raise Exception()
    base, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext in ['.outb']:
        fformat = 'outb'
    elif ext in ['.parquet']:
        fformat = 'parquet'
    elif ext in ['.csv']:
        fformat = 'csv'
    else:
        print('[WARN] defaulting to csv, extension unknown: `{}`'.format(ext))
        fformat = 'csv'

    writeDataFrameToFormat(df, filename, fformat)

def writeFileDataFrames(fileObject, writer, extension='.conv', filename=None, **kwargs):
    """ 
    From a fileObejct, extract dataframes and write them to disk.

    - fileObject: object inheriting from weio.File with at least
                   - the attributes .filename
                   - the method     .toDataFrame()
    - writer: function with the interface:   writer ( dataframe, filename, **kwargs )
    """ 
    if filename is None:
        base, _ = os.path.splitext(fileObject.filename)
        filename = base + extension
    else:
        base, ext = os.path.splitext(filename)
        if len(ext)!=0:
            extension = ext
    if filename == fileObject.filename:
        raise Exception('Not overwritting {}. Specify a filename or an extension.'.format(filename))
        
    dfs = fileObject.toDataFrame()
    if isinstance(dfs, dict):
        for name,df in dfs.items():
            filename = base + name + extension
            if filename == fileObject.filename:
                raise Exception('Not overwritting {}. Specify a filename or an extension.'.format(filename))
            writeDataFrame(df=df, writer=writer, filename=filename, **kwargs)
    else:
        writeDataFrame(df=dfs, writer=writer, filename=filename, **kwargs)

def writeDataFrame(df, writer, filename, **kwargs):
    """ 
    Write a dataframe to disk based on a "writer" function. 
    - df: pandas dataframe
    - writer: function with the interface:   writer ( dataframe, filename, **kwargs )
    - filename: filename 
    """
    writer(df, filename, **kwargs)

# --- Low level writers
def dataFrameToCSV(df, filename, sep=',', index=False, **kwargs):
    df.to_csv(filename, sep=sep, index=index, **kwargs)

def dataFrameToOUTB(df, filename, **kwargs):
    from .fast_output_file import writeDataFrame as writeDataFrameToOUTB
    writeDataFrameToOUTB(df, filename, binary=True)

def dataFrameToParquet(df, filename, **kwargs):
    df.to_parquet(path=filename, **kwargs)

