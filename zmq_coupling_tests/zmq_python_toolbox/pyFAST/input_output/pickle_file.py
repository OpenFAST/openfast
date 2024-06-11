""" 
Input/output class for the pickle fileformats
"""
import numpy as np
import pandas as pd
import os
import pickle
import builtins

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File=dict
    EmptyFileError    = type('EmptyFileError', (Exception,),{})
    WrongFormatError  = type('WrongFormatError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})

class PickleFile(File):
    """ 
    Read/write a pickle file. The object behaves as a dictionary.
    
    Main methods
    ------------
    - read, write, toDataFrame, keys
    
    Examples
    --------
        f = PickleFile('file.pkl')
        print(f.keys())
        print(f.toDataFrame().columns)  
    
    """

    @staticmethod
    def defaultExtensions():
        """ List of file extensions expected for this fileformat"""
        return ['.pkl']

    @staticmethod
    def formatName():
        """ Short string (~100 char) identifying the file format"""
        return 'Pickle file'

    @staticmethod
    def priority(): return 60 # Priority in weio.read fileformat list between 0=high and 100:low


    def __init__(self, filename=None, data=None, **kwargs):
        """ Class constructor. If a `filename` is given, the file is read. """
        self.filename = filename
        if filename and not data:
            self.read(**kwargs)
        if data:
            self._setData(data)
            if filename:
                self.write()

    def _setData(self, data):
        if isinstance(data, dict):
            for k,v in data.items():
                self[k] = v
        else:
            if hasattr(data, '__dict__'):
                self.update(data.__dict__)
            else:
                self['data'] = data

    def addDict(self, data):
        self._setData(data)

    def additem(self, key, data):
        self[key]=data

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
        # Reads self.filename and stores data into self. Self is (or behaves like) a dictionary
        # If pickle data is a dict we store its keys in self, otherwise with store the pickle in the "data" key
        d = pickle.load(open(self.filename, 'rb'))
        self._setData(d)

    def write(self, filename=None):
        """ Rewrite object to file, or write object to `filename` if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        with open(self.filename, 'wb') as fid:
            pickle.dump(dict(self), fid)

    def toDataFrame(self):
        """ Returns object into one DataFrame, or a dictionary of DataFrames"""
        dfs={}
        for k,v in self.items():
            if isinstance(v, pd.DataFrame):
                dfs[k] = v
            elif isinstance(v, np.ndarray):
                if len(v.shape)==2:
                    dfs[k] = pd.DataFrame(data=v, columns=['C{}'.format(i) for i in range(v.shape[1])])
                elif len(v.shape)==1:
                    dfs[k] = pd.DataFrame(data=v, columns=[k])
        if len(dfs)==1:
            dfs=dfs[list(dfs.keys())[0]]
        return dfs

    # --- Optional functions
    def __repr__(self):
        """ String that is written to screen when the user calls `print()` on the object. 
        Provide short and relevant information to save time for the user. 
        """
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|Main attributes:\n'
        s+='| - filename: {}\n'.format(self.filename)
        # --- Example printing some relevant information for user
        s+='|Main keys:\n'
        for k,v in self.items():
            try:
                s+='| - {}: type:{} shape:{}\n'.format(k,type(v),v.shape)
            except:
                try:
                    s+='| - {}: type:{} len:{}\n'.format(k,type(v), len(v))
                except:
                    s+='| - {}: type:{}\n'.format(k,type(v))
        s+='|Main methods:\n'
        s+='| - read, write, toDataFrame, keys'
        return s
    

    # --- Functions speficic to filetype
    def toGlobal(self, namespace=None, overwrite=True, verbose=False, force=False) :
    #def toGlobal(self, **kwargs):
        """ 
        NOTE: very dangerous, mostly works for global, but then might infect everything

        Inject variables (keys of read dict) into namespace (e.g. globals()). 
        By default, the namespace of the caller is used
        To use the global namespace, use namespace=globals()
        """
        import inspect
        st = inspect.stack()
        if len(st)>2:
            if not force:
                raise Exception('toGlobal is very dangerous, only use in isolated script. use `force=True` if you really know what you are doing')
            else:
                print('[WARN] toGlobal is very dangerous, only use in isolated script')
        if namespace is None:
            # Using parent local namespace
            namespace = inspect.currentframe().f_back.f_globals
            #namespace = inspect.currentframe().f_back.f_locals # could use f_globals
            # Using global (difficult, overwriting won't work) It's best if the user sets namespace=globals()
            # import builtins as _builtins
            # namespace = _builtins # NOTE: globals() is for the package globals only, we need "builtins"

        gl_keys = list(namespace.keys())
        for k,v in self.items():
            if k in gl_keys: 
                if not overwrite:
                    print('[INFO] not overwritting variable {}, already present in global namespace'.format(k))
                    continue
                else:
                    print('[WARN] overwritting variable {}, already present in global namespace'.format(k))

            if verbose:
                print('[INFO] inserting in namespace: {}'.format(k))
            namespace[k] = v # OR do: builtins.__setattr__(k,v)

