""" 
Hawc2 AE file
"""
import os
import numpy as np
import pandas as pd

try:
    from .file import File, WrongFormatError, EmptyFileError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    WrongFormatError = type('WrongFormatError', (Exception,),{})

from .wetb.hawc2.ae_file import AEFile

class HAWC2AEFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.dat','.ae','.txt']

    @staticmethod
    def formatName():
        return 'HAWC2 AE file'

    def __init__(self,filename=None,**kwargs):
        if filename:
            self.filename = filename
            self.read(**kwargs)
        else:
            self.filename = None
            self.data = AEFile()

    def read(self, filename=None, **kwargs):
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)
        # ---
        try:
            self.data = AEFile(self.filename)
        except Exception as e:    
            raise WrongFormatError('AE File {}: '.format(self.filename)+e.args[0])

    def write(self, filename=None):
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        # ---
        self.data.save(self.filename)

    def toDataFrame(self):
        cols=['radius_[m]','chord_[m]','thickness_[%]','pc_set_[#]']
        nset = len(self.data.ae_sets)
        if nset == 1:
            return pd.DataFrame(data=self.data.ae_sets[1], columns=cols)
        else:
            dfs = {}
            for iset,aeset in enumerate(self.data.ae_sets):
                name='ae_set_{}'.format(iset+1)
                dfs[name] = pd.DataFrame(data=self.data.ae_sets[iset+1], columns=cols)
            return dfs

    @property
    def sets(self):
        # Returns a list of ae_sets, otherwise not easy to iterate
        sets=[]
        for iset,aeset in enumerate(self.data.ae_sets):
            sets.append(self.data.ae_sets[iset+1])
        return sets

    # --- Convenient utils
    def add_set(self, **kwargs):
        self.data.add_set(**kwargs)

    def __repr__(self):
        cols=['radius_[m]','chord_[m]','thickness_[%]','pc_set_[#]']
        nRows =  [np.asarray(s).shape[0] for s in self.sets]
        s='<{} object>\n'.format(type(self).__name__)
        s+='| Attributes:\n'
        s+='| - filename: {}\n'.format(self.filename)
        s+='| - data: AEFile, with attributes `ae_sets`\n'
        s+='| Derived attributes:\n'
        s+='| * sets: list of {} arrays, length: {}, 4 columns: {}\n'.format(len(self.data.ae_sets), nRows, cols)
        s+='| Methods: add_set, toDataFrame\n'
        return s
