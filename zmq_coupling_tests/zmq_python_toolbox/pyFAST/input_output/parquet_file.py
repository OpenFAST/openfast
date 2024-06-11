import pandas as pd

from .file import File


class ParquetFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.parquet']

    @staticmethod
    def formatName():
        return 'Parquet file'

    def __init__(self,filename=None,**kwargs):
        self.filename = filename
        if filename:
            self.read(**kwargs)


    def _read(self):
        """ use pandas read_parquet function to read parquet file"""
        self.data=pd.read_parquet(self.filename)

    def _write(self):
        """ use pandas DataFrame.to_parquet method to write parquet file """
        self.data.to_parquet(path=self.filename)

    def toDataFrame(self):
        #already stored as a data frame in self.data
        #just return self.data
        return self.data

    def fromDataFrame(self, df):
        #data already in dataframe
        self.data = df

    def toString(self):
        """ use pandas DataFrame.to_string method to convert to a string """
        s=self.data.to_string()
        return s

    def __repr__(self):
        s ='Class Parquet (attributes: data)\n'
        return s


