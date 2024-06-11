import pandas as pd

try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

class NetCDFFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.nc']

    @staticmethod
    def formatName():
        return 'NetCDF file (<=2D)'

    def _read(self):
        try:
            import xarray as xr
        except:
            raise Exception('Python module `xarray` not installed')

        self.data=xr.open_dataset(self.filename)

    def _write(self):
        self.data.to_netcdf(self.filename)

    def _toDataFrame(self):
        dfs={}
        for k in self.data.keys():
            # Not pretty...
            if len(self.data[k].shape)==2:
                dfs[k]=pd.DataFrame(data=self.data[k].values)
            elif len(self.data[k].shape)==1:
                dfs[k]=pd.DataFrame(data=self.data[k].values)
        return dfs

