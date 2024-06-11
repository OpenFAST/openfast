"""Read AMR-Wind NETCDF file

"""
import xarray as xr
import numpy as np

class AMRWindFile(dict):
    """ 
    Read a AMR-Wind output file (.nc)
    """

    @staticmethod
    def defaultExtensions():
        """ List of file extensions expected for this fileformat"""
        return ['.nc']

    @staticmethod
    def formatName():
        """ Short string (~100 char) identifying the file format"""
        return 'NetCDF plane sampling file from AMRWind'

    @staticmethod
    def priority(): return 60 # Priority in weio.read fileformat list between 0=high and 100:low

    def __init__(self, filename=None, timestep=None, output_frequency=None, **kwargs):
        self.filename   = filename
        self.amrwind_dt = timestep
        self.output_dt  = timestep * output_frequency 
            
        if filename:
            self.read(**kwargs)

    def read(self, group_name):        
        """
        Parameters
        ----------
        
        group_name : str,
            group name inside netcdf file that you want to read, e.g. p_slice

        TODO: see if group_name can be avoided, and add a read_group function
        """   
        # --- Standard tests and exceptions (generic code)
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise Exception('File is empty:',self.filename)
        
        
        ds = xr.open_dataset(self.filename,group=group_name)    
    
        coordinates = {"x":(0,"axial"), "y":(1,"lateral"),"z":(2,"vertical")}
        c           = {}
        for coordinate,(i,desc) in coordinates.items():
            c[coordinate] = xr.IndexVariable( 
                                             dims=[coordinate],
                                             data=np.sort(np.unique(ds['coordinates'].isel(ndim=i))), 
                                             attrs={"description":"{0} coordinate".format(desc),"units":"m"}
                                            )
        c["t"]                    = xr.IndexVariable( 
                                                     dims=["t"],
                                                     data=ds.num_time_steps*self.output_dt,
                                                     attrs={"description":"time from start of simulation","units":"s"}
                                                    )    

        self.nt = len(c["t"])
        self.nx = len(c["x"])
        self.ny = len(c["y"])
        self.nz = len(c["z"])

        coordinates = {"x":(0,"axial","u"), "y":(1,"lateral","v"),"z":(2,"vertical","w")}    
        v           = {}    
        for coordinate,(i,desc,u) in coordinates.items():        
            v[u] = xr.DataArray(np.reshape(getattr(ds,"velocity{0}".format(coordinate)).values,(self.nt,self.nx,self.ny,self.nz)), 
                                 coords=c, 
                                 dims=["t","x","y","z"],
                                 name="{0} velocity".format(desc), 
                                 attrs={"description":"velocity along {0}".format(coordinate),"units":"m/s"})

        ds = xr.Dataset(data_vars=v, coords=v[u].coords)           
        ds.attrs = {"original file":self.filename}
        
        self.data = ds


    def write(self, filename=None):
        """ Rewrite object to file, or write object to `filename` if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        raise NotImplementedError()

    def toDataFrame(self):
        """ Returns object into one DataFrame, or a dictionary of DataFrames"""
        # --- Example (returning one DataFrame):
        #  return pd.DataFrame(data=np.zeros((10,2)),columns=['Col1','Col2'])
        # --- Example (returning dict of DataFrames):
        #dfs={}
        #cols=['Alpha_[deg]','Cl_[-]','Cd_[-]','Cm_[-]']
        #dfs['Polar1'] = pd.DataFrame(data=..., columns=cols)
        #dfs['Polar1'] = pd.DataFrame(data=..., columns=cols)
        # return dfs
        raise NotImplementedError()

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
        s+='| - read, write, toDataFrame, keys'
        return s

