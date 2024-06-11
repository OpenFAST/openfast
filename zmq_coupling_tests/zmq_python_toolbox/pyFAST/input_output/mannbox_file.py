""" 
Read/Write Mann Box file

Part of weio library: https://github.com/ebranlard/weio
"""
import pandas as pd
import numpy as np
import os
import re
import struct

try:
    from .file import File, EmptyFileError, BrokenFormatError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})
    File=dict

class MannBoxFile(File):
    """
    Read/Write Mann Box file

    Main keys
    ---------
    - 'field': velocity field, shape (nx x ny x nz)

    Main properties
    ---------
    - 'y', 'z': space coordinates 

    Main methods
    ------------

    - read, write, toDataFrame, keys, valuesAt, vertProfile, fromTurbSim


    Examples
    --------

        # Open a box
        mb = MannBoxFile('Turb.u', N = (1024, 16, 16))
            OR 
        mb = MannBoxFile('Turb_1024x16x16.u')

        # Show info
        print(mb)
        print(mb['field'].shape)  

        # Use methods to extract relevant values
        u = mb.valuesAt(y=10.5, z=90)
        z, means, stds = mb.vertProfile

        # Write to a new file
        mb.write('Output_1024x16x16.u')

    Notes
    --------

    Mann box: 
      - z is the fast index, then y, then x
      - y from from  ly/2 to -ly/2  !<<<<<<< IMPORTANT, we will flip it 
      - z from from -lz/2 to  lz/2 
      - ix==1..nx corresponds to it = nt..1

    The field stored in this class has the following properties:
      - shape: nx x ny x nz
      - y: goes from -ly/2 to ly/2  ! <<<<< IMPORTANT subtlety, it has been flipped
      - z: goes from -lz/2 to lz/2  
      - ix==1..nx corresponds to it = nt..1 (it has not been flipped)
    """


    @staticmethod
    def defaultExtensions():
        return ['.u','.v','.w','.bin']

    @staticmethod
    def formatName():
        return 'HAWC2 Turbulence box'

    def __init__(self,filename=None,**kwargs):
        """ Initialize the class, if a filename is provided, the box is read. 
        See the method `read` for  keywords arguments.
        """
        self.filename = None
        if filename:
            self.read(filename=filename,**kwargs)

    def read(self, filename=None, N=None, dy=1, dz=1, y0=None, z0=0, zMid=None):
        """ read MannBox
        INPUTS (all optional):
        - filename: name of input file to be read
        - N: tuple (nx x ny x nz) of box dimension. If None, dimensions are inferred from filename
        - y0: minimum value of the y vector (default is -ly/2 where ly = ny x dy)
        - z0: minimum value of the z vector (default is 0)
        - zMid: mid value of the z vector (default it lz/2 where lz= nz x dz )

        SET:
         - the keys 'field', array of shape (nx x ny x nz)
        NOTE: y-coord in Mann Box goes from Ly/2 -> -Ly/2 but we flip this to -Ly/2 -> Ly/2
        """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        if N is None:
            # try to infer N's from filename with format 'stringN1xN2xN3'
            basename = os.path.splitext(os.path.basename(self.filename))[0]
            splits = basename.split('x')
            temp = re.findall(r'\d+', basename)
            res = list(map(int, temp))
            if len(res)>=3:
                N=res[-3:]
            else:
                raise BrokenFormatError('Reading a Mann box requires the knowledge of the dimensions. The dimensions can be inferred from the filename, for instance: `filebase_1024x32x32.u`. Try renaming your file such that the three last digits are the dimensions in x, y and z.')
        nx,ny,nz=N

        def _read_buffered():
            data=np.zeros((nx,ny,nz),dtype=np.float32)
            with open(self.filename, mode='rb') as f:            
                for ix in range(nx):
                    Buffer = np.frombuffer(f.read(4*ny*nz), dtype=np.float32) # 4-bytes
                    data[ix,:,:] =  np.flip(Buffer.reshape((ny,nz)),0)
            return data

        def _read_nonbuffered():
            data = np.fromfile(self.filename, np.dtype('<f'), -1)
            assert len(data) == nx*ny*nz, "Size of turbulence box (%d) does not match nx x ny x nz (%d)" % (len(data),nx*ny*nz)
            # Fortran order z the fastest, then y then x. We then flip that back to nx, ny, nz
            data = np.transpose(data.reshape(nz, ny, nx, order='F'), axes=(2,1,0))
            # The issue is the y-coordinate in Mann Boxes go from Ly/2 -> -Ly/2
            # So we flip the y-axis, so that the field is consistent with typical y values
            return np.flip(data, 1) # i.e. data=data[:,::-1,:]

#         self['field']= _read_nonbuffered()
        self['field']= _read_buffered()
        self['dy']=dy
        self['dz']=dz
        self['y0']=y0
        self['z0']=z0
        self['zMid']=zMid

    def write(self, filename=None):
        """ Write mann box """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        nx,ny,nz = self['field'].shape
        sfmt='<{:d}f'.format(ny*nz)
        with open(self.filename, mode='wb') as f:            
            for ix in np.arange(nx):
                data = np.flip(self['field'][ix,:,:],0).ravel() # We have to flip the y axis again
                f.write(struct.pack(sfmt, *data))

    
    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        s+='| - filename: {}\n'.format(self.filename)
        s+='| - field:  shape {}x{}x{}\n'.format(self['field'].shape[0],self['field'].shape[1],self['field'].shape[2])
        s+='|   min: {}, max: {}, mean: {} \n'.format(np.min(self['field']), np.max(self['field']), np.mean(self['field']))
        s+='| - dy, dz:  {}, {}\n'.format(self['dy'], self['dz'])
        s+='| - y0, z0 zMid:  {}, {}, {}\n'.format(self['y0'], self['z0'], self['zMid'])
        z=self.z
        y=self.y
        s+='| * y: [{} ... {}],  dy: {}, n: {} \n'.format(y[0],y[-1],self['dy'],len(y))
        s+='| * z: [{} ... {}],  dz: {}, n: {} \n'.format(z[0],z[-1],self['dz'],len(z))
        s+='|useful functions:\n'
        s+='| - t(dx, U)\n'
        s+='| - valuesAt(y,z), vertProfile, fromTurbSim(*), _iMid()\n'
        return s


    @property
    def z(self):
        zmax = self['z0'] + (self['field'].shape[2]-1+0.1)*self['dz']
        z = np.arange(self['z0'], zmax, self['dz'])
        if self['zMid'] is not None:
            z+= self['zMid']-np.mean(z)
        return z

    @property
    def y(self):
        if self['y0'] is not None:
            ymax = self['y0'] + (self['field'].shape[1]-1+0.1)*self['dy']
            y = np.arange(self['y0'], ymax, self['dy'])
        else:
            ymax = (self['field'].shape[1]-1+0.1)*self['dy']
            y = np.arange(0, ymax, self['dy'])
            y -= np.mean(y)
        return y

    def t(self, dx, U):
        # 1.5939838     dx          - distance (in meters) between points in the x direction    (m)
        # 99.5          RefHt_Hawc  - reference height; the height (in meters) of the vertical center of the grid (m)
        # 6.26          URef        - Mean u-component wind speed at the reference height (m/s)
        dt = dx/U
        nt = self['field'].shape[0]
        return np.arange(0, dt*(nt-0.5), dt)

    # --------------------------------------------------------------------------------}
    # --- Extracting relevant data 
    # --------------------------------------------------------------------------------{
    def valuesAt(self, y, z, method='nearest'):
        """ return wind speed time series at a point """
        if method == 'nearest':
            iy, iz = self.closestPoint(y, z)
            u = self['field'][:,iy,iz]
        else:
            raise NotImplementedError()
        return u

    def closestPoint(self, y, z):
        iy = np.argmin(np.abs(self.y-y))
        iz = np.argmin(np.abs(self.z-z))
        return iy,iz

    def _iMid(self):
        _, ny, nz = self['field'].shape
        return int(ny/2), int(nz/2)

    @property
    def vertProfile(self):
        iy, iz = self._iMid()
        m = np.mean(self['field'][:,iy,:], axis=0)
        s = np.std (self['field'][:,iy,:], axis=0)
        return self.z,m,s


    def toDataFrame(self):
        dfs={}
        ny = len(self.y)
        nz = len(self.z)
        # Index at mid box
        iy,iz = self._iMid()

        # Mean vertical profile
        z, m, s = self.vertProfile
        ti = s/m*100
        cols=['z_[m]','vel_[m/s]','sigma_[m/s]','TI_[%]']
        data = np.column_stack((z,m[:],s[:],ti[:]))
        dfs['VertProfile'] = pd.DataFrame(data = data ,columns = cols)

        # Mid time series
        u = self['field'][:,iy,iz]
        cols=['t/T_[-]','vel_[m/s]']
        fake_t = np.linspace(0, 1, len(u))
        data = np.column_stack((fake_t,u[:]))
        dfs['ZMidLine'] = pd.DataFrame(data = data ,columns = cols)


        # ZMin YEnd time series
        u = self['field'][:,-1,iz]
        cols=['t/T_[-]','vel_[m/s]']
        fake_t = np.linspace(0, 1, len(u))
        data = np.column_stack((fake_t,u[:]))
        dfs['ZMidYEndLine'] = pd.DataFrame(data = data ,columns = cols)

        # ZMin YStart time series
        u = self['field'][:,0,iz]
        cols=['t/T_[-]','vel_[m/s]']
        fake_t = np.linspace(0, 1, len(u))
        data = np.column_stack((fake_t,u[:]))
        dfs['ZMidYStartLine'] = pd.DataFrame(data = data ,columns = cols)



#         # Mid crosscorr y
#         y, rho_uu_y, rho_vv_y, rho_ww_y = self.crosscorr_y()
#         cols = ['y_[m]', 'rho_uu_[-]','rho_vv_[-]','rho_ww_[-]']
#         data = np.column_stack((y, rho_uu_y, rho_vv_y, rho_ww_y))
#         dfs['Mid_xcorr_y'] = pd.DataFrame(data = data ,columns = cols)
# 
#         # Mid crosscorr z
#         z, rho_uu_z, rho_vv_z, rho_ww_z = self.crosscorr_z()
#         cols = ['z_[m]', 'rho_uu_[-]','rho_vv_[-]','rho_ww_[-]']
#         data = np.column_stack((z, rho_uu_z, rho_vv_z, rho_ww_z))
#         dfs['Mid_xcorr_z'] = pd.DataFrame(data = data ,columns = cols)
# 
#         # Mid csd
#         fc, chi_uu, chi_vv, chi_ww = self.csd_longi()
#         cols = ['f_[Hz]','chi_uu_[-]', 'chi_vv_[-]','chi_ww_[-]']
#         data = np.column_stack((fc, chi_uu, chi_vv, chi_ww))
#         dfs['Mid_csd_longi'] = pd.DataFrame(data = data ,columns = cols)
# 
#         # Mid csd
#         fc, chi_uu, chi_vv, chi_ww = self.csd_lat()
#         cols = ['f_[Hz]','chi_uu_[-]', 'chi_vv_[-]','chi_ww_[-]']
#         data = np.column_stack((fc, chi_uu, chi_vv, chi_ww))
#         dfs['Mid_csd_lat'] = pd.DataFrame(data = data ,columns = cols)
# 
#         # Mid csd
#         fc, chi_uu, chi_vv, chi_ww = self.csd_vert()
#         cols = ['f_[Hz]','chi_uu_[-]', 'chi_vv_[-]','chi_ww_[-]']
#         data = np.column_stack((fc, chi_uu, chi_vv, chi_ww))
#         dfs['Mid_csd_vert'] = pd.DataFrame(data = data ,columns = cols)
        return dfs



    # Useful converters
    def fromTurbSim(self, u, icomp=0, removeConstant=None, removeAllMean=False):
        """ 
        Assumes: 
             u (3 x nt x ny x nz)
        Removes the mean of the turbsim file for the "u" component.
        """
        if icomp==0:
            if removeAllMean is True:
                self['field'] = u[icomp, :, : ,: ]-np.mean(u[icomp,:,:,:],axis=0)
            elif removeConstant is not None:
                self['field'] = u[icomp, :, : ,: ]-removeConstant
            else:
                self['field'] = u[icomp, :, : ,: ]
        else:
            self['field'] = u[icomp, :, : ,: ]
        return self

if __name__=='__main__':
    mb = MannBoxFile('mini-u_1024x32x32.bin')
#     mb = MannBoxFile('mann_bin/mini-u.bin', N=(2,4,8))
#     F1=mb['field'].ravel()
#     mb.write('mann_bin/mini-u-out.bin')
# 
#     mb2= MannBoxFile('mann_bin/mini-u-out.bin', N=(2,4,8))
#     F2=mb2['field'].ravel()
#     print(F1-F2)
