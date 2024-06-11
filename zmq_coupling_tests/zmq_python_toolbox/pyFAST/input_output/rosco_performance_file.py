""" 
Input/output class for the ROSCO performance (Cp,Ct,Cq) fileformat 
"""
import numpy as np
import pandas as pd
import os

try:
    from .file import File, EmptyFileError, WrongFormatError, BrokenFormatError
except:
    EmptyFileError    = type('EmptyFileError', (Exception,),{})
    WrongFormatError  = type('WrongFormatError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})
    File=dict

class ROSCOPerformanceFile(File):
    """ 
    Read/write a ROSCO performance file. The object behaves as a dictionary.
    
    Main methods
    ------------
    - read, write, toDataFrame, keys
    
    Examples
    --------
        f = ROSCOPerformanceFile('Cp_Ct_Cq.txt')
        print(f.keys())
        print(f.toDataFrame().columns)  
        fig = f.plotCP3D()
        CPmax, tsr_max, pitch_max = f.CPmax()
        CP = fCP([0, 1], [5, 5])
        CT = fCT([0, 1], [5, 5])
    
    """

    @staticmethod
    def defaultExtensions():
        """ List of file extensions expected for this fileformat"""
        return ['.txt']

    @staticmethod
    def formatName():
        """ Short string (~100 char) identifying the file format"""
        return 'ROSCO Performance file'

    def __init__(self, filename=None, pitch=None, tsr=None, WS=None, CP=None, CT=None, CQ=None, name='',**kwargs):
        """ Class constructor. If a `filename` is given, the file is read. 
        Otherwise values may be provided directly

        INPUTS:
          - filename: input file for ROSCO Performance file
        OR
          - pitch: pitch angle [deg], array of length nPitch
          - tsr: tip-speed ratio [-], array of length nTSR
          - CP,CT,CQ: aerodynamic coefficients, arrays of shape nTSR x nPitch
                      CQ is optional since CP = tsr*CQ
          - name: wind turbine name
        """
        self.filename = filename
        self.name     = name # Turbine name
        self['pitch'] = pitch
        self['TSR']   = tsr
        self['WS']    = WS
        self['CP']    = CP
        self['CT']    = CT
        self['CQ']    = CQ

        if filename:
            self.read(**kwargs)

        if self['pitch'] is not None:
            self.checkConsistency()
        
    def read(self, filename=None, **kwargs):
        """ Reads the file self.filename, or `filename` if provided 
        stores data into self. 
        self is (or behaves like) a dictionary"""
        # --- Standard tests and exceptions (generic code)
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        pitch, TSR, WS, CP, CT, CQ = load_from_txt(self.filename)
        self['pitch'] = pitch
        self['TSR']   = TSR
        self['WS']    = WS
        self['CP']    = CP
        self['CT']    = CT
        self['CQ']    = CQ

    def write(self, filename=None):
        """ Rewrite object to file, or write object to `filename` if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        # Sanity
        self.checkConsistency()
        # Write
        write_rotor_performance(self.filename, self['pitch'], self['TSR'], self['CP'],self['CT'], self['CQ'], self['WS'], TurbineName=self.name)

    def checkConsistency(self):
        """ 
        Check that data makes sense.
         in particular, check if CP=lambda CQ
        """
        if self['WS'] is not None:
            if not hasattr(self['WS'],'__len__' ):
                self['WS'] = np.array([self['WS']]).astype(float)

        CQ  = self['CQ']
        CP  = self['CP']
        tsr = np.asarray(self['TSR'])
        TSR = np.tile(tsr.flatten(), (len(self['pitch']),1)).T
        if CQ is None and CP is not None:
            CQ = CP/TSR
            print('[INFO] Computing CQ from CP')
        elif CQ is not None and CP is None:
            CP = CQ*TSR
            print('[INFO] Computing CP from CQ')
        elif CQ is not None and CP is not None:
            pass
        else:
            raise Exception('CP and CQ cannot be None')
        # Check consistency
        CP2 = CQ*TSR
        deltaCP = np.abs(CP-CP2)/0.5*100 # relative difference in %, for a mean CP of 0.5
        if np.max(deltaCP)>5: # more than 5%
            raise Exception('Inconsitency between power coefficient and torque coefficient')
        self['CP'] = CP
        self['CQ'] = CQ

    def toDataFrame(self):
        """ Returns object into dictionary of DataFrames"""
        dfs={}
        columns = ['TSR_[-]']+['Pitch_{:.2f}_[deg]'.format(p) for p in self['pitch']]
        dfs['CP'] = pd.DataFrame(np.column_stack((self['TSR'], self['CP'])), columns=columns)
        dfs['CT'] = pd.DataFrame(np.column_stack((self['TSR'], self['CT'])), columns=columns)
        dfs['CQ'] = pd.DataFrame(np.column_stack((self['TSR'], self['CQ'])), columns=columns)
        return dfs

    # --- Optional functions
    def toAeroDisc(self, filename, R, csv=False, WS=None, omegaRPM=10):
        """ Convert to AeroDisc Format 
        INPUTS:
         - filename: filename to be written 
         - R: rotor radius [m]
         - csv: if True write to CSV format, else, use OpenFAST .dat format
        either:
         - WS: wind speed [m/s]
          or
         - omegaRPM: rotational speed [rpm]

        Logic to determine wind speed or rotational speed:
         - If user provide a wind speed, we use it. Omega is determined from TSR and WS
         - If user provide a rotational speed, we use it. WS is determined from TSR and omega
         - If ROSCO file contain one wind speed, we use it.
         - Otherwise, we don't know what to do so we raise an exception
        """
        # --- Logic to determine wind speed or rotational speed
        if WS is not None:
            WS = WS
        elif omegaRPM is not None:
            WS = None
            omega = omegaRPM*(2*np.pi)/60
        elif self['WS'] is not None and len(self['WS'])==1:
            WS = self['WS'][0]
        else:
            raise Exception('Provide either a wind speed (`WS`) or a rotational speed (`omegaRPM`)')

        with open(filename,'w') as fid:
            # Header
            if csv:
                fid.write("TSR_(-),  RtSpd_(rpm) ,   VRel_(m/s)  ,  Skew_(deg) ,   Pitch_(deg) ,  C_Fx_(-)  ,   C_Fy_(-)   ,  C_Fz_(-)   ,  C_Mx_(-)   ,  C_My_(-)  ,   C_Mz_(-)\n")
            else:
                fid.write('     TSR       RtSpd      VRel       Skew       Pitch      C_Fx       C_Fy       C_Fz       C_Mx       C_My       C_Mz\n')
                fid.write('     (-)       (rpm)      (m/s)      (deg)      (deg)       (-)        (-)        (-)        (-)        (-)        (-)\n')
            if csv:
                FMT='{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f},{:10.4f}\n'
            else:
                FMT='{:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}\n'
            # Loop on oper
            for j,tsr in enumerate(self['TSR']):
                if WS is None:
                    U0 = omega*R/tsr
                else:
                    U0 = WS
                    omega = tsr*U0/R
                    omegaRPM = omega*60/(2*np.pi)
                for i,p in enumerate(self['pitch']):
                    CP=self['CP'][j,i]
                    CT=self['CT'][j,i]
                    CQ=self['CQ'][j,i]
                    skew=0
                    cfx=CT
                    cfy=0
                    cfz=0
                    cmx=CQ
                    cmy=0
                    cmz=0
                    fid.write(FMT.format(tsr, omegaRPM, U0, skew, p, cfx,cfy,cfz,cmx,cmy,cmz))

    def computeWeights(self):
        """ Compute interpolant weights for fast evaluation of CP and CT at intermediate values"""
        CP = self['CP'].copy()
        CT = self['CT'].copy()
        CP = CP[CP<0]=0
        CT = CT[CT<0]=0
        self._fCP = interp2d_pairs(self['pitch'], self['TSR'], CP, kind='cubic')
        self._fCT = interp2d_pairs(self['pitch'], self['TSR'], CT, kind='cubic')

    def fCP(self, pitch, tsr):
        """ Compute CP for given pitch and tsr, where inputs can be scalar, arrays or matrices"""
        if self._fCP is None:
            self.computeWeights()
        return self.fCP(pitch, tsr)

    def fCT(self, pitch, tsr):
        """ Compute CT for given pitch and tsr, where inputs can be scalar, arrays or matrices"""
        if self._fCT is None:
            self.computeWeights()
        return self.fCT(pitch, tsr)

    def CPmax(self):
        """ return values at CPmax
        TODO: interpolation instead of nearest value..
        """
        CP = self['CP']
        i,j = np.unravel_index(CP.argmax(), CP.shape)
        CPmax, tsr_max, pitch_max = CP[i,j], self['TSR'][i], self['pitch'][j]

        return CPmax, tsr_max, pitch_max

    def plotCP3D(self, plotMax=True, trajectory=None):
        """
        Plot 3D surface of CP
        Optionally plot the maximum and a controller trajectory
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        # Data
        LAMBDA, PITCH = np.meshgrid(self['TSR'], self['pitch'])
        CP = self['CP'].copy()
        CP[CP<0]=0 # 
        CP_max, tsr_max, pitch_max = self.CPmax()
        # plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(LAMBDA, PITCH, np.transpose(CP), cmap=cm.coolwarm, linewidth=0, antialiased=True,alpha=0.8, label='$C_p$')
        if plotMax:
            ax.scatter(tsr_max, pitch_max, CP_max, c='k', marker='o', s=50, label=r'$C_{p,max}$')
        if trajectory is not None:
            if len(trajectory)==3:
                tsr_, pitch_, CP_ = trajectory
            else:
                tsr_, pitch_ = trajectory
                CP_ = self.fCP(tsr_, pitch_)
            ax.plot_surface(tsr_, pitch_, CP_, 'k-', linewidth=1 )
        #fig.tight_layout()
        #fig.colorbar(surf, shrink=0.5, aspect=15)
        ax.view_init(elev=20., azim=26)
        ax.set_xlabel('TSR [-]')
        ax.set_ylabel('Pitch [deg]')
        ax.set_zlabel(r'Power coefficient [-]')
        return fig

    def __repr__(self):
        """ String that is written to screen when the user calls `print()` on the object. 
        Provide short and relevant information to save time for the user. 
        """
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|Main attributes:\n'
        s+='| - filename: {}\n'.format(self.filename)
        # --- Example printing some relevant information for user
        s+='|Main keys:\n'
        s+='| - pitch: {} values: {}\n'.format(len(self['pitch']) if self['pitch'] is not None else 0, self['pitch'])
        s+='| - TSR:   {} values: {}\n'.format(len(self['TSR']  ) if self['TSR']   is not None else 0, self['TSR']  )
        s+='| - WS:    {} values: {}\n'.format(len(self['WS']   ) if self['WS']    is not None else 0, self['WS']   )
        if self['CP'] is not None:
            s+='| - CP,CT,CQ : shape {}\n'.format(self['CP'].shape)
        s+='|Main methods:\n'
        s+='| - read, write, toDataFrame, keys\n'
        s+='| - CPmax, plotCP3d, fCP, fCT, toAeroDisc'
        return s
    



def load_from_txt(txt_filename):
    '''
    Adapted from ROSCO_toolbox/utitities.py by Nikhar Abbas
        https://github.com/NREL/ROSCO
        Apache 2.0 License

    Load rotor performance data from a *.txt file. 
    Parameters:
    -----------
        txt_filename: str
                        Filename of the text containing the Cp, Ct, and Cq data. This should be in the format printed by the write_rotorperformance function
    '''

    pitch = None
    TSR   = None
    WS    = None

    with open(txt_filename) as pfile:
        for iline, line in enumerate(pfile):
            # Read Blade Pitch Angles (degrees)
            if 'Pitch angle' in line:
                pitch = np.array([float(x) for x in pfile.readline().strip().split()])

            # Read Tip Speed Ratios (rad)
            elif 'TSR' in line:
                TSR = np.array([float(x) for x in pfile.readline().strip().split()])

            #Read WS
            elif 'Wind speed' in line:
                WS = np.array([float(x) for x in pfile.readline().strip().split()])
            
            # Read Power Coefficients
            elif 'Power' in line:
                pfile.readline()
                Cp = np.empty((len(TSR),len(pitch)))
                for tsr_i in range(len(TSR)):
                    Cp[tsr_i] = np.array([float(x) for x in pfile.readline().strip().split()])
            
            # Read Thrust Coefficients
            elif 'Thrust' in line:
                pfile.readline()
                Ct = np.empty((len(TSR),len(pitch)))
                for tsr_i in range(len(TSR)):
                    Ct[tsr_i] = np.array([float(x) for x in pfile.readline().strip().split()])

            # Read Torque Coefficients
            elif 'Torque' in line:
                pfile.readline()
                Cq = np.empty((len(TSR),len(pitch)))
                for tsr_i in range(len(TSR)):
                    Cq[tsr_i] = np.array([float(x) for x in pfile.readline().strip().split()])

            if pitch is None and iline>10:
                raise WrongFormatError('This does not appear to be a ROSCO performance file, Pitch vector not found')

        return pitch, TSR, WS, Cp, Ct, Cq


def write_rotor_performance(txt_filename, pitch, TSR, CP, CT, CQ, WS=None, TurbineName=''):
    '''
    Adapted from ROSCO_toolbox/utitities.py by Nikhar Abbas
        https://github.com/NREL/ROSCO
        Apache 2.0 License

    Write text file containing rotor performance data
    Parameters:
    ------------
        txt_filename: str, optional
                        Desired output filename to print rotor performance data. Default is Cp_Ct_Cq.txt
    '''
    file = open(txt_filename,'w')
    # Headerlines
    file.write('# ----- Rotor performance tables for the wind turbine: {} ----- \n'.format(TurbineName))
    file.write('# ------------ Written using weio\n\n')

    # Pitch angles, TSR, and wind speed
    file.write('# Pitch angle vector, {} entries - x axis (matrix columns) (deg)\n'.format(len(pitch)))
    for i in range(len(pitch)):
        file.write('{:0.4}   '.format(pitch[i]))
    file.write('\n# TSR vector, {} entries - y axis (matrix rows) (-)\n'.format(len(TSR)))
    for i in range(len(TSR)):
        file.write('{:0.4}    '.format(TSR[i]))
    if WS is not None:
        file.write('\n# Wind speed vector - z axis (m/s)\n')
        for i in range(len(WS)):
            file.write('{:0.4f}    '.format(WS[i]))
        file.write('\n')
    
    # Cp
    file.write('\n# Power coefficient\n\n')
    for i in range(len(TSR)):
        for j in range(len(pitch)):
            file.write('{0:.6f}   '.format(CP[i,j]))
        file.write('\n')
    file.write('\n')
    
    # Ct
    file.write('\n#  Thrust coefficient\n\n')
    for i in range(len(TSR)):
        for j in range(len(pitch)):
            file.write('{0:.6f}   '.format(CT[i,j]))
        file.write('\n')
    file.write('\n')
    
    # Cq
    file.write('\n# Torque coefficient\n\n')
    for i in range(len(TSR)):
        for j in range(len(pitch)):
            file.write('{0:.6f}   '.format(CQ[i,j]))
        file.write('\n')
    file.write('\n')
    file.close()


def interp2d_pairs(*args,**kwargs):
    """ Same interface as interp2d but the returned interpolant will evaluate its inputs as pairs of values.
    Inputs can therefore be arrays

    example:
       f = interp2d_pairs(vx, vy, M, kind='cubic')

    vx: array of length nx
    vy: array of length ny
    M : array of shape nx x ny
    f : interpolant function
          v = f(x,y) : if x,y are array of length n, v is of length n
                       with  v_i = f(x_i, y_i)
    author: E. Branlard
    """
    import scipy.interpolate as si
    # Internal function, that evaluates pairs of values, output has the same shape as input
    def interpolant(x,y,f):
        x,y = np.asarray(x), np.asarray(y)
        return (si.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    # Wrapping the scipy interp2 function to call out interpolant instead
    return lambda x,y: interpolant(x,y,si.interp2d(*args,**kwargs))


if __name__ == '__main__':
    f = ROSCOPerformanceFile('./tests/example_files/RoscoPerformance_CpCtCq.txt')
    print(f)
    dfs = f.toDataFrame()
    print(dfs['CP'])

