import os
import numpy as np
import re
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class BrokenFormatError(Exception): pass

class SlowReaderNeededError(Exception):
    pass


_lin_vec = ['x','xd','xdot','u','y','z','header']
_lin_mat = ['A','B','C','D','dUdu','dUdy', 'StateRotation', 'M']
_lin_dict = ['x_info','xdot_info','u_info','y_info']

class FASTLinearizationFile(File):
    """ 
    Read/write an OpenFAST linearization file. The object behaves like a dictionary.

    Main keys
    ---------
    - 'x', 'xdot', 'xd', 'u', 'y', 'z', 'A', 'B', 'C', 'D'

    Main methods
    ------------
    - read, write, toDataFrame, keys, xdescr, ydescr, udescr

    Examples
    --------

        f = FASTLinearizationFile('5MW.1.lin')
        print(f.keys())
        print(f['u'])     # input operating point
        print(f.udescr()) # description of inputs

        # use a dataframe with "named" columns and rows
        df = f.toDataFrame()
        print(df['A'].columns)
        print(df['A'])

    """
    @staticmethod
    def defaultExtensions():
        return ['.lin']

    @staticmethod
    def formatName():
        return 'FAST linearization output'

    def __init__(self, filename=None, **kwargs):
        """ Class constructor. If a `filename` is given, the file is read. """
        self.filename = filename
        if filename:
            self.read(**kwargs)

    def read(self, filename=None, starSub=None, removeStatesPattern=None):
        """ Reads the file self.filename, or `filename` if provided

        - starSub: if None, raise an error if `****` are present
                           otherwise replace *** with `starSub` (e.g. 0)
        - removeStatesPattern: if None, do nothing
                           otherwise search for states matching a pattern and remove them
                           e.g:  'tower|Drivetrain'  or '^AD'
                           see removeStates in this file.
        """
        
        # --- Standard tests and exceptions (generic code)
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        # --- Main Data
        self['header']=[]

        # --- StarValues replacement `*****` -> inf
        starPattern = re.compile(r"[\*]+")
        starSubStr = ' inf '
        starSubFn  = lambda si: starPattern.sub(starSubStr, si)

        # Reading function, with slow or fast reader. See sub functions at end of this file
        def doRead(slowReader=False):
            with open(self.filename, 'r', errors="surrogateescape") as f:
                # --- Reader header
                self['header'], lastLine=readToMarker(f, 'Jacobians included', 30)
                self['header'].append(lastLine)
                nx   = extractVal(self['header'],'Number of continuous states:'    , dtype=int, NA=np.nan, missing=None)
                nxd  = extractVal(self['header'],'Number of discrete states:'      , dtype=int, NA=np.nan, missing=None)
                nz   = extractVal(self['header'],'Number of constraint states:'    , dtype=int, NA=np.nan, missing=None)
                nu   = extractVal(self['header'],'Number of inputs:'               , dtype=int, NA=np.nan, missing=None)
                ny   = extractVal(self['header'],'Number of outputs:'              , dtype=int, NA=np.nan, missing=None)
                bJac = extractVal(self['header'],'Jacobians included in this file?', dtype=bool, NA=False, missing=None)
                self['Azimuth']   = extractVal(self['header'], 'Azimuth:'    , dtype=float, NA=np.nan, missing=None)
                self['RotSpeed']  = extractVal(self['header'], 'Rotor Speed:', dtype=float, NA=np.nan, missing=None) # rad/s
                self['WindSpeed'] = extractVal(self['header'], 'Wind Speed:' , dtype=float, NA=np.nan, missing=None)
                self['t']  = extractVal(self['header'],'Simulation time:'    , dtype=float, NA=np.nan, missing=None)
                for i, line in enumerate(f):
                    line = line.strip()
                    if line.find('Order of continuous states:')>=0:
                        self['x'], self['x_info'] = readOP(f, nx, 'x', defaultDerivOrder=1, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('Order of continuous state derivatives:')>=0:
                        self['xdot'], self['xdot_info'] = readOP(f, nx, 'xdot', defaultDerivOrder=2, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('Order of discrete states:')>=0:
                        self['xd'], self['xd_info'] = readOP(f, nxd, 'xd', defaultDerivOrder=2, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('Order of inputs')>=0:
                        self['u'], self['u_info'] = readOP(f, nu, 'u', defaultDerivOrder=0, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('Order of outputs')>=0:
                        self['y'], self['y_info'] = readOP(f, ny, 'y', defaultDerivOrder=0, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('Order of constraint states:')>=0:
                        self['z'], self['z_info'] = readOP(f, nz, 'z', defaultDerivOrder=0, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('A:')>=0:
                        self['A'] = readMat(f, nx, nx, 'A', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('B:')>=0:
                        self['B'] = readMat(f, nx, nu, 'B', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('C:')>=0:
                        self['C'] = readMat(f, ny, nx, 'C', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('D:')>=0:
                        self['D'] = readMat(f, ny, nu, 'D', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('dUdu:')>=0:
                        self['dUdu'] = readMat(f, nu, nu,'dUdu', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('dUdy:')>=0:
                        self['dUdy'] = readMat(f, nu, ny,'dUdy', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
                    elif line.find('StateRotation:')>=0:
                        pass
                        # TODO
                        #StateRotation:
                    elif line.find('ED M:')>=0:
                        self['EDDOF'] = line[5:].split()
                        self['M']     = readMat(f, 24, 24,'M', slowReader=slowReader, filename=self.filename, starSubFn=starSubFn, starSub=starSub)
        try:
            doRead(slowReader=False)
        except SlowReaderNeededError:
            doRead(slowReader=True)

        if removeStatesPattern is not None:
            self.removeStates(pattern=removeStatesPattern)

    def toString(self):
        s=''
        return s

    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString())

    @property
    def nx(self):
        if 'x' in self.keys():
            return len(self['x'])
        return 0

    @property
    def nxd(self):
        if 'xd' in self.keys():
            return len(self['xd'])
        return 0

    @property
    def nu(self):
        if 'u' in self.keys():
            return len(self['u'])
        return 0

    @property
    def ny(self):
        if 'y' in self.keys():
            return len(self['y'])
        return 0

    @property
    def nz(self):
        if 'z' in self.keys():
            return len(self['z'])
        return 0

    @property
    def u_descr(self):
        if self.nu>0:
            return self['u_info']['Description']
        else:
            return []

    @property
    def x_descr(self):
        if self.nx>0:
            return self['x_info']['Description']
        else:
            return []

    @property
    def xd_descr(self): # Discrete states not derivative!
        if self.nxd>0:
            return self['xd_info']['Description']
        else:
            return []

    @property
    def xdot_descr(self):
        if self.nx>0:
            return self['xdot_info']['Description']
        else:
            return []

    @property
    def y_descr(self):
        if self.ny>0:
            return self['y_info']['Description']
        else:
            return []

    @property
    def z_descr(self):
        if self.nz>0:
            return self['z_info']['Description']
        else:
            return []

    def __repr__(self):
        s='<{} object> with attributes:\n'.format(type(self).__name__)
        s+=' - filename: {}\n'.format(self.filename)
        s+=' * nx      : {}\n'.format(self.nx)
        s+=' * nxd     : {}\n'.format(self.nxd)
        s+=' * nu      : {}\n'.format(self.nu)
        s+=' * ny      : {}\n'.format(self.ny)
        s+=' * nz      : {}\n'.format(self.nz)
        s+='keys:\n'
        for k,v in self.items():
            if k in _lin_vec:
                s+=' - {:15s}: shape: ({}) \n'.format(k,len(v))
            elif k in _lin_mat:
                s+=' - {:15s}: shape: ({} x {})\n'.format(k,v.shape[0], v.shape[1])
            elif k in _lin_dict:
                s+=' - {:15s}: dict with keys: {} \n'.format(k,list(v.keys()))
            else:
                s+=' - {:15s}: {}\n'.format(k,v)
        s+='methods:\n'
        s+=' - toDataFrame: convert A,B,C,D to dataframes\n'
        s+=' - removeStates: remove states\n'
        s+=' - eva: eigenvalue analysis\n'

        return s

    def toDataFrame(self):
        import pandas as pd
        dfs={}

        xdescr_short    = short_descr(self.x_descr)
        xddescr_short   = short_descr(self.xd_descr)
        xdotdescr_short = short_descr(self.xdot_descr)
        udescr_short    = short_descr(self.u_descr)
        ydescr_short    = short_descr(self.y_descr)
        zdescr_short    = short_descr(self.z_descr)

        if 'A' in self.keys():
            dfs['A'] = pd.DataFrame(data = self['A'], index=xdescr_short, columns=xdescr_short)
        if 'B' in self.keys():
            dfs['B'] = pd.DataFrame(data = self['B'], index=xdescr_short, columns=udescr_short)
        if 'C' in self.keys():
            dfs['C'] = pd.DataFrame(data = self['C'], index=ydescr_short, columns=xdescr_short)
        if 'D' in self.keys():
            dfs['D'] = pd.DataFrame(data = self['D'], index=ydescr_short, columns=udescr_short)
        if 'x' in self.keys():
            dfs['x'] = pd.DataFrame(data = np.asarray(self['x']).reshape((1,-1)), columns=xdescr_short)
        if 'xd' in self.keys():
            dfs['xd'] = pd.DataFrame(data = np.asarray(self['xd']).reshape((1,-1)))
        if 'xdot' in self.keys():
            dfs['xdot'] = pd.DataFrame(data = np.asarray(self['xdot']).reshape((1,-1)), columns=xdotdescr_short)
        if 'u' in self.keys():
            dfs['u'] = pd.DataFrame(data = np.asarray(self['u']).reshape((1,-1)), columns=udescr_short)
        if 'y' in self.keys():
            dfs['y'] = pd.DataFrame(data = np.asarray(self['y']).reshape((1,-1)), columns=ydescr_short)
        if 'z' in self.keys():
            dfs['z'] = pd.DataFrame(data = np.asarray(self['z']).reshape((1,-1)), columns=zdescr_short)
        if 'M' in self.keys():
            dfs['M'] = pd.DataFrame(data = self['M'], index=self['EDDOF'], columns=self['EDDOF'])
        if 'dUdu' in self.keys():
            dfs['dUdu'] = pd.DataFrame(data = self['dUdu'], index=udescr_short, columns=udescr_short)
        if 'dUdy' in self.keys():
            dfs['dUdy'] = pd.DataFrame(data = self['dUdy'], index=udescr_short, columns=ydescr_short)

        return dfs

    def removeStates(self, pattern=None, Irm=None, verbose=True):
        """ 
        remove states based on pattern or index

        - pattern: e.g:  'tower|Drivetrain'  or '^AD'
        """
        if self.nx==0:
            return
        desc = self['x_info']['Description']
        Iall = set(range(len(desc)))
        sInfo=''
        if pattern is not None:
            Irm = [i for i, s in enumerate(desc) if re.search(pattern, s)]
            sInfo=' with pattern `{}`'.format(pattern)
        if verbose:
            print('[INFO] removing {}/{} states{}'.format(len(Irm), len(Iall), sInfo))
        Ikeep = list(Iall.difference(Irm))
        Ikeep.sort() # safety
        if len(Ikeep)==0:
            raise Exception('All states have been removed{}!'.format(sInfo))
        # Remove states and info in vectors
        self['x']      = self['x'][Ikeep]
        self['xdot']   = self['xdot'][Ikeep]
        for k in self['x_info'].keys():
            self['x_info'][k] = self['x_info'][k][Ikeep]
            self['xdot_info'][k] = self['xdot_info'][k][Ikeep]
        # Remove states in matrices
        if 'A' in self.keys():
            self['A'] = self['A'][np.ix_(Ikeep,Ikeep)]
        if 'B' in self.keys():
            self['B'] = self['B'][Ikeep,:]
        if 'C' in self.keys():
            self['C'] = self['C'][:, Ikeep]


    def eva(self, normQ=None, sort=True, discardIm=True):
        """ Perform eigenvalue analysis of A matrix and return frequencies and damping """
        # --- Excerpt from welib.tools.eva.eigA
        A = self['A'] 
        n,m = A.shape
        if m!=n:
            raise Exception('Matrix needs to be square')
        # Basic EVA
        D,Q    = np.linalg.eig(A)
        Lambda = np.diag(D)
        v      = np.diag(Lambda)

        # Selecting eigenvalues with positive imaginary part (frequency)
        if discardIm:
            Ipos = np.imag(v)>0
            Q = Q[:,Ipos]
            v = v[Ipos]

        # Frequencies and damping based on compled eigenvalues
        omega_0 = np.abs(v)              # natural cylic frequency [rad/s]
        freq_d  = np.imag(v)/(2*np.pi)   # damped frequency [Hz]
        zeta    = - np.real(v)/omega_0   # damping ratio
        freq_0  = omega_0/(2*np.pi)      # natural frequency [Hz]
        # Sorting
        if sort:
            I = np.argsort(freq_0)
            freq_d = freq_d[I]
            freq_0 = freq_0[I]
            zeta   = zeta[I]
            Q      = Q[:,I]

        # Normalize Q
        if normQ=='byMax':
            for j in range(Q.shape[1]):
                q_j = Q[:,j]
                scale = np.max(np.abs(q_j))
                Q[:,j]= Q[:,j]/scale
        return freq_d, zeta, Q, freq_0 


def short_descr(slist):
    """ Shorten and "unify" the description from lin file """
    def shortname(s):
        s=s.strip()
        s = s.replace('(m/s)'   , '_[m/s]'  );
        s = s.replace('(kW)'    , '_[kW]'   );
        s = s.replace('(deg)'   , '_[deg]'  );
        s = s.replace('(N)'     , '_[N]'    );
        s = s.replace('(kN-m)'  , '_[kNm]' );
        s = s.replace('(N-m)'  , '_[Nm]' );
        s = s.replace('(kN)'  , '_[kN]' );
        s = s.replace('(rpm)'   , '_[rpm]'  );
        s = s.replace('(rad)'   , '_[rad]'  );
        s = s.replace('(rad/s)' , '_[rad/s]'  );
        s = s.replace('(rad/s^2)', '_[rad/s^2]'  );
        s = s.replace('(m/s^2)' , '_[m/s^2]');
        s = s.replace('(deg/s^2)','_[deg/s^2]');
        s = s.replace('(m)'     , '_[m]'    );
        s = s.replace(', m/s/s','_[m/s^2]');
        s = s.replace(', m/s^2','_[m/s^2]');
        s = s.replace(', m/s','_[m/s]');
        s = s.replace(', m','_[m]');
        s = s.replace(', rad/s/s','_[rad/s^2]');
        s = s.replace(', rad/s^2','_[rad/s^2]');
        s = s.replace(', rad/s','_[rad/s]');
        s = s.replace(', rad','_[rad]');
        s = s.replace(', -','_[-]');
        s = s.replace(', Nm/m','_[Nm/m]');
        s = s.replace(', Nm','_[Nm]');
        s = s.replace(', N/m','_[N/m]');
        s = s.replace(', N','_[N]');
        s = s.replace('(1)','1')
        s = s.replace('(2)','2')
        s = s.replace('(3)','3')
        s= re.sub(r'\([^)]*\)','', s) # remove parenthesis
        s = s.replace('ED ','');
        s = s.replace('BD_','BD_B');
        s = s.replace('IfW ','');
        s = s.replace('Extended input: ','')
        s = s.replace('1st tower ','qt1');
        s = s.replace('2nd tower ','qt2');
        nd = s.count('First time derivative of ')
        if nd>=0:
            s = s.replace('First time derivative of '     ,'');
            if nd==1:
                s = 'd_'+s.strip()
            elif nd==2:
                s = 'dd_'+s.strip()
        s = s.replace('Variable speed generator DOF ','psi_rot'); # NOTE: internally in FAST this is the azimuth of the rotor
        s = s.replace('fore-aft bending mode DOF '    ,'FA'     );
        s = s.replace('side-to-side bending mode DOF','SS'     );
        s = s.replace('bending-mode DOF of blade '    ,''     );
        s = s.replace(' rotational-flexibility DOF, rad','-ROT'   );
        s = s.replace('rotational displacement in ','rot'   );
        s = s.replace('Drivetrain','DT'   );
        s = s.replace('translational displacement in ','trans'   );
        s = s.replace('finite element node ','N'   );
        s = s.replace('-component position of node ','posN')
        s = s.replace('-component inflow on tower node','TwrN')
        s = s.replace('-component inflow on blade 1, node','Bld1N')
        s = s.replace('-component inflow on blade 2, node','Bld2N')
        s = s.replace('-component inflow on blade 3, node','Bld3N')
        s = s.replace('-component inflow velocity at node','N')
        s = s.replace('X translation displacement, node','TxN')
        s = s.replace('Y translation displacement, node','TyN')
        s = s.replace('Z translation displacement, node','TzN')
        s = s.replace('X translation velocity, node','TVxN')
        s = s.replace('Y translation velocity, node','TVyN')
        s = s.replace('Z translation velocity, node','TVzN')
        s = s.replace('X translation acceleration, node','TAxN')
        s = s.replace('Y translation acceleration, node','TAyN')
        s = s.replace('Z translation acceleration, node','TAzN')
        s = s.replace('X orientation angle, node'  ,'RxN')
        s = s.replace('Y orientation angle, node'  ,'RyN')
        s = s.replace('Z orientation angle, node'  ,'RzN')
        s = s.replace('X rotation velocity, node'  ,'RVxN')
        s = s.replace('Y rotation velocity, node'  ,'RVyN')
        s = s.replace('Z rotation velocity, node'  ,'RVzN')
        s = s.replace('X rotation acceleration, node'  ,'RAxN')
        s = s.replace('Y rotation acceleration, node'  ,'RAyN')
        s = s.replace('Z rotation acceleration, node'  ,'RAzN')
        s = s.replace('X force, node','FxN')
        s = s.replace('Y force, node','FyN')
        s = s.replace('Z force, node','FzN')
        s = s.replace('X moment, node','MxN')
        s = s.replace('Y moment, node','MyN')
        s = s.replace('Z moment, node','MzN')
        s = s.replace('FX', 'Fx')
        s = s.replace('FY', 'Fy')
        s = s.replace('FZ', 'Fz')
        s = s.replace('MX', 'Mx')
        s = s.replace('MY', 'My')
        s = s.replace('MZ', 'Mz')
        s = s.replace('FKX', 'FKx')
        s = s.replace('FKY', 'FKy')
        s = s.replace('FKZ', 'FKz')
        s = s.replace('MKX', 'MKx')
        s = s.replace('MKY', 'MKy')
        s = s.replace('MKZ', 'MKz')
        s = s.replace('Nodes motion','')
        s = s.replace('cosine','cos'   );
        s = s.replace('sine','sin'   );
        s = s.replace('collective','coll.');
        s = s.replace('Blade','Bld');
        s = s.replace('rotZ','TORS-R');
        s = s.replace('transX','FLAP-D');
        s = s.replace('transY','EDGE-D');
        s = s.replace('rotX','EDGE-R');
        s = s.replace('rotY','FLAP-R');
        s = s.replace('flapwise','FLAP');
        s = s.replace('edgewise','EDGE');
        s = s.replace('horizontal surge translation DOF','Surge');
        s = s.replace('horizontal sway translation DOF','Sway');
        s = s.replace('vertical heave translation DOF','Heave');
        s = s.replace('roll tilt rotation DOF','Roll');
        s = s.replace('pitch tilt rotation DOF','Pitch');
        s = s.replace('yaw rotation DOF','Yaw');
        s = s.replace('vertical power-law shear exponent','alpha')
        s = s.replace('horizontal wind speed ','WS')
        s = s.replace('propagation direction','WD')
        s = s.replace(' pitch command','pitch')
        s = s.replace('HSS_','HSS')
        s = s.replace('Bld','B')
        s = s.replace('tower','Twr')
        s = s.replace('Tower','Twr')
        s = s.replace('Nacelle','Nac')
        s = s.replace('Platform','Ptfm')
        s = s.replace('SrvD','SvD')
        s = s.replace('Generator torque','Qgen')
        s = s.replace('coll. blade-pitch command','PitchColl')
        s = s.replace('wave elevation at platform ref point','WaveElevRefPoint')
        s = s.replace('1)','1');
        s = s.replace('2)','2');
        s = s.replace('3)','3');
        s = s.replace(',','');
        s = s.replace(' ','');
        s=s.strip()
        return s
    return [shortname(s) for s in slist]



def extractVal(lines, key, NA=None, missing=None, dtype=float):
    for l in lines:
        if l.find(key)>=0:
            #l = starPattern.sub(starSubStr, l)
            try:
                return dtype(l.split(key)[1].split()[0])
            except:
                return NA
    return missing

def readToMarker(fid, marker, nMax):
    lines=[]
    for i, line in enumerate(fid):
        if i>nMax:
            raise BrokenFormatError('`{}` not found in file'.format(marker))
        if line.find(marker)>=0:
            break
        lines.append(line.strip())
    return lines, line

def readOP(fid, n, name='', defaultDerivOrder=1, filename='', starSubFn=None, starSub=None):
    OP=[]
    Var = {'RotatingFrame': [], 'DerivativeOrder': [], 'Description': []}
    colNames=fid.readline().strip()
    dummy=   fid.readline().strip()
    bHasDeriv= colNames.find('Derivative Order')>=0
    for i, line in enumerate(fid):
        line = line.strip()
        line = starSubFn(line)
        sp   = line.split()
        if sp[1].find(',')>=0:
            #  Most likely this OP has three values (e.g. orientation angles)
            # For now we discard the two other values
            OP.append(float(sp[1][:-1]))
            iRot=4
        else:
            OP.append(float(sp[1]))
            iRot=2
        Var['RotatingFrame'].append(sp[iRot])
        if bHasDeriv:
            Var['DerivativeOrder'].append(int(sp[iRot+1]))
            Var['Description'].append(' '.join(sp[iRot+2:]).strip())
        else:
            Var['DerivativeOrder'].append(defaultDerivOrder)
            Var['Description'].append(' '.join(sp[iRot+1:]).strip())
        if i>=n-1:
            break
    OP = np.asarray(OP)
    nInf = np.sum(np.isinf(OP))
    if nInf>0:
        sErr = 'Some ill-formated/infinite values (e.g. `*******`) were found in the vector `{}`\n\tin linflile: {}'.format(name, filename)
        if starSub is None:
            raise Exception(sErr)
        else:
            print('[WARN] '+sErr)
            OP[np.isinf(OP)] = starSub
            
    Var['RotatingFrame']   = np.asarray(Var['RotatingFrame'])
    Var['DerivativeOrder'] = np.asarray(Var['DerivativeOrder'])
    Var['Description']     = np.asarray(Var['Description'])
    return OP, Var



def readMat(fid, n, m, name='', slowReader=False, filename='', starSubFn=None, starSub=None):
    if not slowReader:
        try:
            return np.array([fid.readline().strip().split() for i in np.arange(n)],dtype=float)
        except:
            print('[INFO] Failed to read some value in matrix {}, trying slower reader'.format(name))
            raise SlowReaderNeededError()
    else:
        #vals = vals.ravel()
        #vals = np.array(list(map(starSubFn, vals))).reshape(n,m)
        vals=np.array([starSubFn( fid.readline().strip() ).split() for i in np.arange(n)], dtype=str)
        try:
            vals = vals.astype(float) # This could potentially fail
        except:
            raise Exception('Failed to convert into an array of float the matrix `{}`\n\tin linfile: {}'.format(name, filename))
        if vals.shape[0]!=n or vals.shape[1]!=m:
            shape1 = vals.shape
            shape2 = (n,m)
            raise Exception('Shape of matrix `{}` has wrong dimension ({} instead of {})\n\tin linfile: {}'.format(name, shape1, shape2, name, filename))

        nNaN = np.sum(np.isnan(vals.ravel()))
        nInf = np.sum(np.isinf(vals.ravel()))
        if nInf>0:
            sErr = 'Some ill-formated/infinite values (e.g. `*******`) were found in the matrix `{}`\n\tin linflile: {}'.format(name, filename)
            if starSub is None:
                raise Exception(sErr)
            else:
                print('[WARN] '+sErr)
                vals[np.isinf(vals)] = starSub
        if nNaN>0:
            raise Exception('Some NaN values were found in the matrix `{}`\n\tin linfile: `{}`.'.format(name, filename))
        return vals

if __name__ == '__main__':
    f = FASTLinearizationFile('../../data/example_files/StandstillSemi_ForID_EDHD.1.lin')
    print(f)
    _, zeta1, _, freq1 = f.eva()
    f.removeStates(pattern=r'^AD')
    print(f)
    dfs = f.toDataFrame()
    _, zeta2, _, freq2 = f.eva()
    print('f',freq1)
    print('f',freq2)
    print('d',zeta1)
    print('d',zeta2)

