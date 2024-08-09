""" 
Taken from python_toolbox https://github.com/openfast/python_toolbox
"""
import os
import numpy as np
import re
class BrokenFormatError(Exception): pass

class FASTLinearizationFile(dict):
    """ 
    Read/write an OpenFAST linearization file. The object behaves like a dictionary.

    Main keys
    ---------
    - 'x', 'xdot' 'u', 'y', 'A', 'B', 'C', 'D'

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
        # --- Calling (children) function to read
        self._read(**kwargs)

    def _read(self, *args, **kwargs):
        self['header']=[]

        def extractVal(lines, key):
            for l in lines:
                if l.find(key)>=0:
                    return l.split(key)[1].split()[0]
            return None

        def readToMarker(fid, marker, nMax):
            lines=[]
            for i, line in enumerate(fid):
                if i>nMax:
                    raise BrokenFormatError('`{}` not found in file'.format(marker))
                if line.find(marker)>=0:
                    break
                lines.append(line.strip())
            return lines, line
        
        def readOP(fid, n, name=''):
            OP=[]
            Var = {'RotatingFrame': [], 'DerivativeOrder': [], 'Description': []}
            colNames=fid.readline().strip()
            fid.readline().strip()
            bHasDeriv= colNames.find('Derivative Order')>=0
            for i, line in enumerate(fid):
                sp = line.strip().split()
                if sp[1].find(',') >= 0:
                    #  Most likely this OP has three values (e.g. orientation angles)
                    # For now we discard the two other values
                    OP.append(np.float32(sp[1][:-1]))
                    iRot=4
                else:
                    OP.append(np.float32(sp[1]))
                    iRot=2
                Var['RotatingFrame'].append(sp[iRot])
                if bHasDeriv:
                    Var['DerivativeOrder'].append(int(sp[iRot+1]))
                    Var['Description'].append(' '.join(sp[iRot+2:]).strip())
                else:
                    Var['DerivativeOrder'].append(-1)
                    Var['Description'].append(' '.join(sp[iRot+1:]).strip())
                if i>=n-1:
                    break
            OP=np.asarray(OP)
            return OP, Var

        def readMat(fid, n, m, name=''):

            # Read rows from file, raise exception on failure
            try:
                vals = np.genfromtxt(fid, dtype=np.float64, max_rows=n)
            except:
                raise Exception('Failed to convert into an array of float the matrix `{}`\n\tin linfile: {}'.format(name, self.filename))
            
            # Raise exception if actual matrix shape does not match expected shape
            if vals.shape[0]!=n or vals.shape[1]!=m:
                shape1 = vals.shape
                shape2 = (n,m)
                raise Exception('Shape of matrix `{}` has wrong dimension ({} instead of {})\n\tin linfile: {}'.format(name, shape1, shape2, name, self.filename))

            # Raise exceptions if any elements are NaN or infinity
            if np.any(np.isnan(vals.ravel())):
                raise Exception('Some ill-formated/infinite values (e.g. `*******`) were found in the matrix `{}`\n\tin linflile: {}'.format(name, self.filename))
            if np.any(np.isinf(vals.ravel())):
                raise Exception('Some NaN values were found in the matrix `{}`\n\tin linfile: `{}`.'.format(name, self.filename))
            return vals


        # Reading 
        with open(self.filename, 'r', errors="surrogateescape") as f:
            # --- Reader header
            self['header'], lastLine=readToMarker(f, 'Jacobians included', 30)
            self['header'].append(lastLine)
            nx  = int(extractVal(self['header'],'Number of continuous states:'))
            nxd = int(extractVal(self['header'],'Number of discrete states:'  ))
            nz  = int(extractVal(self['header'],'Number of constraint states:'))
            nu  = int(extractVal(self['header'],'Number of inputs:'           ))
            ny  = int(extractVal(self['header'],'Number of outputs:'          ))
            bJac = extractVal(self['header'],'Jacobians included in this file?')
            try:
                self['Azimuth'] = np.float32(extractVal(self['header'],'Azimuth:'))
            except:
                self['Azimuth'] = None
            try:
                self['RotSpeed'] = np.float32(extractVal(self['header'],'Rotor Speed:')) # rad/s
            except:
                self['RotSpeed'] = None
            try:
                self['WindSpeed'] = np.float32(extractVal(self['header'],'Wind Speed:'))
            except:
                self['WindSpeed'] = None

            for i, line in enumerate(f):
                line = line.strip()
                if line.find('Order of continuous states:')>=0:
                    self['x'], self['x_info'] = readOP(f, nx, 'x')
                elif line.find('Order of continuous state derivatives:')>=0:
                    self['xdot'], self['xdot_info'] = readOP(f, nx, 'xdot')
                elif line.find('Order of inputs')>=0:
                    self['u'], self['u_info'] = readOP(f, nu, 'u')
                elif line.find('Order of outputs')>=0:
                    self['y'], self['y_info'] = readOP(f, ny, 'y')
                elif line.find('A:')>=0:
                    self['A'] = readMat(f, nx, nx, 'A')
                elif line.find('B:')>=0:
                    self['B'] = readMat(f, nx, nu, 'B')
                elif line.find('C:')>=0:
                    self['C'] = readMat(f, ny, nx, 'C')
                elif line.find('D:')>=0:
                    self['D'] = readMat(f, ny, nu, 'D')
                elif line.find('dUdu:')>=0:
                    self['dUdu'] = readMat(f, nu, nu,'dUdu')
                elif line.find('dUdy:')>=0:
                    self['dUdy'] = readMat(f, nu, ny,'dUdy')
                elif line.find('StateRotation:')>=0:
                    pass
                    # TODO
                    #StateRotation:
                elif line.find('ED M:')>=0:
                    self['EDDOF'] = line[5:].split()
                    self['M']     = readMat(f, 24, 24,'M')

    def toString(self):
        s=''
        return s

    def _write(self):
        with open(self.filename,'w') as f:
            f.write(self.toString())

    def short_descr(self,slist):
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

    def xdescr(self):
        if 'x_info' in self.keys():
            return self.short_descr(self['x_info']['Description'])
        else:
            return []

    def xdotdescr(self):
        if 'xdot_info' in self.keys():
            return self.short_descr(self['xdot_info']['Description'])
        else:
            return []

    def ydescr(self):
        if 'y_info' in self.keys():
            return self.short_descr(self['y_info']['Description'])
        else:
            return []
    def udescr(self):
        if 'u_info' in self.keys():
            return self.short_descr(self['u_info']['Description'])
        else:
            return []

    def toDataFrame(self):
        import pandas as pd
        dfs={}

        xdescr_short    = self.xdescr()
        xdotdescr_short = self.xdotdescr()
        ydescr_short    = self.ydescr()
        udescr_short    = self.udescr()

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
        if 'xdot' in self.keys():
            dfs['xdot'] = pd.DataFrame(data = np.asarray(self['xdot']).reshape((1,-1)), columns=xdotdescr_short)
        if 'u' in self.keys():
            dfs['u'] = pd.DataFrame(data = np.asarray(self['u']).reshape((1,-1)), columns=udescr_short)
        if 'y' in self.keys():
            dfs['y'] = pd.DataFrame(data = np.asarray(self['y']).reshape((1,-1)), columns=ydescr_short)
        if 'M' in self.keys():
            dfs['M'] = pd.DataFrame(data = self['M'], index=self['EDDOF'], columns=self['EDDOF'])
        if 'dUdu' in self.keys():
            dfs['dUdu'] = pd.DataFrame(data = self['dUdu'], index=udescr_short, columns=udescr_short)
        if 'dUdy' in self.keys():
            dfs['dUdy'] = pd.DataFrame(data = self['dUdy'], index=udescr_short, columns=ydescr_short)

        return dfs


