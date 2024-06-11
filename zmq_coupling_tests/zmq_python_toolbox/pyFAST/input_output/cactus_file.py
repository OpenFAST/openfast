import numpy as np
import pandas as pd
import os

try:
    from .file import File, WrongFormatError, BrokenFormatError, EmptyFileError
except:
    File=dict
    EmptyFileError    = type('EmptyFileError', (Exception,),{})
    WrongFormatError  = type('WrongFormatError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})

class CactusFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.in']

    @staticmethod
    def formatName():
        return 'CACTUS file'

    def __init__(self,filename=None,**kwargs):
        self.filename = filename
        if filename:
            self.read(**kwargs)

    def read(self, filename=None, **kwargs):
        """ read self, or read filename if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)
        # Calling children function
        self._read(**kwargs)

    def write(self, filename=None):
        """ write self, or to filename if provided """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        # Calling children function
        self._write()

    def _read(self):
        """ """
        import f90nml
        from .csv_file import CSVFile

        filepath  = self.filename
        basepath  = os.path.splitext(filepath)[0]
        basename  = os.path.basename(basepath)
        parentdir = os.path.dirname(basepath)

        # --- Read main input file
        nml     = f90nml.read(filepath)
        for k in ['configinputs','caseinputs']:
            self[k] = nml[k]

        # --- Try to read geometry file
        arfoilfile = os.path.join(parentdir, nml['caseinputs']['afdpath'])
        geomfile = os.path.join(parentdir, nml['caseinputs']['geomfilepath'])
        if os.path.exists(geomfile):
            with open(geomfile, 'r', errors="surrogateescape") as fid:
                geom=dict()
                nMax=10
                for i, line in enumerate(fid):
                    # remove comment
                    line = line.strip().split('!')[0]
                    sp   = line.strip().split(':')
                    key  = sp[0].strip().lower()
                    if len(key)>0:
                        strvalue = sp[1]
                        try:
                            value    = np.asarray(strvalue.split()).astype(float)
                            if len(value)==1:
                                value = value[0]
                        except:
                            value    = strvalue
                        geom[key]=value
                    if i==nMax:
                        break
            self['geom']=geom
            self['geom']['file']=geomfile
        else:
            print('[FAIL] Geom file not found (quantites will be pooorly scaled):',geomfile)
            self['geom']={'nblade':1, 'refr':3.28084, 'refar':2, 'file':None}

        # --- Try to read element time data file
        timefile = os.path.join(parentdir, 'output', basename+'_TimeData.csv')
        if os.path.exists(timefile):
            df = CSVFile(timefile).toDataFrame()
            nBlades = list(df.columns).count('Blade Fx Coeff. (-)')
            self['geom']['nblade'] = list(df.columns).count('Blade Fx Coeff. (-)')
            cols=list(df.columns[:8])
            bldCols=['Blade{:d} Fx Coeff. (-)', 'Blade{:d} Fy Coeff. (-)', 'Blade{:d} Fz Coeff. (-)', 'Blade{:d} Torque Coeff. (-)']
            for ib in range(self['geom']['nblade']):
                cols+=[b.format(ib+1) for b in bldCols]
            df.columns=cols
            self['dfTime']=df

        else:
            self['dfTime']=None
            print('TimeData file not found:',timefile)

        # --- Try to read element data file
        elemfile = os.path.join(parentdir, 'output', basename+'_ElementData.csv')
        if os.path.exists(elemfile):
            dfElem = CSVFile(elemfile).toDataFrame()
            self['dfElem'] = dfElem
        else:
            self['dfElem'] = None
            print('ElementData file not found:',elemfile)


        # --- Read DS file
        dsfile = os.path.join(parentdir, 'output', basename+'_DSData.csv')
        try:
            dfDS   =CSVFile(dsfile).toDataFrame()
            self['dfDS'] = dfDS
        except (FileNotFoundError, EmptyFileError):
            self['dfDS'] = None
            print('DSData file not found or empty:',dsfile)


    @property
    def omega(self):
        return self['caseinputs']['rpm']*2*np.pi/60

    @property
    def TSR(self):
        return self['caseinputs']['ut']

    @property
    def RPM(self):
        return self['caseinputs']['rpm']

    @property
    def dt(self):
        nRot      = self['configinputs']['nr']
        nPerRot   = self['configinputs']['nti']
        T         = 2*np.pi/(self.omega)
        return   T/nPerRot

    @property
    def R(self):
        if self['geom']['file'] is not None:
            R = self['geom']['refr']/3.28084 # feet to m
        else:
            R=1
        return R

    @property
    def A(self):
        # NOTE: Turbine reference area (for force/torque/power normalization) divided by reference radius squared.
        if self['geom']['refar'] is not None:
            #A = self['geom']['refar']/(3.28084**2) # feet^2 to m^2
            A = self['geom']['refar']*self['geom']['refr']**2
            A /=(3.28084**2) # feet^2 to m^2
        else:
            A    = (2*self.R)**2   # D^2
        return A

    @property
    def U(self):
        return self.omega*self.R/self.TSR


    @property
    def time(self):
        nRot      = self['configinputs']['nr']
        nPerRot   = self['configinputs']['nti']
        timeSteps = np.arange(0,nRot*nPerRot)
        T         = 2*np.pi/(self.omega)
        return timeSteps*self.dt


    def timeDataToOpenFAST(self, df):
        """ Convert to similar labels as OpenFAST"""
        if df is None:
            return None
        nRot      = self['configinputs']['nr']
        nPerRot   = self['configinputs']['nti']
        TSR       = self.TSR
        CTExcrM   = self['caseinputs']['ctexcrm']
        rho       = self['caseinputs']['rho']*1.2/0.0023280000

        time =  self.time
        if df.shape[0]<len(time):
            print('[WARN] Wrong number of time steps ({} instead of {}), assuming the simulation was stopped early'.format(df.shape[0],len(time)))
            time = time[:df.shape[0]]

        theta     = time*self.omega

        #time_norm = df['Normalized Time (-)'].values
        #dt_norm   = time_norm[1]-time_norm[0]

        R    = self.R
        U    = self.omega*self.R/TSR
        Utip = self.omega*self.R
        #D    = 2*R
        #A    = D**2                  # <<<<<<<<<<<<<< NOTE:
        A    = self.A
        #print('R=',R,'AR',self['geom']['refar'],'D^2=',D**2,'A=',self.A)
        P0  = 0.5 *rho * A * U**3 
        F0  = 0.5 *rho * A * U**2
        Q0  = 0.5 *rho * A * U**2 * R
        ones = np.ones(time.shape)

        print('TSR= ',TSR,'omega = ',self.omega,'RPM = ',self.RPM)
        print('U  = ',U,'Utip = ',Utip,'dt = ',self.dt,'R=',R)

        # --- General 
        c=0
        df.insert(c , 'Time_[s]'        , time)                              ; c+=1
        df.insert(c , 'Azimuth_[deg]'   , np.mod(theta*180/np.pi , 360))          ; c+=1
        df.insert(c , 'HWindSpeedX_[m/s]' , U*ones)                               ; c+=1
        df.insert(c , 'RtAeroCp_[-]'   , df['Power Coeff. (-)'])                  ; c+=1
        df.insert(c , 'RtAeroCq_[-]'   , df['Torque Coeff. (-)'])                 ; c+=1
        df.insert(c , 'ZAeroFxg_[N]'   , df['Fx Coeff. (-)']*F0)                  ; c+=1
        df.insert(c , 'ZAeroFyg_[N]'   ,-df['Fz Coeff. (-)']*F0)                  ; c+=1
        df.insert(c , 'ZAeroFzg_[N]'   , df['Fy Coeff. (-)']*F0)                  ; c+=1
        df.insert(c , 'ZAeroMzg_[N-m]' , df['Torque Coeff. (-)']*Q0)              ; c+=1

        # TODO this need cos/sin azimuth
        df.insert(c , 'RtAeroFxh_[N]'   , df['Fy Coeff. (-)']*F0)                  ; c+=1

        df.insert(c , 'RtAeroFyh_[N]'   , (np.cos(theta)*df['Fx Coeff. (-)'] - np.sin(theta)*df['Fz Coeff. (-)'])*F0  ); c+=1
        df.insert(c , 'RtAeroFzh_[N]'   ,-(np.cos(theta)*df['Fz Coeff. (-)'] + np.sin(theta)*df['Fx Coeff. (-)'])*F0  ); c+=1
        df.insert(c , 'RtAeroMxh_[N-m]' , df['Torque Coeff. (-)']*Q0)              ; c+=1
        df.insert(c , 'RtAeroPwr_[W]'   , df['Power Coeff. (-)']*P0)               ; c+=1
        df['RtAeroPwr_[W]'].values[0]=df['RtAeroPwr_[W]'].values[1]

        df.insert(c , 'RotSpeed_[rpm]'   , self.RPM*ones)                               ; c+=1

        # --- 
        nBld = self['geom']['nblade']
        if nBld<4:
            azimuth = np.linspace(0,2*np.pi, nBld+1)[:-1]
        else:
            #raise Exception('Hacky azimuth and blade number')
            azimuth = np.linspace(0,2*np.pi, nBld+1)[:-1]

        for ib, azi in enumerate(azimuth):
            psi = theta - azi
            sb =str(ib+1)
            df.insert(c , 'Y'+sb+'AeroFxb_[N]'   , -(np.cos(psi)*df['Blade'+sb+' Fz Coeff. (-)'] + np.sin(psi)*df['Blade'+sb+' Fx Coeff. (-)'])*F0)           ; c+=1
            df.insert(c , 'Y'+sb+'AeroFyb_[N]'   ,  (np.cos(psi)*df['Blade'+sb+' Fx Coeff. (-)'] - np.sin(psi)*df['Blade'+sb+' Fz Coeff. (-)'])*F0)           ; c+=1
            df.insert(c , 'Y'+sb+'AeroFzb_[N]'   ,-df['Blade'+sb+' Fy Coeff. (-)']*F0)           ; c+=1
            df.insert(c , 'Y'+sb+'AeroMxb_[N-m]' , df['Blade'+sb+' Torque Coeff. (-)']*Q0)       ; c+=1
            df.insert(c , 'Y'+sb+'AeroPwr_[W]'   , df['Blade'+sb+' Torque Coeff. (-)']*Q0*self.omega) ; c+=1

            df.insert(c , 'Z'+sb+'AeroFxg_[N]'   , df['Blade'+sb+' Fx Coeff. (-)']*F0)           ; c+=1
            df.insert(c , 'Z'+sb+'AeroFyg_[N]'   ,-df['Blade'+sb+' Fz Coeff. (-)']*F0)           ; c+=1
            df.insert(c , 'Z'+sb+'AeroFzg_[N]'   , df['Blade'+sb+' Fy Coeff. (-)']*F0)           ; c+=1
            df.insert(c , 'Z'+sb+'AeroMzg_[N-m]' , df['Blade'+sb+' Torque Coeff. (-)']*Q0)       ; c+=1

            df['Y'+sb+'AeroPwr_[W]'].values[0] = df['Y'+sb+'AeroPwr_[W]'].values[1]

        df.insert(c , 'Wind1VelX_[m/s]' , U*ones)                                  ; c+=1

        return df, c


    def elemDataToOpenFAST(self, dfElem,  df, c=0, dfDS=None, alphaSign=-1):
        """ Convert to similar lables as openfast """
        if dfElem is None:
            return df

        U=self.U
        R=self.R

        nBld = self['geom']['nblade']
        if nBld<4:
            azimuth = np.linspace(0,2*np.pi, nBld+1)[:-1]
        else:
            #raise Exception('Hacky azimuth and blade number')
            azimuth = np.linspace(0,2*np.pi, nBld+1)[:-1]

        time=self.time
        if df.shape[0]<len(time):
            #print('[WARN] Wrong number of time steps ({} instead of {}), assuming the simulation was stopped early'.format(df.shape[0],len(time)))
            time = time[:df.shape[0]]
        theta  = time*self.omega


        IBld  = np.unique(dfElem['Blade']  ).astype(int)
        for iB in IBld:
            dfBld=dfElem[dfElem['Blade']==iB]
            if dfDS is not None:
                dfBld_DS = dfDS  [dfDS  ['Blade'] == iB]
            IElem = np.unique(dfBld['Element']).astype(int)
            
            psi = theta - azimuth[iB-1]

            for ie in IElem:
                dfSec = dfBld[dfBld['Element']==ie]
                if dfSec.shape[0]!=df.shape[0]:
                    print('>>> Inconsistent shape, ',iB,ie)
                else:
                    # TODO x/y/R
                    uix_g=dfSec['IndU (-)'].values*U
                    uiy_g=dfSec['IndV (-)'].values*U
                    uiz_g=dfSec['IndW (-)'].values*U

                    uix=-(np.cos(psi)*uiz_g + np.sin(psi)*uix_g)
                    uiy= (np.cos(psi)*uix_g - np.sin(psi)*uiz_g)

                    df.insert(c , 'AB{:d}N{:03d}Alpha_[deg]'.format(iB,ie)   , alphaSign*dfSec['AOA25 (deg)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Alpha50_[deg]'.format(iB,ie) , alphaSign*dfSec['AOA50 (deg)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Alpha75_[deg]'.format(iB,ie) , alphaSign*dfSec['AOA75 (deg)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Cl_[-]'     .format(iB,ie)   , alphaSign*dfSec['CL (-)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Cd_[-]'     .format(iB,ie)   ,           dfSec['CD (-)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Cm_[-]'     .format(iB,ie)   , alphaSign*dfSec['CM25 (-)'].values); c+=1

                    #BladeElemOutData(BladeElemOutRow,24)=CN                       ! Element normal force coefficient (per span) based on local chord and flow velocity
                    #BladeElemOutData(BladeElemOutRow,25)=CT                       ! Element tangential force coefficient (per span) based on local chord and flow velocity
                    df.insert(c , 'AB{:d}N{:03d}Cn_[-]'     .format(iB,ie)   , alphaSign*dfSec['CN (-)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Ct_[-]'     .format(iB,ie)   , alphaSign*dfSec['CT (-)'].values); c+=1 
                    #   CT=-CL5*sin(alpha25)+CD5*cos(alpha25)
                    #   CT=-CL5*sin(alpha50)+CD5*cos(alpha50)
                    Cl= dfSec['CL (-)'].values; Cd= dfSec['CD (-)'].values; alpha= dfSec['AOA25 (deg)'].values*np.pi/180
                    df.insert(c , 'AB{:d}N{:03d}Ct2_[-]'    .format(iB,ie)   , alphaSign*(-Cl*np.sin(alpha) + Cd*np.cos(alpha))); c+=1 

                    df.insert(c , 'AB{:d}N{:03d}Cxg_[-]'     .format(iB,ie)  ,           dfSec['Fx (-)'].values); c+=1 # TODO, this is likely coefficients related to global coords
                    df.insert(c , 'AB{:d}N{:03d}Cyg_[-]'     .format(iB,ie)  , -         dfSec['Fz (-)'].values); c+=1 # TODO
                    df.insert(c , 'AB{:d}N{:03d}Czg_[-]'     .format(iB,ie)  ,           dfSec['Fy (-)'].values); c+=1 # TODO
                    df.insert(c , 'AB{:d}N{:03d}ClC_[-]'    .format(iB,ie)   , alphaSign*dfSec['CLCirc (-)'].values); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Re_[-]'     .format(iB,ie)   ,           dfSec['Re (-)'].values/1e6); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Gam_[m^2/s]'.format(iB,ie)   , alphaSign*dfSec['GB (?)'].values*U*R); c+=1 # TODO
                    df.insert(c , 'AB{:d}N{:03d}Vrel_[m/s]' .format(iB,ie)   ,           dfSec['Ur (-)'].values*U); c+=1
                    df.insert(c , 'AB{:d}N{:03d}Vindx_[m/s]'.format(iB,ie)   ,        uix   ); c+=1 # TODO
                    df.insert(c , 'AB{:d}N{:03d}Vindy_[m/s]'.format(iB,ie)   ,        uiy   ); c+=1 # TODO

                    if dfDS is not None:
                        dfSecDS = dfBld_DS[dfBld_DS['Element']==ie]
                        df.insert(c , 'AB{:d}N{:03d}alpha_34_[deg]'.format(iB,ie) , alphaSign*dfSecDS['alpha (deg)'].values); c+=1
                        df.insert(c , 'AB{:d}N{:03d}alphaE_[deg]'.format(iB,ie)   , alphaSign*dfSecDS['alrefL (deg)'].values); c+=1 
                        df.insert(c , 'AB{:d}N{:03d}alphaED_[deg]'.format(iB,ie)  , alphaSign*dfSecDS['alrefD (deg)'].values); c+=1 
                        df.insert(c , 'AB{:d}N{:03d}adotnorm_[-]'.format(iB,ie)   , alphaSign*dfSecDS['adotnorm (-)'].values); c+=1 
                        #df.insert(c , 'AB{:d}N{:03d}AlphaDot_[-]'.format(iB,ie)  , alphaSign*dfSec['AdotNorm (-)'].values); c+=1 # TODO
                        try:
                            df.insert(c , 'AB{:d}N{:03d}activeL_[-]'.format(iB,ie)    ,           dfSecDS['DynamicFlagL'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}activeD_[-]'.format(iB,ie)    ,           dfSecDS['DynamicFlagD'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}alphaLagD_[deg]'.format(iB,ie), alphaSign*dfSecDS['alphaLagD (deg)'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}delP_[-]'.format(iB,ie)       ,           dfSecDS['delN'].values); c+=1   # NOTE SWAPPING N AND P!!!!
                            df.insert(c , 'AB{:d}N{:03d}delN_[-]'.format(iB,ie)       ,           dfSecDS['delP'].values); c+=1   # NOTE SWAPPING N AND P!!!!
                            df.insert(c , 'AB{:d}N{:03d}transA_[-]'.format(iB,ie)     ,           dfSecDS['transA'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}gammaL_[-]'.format(iB,ie)     ,           dfSecDS['gammaL'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}gammaD_[-]'.format(iB,ie)     ,           dfSecDS['gammaD'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}dalphaL_[deg]'.format(iB,ie)  , alphaSign*dfSecDS['dalphaL (deg)'].values); c+=1 
                            df.insert(c , 'AB{:d}N{:03d}dalphaD_[deg]'.format(iB,ie)  , alphaSign*(dfSecDS['alpha (deg)'].values -dfSecDS['alrefD (deg)'].values))
                            df.insert(c , 'AB{:d}N{:03d}Tu_[s]'.format(iB,ie)         ,           dfSecDS['Tu'].values); c+=1  # TODO TODO WRONG
                            df.insert(c , 'AB{:d}N{:03d}alphaDot_[rad/s]'.format(iB,ie),alphaSign*dfSecDS['alphadot'].values); c+=1  
                            df.insert(c , 'AB{:d}N{:03d}alphaDot2_[rad/s]'.format(iB,ie),np.concatenate(([0],np.diff(-dfSecDS['alpha (deg)'].values))))*np.pi/180; c+=1  # TODO TODO WRONG
                        except:
                            pass
        return df,c


    def _write(self):
        """ """
        with open(self.filename,'w') as f:
            f.write(self.toString)

    def toDataFrame(self, format='OpenFAST', alphaSign=-1):
        # ---
        df,c = self.timeDataToOpenFAST(df = self['dfTime'])
        df,c = self.elemDataToOpenFAST(dfElem=self['dfElem'],  df=df, c=c, dfDS=self['dfDS'], alphaSign=alphaSign)
        return df


    def toString(self):
        s=''
        return s

    def __repr__(self):
        s ='Class XXXX (attributes: data)\n'
        return s


