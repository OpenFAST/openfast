# -*- coding: utf-8 -*-
"""
  - The x-y extents of the box are large enough to accomodate all turbines
  - The z extent is large enough to include the rotor, start from the user specified `zbot`,
      and accomodates for the meandering in the vertical direction .
  - The dy and dz resolution is set of the maximum chord of the turbine
  - The dt resolution is set based on the maximum frequency expected to be relevant for dynamics

@author: kshaler
"""
import os, glob, struct
import numpy as np

class TSCaseCreation:
    
    def __init__(self, D, HubHt, Vhub, TI, PLexp, x, y, z=None, zbot=1.0, cmax=5.0, fmax=5.0,
                 Cmeander=1.9, boxType='highres', high_ext=1.2, low_ext=None, ds_low=None):
        """
        Instantiate the object. 
        
        Parameters
        ----------
        D        :   float,
                    rotor diameter (m)
        HubHt    :   float,
                    turbine hub height (m)
        Vhub	 :   float,
                    mean wind speed at hub height (m/s)
        TI	     :   float,
                   turbulence intensity at hub height
        PLexp    :   float,
                    power law exponent for shear (-)
        x, y, z  :   float,
                    x-, y-, and z-location of turbine, respectively
                    if z is None, z is set to 0
        cmax     :   float,
                    maximum blade chord (m). If not specified, set to NREL 5MW value.
        fmax     :   float,
                    maximum excitation frequency (Hz). If not specified set to NREL 5MW tower value.
        boxType  :   str,
                    box type, either 'lowres' or 'highres'. Sets the appropriate dt, dy, and dz for discretization
                    Defaults to `highres` for backward compatibility
        high_ext :   float
                    extent of the high-res box around individual turbines (in D). This is the total length
        low_ext  :  list of floats [xmin, xmax, ymin, ymax, zabovehub]
                    extents for the low-res box. All values should be positive  If not specified, resorts to
                    computations by the manual
        """

        # Perform some checks on the input
        if low_ext is not None and len(low_ext) != 5:
            raise ValueError('low_ext not defined properly. It should be [xmin, xmax, ymin, ymax, zabovehub]')

        if low_ext is None:
            manual_mode = False
        else:
            manual_mode = True

        if ds_low is None:
            manual_ds_low = False
        else:
            manual_ds_low = True

        # Set parameters for convenience
        self.Cmeander = Cmeander
        self.boxType  = boxType
        self.high_ext = high_ext
        self.low_ext  = low_ext
        self.ds_low   = ds_low
        # Turbine parameters
        self.Turb(D, HubHt, cmax, fmax)
        # Turbine location
        self.turbLocs(x,y,z)
        # Discretization
        self.discretization(Vhub, TI, PLexp, manual_ds_low)
        # Setup domain size
        self.domainSize(zbot=zbot, manual_mode=manual_mode)
        # Determine origin
        # self.originLoc()

    def Turb(self, D, HubHt, cmax=5.0, fmax=5.0):
        """
        Define turbine parameters
        
        Parameters
        __________
        D       :   float,
                   rotor diameter (m)
        HubHt   :   float,
                   turbine hub height (m)
        tpath   :   string,
                   path to base turbine location (.fst)
        cmax    :   float,
                   maximum blade chord (m). If not specified, set to NREL 5MW value.
        fmax    :   float,
                   maximum excitation frequency (Hz). If not specified set to NREL 5MW tower value.
        """
        
        self.D     = D
        self.RefHt = HubHt
        self.cmax  = cmax
        self.fmax  = fmax
        
    def turbLocs(self,x,y,z=None):
        """
        Specify turbine locations
        
        Parameters
        ----------
        x, y, z   :   float,
               x-, y-, and z-location of turbine, respectively
        """
        self.x     = np.asarray(x)
        self.y     = np.asarray(y)
        if z is None:
            self.z     = np.asarray(y)*0
        else:
            self.z     = np.asarray(z)

    def discretization(self, Vhub, TI, Shear, manual_ds_low=False):
        '''
        Specify discretization for both the high-res and low-res boxes. Follows guidelines present at
        https://openfast.readthedocs.io/en/main/source/user/fast.farm/ModelGuidance.html#low-resolution-domain

        '''
        
        self.URef = Vhub
        self.TI = TI
        self.PLexp = Shear
        
        # Derived properties
        if self.boxType == 'lowres':
            self.dt = self.Cmeander*self.D/(10*Vhub)
            ds_low  = self.Cmeander*self.D*Vhub/150
            if manual_ds_low:
                ds_low = self.ds_low
            ds_high = self.cmax
            self.dy = np.floor(ds_low/ds_high)*ds_high
            self.dz = np.floor(ds_low/ds_high)*ds_high
            #self.dt = 1.0/(2.0*self.fmax)
            #self.dy = self.cmax
            #self.dz = self.cmax

        elif self.boxType == 'highres':
            self.dt = 1.0/(2.0*self.fmax)
            self.dy = self.cmax
            self.dz = self.cmax

        else:
            raise ValueError("boxType can only be 'lowres' or 'highres'. Stopping.")

    def domainSize(self, zbot, manual_mode=False):
    
        # Set default
        self.ymin = None
        self.ymax = None

        if self.boxType == 'lowres':
            if manual_mode:
                self.ymin = min(self.y) - self.low_ext[2]*self.D
                self.ymax = max(self.y) + self.low_ext[3]*self.D
                Zdist_Low = self.RefHt + self.low_ext[4]*self.D
            else:
                self.ymin = min(self.y)-2.23313*self.Cmeander*self.D/2  # JJ: I don't recall where these recommendations came from. I can't find them on the modelling guidance document
                self.ymax = max(self.y)+2.23313*self.Cmeander*self.D/2  # JJ: I only see the y0_Low <= WT_Y_min -3*D recommendation
                Zdist_Low = self.RefHt + self.D/2 + 2.23313*self.Cmeander*self.D/2 # JJ: ditto
            
            Ydist_Low = self.ymax - self.ymin

            self.ny = np.ceil(Ydist_Low/self.dy)+1
            self.nz = np.ceil(Zdist_Low/self.dz)+1

            # We need to make sure the number of points is odd.
            if self.ny%2 == 0:
                self.ny += 1
            if self.nz%2 == 0:
                self.nz += 1
               
            
            self.Width  = self.dy*(self.ny-1)
            self.Height = self.dz*(self.nz-1)
            
            Dgrid=min(self.Height,self.Width)

            # Set the hub height using half of the total grid height 
            self.HubHt_for_TS = zbot - 0.5*Dgrid + self.Height

        elif self.boxType=='highres':
            Ydist_high = self.high_ext*self.D
            Zdist_high = self.RefHt + self.high_ext*self.D/2 - zbot
           
            self.ny = np.ceil(Ydist_high/self.dy)+1
            self.nz = np.ceil(Zdist_high/self.dz)+1
           
            # We need to make sure the number of points is odd.
            if self.ny%2 == 0:
                self.ny += 1
            if self.nz%2 == 0:
                self.nz += 1

            self.Width  = self.dy*(self.ny-1)
            self.Height = self.dz*(self.nz-1)
           
            Dgrid = min(self.Height,self.Width)

            # Set the hub height using half of the total grid height 
            self.HubHt_for_TS = zbot - 0.5*Dgrid + self.Height

        else:
            raise ValueError("boxType can only be 'lowres' or 'highres'. Stopping.")


    def originLoc(self):
        raise NotImplementedError
        


    def plotSetup(self, fig=None, ax=None):
        """ Plot a figure showing the turbine locations and the extent of the turbulence box"""
        if fig is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(6,5))
            ax  = fig.add_subplot(111,aspect="equal")
        xmin = min(self.x)-self.D
        xmax = max(self.x)+self.D

        # high-res boxes
        for wt in range(len(self.x)):
            ax.plot(self.x[wt],self.y[wt],'x',ms=8,mew=2,label=f"WT{wt+1}")

        # low-res box
        #         ax.plot([xmin,xmax,xmax,xmin,xmin],
        #                 [ymin,ymin,ymax,ymax,ymin],'--k',lw=2,label='Low')
        if self.ymin is not None:
            ax.axhline(self.ymin, ls='--', c='k', lw=2, label='Low')
        if self.ymax is not None:
            ax.axhline(self.ymax, ls='--', c='k', lw=2)

        ax.legend(bbox_to_anchor=(1.05,1.015),frameon=False)
        ax.set_xlabel("x-location [m]")
        ax.set_ylabel("y-location [m]")
        fig.tight_layout
        return fig, ax

    def writeTSFile(self, fileIn, fileOut, NewFile=True, tpath=None, tmax=50, turb=None, verbose=0):
        """ Write a TurbSim primary input file, 
        See WriteTSFile below.
        """
        WriteTSFile(fileIn, fileOut, self, NewFile=NewFile, tpath=tpath, tmax=tmax, turb=turb, verbose=verbose)


        
def WriteTSFile(fileIn, fileOut, params, NewFile=True, tpath=None, tmax=50, turb=None, verbose=0):
    """ Write a TurbSim primary input file, 

        tpath:     string,
                   path to base turbine location (.fst)
                   only used if NewFile is False
        boxType:   string,
                   Box type, either 'lowres' or 'highres'. Writes the proper `TurbModel`
                   if boxType=='highres', `turb` needs to be specified
        turb:      int,
                   turbine number to be printed on the time series file. Only needed
                   if boxType='highres'

    """

    if params.boxType=='highres' and not isinstance(turb, int):
        raise ValueError("turb needs to be an integer when boxType is 'highres'")
    if params.boxType=='lowres' and turb is not None:
        print("WARNING: `turb` is not used when boxType is 'lowres'. Remove `turb` to dismiss this warning.")

    if NewFile == True:
        if verbose>1:  print(f'Writing a new {fileOut} file from scratch')
        # --- Writing FFarm input file from scratch
        with open(fileOut, 'w') as f:
            f.write(f'--------TurbSim v2.00.* Input File------------------------\n')
            f.write(f'for Certification Test #1 (Kaimal Spectrum, formatted FF files).\n')
            f.write(f'---------Runtime Options-----------------------------------\n')
            f.write(f'False\tEcho\t\t- Echo input data to <RootName>.ech (flag)\n')
            f.write(f'123456\tRandSeed1\t\t- First random seed  (-2147483648 to 2147483647)\n')
            f.write(f'RanLux\tRandSeed2\t\t- Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n')
            f.write(f'False\tWrBHHTP\t\t- Output hub-height turbulence parameters in binary form?  (Generates RootName.bin)\n')
            f.write(f'False\tWrFHHTP\t\t- Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)\n')
            f.write(f'False\tWrADHH\t\t- Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)\n')
            f.write(f'True\tWrADFF\t\t- Output full-field time-series data in TurbSim/AeroDyn form? (Generates RootName.bts)\n')
            f.write(f'False\tWrBLFF\t\t- Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)\n')
            f.write(f'False\tWrADTWR\t\t- Output tower time-series data? (Generates RootName.twr)\n')
            f.write(f'False\tWrFMTFF\t\t- Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)\n')
            f.write(f'False\tWrACT\t\t- Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)\n')
            f.write(f'True\tClockwise\t\t- Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)\n')
            f.write(f'0\tScaleIEC\t\t- Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]\n')
            f.write(f'\n')
            f.write(f'--------Turbine/Model Specifications-----------------------\n')
            f.write(f'{params.nz:.0f}\tNumGrid_Z\t\t- Vertical grid-point matrix dimension\n')
            f.write(f'{params.ny:.0f}\tNumGrid_Y\t\t- Horizontal grid-point matrix dimension\n')
            f.write(f'{params.dt:.6f}\tTimeStep\t\t- Time step [seconds]\n')
            f.write(f'{tmax:.4f}\tAnalysisTime\t\t- Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n')
            f.write(f'"ALL"\tUsableTime\t\t- Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n')
            f.write(f'{params.HubHt_for_TS:.3f}\tHubHt\t\t- Hub height [m] (should be > 0.5*GridHeight)\n')
            f.write(f'{params.Height:.3f}\tGridHeight\t\t- Grid height [m]\n')
            f.write(f'{params.Width:.3f}\tGridWidth\t\t- Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n')
            f.write(f'0\tVFlowAng\t\t- Vertical mean flow (uptilt) angle [degrees]\n')
            f.write(f'0\tHFlowAng\t\t- Horizontal mean flow (skew) angle [degrees]\n')
            f.write(f'\n')
            f.write(f'--------Meteorological Boundary Conditions-------------------\n')
            if params.boxType=='lowres':
                f.write(f'"IECKAI"\tTurbModel\t\t- Turbulence model ("IECKAI","IECVKM","GP_LLJ","NWTCUP","SMOOTH","WF_UPW","WF_07D","WF_14D","TIDAL","API","IECKAI","TIMESR", or "NONE")\n')
                f.write(f'"unused"\tUserFile\t\t- Name of the file that contains inputs for user-defined spectra or time series inputs (used only for "IECKAI" and "TIMESR" models)\n')
            elif params.boxType=='highres':
                f.write(f'"TIMESR"\tTurbModel\t\t- Turbulence model ("IECKAI","IECVKM","GP_LLJ","NWTCUP","SMOOTH","WF_UPW","WF_07D","WF_14D","TIDAL","API","USRINP","TIMESR", or "NONE")\n')
                f.write(f'"USRTimeSeries_T{turb}.txt"\tUserFile\t\t- Name of the file that contains inputs for user-defined spectra or time series inputs (used only for "IECKAI" and "TIMESR" models)\n')
            else:
                raise ValueError("boxType can only be 'lowres' or 'highres'. Stopping.")
            f.write(f'1\tIECstandard\t\t- Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )\n')
            f.write(f'"{params.TI:.3f}\t"\tIECturbc\t\t- IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP model, not used for other models)\n')
            f.write(f'"NTM"\tIEC_WindType\t\t- IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n')
            f.write(f'"default"\tETMc\t\t- IEC Extreme Turbulence Model "c" parameter [m/s]\n')
            f.write(f'"PL"\tWindProfileType\t\t- Velocity profile type ("LOG";"PL"=power law;"JET";"H2L"=Log law for TIDAL model;"API";"PL";"TS";"IEC"=PL on rotor disk, LOG elsewhere; or "default")\n')
            f.write(f'"unused"\tProfileFile\t\t- Name of the file that contains input profiles for WindProfileType="USR" and/or TurbModel="USRVKM" [-]\n')
            f.write(f'{params.RefHt:.3f}\tRefHt\t\t- Height of the reference velocity (URef) [m]\n')
            f.write(f'{params.URef:.3f}\tURef\t\t- Mean (total) velocity at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n')
            f.write(f'350\tZJetMax\t\t- Jet height [m] (used only for JET velocity profile, valid 70-490 m)\n')
            f.write(f'"{params.PLexp:.3f}"\tPLExp\t\t- Power law exponent [-] (or "default")\n')
            f.write(f'"default"\tZ0\t\t- Surface roughness length [m] (or "default")\n')
            f.write(f'\n')
            f.write(f'--------Non-IEC Meteorological Boundary Conditions------------\n')
            f.write(f'"default"\tLatitude\t\t- Site latitude [degrees] (or "default")\n')
            f.write(f'0.05\tRICH_NO\t\t- Gradient Richardson number [-]\n')
            f.write(f'"default"\tUStar\t\t- Friction or shear velocity [m/s] (or "default")\n')
            f.write(f'"default"\tZI\t\t- Mixing layer depth [m] (or "default")\n')
            f.write(f'"default"\tPC_UW\t\t- Hub mean u\'w\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write(f'"default"\tPC_UV\t\t- Hub mean u\'v\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write(f'"default"\tPC_VW\t\t- Hub mean v\'w\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write(f'\n')
            f.write(f'--------Spatial Coherence Parameters----------------------------\n')
            f.write(f'"IEC"\tSCMod1\t\t- u-component coherence model ("GENERAL","IEC","API","NONE", or "default")\n')
            f.write(f'"IEC"\tSCMod2\t\t- v-component coherence model ("GENERAL","IEC","NONE", or "default")\n')
            f.write(f'"IEC"\tSCMod3\t\t- w-component coherence model ("GENERAL","IEC","NONE", or "default")\n')
            f.write(f'"12.0 {0.12/(8.1*42):.8f}"\tInCDec1\t- u-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write(f'"12.0 {0.12/(8.1*42):.8f}"\tInCDec2\t- v-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write(f'"12.0 {0.12/(8.1*42):.8f}"\tInCDec3\t- w-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write(f'"0.0"\tCohExp\t\t- Coherence exponent for general model [-] (or "default")\n')
            f.write(f'\n')
            f.write(f'--------Coherent Turbulence Scaling Parameters-------------------\n')
            f.write(f'".\\EventData"\tCTEventPath\t\t- Name of the path where event data files are located\n')
            f.write(f'"random"\tCTEventFile\t\t- Type of event files ("LES", "DNS", or "RANDOM")\n')
            f.write(f'true\tRandomize\t\t- Randomize the disturbance scale and locations? (true/false)\n')
            f.write(f'1\tDistScl\t\t- Disturbance scale [-] (ratio of event dataset height to rotor disk). (Ignored when Randomize = true.)\n')
            f.write(f'0.5\tCTLy\t\t- Fractional location of tower centerline from right [-] (looking downwind) to left side of the dataset. (Ignored when Randomize = true.)\n')
            f.write(f'0.5\tCTLz\t\t- Fractional location of hub height from the bottom of the dataset. [-] (Ignored when Randomize = true.)\n')
            f.write(f'30\tCTStartTime\t\t- Minimum start time for coherent structures in RootName.cts [seconds]\n')
            f.write(f'\n')
            f.write(f'====================================================\n')
            f.write(f'! NOTE: Do not add or remove any lines in this file!\n')
            f.write(f'====================================================\n')

    else:
        print(f'Modifying {fileIn} to be {fileOut}')

        NewPars = [int(params.nz), int(params.ny), int(params.dt), format(params.HubHt,'.2f'), format(params.Height,'.2f'), format(params.Width,'.2f'), format(params.TI,'.2f'), format(params.RefHt,'.2f'), format(params.URef,'.2f'), int(params.PLexp)]
        ModVars = ['NumGrid_Z','NumGrid_Y','TimeStep','HubHt','GridHeight','GridWidth','IECturb','RefHt','URef','PLExp']
        wt=0
        with open(fileOut, 'w+') as new_file:
            with open(fileIn) as old_file:
                for line in old_file.readlines():
                    newline = line
                    for index,tmpVar in enumerate(ModVars):
                        if tmpVar in line:
                            newline = str(NewPars[index])+'\t!!Orig is:  '+line
                    if '.fst' in line:
                        newline =str('{params.x[wt]:.3f}\t\t{params.y[wt]:.3f}\t\t{params.z[wt]:.3f}\t\t{tpath}_WT{wt+1:d}.fst"\t{params.X0_High[wt]:.3f}\t\t{params.Y0_High[wt]:.3f}\t\t{params.Z0_High:.3f}\t\t{params.dX_High:.3f}\t\t{params.dY_High:.3f}\t\t{params.dZ_High:.3f}\n')
                        wt+=1
                    new_file.write(newline)


def writeTimeSeriesFile(fileOut,yloc,zloc,u,v,w,time):
    """ Write a TurbSim primary input file, 

    """

    print(f'Writing {fileOut}')
    # --- Writing TurbSim user-defined time series file
    with open(fileOut, 'w') as f:
        f.write( '--------------TurbSim v2.00.* User Time Series Input File-----------------------\n')
        f.write( '     Time series input from low-res turbsim run\n')
        f.write( '--------------------------------------------------------------------------------\n')
        f.write( '          3 nComp - Number of velocity components in the file\n')
        f.write( '          1 nPoints - Number of time series points contained in this file (-)\n')
        f.write( '          1 RefPtID - Index of the reference point (1-nPoints)\n')
        f.write( '     Pointyi Pointzi ! nPoints listed in order of increasing height\n')
        f.write( '       (m)     (m)\n')
        f.write(f'       {yloc:.5f}   {zloc:.5f}\n')
        f.write( '--------Time Series-------------------------------------------------------------\n')
        f.write( 'Elapsed Time         Point01u             Point01v           Point01w\n')
        f.write( '         (s)            (m/s)                (m/s)              (m/s)\n')
        for i in range(len(time)):
            f.write(f'\t{time[i]:.2f}\t\t\t  {u[i]:.5f}\t\t\t  {v[i]:.5f}\t\t\t {w[i]:.5f}\n')

