import os
import shutil
import subprocess

class TurbsimReader(object):

    def read_input_file(self, input_file_name):
        inpf = open(input_file_name, 'r')

        # Runtime Options
        inpf.readline()
        inpf.readline()
        inpf.readline()
        self.Echo = inpf.readline().split()[0]
        self.RandSeed1 = inpf.readline().split()[0]
        self.RandSeed2 = inpf.readline().split()[0]
        self.WrBHHTP = inpf.readline().split()[0]
        self.WrFHHTP = inpf.readline().split()[0]
        self.WrADHH = inpf.readline().split()[0]
        self.WrADFF = inpf.readline().split()[0]
        self.WrBLFF = inpf.readline().split()[0]
        self.WrADTWR = inpf.readline().split()[0]
        self.WrFMTFF = inpf.readline().split()[0]
        self.WrACT = inpf.readline().split()[0]
        self.Clockwise = inpf.readline().split()[0]
        self.ScaleIEC = inpf.readline().split()[0]

        # Turbine/Model Specifications
        inpf.readline()
        inpf.readline()
        self.NumGrid_Z = inpf.readline().split()[0]
        self.NumGrid_Y = inpf.readline().split()[0]
        self.TimeStep = inpf.readline().split()[0]
        self.AnalysisTime = inpf.readline().split()[0]
        self.UsableTime = inpf.readline().split()[0]
        self.HubHt = inpf.readline().split()[0]
        self.GridHeight = inpf.readline().split()[0]
        self.GridWidth = inpf.readline().split()[0]
        self.VFlowAng = inpf.readline().split()[0]
        self.HFlowAng = inpf.readline().split()[0]

        # Meteorological Boundary Conditions 
        inpf.readline()
        inpf.readline()
        self.TurbModel = inpf.readline().split()[0]
        self.UserFile = inpf.readline().split()[0].replace("'","").replace('"','')
        self.IECstandard = inpf.readline().split()[0]
        self.IECturbc = inpf.readline().split()[0]
        self.IEC_WindType = inpf.readline().split()[0]
        self.ETMc = inpf.readline().split()[0]
        self.WindProfileType = inpf.readline().split()[0]
        self.ProfileFile = inpf.readline().split()[0].replace("'","").replace('"','')
        self.RefHt = inpf.readline().split()[0]
        self.URef = inpf.readline().split()[0]
        self.ZJetMax = inpf.readline().split()[0]
        self.PLExp = inpf.readline().split()[0]
        self.Z0 = inpf.readline().split()[0]


        # Meteorological Boundary Conditions 
        inpf.readline()
        inpf.readline()
        self.Latitude = inpf.readline().split()[0]
        self.RICH_NO = inpf.readline().split()[0]
        self.UStar = inpf.readline().split()[0]
        self.ZI = inpf.readline().split()[0]
        self.PC_UW = inpf.readline().split()[0]
        self.PC_UV = inpf.readline().split()[0]
        self.PC_VW = inpf.readline().split()[0]

        # Spatial Coherence Parameters
        inpf.readline()
        inpf.readline()
        self.SCMod1 = inpf.readline().split()[0]
        self.SCMod2 = inpf.readline().split()[0]
        self.SCMod3 = inpf.readline().split()[0]
        self.InCDec1 = inpf.readline().split()[0]
        self.InCDec2 = inpf.readline().split()[0]
        self.InCDec3 = inpf.readline().split()[0]
        self.CohExp = inpf.readline().split()[0]

        # Spatial Coherence Parameters
        inpf.readline()
        inpf.readline()
        self.CTEventPath = inpf.readline().split()[0]
        self.CTEventFile = inpf.readline().split()[0]
        self.Randomize = inpf.readline().split()[0]
        self.DistScl = inpf.readline().split()[0]
        self.CTLy = inpf.readline().split()[0]
        self.CTLz = inpf.readline().split()[0]
        self.CTStartTime = inpf.readline().split()[0]


class TurbsimWriter(object):

    def __init__(self, turbsiminputs):

        self.turbsiminputs = turbsiminputs
        
    def execute(self, turbsim_input_file):

         tsinp = open(turbsim_input_file, 'w')
         tsinp.write("---------TurbSim v2.00.* Input File------------------------\n")
         tsinp.write(" Turbsim input file\n")
         tsinp.write("---------Runtime Options-----------------------------------\n")

         # runtime options
         tsinp.write('{!s:<12}   Echo            - Echo input data to <RootName>.ech (flag)\n'.format(self.turbsiminputs.Echo))
         tsinp.write('{!s:<12}   RandSeed1       - First random seed  (-2147483648 to 2147483647)\n'.format(int(self.turbsiminputs.RandSeed1)))
         tsinp.write('{!s:<12}   RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'.format(self.turbsiminputs.RandSeed2))
         tsinp.write('{!s:<12}   WrBHHTP         - Output hub-height turbulence parameters in binary form?  (Generates RootName.bin)\n'.format(self.turbsiminputs.WrBHHTP))
         tsinp.write('{!s:<12}   WrFHHTP         - Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)\n'.format(self.turbsiminputs.WrFHHTP))
         tsinp.write('{!s:<12}   WrADHH          - Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)\n'.format(self.turbsiminputs.WrADHH))
         tsinp.write('{!s:<12}   WrADFF          - Output full-field time-series data in TurbSim/AeroDyn form? (Generates RootName.bts)\n'.format(self.turbsiminputs.WrADFF))
         tsinp.write('{!s:<12}   WrBLFF          - Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)\n'.format(self.turbsiminputs.WrBLFF))
         tsinp.write('{!s:<12}   WrADTWR         - Output tower time-series data? (Generates RootName.twr)\n'.format(self.turbsiminputs.WrADTWR))
         tsinp.write('{!s:<12}   WrFMTFF         - Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)\n'.format(self.turbsiminputs.WrFMTFF))
         tsinp.write('{!s:<12}   WrACT           - Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)\n'.format(self.turbsiminputs.WrACT))
         tsinp.write('{!s:<12}   Clockwise       - Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)\n'.format(self.turbsiminputs.Clockwise))
         tsinp.write('{!s:<12}   ScaleIEC        - Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]\n'.format(self.turbsiminputs.ScaleIEC))

         # Turbine/Model Specifications
         tsinp.write("\n")
         tsinp.write("--------Turbine/Model Specifications-----------------------\n")
         tsinp.write('{!s:<12}   NumGrid_Z       - Vertical grid-point matrix dimension\n'.format(self.turbsiminputs.NumGrid_Z))
         tsinp.write('{!s:<12}   NumGrid_Y       - Horizontal grid-point matrix dimension\n'.format(self.turbsiminputs.NumGrid_Y))
         tsinp.write('{!s:<12}   TimeStep        - Time step [seconds]\n'.format(self.turbsiminputs.TimeStep))
         tsinp.write('{!s:<12}   AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'.format(self.turbsiminputs.AnalysisTime))
         tsinp.write('{!s:<12}   UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n'.format(self.turbsiminputs.UsableTime))
         tsinp.write('{!s:<12}   HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n'.format(self.turbsiminputs.HubHt))
         tsinp.write('{!s:<12}   GridHeight      - Grid height [m]\n'.format(self.turbsiminputs.GridHeight))
         tsinp.write('{!s:<12}   GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'.format(self.turbsiminputs.GridWidth))
         tsinp.write('{!s:<12}   VFlowAng        - Vertical mean flow (uptilt) angle [degrees]\n'.format(self.turbsiminputs.VFlowAng))
         tsinp.write('{!s:<12}   HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'.format(self.turbsiminputs.HFlowAng))

         # Meteorological Boundary Conditions
         tsinp.write("\n")
         tsinp.write("--------Meteorological Boundary Conditions-------------------\n")
         tsinp.write('{!s:<12}   TurbModel       - Turbulence model ("IECKAI","IECVKM","GP_LLJ","NWTCUP","SMOOTH","WF_UPW","WF_07D","WF_14D","TIDAL","API","USRINP","TIMESR", or "NONE")\n'.format(self.turbsiminputs.TurbModel))
         tsinp.write('{!s:<12}   UserFile        - Name of the file that contains inputs for user-defined spectra or time series inputs (used only for "USRINP" and "TIMESR" models)\n'.format(self.turbsiminputs.UserFile))
         tsinp.write('{!s:<12}   IECstandard     - Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )\n'.format(self.turbsiminputs.IECstandard))
         tsinp.write('{!s:<12}   IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP model, not used for other models)\n'.format(self.turbsiminputs.IECturbc))
         tsinp.write('{!s:<12}   IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n'.format(self.turbsiminputs.IEC_WindType))
         tsinp.write('{!s:<12}   ETMc            - IEC Extreme Turbulence Model "c" parameter [m/s]\n'.format(self.turbsiminputs.ETMc))
         tsinp.write('{!s:<12}   WindProfileType - Velocity profile type ("LOG";"PL"=power law;"JET";"H2L"=Log law for TIDAL model;"API";"USR";"TS";"IEC"=PL on rotor disk, LOG elsewhere; or "default")\n'.format(self.turbsiminputs.WindProfileType))
         tsinp.write('{!s:<12}   ProfileFile     - Name of the file that contains input profiles for WindProfileType="USR" and/or TurbModel="USRVKM" [-]\n'.format(self.turbsiminputs.ProfileFile))
         tsinp.write('{!s:<12}   RefHt           - Height of the reference velocity (URef) [m]\n'.format(self.turbsiminputs.RefHt))
         tsinp.write('{!s:<12}   URef            - Mean (total) velocity at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'.format(self.turbsiminputs.URef))
         tsinp.write('{!s:<12}   ZJetMax         - Jet height [m] (used only for JET velocity profile, valid 70-490 m)\n'.format(self.turbsiminputs.ZJetMax))
         tsinp.write('{!s:<12}   PLExp           - Power law exponent [-] (or "default")\n'.format(self.turbsiminputs.PLExp))
         tsinp.write('{!s:<12}   Z0              - Surface roughness length [m] (or "default")\n'.format(self.turbsiminputs.Z0))

         # Non-IEC Meteorological Boundary Conditions
         tsinp.write("\n")
         tsinp.write("--------Non-IEC Meteorological Boundary Conditions------------\n")
         tsinp.write('{!s:<12}   Latitude        - Site latitude [degrees] (or "default")\n'.format(self.turbsiminputs.Latitude))
         tsinp.write('{!s:<12}   RICH_NO         - Gradient Richardson number [-]\n'.format(self.turbsiminputs.RICH_NO))
         tsinp.write('{!s:<12}   UStar           - Friction or shear velocity [m/s] (or "default")\n'.format(self.turbsiminputs.UStar))
         tsinp.write('{!s:<12}   ZI              - Mixing layer depth [m] (or "default")\n'.format(self.turbsiminputs.ZI))
         tsinp.write('{!s:<12}   PC_UW           - Hub mean uw Reynolds stress [m^2/s^2] (or "default" or "none")\n'.format(self.turbsiminputs.PC_UW))
         tsinp.write('{!s:<12}   PC_UV           - Hub mean uv Reynolds stress [m^2/s^2] (or "default" or "none")\n'.format(self.turbsiminputs.PC_UV))
         tsinp.write('{!s:<12}   PC_VW           - Hub mean vw Reynolds stress [m^2/s^2] (or "default" or "none")\n'.format(self.turbsiminputs.PC_VW))

         # Spatial Coherence Parameters
         tsinp.write('\n')
         tsinp.write(
             '--------Spatial Coherence Parameters----------------------------\n')
         tsinp.write('{!s:<12}   SCMod1           - u-component coherence model ("GENERAL", "IEC", "API", "NONE", or "default")\n'.format(
             self.turbsiminputs.SCMod1))
         tsinp.write('{!s:<12}   SCMod2           - v-component coherence model ("GENERAL", "IEC", "NONE", or "default")\n'.format(
             self.turbsiminputs.SCMod2))
         tsinp.write('{!s:<12}   SCMod3           - w-component coherence model ("GENERAL", "IEC", "NONE", or "default")\n'.format(
             self.turbsiminputs.SCMod3))
         if not type(self.turbsiminputs.InCDec1) is str:
            tsinp.write('{:<5.2f}  {:<5.2f}   InCDec1        - u-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                float(self.turbsiminputs.InCDec1[0]), float(self.turbsiminputs.InCDec1[1])))
         else:
            tsinp.write('{!s:<12}   InCDec1        - u-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                self.turbsiminputs.InCDec1))
         if not type(self.turbsiminputs.InCDec2) is str:
            tsinp.write('{:<5.2f}  {:<5.2f}   InCDec2        - v-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                float(self.turbsiminputs.InCDec2[0]), float(self.turbsiminputs.InCDec2[1])))
         else:
            tsinp.write('{!s:<12}   InCDec2        - v-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                self.turbsiminputs.InCDec2))
         if not type(self.turbsiminputs.InCDec3) is str:
            tsinp.write('{:<5.2f}  {:<5.2f}   InCDec3        - w-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                float(self.turbsiminputs.InCDec3[0]), float(self.turbsiminputs.InCDec3[1])))
         else:
            tsinp.write('{!s:<12}   InCDec3        - w-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'.format(
                self.turbsiminputs.InCDec3))
         tsinp.write('{!s:<12}   CohExp           - Coherence exponent for general model [-] (or "default")\n'.format(
             self.turbsiminputs.CohExp))

         # Coherent Turbulence Scaling Parameters
         tsinp.write('\n')
         tsinp.write('--------Coherent Turbulence Scaling Parameters-------------------\n')
         tsinp.write('{!s:<12}   CTEventPath     - Name of the path where event data files are located\n'.format(self.turbsiminputs.CTEventPath))
         tsinp.write('{!s:<12}   CTEventFile     - Type of event files ("LES", "DNS", or "RANDOM")\n'.format(self.turbsiminputs.CTEventFile))
         tsinp.write('{!s:<12}   Randomize       - Randomize the disturbance scale and locations? (true/false)\n'.format(self.turbsiminputs.Randomize))
         tsinp.write('{!s:<12}   DistScl         - Disturbance scale [-] (ratio of event dataset height to rotor disk). (Ignored when Randomize = true.)\n'.format(self.turbsiminputs.DistScl))
         tsinp.write('{!s:<12}   CTLy            - Fractional location of tower centerline from right [-] (looking downwind) to left side of the dataset. (Ignored when Randomize = true.)\n'.format(self.turbsiminputs.CTLy))
         tsinp.write('{!s:<12}   CTLz            - Fractional location of hub height from the bottom of the dataset. [-] (Ignored when Randomize = true.)\n'.format(self.turbsiminputs.CTLz))
         tsinp.write('{!s:<12}   CTStartTime     - Minimum start time for coherent structures in RootName.cts [seconds]\n'.format(self.turbsiminputs.CTStartTime))

         tsinp.close()
