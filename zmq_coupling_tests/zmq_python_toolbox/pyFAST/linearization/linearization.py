"""
Tools for OpenFAST linearization, that rely on pyFAST

Generic tools are found in campbell.py and mbc


Main functions
--------------
def campbell(caseFile, mainFst, tStart, nPerPeriod, workDir, toolboxDir, matlabExe, fastExe,  
           baseDict=None, generateInputs=True, runFast=True, runMBC=True, prefix='',sortedSuffix=None, ylim=None)
  Wrapper function to perform a Campbell diagram study
 
def writeLinearizationFiles(main_fst, workDir, operatingPointsFile, 
        nPerPeriod=36, baseDict=None, tStart=100,
        LinInputs=0, LinOutputs=0)
  Write FAST inputs files for linearization, to a given directory `workDir`.




Sub functions
-------------
def readOperatingPoints(OP_file):
    Reads an "Operating Point" delimited file to a pandas DataFrame

def defaultFilenames(OP, rpmSweep=None):
    Generate default filenames for linearization based on RotorSpeed (for rpmSweep) or WindSpeed 
 
"""
import os, glob
import pandas as pd
import numpy as np
from pyFAST.linearization.campbell import postproCampbell, plotCampbell

# TODO Alternative, Aeroelastic SE
#--- Used for functions campbell, 
import pyFAST.case_generation.runner as runner
from pyFAST.input_output import FASTInputFile
from pyFAST.input_output import FASTInputDeck
from pyFAST.case_generation.case_gen import templateReplace

def campbell(templateFstFile, operatingPointsFile, workDir, toolboxDir, fastExe,  
             nPerPeriod=36, baseDict=None, tStart=5400, trim=True, viz=False,
             trimGainPitch = 0.001, trimGainGenTorque = 300,
             maxTrq= None, 
             generateInputs=True, runFast=True, runMBC=True, prefix='',sortedSuffix=None, ylim=None, 
             removeTwrAzimuth=False, starSub=None, removeStatesPattern=None # Options for A matrices
             ):
    """ 
    Wrapper function to perform a Campbell diagram study
       see: writeLinearizationFiles, postproLinearization and postproMBC for more description of inputs.

    INPUTS:
      - templateFstFile  Main file, used as a template
      - operatingPointsFile  input file with WS, RPM, Pitch, and potentially gen torque and tower top FA
      - tStart       Time after which linearization is done (need to reach periodic steady state)
      - nPerPeriod   Number of linearizations per revolution
      - workDir      Output folder for FAST input files and linearization (will be created)
      - fastExe      Full path to a FAST exe (and dll)
      - toolboxDir   path to matlab-toolbox
      - matlabExe    path the matlab or octave exe
      - prefix:     strings such that the output files will looked like: [folder prefix ]
      - sortedSuffix use a separate file where IDs have been sorted
      - runFast      Logical to specify whether to run the simulations or not
      - removeTwrAzimuth: if False do nothing
                 otherwise discard lin files where azimuth in [60, 180, 300]+/-4deg (close to tower). 
      - starSub: if None, raise an error if `****` are present
                 otherwise replace *** with `starSub` (e.g. 0)
                 see FASTLinearizationFile. 
      - removeStatesPattern: remove states matching a giving description pattern.
                e.g:  'tower|Drivetrain'  or '^AD'
                see FASTLinearizationFile. 
    """
    Cases=pd.read_csv(caseFile); Cases.rename(columns=lambda x: x.strip(), inplace=True)

    # --- Generating input files
    if generateInputs:
        #GenTorq = Cases['GenTrq_[kNm]'] # TODO
        # Generate input files
        fastfiles= writeLinearizationFiles(mainFst, workDir, operatingPointsFile, 
                    nPerPeriod=nPerPeriod, baseDict=baseDict, tStart=tStart, trim=trim, viz=viz,
                    trimGainPitch=trimGainPitch, trimGainGenTorque=trimGainGenTorque,
                    prefix=prefix, maxTrq=maxTrq)
        # Create a batch script (optional)
        runner.writeBatch(os.path.join(workDir,'_RUN_ALL.bat'),fastfiles,fastExe=fastExe)

    # --- Run the simulations
    if runFast:
        runner.run_fastfiles(fastfiles, fastExe=fastExe, parallel=True, showOutputs=True, nCores=3)

    # --- Postprocess linearization outputs (MBC + modes ID)
    if runMBC:
        OP, Freq, Damp, _, _, modeID_file = postproCampbell(FSTfilenames, removeTwrAzimuth=removeTwrAzimuth, starSub=starSub, removeStatesPattern=removeStatesPattern, suffix=suffix)

    # ---  Plot Campbell
    fig, axes = plotCampbell(OP, Freq, Damp, sx='WS_[m/s]', UnMapped=UnMapped, ylim=ylim)
    #fig, axes = plotCampbell(Freq, Damp, sx='RPM_[rpm]', UnMapped=UnMapped)
    #  fig.savefig('{:s}CampbellWS.png'.format(suffix))
    return OP, Freq, Damp, fig


def readOperatingPoints(OP_file):
    """ 
    Reads an "Operating Point" delimited file to a pandas DataFrame

    The standard column names are (in any order):
      - RotorSpeed_[rpm], WindSpeed_[m/s], PitchAngle_[deg], GeneratorTorque_[Nm], Filename_[-]

    """
    OP=pd.read_csv(OP_file);
    for c in OP.columns:
        if c.lower().find('kn-m')>1 or c.lower().find('knm')>1:
            OP[c]*=1000
    OP.rename(columns=lambda x: x.strip().lower().replace(' ','').replace('_','').replace('(','[').split('[')[0], inplace=True)
    # Perform column replacements (tolerating small variations)
    #   format:   old:new
    OP.rename(columns={'wind':'windspeed','ws': 'windspeed'}, inplace=True)
    OP.rename(columns={'rotorspeed': 'rotorspeed', 'rpm': 'rotorspeed','omega':'rotorspeed'}, inplace=True)
    OP.rename(columns={'file': 'filename'}, inplace=True)
    OP.rename(columns={'pitch': 'pitchangle', 'bldpitch':'pitchangle'}, inplace=True)
    OP.rename(columns={'gentrq': 'generatortorque', 'gentorque':'generatortorque'}, inplace=True)
    OP.rename(columns={'ttdspfa': 'ttdspfa'}, inplace=True)
    OP.rename(columns={'oopdefl': 'oopdefl'}, inplace=True)
    OP.rename(columns={'amean': 'a_bar', 'ameanfromct': 'a_bar'}, inplace=True)

    # Standardizing column names
    OP.rename(columns={
         'rotorspeed': 'RotorSpeed_[rpm]', 
         'rotspeed': 'RotorSpeed_[rpm]', 
         'windspeed' : 'WindSpeed_[m/s]',
         'pitchangle': 'PitchAngle_[deg]',
         'generatortorque': 'GeneratorTorque_[Nm]',
         'filename': 'Filename_[-]',
         'ttdspfa': 'TTDspFA_[m]',
         'oopdefl': 'OoPDefl_[m]',
         'a_bar': 'a_bar_[-]',
         }, inplace=True)
    if 'FileName_[-]' not in OP.columns:
        OP['Filename_[-]']=defaultFilenames(OP)
    return OP

def defaultFilenames(OP, rpmSweep=None):
    """
    Generate default filenames for linearization based on RotorSpeed (for rpmSweep) or WindSpeed 
    If RotorSpeed is a field, windspeed is preferred for filenames, unless, rpmSweep is set to true as input

    INPUTS:
      - OP: structure with field RotorSpeed, and optionally WindSpeed
    OPTIONAL INPUTS:
      - rpmSweep : if present, overrides the logic: if true, rpm is used, otherwise ws
    OUTPUTS:
      - filenames: list of strings with default filenames for linearization simulations
    """

    if rpmSweep is None:
        if 'WindSpeed_[m/s]' not in OP.columns:
            rpmSweep=True
        elif len(np.unique(OP['WindSpeed_[m/s]']))==1:
            rpmSweep=True
        else:
            rpmSweep=False
    nOP=len(OP['RotorSpeed_[rpm]']);
    filenames=['']*nOP;
    for iOP, line in OP.iterrows():
        if rpmSweep:
            filenames[iOP]='rpm{:05.2f}.fst'.format(line['RotorSpeed_[rpm]'])
        else:
            filenames[iOP]='ws{:04.1f}.fst'.format(line['WindSpeed_[m/s]'])
    return filenames

def writeLinearizationFiles(main_fst, workDir, operatingPointsFile, 
        nPerPeriod=36, baseDict=None, tStart=5400, trim=True, viz=False,
        trimGainPitch = 0.001, trimGainGenTorque = 300,
        maxTrq= None, 
        LinInputs=0, LinOutputs=0):
    """
    Write FAST inputs files for linearization, to a given directory `workDir`.


    INPUTS:
      - main_fst: path to an existing  .fst file
                  This file (and the ones it refers to) will be used as templates.
                  Values of the templates can be modified using `baseDict`.
                  The parent directory of the main fst file will be copied to `workDir`.
      - workDir: directory (will be created) where the simulation files will be generated

      - operatingPointsFile: csv file containing operating conditions
      - nPerPeriod : number of linearization points per rotation (usually 12 or 36)
      - baseDict : a dictionary of inputs files keys to be applied to all simulations
                   Ignored if not provided.
                   e.g. baseDict={'DT':0.01, 'EDFile|ShftTilt':-5, 'InflowFile|PLexp':0.0}
                   see templateReplaceGeneral.
      - tStart: time at which the linearization will start. 
                When triming option is not available, this needs to be sufficiently large for
                the rotor to reach an equilibrium
      - trim: flag to use the trim option (more TODO)
      - viz : flag to setup the files to use VTK vizualization of modes (more TODO)
      - trimGainPitch     : [only for OpenFAST>2.3] Gain for Pitch trim (done around Max Torque)
      - trimGainGenTorque : [only for OpenFAST>2.3] Gain for GenTrq trim (done below Max Torque)
      - maxTrq: maximum/rated generator torque in Nm, we'll use this to determine which trim to do 
      - LinInputs:  linearize wrt. inputs (see OpenFAST documentation). {0,1,2, default:1}
      - LinOutputs: linearize wrt. outputs (see OpenFAST documentation). {0,1, default:0}

    OUTPUTS:
       - list of fst files created
    """

    # --- Optional values
    if baseDict is None:
        baseDict=dict()

    # --- Checking main fst file
    deck = FASTInputDeck(main_fst, ['ED','AD'])
    fst = deck.fst
    ED  = deck.fst_vt['ElastoDyn'] # TODO update using basedict
    AD  = deck.AD                  # TODO update using basedict

    # Update fst file with basedict if needed
    fst_keys = fst.keys()
    for k,v in baseDict.items():
        if k in fst_keys:
            print('Updating fst file key {} from {} to {}'.format(k,fst[k],v))
            fst[k] = v

    # Sanity checks
    if fst['CompAero']>0:
        if AD is None:
            raise Exception('Unable to read AeroDyn file but the file will be needed to generate cases.')
    if ED is None:
        raise Exception('Unable to infer BladeLen, ElastoDyn file not found.')
    BladeLen = ED['TipRad'] - ED['HubRad']
    hasTrim = 'TrimCase' in fst.keys()


    # --- Reading operating points
    OP = readOperatingPoints(operatingPointsFile)

    if fst['CompServo']==1:
        if 'GeneratorTorque_[Nm]' not in OP:
            raise Exception('`GeneratorTorque_[Nm]` not found in operating point file. It needs to be provided when the controller is active (`CompServo>0`).')
        if maxTrq is None:
            maxTrq = np.max(OP['GeneratorTorque_[Nm]'])
        # TODO, safety checks
        #         if SD['PCMode'),[0])
        #             error('When CompServo is 1, the PCMode should be 0 or 1');
        #         if ~ismember(GetFASTPar(paramSD,'VSContrl'),[0,1])
        #             error('When CompServo is 1, the VSContrl should be 0 or 1');

    if trim and not hasTrim:
        trim=False
        print('[WARN] Deactivating trim since not available in this version of OpenFAST')
    if viz and not hasTrim:
        viz=False
        print('[WARN] Deactivating VTK vizualization since not available in this version of OpenFAST')

    if AD is not None:
        if AD['WakeMod']==2 and AD['DBEMT_Mod'] in [1,3]:
            if 'a_bar_[-]' not in OP.keys():
                print('[WARN] Axial induction `a` not present in Operating point file, but DBEMT needs `tau1_constant`. Provide this column, or make sure your value of `tau1_const` is valid')


    # --- Generating list of parameters that vary based on the operating conditions provided
    PARAMS     = []
    for i, op in OP.iterrows():
        # Extract operating conditions (TODO, handling of missing fields)
        ws       = op['WindSpeed_[m/s]']
        rpm      = op['RotorSpeed_[rpm]']
        pitch    = op['PitchAngle_[deg]']
        filename = op['Filename_[-]']
        if 'TTDspFA_[m]' in op.keys():
            tt= op['TTDspFA_[m]']
        else:
            tt=None
        if 'OoPDefl_[m]' in op.keys():
            oop= op['OoPDefl_[m]']
        else:
            oop=None

        # Main Flags
        noAero=abs(ws)<0.001 

        nLinTimes = nPerPeriod
        if abs(rpm)<0.001:
            nLinTimes=1
            ws=1e-4

        # Determine linearization times based on RPM and nPerPeriod
        Omega = rpm/60*2*np.pi
        if trim:
            LinTimes=[9999]*nLinTimes
            Tmax   = tStart
        else:
            if abs(Omega)<0.001:
                LinTimes = [tStart]
                Tmax     = tStart+1
            else:
                T = 2*np.pi/Omega
                LinTimes = np.linspace(tStart,tStart+T,nLinTimes+1)[:-1]
                Tmax       = tStart+1.01*T
        # --- Creating "linDict", dictionary of changes to fast input files for linearization
        linDict=dict()
        linDict['__name__']     = os.path.splitext(filename)[0]
        # --- Main fst options
        linDict['TMax']         = Tmax
        linDict['TStart']       = 0

        if noAero:
            if not viz: # viz require aerodyn..
                linDict['CompAero']    = 0
                linDict['CompInflow']  = 0
        # --- Linearization options
        linDict['Linearize']    = True
        if trim:
            linDict['CalcSteady'] = True
            if fst['CompServo']==1:
                GenTrq= op['GeneratorTorque_[Nm]']
                trq_rat = np.abs((GenTrq -maxTrq))/maxTrq*100 
                print('ws {:5.2f} - GenTrq {:9.1f} - ratio: {:9.2f} '.format(ws, GenTrq, trq_rat), end='')
                if abs(op['GeneratorTorque_[Nm]']-maxTrq)/maxTrq*100 < 5 and (not noAero):
                    linDict['TrimCase'] = 3 # Adjust Pitch to get desired RPM
                    linDict['TrimGain'] = trimGainPitch; 
                    print('Trim: {} - Gain: {}'.format(3, trimGainPitch))
                    # TrimGain = .1 / (RotSpeed(iOP) * pi/30); %-> convert RotSpeed to rad/s
                    # TrimGain = TrimGain*0.1
                    # NOTE: set rated speed to a small value so that Qgen=Qrated
                    linDict['ServoFile|VS_RtGnSp'] = 0.01 
                else:
                    linDict['TrimCase'] = 2 # Adjust GenTorque to get desired RPM
                    linDict['TrimGain'] = trimGainGenTorque; 
                    print('Trim: {} - Gain: {}'.format(2, trimGainGenTorque))
                    # TrimGain = 3340 / (RotSpeed(iOP) * pi/30); %-> convert RotSpeed to rad/s
            else:
                # NOTE: in that case, trimming will just "wait", trim variable is not relevant
                linDict['TrimCase'] = 3
                linDict['TrimGain'] = 0.001
            # TrimCase - Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only when CalcSteady=True]
            # TrimTol - Tolerance for the rotational speed convergence [>eps] [used only when CalcSteady=True]
            # TrimGain - Proportional gain for the rotational speed error (rad/(rad/s) or Nm/(rad/s)) [>0] [used only when CalcSteady=True]
            # Twr_Kdmp - Damping factor for the tower (N/(m/s)) [>=0] [used only when CalcSteady=True]
            # Bld_Kdmp - Damping factor for the blade (N/(m/s)) [>=0] [used only when CalcSteady=True]
        else:
            linDict['CalcSteady'] = False
        linDict['NLinTimes']    = len(LinTimes)
        linDict['LinTimes']     = list(LinTimes)
        linDict['OutFmt']       = '"ES20.11E3"'  # Important for decent resolution
        linDict['LinInputs']    = LinInputs     # 0: none, 1: standard, 2: to get full linearizations
        linDict['LinOutputs']   = LinOutputs    # 0: none, 1: based on outlist 
        # --- Mode shape vizualization options
        if viz:
            linDict['WrVTK']        = 3
            linDict['VTK_type']     = 1
            linDict['VTK_fields']   = True
            linDict['VTK_fps']      = 30
        else:
            linDict['WrVTK']        = 0
        # --- Aero options
        if fst['CompAero']>0:
            if AD['WakeMod']==2 and AD['DBEMT_Mod'] in [1,3]:
                if 'a_bar_[-]' in OP.keys():
                    a_bar= op['a_bar_[-]']
                    linDict['AeroFile|tau1_const'] = np.around(1.1/(1-1.3*min(a_bar,0.5))*BladeLen/ws, 3)
                    print('>>> setting tau_1 to ', linDict['AeroFile|tau1_const'])
        #    linDict['AeroFile|WakeMod']    = 1 # Needed for linearization
        #    linDict['AeroFile|AFAeroMod']  = 1 # Needed for linearization
        #    linDict['AeroFile|FrozenWake'] = True # Needed for linearization
        # --- Inflow options
        if fst['CompInflow']>0:
            linDict['InflowFile|WindType'] = 1
            linDict['InflowFile|HWindSpeed'] = ws
        # --- ElastoDyn options
        linDict['EDFile|BlPitch(1)'] = pitch
        linDict['EDFile|BlPitch(2)'] = pitch
        linDict['EDFile|BlPitch(3)'] = pitch
        linDict['EDFile|RotSpeed']   = rpm
        if tt is not None:
            linDict['EDFile|TTDspFA']    = tt
        if oop is not None:
            linDict['EDFile|OoPDefl']    = oop
        # --- Servo options

        # --- Merging linDict dictionary with user override inputs
        for k,v in baseDict.items():
            if k in linDict and v != linDict[k]:
                print('Overriding key {} with value {} (previous value {})'.format(k,v,linDict[k]))
            linDict[k]=v

        PARAMS.append(linDict)

    # --- Generating all files in a workDir
    refDir    = os.path.dirname(main_fst)
    main_file = os.path.basename(main_fst)
    fastfiles = templateReplace(PARAMS,refDir,outputDir=workDir,removeRefSubFiles=True,main_file=main_file)

    return fastfiles


