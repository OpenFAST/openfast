import numpy as np
import os

import pyFAST.case_generation.case_gen as case_gen
import pyFAST.case_generation.runner as runner
import pyFAST.input_output.postpro as postpro

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

def PowerCurveParametricExample1():
    """ Example to run a set of FAST simulations to determine a power curve.
    In this example, the WS, RPM and Pitch are set within a for loop.
    If the controller and generator are active, these are just "initial conditions".
    Additional parameters may be set by adjusting the BaseDict.

    This script is based on a reference directory which contains a reference main input file (.fst)
    Everything is copied to a working directory.
    The different fast inputs are generated based on a list of dictionaries, named `PARAMS`.
    For each dictionary:
       - they keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `TMax`.
           These should correspond to whater name of the variable is used in the FAST inputs files.
       - they values are the values corresponding to this parameter
    """
    # --- Parameters for this script
    FAST_EXE  = os.path.join(MyDir, '../../../data/openfast.exe') # Location of a FAST exe (and dll)
    ref_dir   = os.path.join(MyDir, '../../../data/NREL5MW/')     # Folder where the fast input files are located (will be copied)
    main_file = 'Main_Onshore.fst'               # Main file in ref_dir, used as a template
    work_dir  = '_NREL5MW_PowerCurveParametric/'     # Output folder (will be created)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS    = [3,5,7,9 ,11,13,15]
    RPM   = [5,6,7,10,10,10,10] # initial conditions
    PITCH = [0,0,0,0 ,5 ,10,15] # initial conditions
    BaseDict = {'TMax': 100, 'DT': 0.01, 'DT_Out': 0.1}
    #BaseDict = case_gen.paramsNoController(BaseDict)
    #BaseDict = case_gen.paramsStiff(BaseDict)
    #BaseDict = case_gen.paramsNoGen(BaseDict)
    PARAMS=[]
    for wsp,rpm,pitch in zip(WS,RPM,PITCH): # NOTE: same length of WS and RPM otherwise do multiple for loops
        p=BaseDict.copy()
        p['EDFile|RotSpeed']       = rpm
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch
        p['InflowFile|HWindSpeed'] = wsp
        p['InflowFile|WindType']   = 1 # Setting steady wind
        p['__name__']              = 'ws{:04.1f}'.format(p['InflowFile|HWindSpeed'])
        PARAMS.append(p)

    # --- Generating all files in a workdir
    fastFiles=case_gen.templateReplace(PARAMS, ref_dir, work_dir, removeRefSubFiles=True, main_file=main_file)
    print(fastFiles)

    # --- Creating a batch script just in case
    runner.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'), fastFiles,fastExe=FAST_EXE)
    # --- Running the simulations
    print('>>> Running {} simulations in {} ...'.format(len(fastFiles), work_dir))
    runner.run_fastfiles(fastFiles, fastExe=FAST_EXE, parallel=True, showOutputs=False, nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastFiles]

    avg_results = postpro.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print('>>> Average results:')
    print(avg_results)
    avg_results.to_csv('_PowerCurve1.csv',sep='\t',index=False)


def PowerCurveParametricExample2():
    """ Example to run a set of FAST simulations to determine a power curve.
    In this example, the WS, RPM and Pitch are set within a for loop.
    If the controller and generator are active, these are just "initial conditions".
    Additional parameters may be set by adjusting the BaseDict.

    This script is based on a reference directory which contains a reference main input file (.fst)
    Everything is copied to a working directory.
    The different fast inputs are generated based on a list of dictionaries, named `PARAMS`.
    For each dictionary:
       - they keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `TMax`.
           These should correspond to whater name of the variable is used in the FAST inputs files.
       - they values are the values corresponding to this parameter
    """
    # --- Parameters for this script
    FAST_EXE  = os.path.join(MyDir, '../../../data/openfast.exe') # Location of a FAST exe (and dll)
    ref_dir   = os.path.join(MyDir, '../../../data/NREL5MW/')     # Folder where the fast input files are located (will be copied)
    main_file = 'Main_Onshore.fst'                # Main file in ref_dir, used as a template
    work_dir  = '_NREL5MW_PowerCurveParametric2/'     # Output folder (will be created)
    out_Ext   = '.outb' # Output extension

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS    = [3,5,7,9 ,11,13,15]
    RPM   = [5,6,7,10,10,10,10] 
    PITCH = [0,0,0,0 ,5 ,10,15] 
    BaseDict = {'TMax': 10, 'DT': 0.01, 'DT_Out': 0.1}
    PARAMS = case_gen.paramsWS_RPM_Pitch(WS, RPM, PITCH, baseDict=BaseDict, flatInputs=True)

    # --- Generating all files in a workdir
    fastFiles = case_gen.templateReplace(PARAMS, ref_dir, work_dir, removeRefSubFiles=True, removeAllowed=True, main_file=main_file)

    # --- Creating a batch script just in case
    runner.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'), fastFiles,fastExe=FAST_EXE)

    # --- Running the simulations
    runner.run_fastfiles(fastFiles, fastExe=FAST_EXE, parallel=True, showOutputs=False, nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+out_Ext for f in fastFiles]
    avg_results = postpro.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=10, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    print('>>> Average results:')
    print(avg_results)
    avg_results.to_csv('_PowerCurve2.csv',sep='\t',index=False)



if __name__=='__main__':
    PowerCurveParametricExample1()
    PowerCurveParametricExample2()

if __name__=='__test__':
    # Need openfast.exe, doing nothin
    pass

