""" 
Example script to create a Campbell diagram with OpenFAST
This script does not use the "trim" option, which means the user needs to provide a large simulation time (simTime) after which linearization will be done.

NOTE: This script is only an example.
      The example data is suitable for OpenFAST 2.5.

Adapt this script to your need, by calling the different subfunctions presented.

The script should be consistent with the one found in the matlab toolbox

"""

import numpy as np
import pandas as pd
import os
import pyFAST.linearization as lin
import pyFAST.case_generation.runner as runner

import matplotlib.pyplot as plt

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

def campbell_example(writeFSTfiles=True, runFAST=True, postproLin=True):
    """  Main flags
     writeFSTfiles = True # Write OpenFAST input files based on template and operatingPointsFile
     runFAST       = True # Run OpenFAST
     postproLin    = True # Postprocess the linearization outputs (*.lin)
    """

    # --- Parameters to generate linearization input files
    simulationFolder    = os.path.join(scriptDir, '../../../data/NREL5MW/_5MW_Land_Lin_Trim/')  # Output folder for input files and linearization (will be created)
    templateFstFile     = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Templates/Main.fst')  # Main file, used as a template
    operatingPointsFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Templates/LinearizationPoints_NoServo.csv')
    nPerPeriod       = 36   # Number of linearization per revolution

    # --- Parameters to run OpenFAST
    fastExe = os.path.join(scriptDir, '../../../data/openfast.exe') # Path to a FAST exe (and dll) 

    # --- Step 1: Write OpenFAST inputs files for each operating points 
    baseDict = {'DT':0.01} # Example of how inputs can be overriden (see case_gen.py templateReplace)
    fstFiles = lin.writeLinearizationFiles(templateFstFile, simulationFolder, operatingPointsFile, nPerPeriod=nPerPeriod, baseDict=baseDict)

    # Create a batch script (optional)
    runner.writeBatch(os.path.join(simulationFolder,'_RUN_ALL.bat'), fstFiles, fastExe=fastExe)

    # --- Step 2: run OpenFAST 
    if runFAST:
        runner.run_fastfiles(fstFiles, fastExe=fastExe, parallel=True, showOutputs=True, nCores=4)

    # --- Step 3: Run MBC, identify Modes, generate CSV files, and binary modes
    if postproLin:
        OP, Freq, Damp, _, _, modeID_file = lin.postproCampbell(fstFiles, writeModes=True, verbose=True)
        # Edit the modeID file manually to identify the modes
        print('[TODO] Edit this file manually: ',modeID_file)

        # --- Step 4: Campbell diagram plot
        csvFile = os.path.join(simulationFolder, 'Campbell_ModesID.csv') # <<< TODO Change me if manual identification is done
        fig, axes, figName = lin.plotCampbellDataFile(csvFile, ws_or_rpm='ws', ylim=[0,4])
        #  fig.savefig(figName+'.png')



        # --- Step 5: Generate visualization data (for advanced users)

        # --- Step 5a: Write viz files (Only useful if OpenFAST was run with WrVTK=3)
        vizDict = {'VTKLinModes':10, 'VTKLinScale':10}  # Options for .viz file. Default values are:
        vizFiles = lin.writeVizFiles(fstFiles, verbose=True, **vizDict)

        # --- Step 5b: Run FAST with VIZ files to generate VTKs
        simDir = os.path.dirname(fstFiles[0])
        ### Option 1 write a batch file and run it
        # batchfile = runner.writeBatch(os.path.join(simDir,'_RUNViz.bat'), vizFiles, fastExe=fastExe, flags='-VTKLin')
        # runner.runBatch(batchfile)
        ### Option 2: direct calls
        # runner.run_cmds(vizFiles, fastExe, showOutputs=True, flags=['-VTKLin'])

        # --- Step 5c: Convert VTKs to AVI - TODO
        # %       Also, this is experimental and users might need to adapt the inputs and batchfile content
        #     pvPython          = 'pvpython'; % path to paraview-python binary
        #     pythonPlotScript  = 'C:/Work/FAST/matlab-toolbox/Campbell/plotModeShapes.py'; % path to python plot script
        #     paraviewStateFile = 'C:/Work/FAST/matlab-toolbox/Campbell/ED_Surfaces.pvsm';  % path  to paraview State file
        #     writeAVIbatch([simulationFolder '/_RunAVI.bat'], simulationFolder, operatingPointsFile, pvPython, pythonPlotScript, paraviewStateFile);


if __name__=='__main__':
    campbell_example(runFAST=False)
    plt.show()

if __name__=='__test__':
    campbell_example(runFAST=False, postproLin=False)
    pass # this example needs an openfast binary

