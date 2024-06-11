""" 
Setup a FAST.Farm input file based on a TurbSim box.

The extent of the high res and low res domain are setup according to the guidelines:
    https://openfast.readthedocs.io/en/dev/source/user/fast.farm/ModelGuidance.html

NOTE: the box SampleFiles/TestCase.bts is not provided as part of this repository
      Run TurbSim on SampleFiles/TestCase.inp to generate the box before running this example.

"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Local packages
from pyFAST.fastfarm import fastFarmTurbSimExtent, writeFastFarm, plotFastFarmSetup
from pyFAST.input_output.fast_input_file import FASTInputFile

MyDir=os.path.dirname(__file__)

if __name__ == '__main__':
    # --- FAST Farm input files
    templateFSTF = os.path.join(MyDir,'SampleFiles/TestCase.fstf')     # template file used for FastFarm input file, need to exist
    outputFSTF   = os.path.join(MyDir,'SampleFiles/_TestCase_mod.fstf')# new file that will be written
    # --- Parameters for TurbSim Extent
    D              = 77.0                       # Turbine diameter (m)
    hubHeight      = 78.045                     # Hub Height (m)
    extent_X_high  = 1.2                        # x-extent of high res box in diamter around turbine location
    extent_YZ_high = 1.2                         # y-extent of high res box in diamter around turbine location
    chord_max      = 5                          # maximum blade chord (m). Turbine specific.
    Cmeander       = 1.9                        # Meandering constant (-)
    BTSFilename    = os.path.join(MyDir,'SampleFiles/TestCase.bts') # TurbSim Box to be used in FAST.Farm simulation, need to exist.
    # --- Layout
    xWT = [0.0, 265.]  # x positions of turbines
    yWT = [0.0, 50.0]  # y postitions of turbines
    zWT = [0.0, 0.0 ]   # z postitions of turbines
    # --- Output list for turbine 1 (will be replicated for other turbines)
    #OutList_Sel=[
    #     'RtAxsXT1, RtAxsYT1, RtAxsZT1',
    #     'RtPosXT1, RtPosYT1, RtPosZT1',
    #     'RtDiamT1',
    #     'YawErrT1',
    #     "TIAmbT1",
    #     'RtVAmbT1',
    #     'RtVRelT1',
    #     'W1VAmbX, W1VAmbY, W1VAmbZ']
    OutList_Sel=None 

    # --- Get box extents
    FFTS = fastFarmTurbSimExtent(BTSFilename, hubHeight, D, xWT, yWT, Cmeander=Cmeander, chord_max=chord_max, extent_X=extent_X_high, extent_YZ=extent_YZ_high, meanUAtHubHeight=True)

    # --- Write Fast Farm file with layout and Low and High res extent
    writeFastFarm(outputFSTF, templateFSTF, xWT, yWT, zWT, FFTS=FFTS, OutListT1=OutList_Sel)
    print('Output file:',outputFSTF)

    # --- Visualize low&high extent and turbine positions
    fig = plotFastFarmSetup(outputFSTF, grid=True, D=D, hubHeight=hubHeight, plane='XY')
    fig = plotFastFarmSetup(outputFSTF, grid=True, D=D, hubHeight=hubHeight, plane='XZ')
    fig = plotFastFarmSetup(outputFSTF, grid=True, D=D, hubHeight=hubHeight, plane='YZ')

    # --- Finer tuning
    #fst = FASTInputFile(outputFSTF)
    #fst['InflowFile']='"../Inflow/{}_IW.dat"'.format(Case)
    #fst['WrDisDT']=1.0
    #fst['DT']=1.0
    #fst.write(outputFile)
    plt.show()



if __name__ == '__test__':
    # This example cannot be run as a test case because the BTS file is not provided in the repository
    pass
