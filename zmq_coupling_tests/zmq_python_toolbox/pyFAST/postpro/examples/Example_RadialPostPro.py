"""
Plot average OpenFAST outputs as function of the radial position
"""
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import pyFAST.input_output as io 
import pyFAST.postpro as postpro

def main():

    # Get current directory so this script can be called from any location
    scriptDir=os.path.dirname(__file__)

    # --- Step 1: Read an openfast output file
    outFile = os.path.join(scriptDir,'../../../data/example_files/fastout_allnodes.outb')
    df = io.fast_output_file.FASTOutputFile(outFile).toDataFrame()

    # --- Step2 : Average data and extrat the radial stations
    # Averaging here is done over 1 period (avgParam=1, avgMethod='periods')
    # To get the output radial stations, a .fst file is needed
    fstFile = os.path.join(scriptDir,'../../../data/NREL5MW/Main_Onshore.fst')
    out = postpro.spanwisePostPro(FST_In=fstFile, avgMethod='periods', avgParam=1, df=df)
    dfRad_ED=out['ED_bld']; dfRad_AD = out['AD']; dfRad_BD = out['BD']

    # --- Step1&2 at once (when .outb and .fst are next to each other in same folder, with same name)
    # out = postpro.spanwisePostPro(FST_In=fstFile, avgMethod='periods', avgParam=1, out_ext='.outb')

    # --- (Optional, compute time series average)
    # Averaging here is done over the last 100s (avgParam=100, avgMethod='constantwindow')
    #dfAvg = postpro.averageDF(df, avgMethod='constantwindow' ,avgParam=100) 

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(dfRad_AD['r/R_[-]'].values, dfRad_AD['B1Cl_[-]'].values/dfRad_AD['B1Cl_[-]'].max()  , label='Lift coefficient (AeroDyn)')
    ax.plot(dfRad_ED['r/R_[-]'].values, dfRad_ED['B1TDx_[m]'].values/dfRad_ED['B1TDx_[m]'].max(), label='Flap displacement (ElastoDyn)')
    ax.set_xlabel('r/R [-]')
    ax.set_ylabel('Values normalized by max [-]')
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('FAST - Average radial outputs')

if __name__=='__main__':
    main()
    plt.show()

if __name__=='__test__':
    main()

if __name__=='__export__':
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

