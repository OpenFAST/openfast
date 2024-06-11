""" 
This example opens a fast output file, and interpolate a timeseries to a given radial location.
This is convenient when outputs are required at a station different from the ones used in the OpenFAST outputs.

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pyFAST.input_output as io 
import pyFAST.input_output.postpro as postpro

def main():
    # Get current directory so this script can be called from any location
    MyDir=os.path.dirname(__file__)

    # --- Read an openfast output file
    outFile = os.path.join(MyDir,'../../../data/example_files/fastout_allnodes.outb')
    df = io.fast_output_file.FASTOutputFile(outFile).toDataFrame()

    # --- Define output radial stations
    # Option 1 - Get all these locations automatically (recommended)
    fstFile = os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore.fst') # must correspond to the one used to generate outputs
    d = postpro.FASTSpanwiseOutputs(fstFile, df.columns.values)
    r_AD      = d['r_AD']
    r_ED_bld  = d['r_ED_bld']
    r_ED_twr  = d['r_ED_twr']
    r_BD      = d['r_BD']
    IR_AD     = d['IR_AD']
    IR_ED_bld = d['IR_ED_bld']
    IR_ED_twr = d['IR_ED_twr']
    IR_BD     = d['IR_BD']
    TwrLen    = d['TwrLen']
    R         = d['R']
    r_hub     = d['r_hub']
    fst       = d['fst']

    # Option 2 - Get ouputs locations for each module
    #r_ED_gag, IR_ED = ED_BldGag(fstFile)
    #r_AD_gag, IR_AD = AD_BldGag(fstFile)

    # Option 3 - Define them manually..
    #r_AD = [0.,30.,60.]  
    #r_ED = [0.,30.,60.] 

    # --- Interpolate Cl and TDx at desired radial position
    # NOTE: format need to be adjusted if you use AllOuts, or outputs at few nodes
    r    = 60         # Radial location where outputs are to be interpolated
    Cl_interp  = postpro.radialInterpTS(df, r, 'Cl_[-]', r_AD,  bldFmt='AB{:d}', ndFmt='N{:03d}')
    TDx_interp = postpro.radialInterpTS(df, r, 'TDx_[m]', r_ED_bld, bldFmt='B{:d}' , ndFmt='N{:03d}')
    #TDx_interp = postpro.radialInterpTS(df, r, 'TDx_[m]', r_ED_bld, bldFmt='B{:d}' , ndFmt='N{d}')

    # --- Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(df['Time_[s]'].values, df['AB1N017Cl_[-]'].values, label='Section before (r={}m)'.format(r_AD[16]))
    ax.plot(df['Time_[s]'].values, Cl_interp.values          , label='Interpolated (r={}m)'.format(r))
    ax.plot(df['Time_[s]'].values, df['AB1N018Cl_[-]'].values, label='Section after (r={}m)'.format(r_AD[17]))
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Cl [-]')
    ax.set_xlim([7,10])
    ax.set_ylim([0.35,0.48])
    ax.legend()
    ax.tick_params(direction='in')
    ax.set_title('FAST - interpolate radial time series')

if __name__=='__main__':
    main()
    plt.show()

if __name__=='__test__':
    main()

if __name__=='__export__':
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

