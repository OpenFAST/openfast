""" 
Remap a dataframe: change names, and perform operations (e.g. change units, scale, combine signals)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
import pyFAST.input_output as io 
import pyFAST.input_output.postpro as postpro

if __name__ == '__main__':
    ColumnMap={
      'WS_[m/s]'         : '{Wind1VelX_[m/s]}'             , # create a new column from existing one
      'RtTSR_[-]'        : '{RtTSR_[-]} * 2  +  {RtAeroCt_[-]}'    , # change value of column
      'RotSpeed_[rad/s]' : '{RotSpeed_[rpm]} * 2*np.pi/60 ', # new column [rpm] -> [rad/s]
    }

    # Get current directory so this script can be called from any location
    MyDir=os.path.dirname(__file__)
    # --- Step 1: Read an openfast output file
    outFile = os.path.join(MyDir,'../../../data/example_files/fastout_allnodes.outb')
    df = io.fast_output_file.FASTOutputFile(outFile).toDataFrame()
    # Change columns based on formulae, potentially adding new columns
    df = postpro.remap_df(df, ColumnMap)

    #df.to_csv('_Out.csv', index=False)
