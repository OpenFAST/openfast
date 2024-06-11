""" 
Convert HAWC2 aerodynamic data to AeroDyn

NOTE:
 - Position of aerodynamic center assumed to be at c/4 from c2def. TODO
"""
import os
import numpy as np
import pandas as pd

from pyFAST.converters.hawc2ToOpenfast import hawc2toAD
import matplotlib.pyplot as plt

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

def main():
    np.set_printoptions(linewidth=300)

    # ---  Hawc2 aero to AeroDyn
    # See documentation in hawc2ToOpenfast.py

    # --- hawc2toAD
    htcFilename     = os.path.join(MyDir,'../../../data/NREL5MW/hawc2/NREL_5MW_reference_wind_turbine_hs2.htc')     # readonly, hawc2 model file
    ADbldFilename_out = '_NREL5MW_AD_bld.dat' # full path of AeroDyn blade file to be written
    polarFilebase_out = '_Polars/_Polar_'     # base path where polars will be written
    correction3D = True                       # Apply 3D correction to polar data
    tsr  = 9                                  # Tip speed ratio used for 3D correction
    r_AD=None                                 # Radial position for AeroDyn. None = same as HAWC2
    #r_AD=[0.000000 , 2.625000 , 5.250000 , 7.875000 , 10.500000, 13.125000, 15.750000, 18.375000, 21.000000, 23.625000, 26.250000, 28.875000, 31.500000, 34.125000, 36.750000, 39.375000, 42.000000, 44.625000, 47.250000, 49.875000, 52.500000, 55.125000, 57.750000, 60.375000, 63.000000, 63.100000, 63.200000, 63.300000, 63.400000, 63.438000]

    # Convert 
    return hawc2toAD(htcFilename, r_AD=r_AD, ADbldFilename_out=ADbldFilename_out, polarFilebase_out=polarFilebase_out, correction3D=correction3D, tsr=tsr)

if __name__=='__main__':
    aeroNodes, polars, polarFilenames = main() 
    print('Polar files: ',polarFilenames)
    print(aeroNodes)
    plt.show()

if __name__=='__test__':
    main() 
    import shutil
    os.remove('_NREL5MW_AD_bld.dat')
    shutil.rmtree('_Polars')
