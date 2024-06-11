""" 
Script to postprocess multiple lin files from OpenFAST.
 - the multiple lin files are obtained for the same periodic operating point (OP)
 - the different lin files are assumed to be at different azimuthal positions
 - MBC3 is performed to go from rotating frame to fixed frame
 - binary mode file is written so that openFAST can be rerun
 - viz file is written
"""
import os
import glob
import numpy as np
import pyFAST.linearization as lin

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

# --- Script Parameters
simDir      = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Land_Lin_Rotating/') # Simulation directory
fstFile     = os.path.join(simDir,'./Main.fst') # fstFile, lin files will be assumed to be basename.i.lin
vizDict = {'VTKLinModes':15, 'VTKLinScale':10}  # Options for .viz file. Default values are: VTKLinModes=15, VTKLinScale=10, VTKLinTim=1, VTKLinTimes1=True, VTKLinPhase=0, VTKModes=None


# --- Get Campbell Diagram Data for one Operating Point (CDDOP) given an OpenFAST input file
# Perform MBC transformation based on all lin files found next to .fst file
CDDOP, MBCOP = lin.getCampbellDataOP(fstFile, writeModes=True, verbose=True) # alternatively provide a list of lin files

# --- Write Viz file
vizfile = lin.writeVizFile(fstFile, verbose=True, **vizDict)

# --- Outputs to screen
Freq,Damp = lin.printCampbellDataOP(CDDOP, nModesMax=10, nCharMaxDesc=50)

if __name__=='__main__':
    pass

if __name__=='__test__':
    np.testing.assert_almost_equal(Freq[:4],     [0.588,  0.722 , 0.842, 0.937],3)
    np.testing.assert_almost_equal(Damp[:4]*100, [63.106, 52.529, 44.01, 1.634],3)
