""" 
Script to postprocess one linearization file from OpenFAST.
NOTE: this should not be used if the rotor is turning (multiple lin files would be needed).
This script would typically be used for a standstill analysis (no rotation),
or to compute modes of when only isolated degrees of freedom are turned on (and no rotation).

NOTE: should match the script found in the matlab-toolbox
"""
import os
import numpy as np
import pyFAST.linearization as lin

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

# Script Parameters
BladeLen     = 40.04                # Blade length, used to tune relative modal energy [m]
TowerLen     = 55.59                # Tower length, used to tune relative modal energy [m]
lin_file     = os.path.join(scriptDir,'../../../data/example_files/Standstill.1.lin') # Linearization file

# Get Campbell Diagram Data for one Operating Point (CDDOP) given a .lin file
# Performing MBC (NOTE: not stricly necessary without rotation)
CDDOP, MBC = lin.getCampbellDataOP([lin_file], BladeLen=BladeLen, TowerLen=TowerLen)

# Outputs to screen
Freq,Damp = lin.printCampbellDataOP(CDDOP, nModesMax=10, nCharMaxDesc=50)

if __name__=='__main__':
    pass

if __name__=='__test__':
    np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.450, 0.669], 3)
    np.testing.assert_almost_equal(Damp[:3]*np.pi*2*100,[1.9505,2.1309,5.0649], 4)
