"""
Convert a BeamDyn blade to a HAWC2 blade.
The part of the htc file defining the c2def line is written based on a template htc file
A new hawc2 "st" is written for the blade, using either:
  - "Fully Populated Matrix" (FPM) 
  - Timoshenko beam peoperties
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from shutil import copyfile
np.set_printoptions(linewidth=1500)
import pyFAST.converters.beamdyn as bd

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

# --- Parameters
FPM                 = True # Use fully populated matrix or regular st file
ConvertBackAndForth = False # True to check how the model is converted back and forth between hawc2 and beamdyn
# Optional give some values for those, otherwise inferred by minimization...
E       = None
G       = None
A       = None
theta_p = None  # Principal axis direction
# TODO Elastic Center location for FPM

# Derived parameters
suffix='_FPM' if FPM  else ''

# --- BeamDynToHawc2
htc_template   = os.path.join(MyDir,'../../../data/templates/hawc2_template.htc') # readonly
BD_mainfile    = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_BeamDyn.dat')      # readonly
BD_bladefile   = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_BeamDyn_Blade.dat') # readonly
H2_htcfile_new = '_NREL5MW{}.htc'.format(suffix) # will be created
H2_stfile      = '_blade_st{}.st'.format(suffix) # will be created
copyfile(htc_template,  H2_htcfile_new) # Backup template
df_c2, df_st = bd.beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile_new, H2_stfile, 'blade1', A=A, E=E, G=G, theta_p_in = theta_p, FPM=FPM, verbose=True)

if ConvertBackAndForth:
    # NOTE: beamdyn uses stiffness proportional damping
    # Mu = 2 zeta / omega,  with omega the blade mode frequencies, zeta damping ratio
    #   flap        edge      torsion     "edge"      "flap"    "torsion"
    Mu=[0.021955,   0.012818,   0.012818,   0.012818,   0.021955,   0.012818]

    # --- Htc2BeamDyn
    H2_htcfile     = H2_htcfile_new                     # readonly
    BDMainTemplate = BD_mainfile                        # readonly, for template
    BD_mainfile  = '_NREL5MW_BeamDyn_2.dat'      # will be created
    BD_bladefile = '_NREL5MW_BeamDyn_Blade_2.dat' # will be created
    fig = bd.htcToBeamDyn(H2_htcfile_new, 'blade1', BD_bladefile, BD_mainfile, BDMainTemplate, Mu=Mu, poly_exp=[2,3,4,5,6], ref_axis='c2def', bPlot=False)

    # --- BeamDynToHawc2 
    H2_htcfile_new2 = '_NREL5MW_2.htc'  
    H2_stfile2      = '_blade_st_2.st'
    copyfile(H2_htcfile_new, H2_htcfile_new2)
    bd.beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile_new2, H2_stfile2, 'blade1', A=A, E=E, G=G, FPM=False, verbose=True)


if __name__ == '__test__':
    # NOTE: NREL5MW is too simple of a tests since straight
    np.testing.assert_almost_equal(df_c2['z_[m]'].values[-1]      , 61.5              )
    np.testing.assert_almost_equal(df_c2['twist_[deg]'].values[-1], 0                 )
    np.testing.assert_almost_equal(df_c2['twist_[deg]'].values[0] , -13.308           )
    np.testing.assert_almost_equal(df_st['x_cg_[m]'].values[10]   , 0.0               )
    np.testing.assert_almost_equal(df_st['x_e_[m]'].values[10]    , 0.0               )
    np.testing.assert_almost_equal(df_st['pitch_[deg]'].values[10], 0.0               )
    np.testing.assert_almost_equal(df_st['K11'].values[10]        , 403729000.0       )
    np.testing.assert_almost_equal(df_st['K44'].values[10]        , 4936840000.0      )
    np.testing.assert_almost_equal(df_st['K55'].values[10]        , 7009180000.0      )
    np.testing.assert_almost_equal(df_st['K66'].values[10]        , 1002120000.0      )
    np.testing.assert_almost_equal(df_st['K11']                   , df_st['K22'].values)

    os.remove(H2_htcfile_new)
    os.remove(H2_stfile)
