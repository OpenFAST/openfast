""" 
Convert a HAWC2 blade data to BeamDyn
"""
import os
import numpy as np
import pandas as pd

import pyFAST.converters.beamdyn as bd
import matplotlib.pyplot as plt

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

def main():
    np.set_printoptions(linewidth=300)

    # ---  Hawc2 to BeamDyn
    # See documentation in hawc2ToBeamDyn.py
    #   ref_axis: string defining how the main axis of beamdyn will be defined.
    #            'c2def': the reference axis is taken directly as Hawc2 c2def
    #            'c2def-polyfit': the reference axis is Hawc2 c2def, smoothened out using a polyfit (see poly_exp)
    #            'straight': the reference axis is straight (prebend and sweep are still included as offsets) 

    # NOTE: beamdyn uses stiffness proportional damping
    # Mu = 2 zeta / omega,  with omega the blade mode frequencies, zeta damping ratio
    #   flap,       edge,       torsion,    "edge"        "flap"    "torsion"
    Mu=[1.0E-03,   1.0E-03,   1.0E-03,    0.0014,     0.0022,    0.0022]

    # --- Htc2BeamDyn
    H2_htcfile     = os.path.join(MyDir,'../../../data/NREL5MW/hawc2/NREL_5MW_reference_wind_turbine_hs2.htc') # readonly, hawc2 model file
    BDMainTemplate = os.path.join(MyDir,'../../../data/templates/BeamDyn.dat') # readonly, template file to write main BD file
    BD_mainfile  = '_NREL5MW_BeamDyn_Created.dat'       # Name of BeamDyn file to be writen
    BD_bladefile = '_NREL5MW_BeamDyn_Blade_Created.dat' # Name of BeamDyn blade file to be written
    fig = bd.htcToBeamDyn(H2_htcfile, 'blade1', BD_bladefile, BD_mainfile, BDMainTemplate, Mu=Mu, poly_exp=[2,3,4,5,6], ref_axis='c2def', bPlot=True, interpCurvilinear=False)

    return BD_mainfile, BD_bladefile


if __name__=='__main__':
    main() 
    plt.show()

if __name__=='__test__':
    BD, BDbld = main() 
    os.remove(BD)
    os.remove(BDbld)
