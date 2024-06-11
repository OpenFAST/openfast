""" 
Create a TurbSim input file for a FAST.Farm simulation:
  - The x-y extents of the box are large enough to accomodate all turbines
  - The z extent is large enough to include the rotor, start from the user specified `zbot`,
      and accomodates for the meandering in the vertical direction .
  - The dy and dz resolution is set of the maximum chord of the turbine
  - The dt resolution is set based on the maximum frequency expected to be relevant for dynamics
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local
from pyFAST.fastfarm.TurbSimCaseCreation import TSCaseCreation

MyDir=os.path.dirname(__file__)

# --- Define parameters necessary for this script
OldTSFile = os.path.join(MyDir, 'SampleFiles/TestCase.inp'     )   # Template file used for TurbSim, need to exist
NewTSFile = os.path.join(MyDir, 'SampleFiles/_TestCase_mod.inp')   # New file that will be written
D     = 77.0                                                # Turbine diameter (m)
HubHt = 78.045                                              # Hub Height (m)
Vhub  = 6                                                   # mean wind speed at hub height (m/s)
TI    = 10                                                  # turbulence intensity at hub height
PLExp = 0.2                                                 # power law exponent for shear (-)
xlocs = [0.0, 265.643]  # x positions of turbines
ylocs = [0.0, 50.0   ]  # y postitions of turbines

# --- "Optional" inputs
cmax     = 5   # maximum blade chord (m). Turbine specific.
fmax     = 5.0 # maximum excitation frequency (Hz). Turbine specific, 5Hz is satisfactory for modern multi-MW turbine.
zbot     = 1.0 # vertical start of the turbulence box (m). Depend on hub height and expected vertical meandering of wakes.
Cmeander = 1.9 # Meandering constant (-)

# --- Use TurbSim Case Creation class to write a new TurbSim file
Case = TSCaseCreation(D, HubHt, Vhub, TI, PLExp, x=xlocs, y=ylocs, zbot=zbot, cmax=cmax, fmax=fmax, Cmeander=Cmeander)
# Rewrite TurbSim Input File
Case.writeTSFile(OldTSFile, NewTSFile, tmax=5, turb=1)
print('NOTE: run TurbSim to generate this new BTS file.')

# --- Visualize low extent and turbine positions
fig, ax  = Case.plotSetup()




if __name__ == '__main__':
    plt.show()

if __name__ == '__test__':
    pass




