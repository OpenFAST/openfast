""" 
Converts an OpenFAST model to HAWC2 for a horizontal axis wind turbine.
A template "htc" file is used, and modified:
    - geometry of bodies (tower, towertop/nacelle, shaft, hub, blades)
    - orientations
    - structural files (blade only, tower todo)
    - aerodynamic files (pc-ae) generated
    - wind and aero options
    - hawcstab2 inputs
"""

import numpy as np
import pandas as pd
import os
from pyFAST.converters.openfastToHawc2 import FAST2Hawc2

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)


htcTemplate = os.path.join(MyDir,'../../../data/templates/hawc2_template.htc')
fstIn       = os.path.join(MyDir,'../../../data/NREL5MW/Main_Onshore_OF2_BD.fst')
htcOut      = os.path.join(MyDir,'_NREL5MW_from_OF.htc') # will be created
OPfile      = None

# htcTemplate = 'data/templates/hawc2_template.htc'
# fstIn       = 'data/BAR0_OF/BAR0.fst'
# htcOut      = 'data/BAR0_hawc2/BAR0_from_OF.htc' # will be created
# OPfile      = 'data/BAR0_OF/Performance/ccblade.dat' # Optional oper. file

# htcTemplate = 'templates/hawc2_template.htc'
# fstIn       = './IEA-15-240-RWT/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile_BD.fst'
# htcOut      = './IEA-15-240-RWT/IEA_15MW_RWT_from_OF.htc'
# OPfile      = None


FAST2Hawc2(fstIn, htcTemplate, htcOut, OPfile=OPfile, TwrFAFreq=0.1, TwrSSFreq=0.1, SftTorFreq=4, FPM = False)



if __name__ == '__main__':
    pass
if __name__ == '__test__':
    try:
        os.remove(htc_out)
    except:
        pass
