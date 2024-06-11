""" 
Open an AeroDyn input file, change some parameters (air density) and write to a new file.
"""

import os
from pyFAST.input_output import FASTInputFile

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

filename = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/AeroDyn.dat')
f = FASTInputFile(filename)
f['TwrAero'] = True
f['AirDens'] = 1.225
# f.write('_AeroDyn_Changed.dat')
