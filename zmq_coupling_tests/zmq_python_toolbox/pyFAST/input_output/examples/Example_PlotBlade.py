""" 
- Open an AeroDyn blade file
- Plot the chord as function of span
"""

import os
import matplotlib.pyplot as plt
from pyFAST.input_output import FASTInputFile

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

df = FASTInputFile(os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn_blade.dat')).toDataFrame()
print(df.keys())
plt.plot(df['BlSpn_[m]'], df['BlChord_[m]'])
plt.xlabel('Span [m]')
plt.ylabel('Chord [m]')

if __name__ == '__main__':
    plt.show()
