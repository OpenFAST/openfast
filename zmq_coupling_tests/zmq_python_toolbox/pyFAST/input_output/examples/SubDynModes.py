""" 
- Read a SubDyn summary file and extract the Guyan and Craig Bampton Modes
- Convert the SubDyn file to JSON for 3D visualization
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
from pyFAST.input_output.fast_summary_file import FASTSummaryFile

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)


# Read SubDyn summary file
filename = os.path.join(scriptDir, '../../../data/example_files/FASTSum_5MW_OC3Mnpl.SD.sum.yaml')
sds = FASTSummaryFile(filename)

# Optional: Extract Guyan and Craig-Bampton modes
#dispGY, posGY, _, dispCB, posCB, _ = sds.getModes()

# Extract Guyan and Craig-Bampton modes and store them in a DataFrame
df = sds.toDataFrame() #sortDim=2, removeZero=True)
print(df.keys())
plt.plot(df['z_[m]'], df['GuyanMode1x_[m]'], 'o')
plt.plot(df['z_[m]'], df['GuyanMode5x_[m]'], 'o')


# Store as JSON for 3d visualization with viz3danim
sds.toJSON('_OUT.json')
sds.toGraph().toJSON('_OUT2.json')

if __name__ == '__main__':
    plt.show()
if __name__=='__test__':
    np.testing.assert_array_almost_equal(df['z_[m]'].values[::2], [-20.   , -13.333,  -6.667,   0.   ,   6.667])
    np.testing.assert_array_almost_equal(df['GuyanMode1x_[m]'].values[::2], [0.      , 0.446041, 1.270304, 2.170204, 2.843171])
    try:
        os.remove('_OUT.json')
        os.remove('_OUT2.json')
    except:
        pass
