""" 
- Open and OpenFAST binary file
- Convert it to a pandas dataframe
- Compute damage equivalent load for a given Wohler exponent
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pyFAST.input_output import FASTOutputFile
from pyFAST.postpro import equivalent_load

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

# Read an openFAST binary
fastoutFilename = os.path.join(scriptDir, '../../../data/example_files/fastout_allnodes.outb')
df = FASTOutputFile(fastoutFilename).toDataFrame()


# Compute equivalent load for one signal and Wohler slope
m = 1 # Wohler slope
Leq = equivalent_load(df['Time_[s]'], df['RootMyc1_[kN-m]'], m=m)
print('Leq ',Leq)
# Leq = equivalent_load(df['Time_[s]'], df['RootMyc1_[kN-m]'], m=m, method='fatpack') # requires package fatpack


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    np.testing.assert_almost_equal(Leq , 284.30398, 3)
