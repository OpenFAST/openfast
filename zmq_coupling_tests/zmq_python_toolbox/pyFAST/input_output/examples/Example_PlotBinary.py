""" 
- Open and OpenFAST binary file
- Convert it to a pandas dataframe
- Plot a given output channel
"""
import os
import matplotlib.pyplot as plt
from pyFAST.input_output import FASTOutputFile

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

fastoutFilename = os.path.join(scriptDir, '../../../data/example_files/fastout_allnodes.outb')
df = FASTOutputFile(fastoutFilename).toDataFrame()
print(df.keys())
time  = df['Time_[s]']
Omega = df['RotSpeed_[rpm]']
plt.plot(time, Omega)
plt.xlabel('Time [s]')
plt.ylabel('RotSpeed [rpm]')

if __name__ == '__main__':
    plt.show()
