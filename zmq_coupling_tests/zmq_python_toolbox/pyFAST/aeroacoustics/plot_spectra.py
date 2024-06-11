import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from parse import *
import re, os, platform
from pyFAST.input_output.fast_output_file import FASTOutputFile

#########################################################################################################################################
## User inputs
# Save plot and/or data?
save_fig  = False
save_data = False
fig_ext   = '.png'

# Number of revolutions (n) to average spectra
n = 1

#########################################################################################################################################
## Paths to files
if platform.system() == 'Windows':
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
else:
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'openfast' + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
AAfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics_2.out'
OFfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics.out'
output_dir      = os.path.dirname( os.path.realpath(__file__) )
outputfilename  = output_dir + os.sep + "data_output2"

#########################################################################################################################################
## Read in data, manipulate it, and plot it

# Read in file data
AA_2 = FASTOutputFile(AAfilename).toDataFrame()
OF   = FASTOutputFile(OFfilename).toDataFrame()

# Determine number of observers
num_obs = int((AA_2.shape[1]-1)/34)

# Calculate sample time for n revolutions
rpm = OF[["RotSpeed_[rpm]"]].mean()[0]
time_revs = n*60/rpm
tot_time = AA_2["Time_[s]"].max()
if time_revs < tot_time:
    sample_time = tot_time - time_revs
else:
    print("Error: Time for number of revolutions exceeds simulation time. Reduce n.")
    raise SystemExit('')

# Slice AA dataframe for t > sample_time
AA_2 = AA_2[AA_2["Time_[s]"] > sample_time]
AA_2=AA_2.drop("Time_[s]",axis=1)

# Average SPL for each observer
AA_2 = AA_2.mean()

# Manipulate PD dataframes
cols = ['Observer','Frequency (Hz)','SPL (dB)'] 
aa_2 = pd.DataFrame(columns=cols)
for i in AA_2.index:
    nums = re.findall(r"[-+]?\d*\.\d+|\d+",i)
    aa_2.loc[len(aa_2)] = [nums[0],nums[1],AA_2[i]]
AA_2 = aa_2
AA_2["Frequency (Hz)"]=AA_2["Frequency (Hz)"].apply(pd.to_numeric)
AA_2["SPL (dB)"]=AA_2["SPL (dB)"].apply(pd.to_numeric)

# Plot spectra
fs = 10
fig,ax=plt.subplots()
plt.xscale('log')
ax.set_xlabel('Frequency (Hz)', fontsize=fs+2, fontweight='bold')
ax.set_ylabel('SPL (dB)', fontsize=fs+2, fontweight='bold')
for i in range(num_obs):
    plt.plot(AA_2["Frequency (Hz)"][i*34:i*34 + 34], AA_2["SPL (dB)"][i*34:i*34 + 34], label = 'Observer ' + str(i))
plt.grid(color=[0.8,0.8,0.8], linestyle='--')
ax.legend()
if save_fig == True:
    fig_name = 'spectra' + fig_ext
    fig.savefig(output_dir + os.sep + fig_name)
plt.show()

# Export to csv
if save_data == True:
    AA_2.to_csv(r'{}-data.csv'.format(outputfilename))




