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
AAfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics_3.out'
OFfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics.out'
output_dir      = os.path.dirname( os.path.realpath(__file__) )
outputfilename  = output_dir + os.sep + "data_output2"

#########################################################################################################################################
## Read in data, manipulate it, and plot it

# Read in file data
AA_3 = FASTOutputFile(AAfilename).toDataFrame()
OF   = FASTOutputFile(OFfilename).toDataFrame()

# Determine number of observers
num_obs = int((AA_3.shape[1]-1)/(7*34))

# Calculate sample time for n revolutions
rpm = OF[["RotSpeed_[rpm]"]].mean()[0]
time_revs = n*60/rpm
tot_time = AA_3["Time_[s]"].max()
if time_revs < tot_time:
    sample_time = tot_time - time_revs
else:
    print("Error: Time for number of revolutions exceeds simulation time. Reduce n.")
    raise SystemExit('')

# Slice AA dataframe for t > sample_time
AA_3 = AA_3[AA_3["Time_[s]"] > sample_time]
AA_3=AA_3.drop("Time_[s]",axis=1)

# Average SPL for each observer
AA_3 = AA_3.mean()

# Manipulate PD dataframes
# convert to dataframe with appropriate columns
cols = ['Observer','Mechanism','Frequency (Hz)','SPL (dB)']
aa_3 = pd.DataFrame(columns=cols)
for i in AA_3.index:
    nums = re.findall(r"[-+]?\d*\.\d+|\d+",i)
    aa_3.loc[len(aa_3)] = [nums[0],nums[2],nums[1],AA_3[i]]

AA_3 = aa_3

# rename mechanism for legend
for i in range(0,AA_3.last_valid_index()+1):
    if AA_3.loc[i,"Mechanism"]=='1':
        AA_3.loc[i,"Mechanism"]="LBL-VS"
    if AA_3.loc[i,"Mechanism"]=='2':
        AA_3.loc[i,"Mechanism"]="TBL-TE-PS"
    if AA_3.loc[i,"Mechanism"]=='3':
        AA_3.loc[i,"Mechanism"]="TBL-TE-SS"
    if AA_3.loc[i,"Mechanism"]=='4':
        AA_3.loc[i,"Mechanism"]="TBL-TE-AoA"
    if AA_3.loc[i,"Mechanism"]=='5':
        AA_3.loc[i,"Mechanism"]="TE Bluntness"
    if AA_3.loc[i,"Mechanism"]=='6':
        AA_3.loc[i,"Mechanism"]="Tip Vortex"
    if AA_3.loc[i,"Mechanism"]=='7':
        AA_3.loc[i,"Mechanism"]="TI"

AA_3["Observer"]=AA_3["Observer"].apply(pd.to_numeric)
AA_3["Frequency (Hz)"]=AA_3["Frequency (Hz)"].apply(pd.to_numeric)
AA_3["SPL (dB)"]=AA_3["SPL (dB)"].apply(pd.to_numeric)


# Plot spectra
fs = 10
for j in range(num_obs):
    fig,ax=plt.subplots()
    plt.xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=fs+2, fontweight='bold')
    ax.set_ylabel('SPL (dB)', fontsize=fs+2, fontweight='bold')
    for i in range(7):
        plt.plot(AA_3["Frequency (Hz)"][j*34*7 + i : j*34*7 + i + 34 * 7:7], AA_3["SPL (dB)"][j*34*7 + i : j*34*7 + i + 34 * 7:7], label = AA_3.loc[i,"Mechanism"])
    ax.set_title('Observer ' + str(j), fontsize=fs+2, fontweight='bold')
    plt.grid(color=[0.8,0.8,0.8], linestyle='--')
    ax.set_ylim(0,)
    ax.legend()
    if save_fig == True:
        fig_name = 'spectra_Obs' + str(j) + fig_ext
        fig.savefig(output_dir + os.sep + fig_name)
    plt.show()

# Export to csv
if save_data == True:
    AA_3.to_csv(r'{}-data.csv'.format(outputfilename))




