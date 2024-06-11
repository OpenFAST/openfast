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
AAfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics_1.out'
OFfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics.out'
locfilename     = FAST_directory + os.sep + 'AA_ObserverLocations_Map.dat'
output_dir      = os.path.dirname( os.path.realpath(__file__) )
outputfilename  = output_dir + os.sep + "data_output1"

#########################################################################################################################################
## Read in data, manipulate it, and plot it
# reads in file data
AA_1     = FASTOutputFile(AAfilename).toDataFrame()
OF       = FASTOutputFile(OFfilename).toDataFrame()
location = pd.read_csv(locfilename,delimiter='\s+',skiprows=[0,1],names=['x','y','z'])

# determine number of observers
num_obs = AA_1.shape[1]-1

# calculate sample time for n revolutions
rpm = OF[["RotSpeed_[rpm]"]].mean()[0]
yaw = OF[["YawPzn_[deg]"]].mean()[0] / 180. * np.pi
time_revs = n*60/rpm
tot_time = AA_1["Time_[s]"].max()
if time_revs < tot_time:
    sample_time = tot_time - time_revs
else:
    print("Error: Time for number of revolutions exceeds simulation time. Reduce n.")
    raise SystemExit('')

# slice AA dataframe for t > sample_time
AA_1 = AA_1[AA_1["Time_[s]"] > sample_time]
AA_1=AA_1.drop("Time_[s]",axis=1)

# average P over rotor revolution
AA_1 = AA_1.mean()

# merge location info with SPL info
AA_1=AA_1.reset_index()
AA_1=AA_1.drop("index",axis=1)
AA_1=pd.merge(location,AA_1,left_index=True,right_index=True)
AA_1=AA_1.rename(index=str,columns={0:"SPL"})

# contour plot of SPL for each location
if num_obs < 3:
    print("Error: Need at least 3 observers to generate contour.")
else:
    x=AA_1['x'];
    y=AA_1['y'];
    z=AA_1['SPL'];
    fs = 10
    fig,ax=plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]', fontsize=fs+2, fontweight='bold')
    ax.set_ylabel('y [m]', fontsize=fs+2, fontweight='bold')
    tcf=ax.tricontourf(x,y,z, range(58, 84, 1))
    fig.colorbar(tcf,orientation="vertical").set_label(label = 'Overall SPL [dB]', fontsize=fs+2,weight='bold')
    if save_fig == True:
        fig_name = 'directivity_map'
        fig.savefig(output_dir + os.sep + fig_name + fig_ext)
    plt.show()

# export to csv
if save_data == True:
    AA_1.to_csv(r'{}-data.csv'.format(outputfilename))

