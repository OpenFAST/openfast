import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from parse import *
import re, os, platform
import matplotlib.colors
from pyFAST.input_output.fast_output_file import FASTOutputFile

#########################################################################################################################################
## User inputs
# Save plot and/or data?
save_fig  = False
save_data = False
fig_ext   = '.png'
R         = 65.
R_hub     = 2.

# Number of revolutions (n) to average spectra
n = 1

#########################################################################################################################################
## Paths to files
if platform.system() == 'Windows':
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
else:
    FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'openfast' + os.sep + 'reg_tests' + os.sep + 'r-tests' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
AAfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics_4.out'
OFfilename      = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics.out'
locfilename     = FAST_directory + os.sep + 'AA_ObserverLocations.dat'
output_dir      = os.path.dirname( os.path.realpath(__file__) )
outputfilename  = output_dir + os.sep + "data_output4"

#########################################################################################################################################


location = pd.read_csv(locfilename,delimiter='\s+',skiprows=[0,1],names=['x','y','z'])


AA_1 = FASTOutputFile(AAfilename).toDataFrame()
OF  = FASTOutputFile(OFfilename).toDataFrame()


with open(AAfilename, 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    n_obs    = int(f.readline().split()[-1])
    n_blades = int(f.readline().split()[-1])
    n_nodes  = int(f.readline().split()[-1])
f.close()

k        = np.ones(n_obs)
for i in range(n_obs):
    if location['x'][i] < 0:
        k[i] = -1

phi = OF['Azimuth_[deg]'] / 180. * np.pi
phi_interp = np.interp(AA_1['Time_[s]'], OF['Time_[s]'], phi)
index = []
for i in range(1, len(phi_interp)):
    if phi_interp[i] < phi_interp[i-1]:
        index.append(i)

y_b = np.linspace(R_hub, R, n_nodes)
x_b = np.zeros_like(y_b)

n_pts = index[-1] - index[-2]

for j in range(n_obs):
    x = np.zeros((n_pts,n_nodes))
    y = np.zeros((n_pts,n_nodes))

    for i in range(n_pts):
        x[i,:] = x_b * np.cos(k[j]*phi_interp[i + index[-2]]) - y_b * np.sin(k[j]*phi_interp[i + index[-2]])
        y[i,:] = x_b * np.sin(k[j]*phi_interp[i + index[-2]]) + y_b * np.cos(k[j]*phi_interp[i + index[-2]])

    z = np.array(AA_1)[index[-2]:index[-1], 1+j:1 + 30*n_obs + j:n_obs]
    fs = 10
    fig,ax=plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('y [m]', fontsize=fs+2, fontweight='bold')
    ax.set_ylabel('z [m]', fontsize=fs+2, fontweight='bold')
    ax.set_title('Observer ' + str(j), fontsize=fs+2, fontweight='bold')
    tcf=ax.tricontourf(x.flatten(),y.flatten(),z.flatten(), range(20,75))
    fig.colorbar(tcf,orientation="vertical").set_label(label = 'Overall SPL [dB]', fontsize=fs+2,weight='bold')
    if save_fig == True:
        fig_name = 'rotor_map_Obs' + str(j) + fig_ext
        fig.savefig(output_dir + os.sep + fig_name)
    plt.show()



