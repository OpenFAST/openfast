import numpy as np
import os
import matplotlib.pyplot as plt

import pyFAST.case_generation.case_gen as case_gen
import pyFAST.input_output.postpro as postpro

# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)

def CPLambdaExample():
    """ Example to determine the CP-CT Lambda Pitch matrices of a turbine.
    This script uses the function CPCT_LambdaPitch which basically does the same as Parametric Examples given in this folder
    """
    FAST_EXE  = os.path.join(MyDir, '../../../data/openfast.exe') # Location of a FAST exe (and dll)
    ref_dir   = os.path.join(MyDir, '../../../data/NREL5MW/')     # Folder where the fast input files are located (will be copied)
    main_file = 'Main_Onshore.fst'  # Main file in ref_dir, used as a template

    # --- Computing CP and CT matrices for range of lambda and pitches
    Lambda = np.linspace(0.1,10,3)
    Pitch  = np.linspace(-10,10,4)

    CP,CT,Lambda,Pitch,MaxVal,result = case_gen.CPCT_LambdaPitch(ref_dir,main_file,Lambda,Pitch,fastExe=FAST_EXE,showOutputs=False,nCores=4,TMax=10)

    print('CP max',MaxVal)

    # --- Plotting matrix of CP values
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    LAMBDA, PITCH = np.meshgrid(Lambda, Pitch)
    CP[CP<0]=0
    surf = ax.plot_surface(LAMBDA, PITCH, np.transpose(CP), cmap=cm.coolwarm, linewidth=0, antialiased=True,alpha=0.8)
    ax.scatter(MaxVal['lambda_opt'],MaxVal['pitch_opt'],MaxVal['CP_max'],c='k',marker='o',s=20)
    ax.set_xlabel('lambda')
    ax.set_ylabel('pitch')
    ax.set_zlabel('CP')
    fig.colorbar(surf, shrink=0.5, aspect=5)



if __name__=='__main__':
    CPLambdaExample()
    plt.show()
if __name__=='__test__':
    # Need openfast.exe, doing nothing
    pass
