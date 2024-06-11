"""
Example to illustrate how to manipulate different files used within OpenFAST:
    - read data from an OpenFAST model
    - write modified files

The main syntax and methods are as follows:
- object = FASTInputFile(filename): read an input file, returns a dictionary-like
- object.keys(): list the available keys
- object.toDataFrame(): attempts to convert object to a pandas DataFrame
- object.write(filename): write the object to a file (may overwrite)


"""
import os
import numpy as np
from pyFAST.input_output import FASTInputFile
import matplotlib.pyplot as plt


# Get current directory so this script can be called from any location
MyDir=os.path.dirname(__file__)


# --- ElastoDyn file (geometry, inertias, etc.)
# Read some geometrical data, modify them, write a new file
print('------------------- ElastoDyn -------------------------------------')
EDFilename = os.path.join(MyDir,'../../../data/NREL5MW/onshore/ElastoDyn.dat')
ED = FASTInputFile(EDFilename)
print('> Keys:',ED.keys())
print('> Hub radius: ',ED['HubRad'])
print('> Tip radius: ',ED['TipRad'])
print('> Hub mass:   ',ED['HubMass'])
ED['TipRadius'] = 64 # Modifying the data
#ED.write('_NewFile.dat') # write a new file with modified data


# --- AeroDyn blade 
# Read an AeroDyn blade file, convert to a dataframe for convenience
print('------------------- AeroDyn blade ---------------------------------')
bldFilename = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_AeroDyn_blade.dat')
bld = FASTInputFile(bldFilename)
print('> Keys:',bld.keys())
data = bld['BldAeroNodes']   # Accessing the blade nodal data directly (numpy array, can be modified)
bld['BldAeroNodes'][:,5] *=2 # Modify (multiply chord by 2)
df  = bld.toDataFrame()      # Using a dataframe for easy export/manipulation
print('> Dataframe columns:',df.columns.values)
df.plot('BlSpn_[m]','BlChord_[m]')

# --- Profile file 
print('------------------- Profile data ----------------------------------')
PFilename = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/Airfoils/DU21_A17.dat')
P = FASTInputFile(PFilename)
print('> Keys:',P.keys())
data = P['AFCoeff']  # Accessing the polar data directly (numpy array, can be modified)
df = P.toDataFrame() # Using a dataframe for easy export/manipulation
print('> Dataframe columns:',df.columns.values)
df.plot('Alpha_[deg]','Cl_[-]')


# --- ElastoDyn blade file 
print('------------------- ElastoDyn blade -------------------------------')
bldFilename = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_Blade.dat')
bld = FASTInputFile(bldFilename)
print('> Keys:',bld.keys())
data = bld['BldProp']   # Accessing the blade nodal data directly (numpy array, can be modified)
df = bld.toDataFrame()  # Using a dataframe for easy export/manipulation
print('> Dataframe columns:',df.columns.values)
df.plot('BlFract_[-]','BMassDen_[kg/m]')

# --- BeamDyn blade file 
print('------------------- BeamDyn blade file ----------------------------')
BDFilename = os.path.join(MyDir,'../../../data/NREL5MW/5MW_Baseline/NRELOffshrBsline5MW_BeamDyn_Blade.dat')
BD = FASTInputFile(BDFilename)
print('> Keys:',BD.keys())
BP_ref = BD['BeamProperties'] # access data directly (dictionary)
# Modify blade property in for loop, write new files
#for loop:
#    BP     = BP_ref.copy()
#    BP['K']= BP['K']*i # modify some of the dictionary
#    BD['BeamProperties'] =BP
#    BD.write('newfile')





if __name__ == '__main__':
    plt.show()
