#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
    This program executes OpenFAST and a regression test for a single test case.
    The test data is contained in a git submodule, r-test, which must be initialized
    prior to running. See the r-test README or OpenFAST documentation for more info.

    Get usage with: `executeOpenfastLinearRegressionCase.py -h`
"""

import os
import sys
basepath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.sep.join([basepath, "lib"]))
import argparse
import numpy as np
import shutil
import subprocess
from eva import eigA, eig
import rtestlib as rtl
import openfastDrivers
import pass_fail
from errorPlotting import exportCaseSummary
from fast_linearization_file import FASTLinearizationFile
# from weio.fast_linearization_file import FASTLinearizationFile

##### Helper functions
excludeExt=['.out','.outb','.ech','.yaml','.sum','.log','.md']

def file_line_count(filename):
    file_handle = open(filename, 'r')
    for i, _ in enumerate(file_handle):
        pass
    file_handle.close()
    return i + 1

def isclose(a, b, rtol=1e-09, atol=0.0):
    return abs(a-b) <= max(rtol * max(abs(a), abs(b)), atol)

##### Main program

### Store the python executable for future python calls
pythonCommand = sys.executable

### Verify input arguments
parser = argparse.ArgumentParser(description="Executes OpenFAST and a regression test for a single test case.")
parser.add_argument("caseName", metavar="Case-Name", type=str, nargs=1, help="The name of the test case.")
parser.add_argument("executable", metavar="OpenFAST", type=str, nargs=1, help="The path to the OpenFAST executable.")
parser.add_argument("sourceDirectory", metavar="path/to/openfast_repo", type=str, nargs=1, help="The path to the OpenFAST repository.")
parser.add_argument("buildDirectory", metavar="path/to/openfast_repo/build", type=str, nargs=1, help="The path to the OpenFAST repository build directory.")
parser.add_argument("rtol", metavar="Relative-Tolerance", type=float, nargs=1, help="Relative tolerance to allow the solution to deviate; expressed as order of magnitudes less than baseline.")
parser.add_argument("atol", metavar="Absolute-Tolerance", type=float, nargs=1, help="Absolute tolerance to allow small values to pass; expressed as order of magnitudes less than baseline.")
parser.add_argument("-p", "-plot", dest="plot", action='store_true', help="bool to include plots in failed cases")
parser.add_argument("-n", "-no-exec", dest="noExec", action='store_true', help="bool to prevent execution of the test cases")
parser.add_argument("-v", "-verbose", dest="verbose", action='store_true', help="bool to include verbose system output")
parser.add_argument("-highpass", dest='highpass', metavar="LowPass-Filter", type=float, nargs='?', default=0.0, help="high pass filter on linearization frequencies to compare")

args = parser.parse_args()

caseName = args.caseName[0]
executable = args.executable[0]
sourceDirectory = args.sourceDirectory[0]
buildDirectory = args.buildDirectory[0]
rtol = args.rtol[0]
atol = args.atol[0]
plotError = args.plot
noExec = args.noExec
verbose = args.verbose
highpass = args.highpass

# --- Tolerances for matrix comparison
# Outputs of lin matrices have 3 decimal digits leading to minimum error of 0.001  
# Therefore the relative error to detect a change in the third decimal place
# is between 1e-3 and 1e-4. We allow a bit of margin and use rtol=2e-3
# Lin matrices have a lot of small values, so atol is quite important
rtol = 2e-3
atol = 5e-4

# --- Tolerances for frequencies 
# Low frequencies are hard to match, so we use a high atol
rtol_f=1e-2
atol_f=1e-2 

# --- Tolerances for damping
# damping ratio is in [%] so we relax the atol
rtol_d=1e-2
atol_d=1e-1  

# --- Filenames for frequency info
fileNameFreqRef="frequencies_ref.txt"
fileNameFreqNew="frequencies_new.txt"


CasePrefix=' Case: {}: '.format(caseName)
def exitWithError(msg):
    rtl.exitWithError(CasePrefix+msg)
def indent(msg, sindent='\t'):
    return '\n'.join([sindent+s for s in msg.split('\n')])


# validate inputs
rtl.validateExeOrExit(executable)
rtl.validateDirOrExit(sourceDirectory)
if not os.path.isdir(buildDirectory):
    os.makedirs(buildDirectory, exist_ok=True)

### Build the filesystem navigation variables for running openfast on the test case
regtests = os.path.join(sourceDirectory, "reg_tests")
lib = os.path.join(regtests, "lib")
rtest = os.path.join(regtests, "r-test")
moduleDirectory = os.path.join(rtest, "glue-codes", "openfast")
inputsDirectory = os.path.join(moduleDirectory, caseName)
targetOutputDirectory = os.path.join(inputsDirectory)
testBuildDirectory = os.path.join(buildDirectory, caseName)
fNameFreqRef = os.path.join(testBuildDirectory, fileNameFreqRef)
fNameFreqNew = os.path.join(testBuildDirectory, fileNameFreqNew)

# verify all the required directories exist
if not os.path.isdir(rtest):
    exitWithError("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(rtest))
if not os.path.isdir(targetOutputDirectory):
    exitWithError("The test data outputs directory, {}, does not exist. Try running `git submodule update`".format(targetOutputDirectory))
if not os.path.isdir(inputsDirectory):
    exitWithError("The test data inputs directory, {}, does not exist. Verify your local repository is up to date.".format(inputsDirectory))

# create the local output directory if it does not already exist
# and initialize it with input files for all test cases
for data in ["Ideal_Beam", "WP_Baseline"]:
    dataDir = os.path.join(buildDirectory, data)
    if not os.path.isdir(dataDir):
        rtl.copyTree(os.path.join(moduleDirectory, data), dataDir, excludeExt=excludeExt)

# Special copy for the 5MW_Baseline folder because the Windows python-only workflow may have already created data in the subfolder ServoData
dst = os.path.join(buildDirectory, "5MW_Baseline")
src = os.path.join(moduleDirectory, "5MW_Baseline")
if not os.path.isdir(dst):
    rtl.copyTree(src, dst, excludeExt=excludeExt)
else:
    names = os.listdir(src)
    for name in names:
        if name == "ServoData":
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        if os.path.isdir(srcname):
            if not os.path.isdir(dstname):
                rtl.copyTree(srcname, dstname, excludeExt=excludeExt)
        else:
            shutil.copy2(srcname, dstname)
# 
# Copying the actual test directory
# if not os.path.isdir(testBuildDirectory):
rtl.copyTree(inputsDirectory, testBuildDirectory, excludeExt=excludeExt, renameExtDict={'.lin':'.ref_lin'})

### Run openfast on the test case
if not noExec:
    caseInputFile = os.path.join(testBuildDirectory, caseName + ".fst")
    returnCode = openfastDrivers.runOpenfastCase(caseInputFile, executable)
    if returnCode != 0:
        sys.exit(returnCode*10)

### Get a all the .lin files in the baseline directory
baselineOutFiles = [f for f in os.listdir(targetOutputDirectory) if '.lin' in f]
if len(baselineOutFiles)==0:
    exitWithError("No lin files present in baseline.")


# these should all exist in the local outputs directory
localFiles = os.listdir(testBuildDirectory)
localOutFiles = [f for f in localFiles if f in baselineOutFiles]
if len(localOutFiles) != len(baselineOutFiles):
    print("Number of local files: {}".format(len(localOutFiles)))
    print ("   {}".format(localOutFiles))
    print("Number of baseline files: {}".format(len(baselineOutFiles)))
    print ("   {}".format(baselineOutFiles))
    exitWithError("An expected local solution file does not exist:")

### test for regression (compare lin files only)

def compareLin(f,file_freq_ref,file_freq_new):
    Errors     = []
    ElemErrors = []

    local_file = os.path.join(testBuildDirectory, f)
    baseline_file = os.path.join(targetOutputDirectory, f)
    local_file2 = local_file # Ellude the long filename

    # Set a prefix for all errors to identify where it comes from
    basename = os.path.splitext(os.path.basename(local_file))[0]
    ext2 = os.path.splitext(basename)[1][1:] +'.lin' # '.1' or '.AD' '.BD'
    errPrefix = CasePrefix[:-1]+ext2+': '

    if verbose:
        print(errPrefix+'file_ref:', baseline_file)
        print(errPrefix+'file_new:', local_file)

    def newError(msg):
        msg=errPrefix+msg
        if verbose:
            print(msg)
        Errors.append(msg)

    def ApplyHighPass(freq,zeta):
        freqL=np.array([])
        zetaL=np.array([])
        for i in range(len(freq)):
            if freq[i]>highpass:
                freqL = np.append(freqL,freq[i])
                zetaL = np.append(zetaL,zeta[i])
        return freqL,zetaL



    # verify both files have the same number of lines
    local_file_line_count = file_line_count(local_file)
    baseline_file_line_count = file_line_count(baseline_file)
    if local_file_line_count != baseline_file_line_count:
        Err="Local and baseline solutions have different line counts in"
        Err+="\n\tFile1:{}".format(local_file)
        Err+="\n\tFile2:{}\n\n".format(baseline_file)
        newError(Err)
        #raise Exception(Err)

    # open both files
    floc = FASTLinearizationFile(local_file)
    fbas = FASTLinearizationFile(baseline_file)

    # --- Test that they have the same variables
    kloc = floc.keys()
    kbas = fbas.keys()
    try:
        np.testing.assert_equal(kloc, kbas)
    except Exception as e:
        Err = 'Different keys in local linfile.\n'
        Err+= '\tNew:{}\n'.format(kloc)
        Err+= '\tRef:{}\n'.format(kbas)
        Err+= '\tin linfile: {}.\n'.format(local_file2)
        newError(Err)
        #raise Exception(Err)

    # --- Compare 10 first frequencies and damping ratios in 'A' matrix
    if 'A' in fbas.keys(): 
        Abas = fbas['A']
        Aloc = floc['A']
        # Note: we could potentially reorder states like MBC does, but no need for freq/damping
        _, zeta_bas, _, freq_bas = eigA(Abas, nq=None, nq1=None, sort=True, fullEV=True)
        _, zeta_loc, _, freq_loc = eigA(Aloc, nq=None, nq1=None, sort=True, fullEV=True)

        freq_bas, zeta_bas = ApplyHighPass( freq_bas, zeta_bas )
        freq_loc, zeta_loc = ApplyHighPass( freq_loc, zeta_loc )

        if len(freq_bas)==0:
            # We use complex eigenvalues instead of frequencies/damping
            # If this fails often, we should discard this test.
            _, Lambda = eig(Abas, sort=False)
            v_bas = np.diag(Lambda)
            _, Lambda = eig(Aloc, sort=False)
            v_loc = np.diag(Lambda)

            if verbose:
                print(errPrefix+'val_ref:', v_bas[:7])
                print(errPrefix+'val_new:', v_loc[:7])
            try:
                np.testing.assert_allclose(v_bas[:10], v_loc[:10], rtol=rtol_f, atol=atol_f)
            except Exception as e:
                Err='Failed to compare A-matrix frequencies\n\tLinfile: {}.\n\tException: {}'.format(local_file2, indent(e.args[0]))
                newError(Err)
        else:
            #if verbose:
            print('\n'+errPrefix+':')
            print('              Frequency (Hz)                        Damping (%)')
            print('           ----------------------------         ----------------------------')
            print('            Ref              New                 Ref              New')

            #write frequencies to file
            try:
                file_freq_ref.write('\n'+errPrefix+':\n')
                file_freq_ref.write('            Freq (Hz)        Damp (%)\n')
                file_freq_new.write('\n'+errPrefix+':\n')
                file_freq_new.write('            Freq (Hz)        Damp (%)\n')
            except Exception:
                pass    # ignore all writing errors


            for j in range(min(10,max(len(freq_bas),len(freq_loc)))):
                if j<len(freq_bas):
                    str_freq_bas='{:14.9f}'.format(freq_bas[j])
                    str_zeta_bas='{:14.9f}'.format(zeta_bas[j])
                else:
                    str_freq_bas='              '
                    str_zeta_bas='              '
                if j<len(freq_loc):
                    str_freq_loc='{:14.9f}'.format(freq_loc[j])
                    str_zeta_loc='{:14.9f}'.format(zeta_loc[j])
                else:
                    str_freq_loc='              '
                    str_zeta_loc='              '
                print('        '+str_freq_bas+'   '+str_freq_loc+'      '+str_zeta_bas+'   '+str_zeta_loc)

                #write frequencies to file
                try:
                    file_freq_ref.write('        '+str_freq_bas+'   '+str_zeta_bas+'\n')
                    file_freq_new.write('        '+str_freq_loc+'   '+str_zeta_loc+'\n')
                except Exception:
                    pass


            try:
                np.testing.assert_allclose(freq_loc[:10], freq_bas[:10], rtol=rtol_f, atol=atol_f)
            except Exception as e:
                Err = 'Failed to compare A-matrix frequencies\n\tLinfile: {}.\n\tException: {}'.format(local_file2, indent(e.args[0]))
                newError(Err)


            if caseName=='Ideal_Beam_Free_Free_Linear':
                # The free-free case is a bit weird, smae frequencies but dampings are +/- a value
                zeta_loc=np.abs(zeta_loc)
                zeta_bas=np.abs(zeta_bas)


            try:
                # Note: damping ratios in [%]
                np.testing.assert_allclose(zeta_loc[:10]*100, zeta_bas[:10]*100, rtol=rtol_d, atol=atol_d)
            except Exception as e:
                Err = 'Failed to compare A-matrix damping ratios\n\tLinfile: {}.\n\tException: {}'.format(local_file2, indent(e.args[0]))
                newError(Err)

    # --- Compare individual matrices/vectors
    KEYS= ['A','B','C','D','dUdu','dUdy', 'x','y','u','xdot']
    for k,v in fbas.items():
        if k in KEYS and v is not None:
            if verbose:
                print(errPrefix+'key:', k)

            # Arrays
            Mloc = np.atleast_2d(floc[k])
            Mbas = np.atleast_2d(fbas[k])

            # --- Compare dimensions
            try:
                np.testing.assert_equal(Mloc.shape, Mbas.shape)
                dimEqual=True
            except Exception as e:
                Err = f'Different dimensions for variable `{k}`.\n'
                Err+= f'\tNew: {Mloc.shape}\n'
                Err+= f'\tRef: {Mbas.shape}\n'
                Err+= f'\tLinfile: {local_file2}.\n'
                newError(Err)
                dimEqual=False

            if not dimEqual:
                # We don't compare elements if shapes are different
                continue
            
            # Get matrix of which M elements are in/not in tolerance
            M_diff_in_tol = np.isclose(Mloc, Mbas, rtol=rtol, atol=atol)

            # Loop through indices of elements not in tolerance
            for (i,j) in zip(*np.where(M_diff_in_tol == False)):
                try:
                    # Redo compare to get error message
                    np.testing.assert_allclose(Mloc[i,j], Mbas[i,j], rtol=rtol, atol=atol)
                except Exception as e:
                    sElem = f'Element [{i+1},{j+1}], new: {Mloc[i,j]}, baseline: {Mbas[i,j]}'
                    if k in ['dXdx', 'A', 'dXdu', 'B']:
                        sElem += '\n\t\t  row: ' + fbas['x_info']['Description'][i]
                    if k in ['dYdx', 'C', 'dYdu', 'D']:
                        sElem += '\n\t\t  row: ' + fbas['y_info']['Description'][i]
                    if k in ['dUdu', 'dUdy']:
                        sElem += '\n\t\t  row: ' + fbas['u_info']['Description'][i]
                    if k in ['dXdx', 'A', 'dYdx', 'C']:
                        sElem += '\n\t\t  col: ' + fbas['x_info']['Description'][j]
                    if k in ['dXdu', 'B', 'dYdu', 'D', 'dUdu']:
                        sElem += '\n\t\t  col: ' + fbas['u_info']['Description'][j]
                    if k in ['dUdy']:
                        sElem += '\n\t\t  col: ' + fbas['y_info']['Description'][j]
                    Err = errPrefix + f'Failed to compare variable `{k}`, {sElem} \n\tLinfile: {local_file2}.\n\tException: {indent(e.args[0])}'
                    ElemErrors.append(Err)
    return Errors, ElemErrors


# Header for file containing all frequency information
def freqFileHdr():
    Err=[]
    file_freq_ref=[]
    file_freq_new=[]
    try:
        file_freq_ref=open(fNameFreqRef,"w")
    except Exception as e:
        Err = 'Could not open file '+fileNameFreqRef+' for writing'

    try:
        file_freq_new=open(fNameFreqNew,"w")
    except Exception as e:
        Err = 'Could not open file '+fileNameFreqNew+' for writing'

    return Err,file_freq_ref,file_freq_new

# close files with frequency info
def freqFileClose(file_freq_ref,file_freq_new):
    try:
        file_ref.close()
    except Exception:
        pass
    try:
        file_new.close()
    except Exception:
        pass


# Compare
Errors=[]

ErrorsLoc,ff1,ff2 = freqFileHdr()
Errors += ErrorsLoc

for i, f in enumerate(localOutFiles):
    ErrorsLoc, ElemErrorsLoc = compareLin(f,ff1,ff2)
    Errors += ErrorsLoc
    if len(ElemErrorsLoc)>0:
        Errors += ElemErrorsLoc[:3] # Just a couple of them

freqFileClose(ff1,ff2)


if len(Errors)>0:
    exitWithError('See errors below: \n'+'\n'.join(Errors))

# passing case
sys.exit(0)

