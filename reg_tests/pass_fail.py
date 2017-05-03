import sys, os
import numpy as np
from fast_io import load_output

def exitWithError(error):
    print error
    sys.exit(1)

# validate input arguments
nArgsExpected = 4
if len(sys.argv) < nArgsExpected:
    exitWithError("Error: {} arguments given, expected {}\n".format(len(sys.argv), nArgsExpected) +
        "Usage: {} solution1 solution2 tolerance".format(sys.argv[0]))

solutionFile1 = sys.argv[1]
solutionFile2 = sys.argv[2]
solutionTolerance = sys.argv[3]

if not os.path.isfile(solutionFile1):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile1))

if not os.path.isfile(solutionFile2):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile2))

try:
    solutionTolerance = float(solutionTolerance)
except ValueError:
    exitWithError("Error: invalid tolerance given, {}".format(solutionTolerance))

# parse the FAST solution files
try:
    dict1, info1 = load_output(solutionFile1)
    dict2, info2 = load_output(solutionFile2)
except Exception as e:
    exitWithError("Error: {}".format(e))

# calculate the difference in solutions
nColumns = np.size(dict1,1)
variance = np.ones(nColumns)
for j in range(nColumns):
    variance[j] = (dict1[:,j]-dict2[:,j]).var()

if max(variance) < solutionTolerance:
    print "ok"
    sys.exit(0)
else:
    print info1['attribute_names']
    print "var = ", variance
    sys.exit(1)
