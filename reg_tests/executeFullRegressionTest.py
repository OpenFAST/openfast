"""
    This program executes the openfast regression test suite through the use of
    CTest and other custom scripts. The test data is contained in a git submodule,
    r-test, which must be initialized prior to running. r-test can be initialized
    with `git submodule update --init --recursive` or updated with `git submodule update`.

    Required dependencies are:
    - Python 2.7
    - CTest

    Usage: python executeRegressionTestCase.py openfast_executable
    Example: python executeRegressionTestCase.py /path/to/openfast
"""

import os
import sys
import shutil
import subprocess

def exitWithError(error):
    print error
    sys.exit(1)

if len(sys.argv) != 2:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTestCase.py openfast_executable")

os.chdir("..")

executable = sys.argv[1]
steerScript = os.path.join("ctest", "steer.cmake")

# verify that steer script exists
if not os.path.isfile(steerScript):
    exitWithError("The CTest steering script is not found where it was expected: {}".format(steerScript))

ctest_command = " ".join(["ctest", "-S", steerScript, "-V", "-DEXECUTABLE="+executable])
ctest_return_code = subprocess.call(ctest_command, shell=True)

print "CTest finished with return code: {}".format(ctest_return_code)
