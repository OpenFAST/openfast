# openfast/reg_tests

This directory contains the regression test suite for OpenFAST and its modules. Its contents are listed here and further described below.
- [r-test](https://github.com/openfast/r-test), a standalone repository containing the regression test data
- CMake/CTest configuration files
- Module specific regression test execution scripts
- A `lib` subdirectory with lower level python scripts

Dependencies required to run the regression test suite are
- Python 3.7+
- Numpy
- CMake and CTest
- Bokeh 2.4+ (optional)

## Execution
The automated regression test runs CTest and can be executed by running either of the commands `make test` or `ctest` from the build directory. If the entire OpenFAST package is to be built, CMake will configure CTest to find the new binary at `openfast/build/glue-codes/openfast/openfast`. However, if the intention is to build only the test suite, the OpenFAST binary should be specified in the CMake configuration under the `CTEST_OPENFAST_EXECUTABLE` flag. There is also a corresponding `CTEST_[MODULE]_NAME` flag for each module that is included in the regression test.

The regression test can be executed manually with the included driver `manualRegressionTest.py`. Run `manualRegressionTest.py -h` for usage.

In both modes of execution a subdirectory is created in the build directory called `reg_tests` where all of the input files for the test cases are copied and all of the locally generated outputs are stored.

## r-test
This repository serves as a container for regression test data for system level and module level testing of OpenFAST. The repository contains:
- input files for test case execution
- baseline solutions for various machine and compiler combinations
- turbine specific inputs

The baseline solutions serve as "gold standards" for the regression test suite and are updated periodically as OpenFAST and its modules are improved.

r-test is brought into OpenFAST as a git submodule and should be initialized after cloning with `git submodule update --init --recursive` or updated with `git submodule update`.

## CTest/CMake
The configuration files consist of
- CMakeLists.txt
- CTestList.cmake

#### CMakeLists.txt
This is a CMake file which configures the regression test in the CMake build directory. It should be left untouched unless advanced configuration of CMake or CTest is required.

#### CTestList.txt
This is the CTest configuration file which lists the test cases that run in the automated test. The test list can be modified as needed by commenting lines with a `#`, but the full regression test consists of all the tests listed in this file.

## Python Scripts
The included Python scripts are used to execute various parts of the automated regression test, so they should remain in their current location with their current name. Each script can be executed independently. The syntax and options for using the scripts can be found by running each with the `-h` flag.

#### executeOpenfastRegressionCase.py
This program executes OpenFAST and a regression test for a single test case.
The test data is contained in a git submodule, r-test, which must be initialized
prior to running. See the r-test README or OpenFAST documentation for more info.

Get usage with: `executeOpenfastRegressionCase.py -h`

#### executeBeamdynRegressionCase.py
This program executes BeamDyn and a regression test for a single test case.
The test data is contained in a git submodule, r-test, which must be initialized
prior to running. See the r-test README or OpenFAST documentation for more info.

Get usage with: `executeBeamdynRegressionCase.py -h`

#### manualRegressionTest.py
This program executes OpenFAST on all of the CertTest cases. It mimics the
regression test execution through CMake/CTest. All generated data goes into
`openfast/build/reg_tests`.

Get usage with: `manualRegressionTest.py -h`

#### lib/errorPlotting.py
This library provides tools for plotting the output channels over time of a 
given solution attribute for two OpenFAST solutions, with the second solution
assumed to be the baseline for comparison. There are functions for solution
file I/O, plot creation, and html creation for navigating the plots.

#### lib/fast_io.py
This program reads OpenFAST structured output files in binary or ascii format
and returns the data in a Numpy array.
  
#### lib/openfastDrivers.py
This library provides tools for executing cases with drivers contained in the
OpenFAST framework. Any new drivers should have a corresponding public driver
function called `def run[NewDriver]Case()` in this library.

#### lib/pass_fail.py
This library provides tools for comparing a test solution to a baseline solution
for any structured output file generated within the OpenFAST framework.

#### lib/rtestlib.py
This library contains utility functions for the custom python programs making
up the regression test system.
