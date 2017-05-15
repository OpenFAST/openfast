# OpenFAST/reg_tests

This directory contains the regression testing suite for OpenFAST. Its contents are list here and further described below.
- r-test, a standalone repository containing the regression test data
- CMake/CTest configuration files
- Python scripts

The automated regression test runs CTest and can be executed in two ways:
- `make test`
  - Requires OpenFAST to have been built with `make`. Specifically, it is assumed that a binary executable exists at `openfast/build/glue-codes/fast/openfast`. This method creates a subdirectory in the CMake build directory called `reg_tests` which contains the inputs to run the test cases and the locally generated outputs.


- `executeFullRegressionTest.py`
  - Runs CTest independently of CMake using a steering script at ``openfast/ctest/steer.cmake``. This method requires the user to specify an OpenFAST executable in the Python program call. A build directory is created at ``openfast/ctest-build`` which contains the inputs to run the test cases and the locally generated outputs.

Dependencies to run the regression test suite are
- Python 2.7
- Numpy
- CTest distributed through CMake

## r-test
This repository serves as a container for regression test data. The test cases are taken from the FAST V8 CertTests. The repository contains:
- input files for test execution
- outputs for various machine and compiler combinations which serve as "gold standard" solutions for the regression test suite

r-test is brought into OpenFAST as a git submodule and can be initialized with `git submodule update --init --recursive` or updated with `git submodule update`.

## CTest/CMake
The configuration files consist of
- CMakeLists.txt
- CTestList.cmake

#### CMakeLists.txt
This is a CMake file which configures the regression test in the CMake build directory. It should be left untouched unless advanced configuration of CMake or CTest is required.

#### CTestList.txt
This is the CTest configuration file which lists the test cases that run in the automated test. The test list can be modified as needed by commenting lines with a `#`, but the full regression test consists of all the tests listed in this file.

## Python Scripts
The included Python scripts are used to execute various parts of the automated regression test, so they should remain in their current location with their current name. Each script can be executed independently.

#### executeFullRegressionTest.py
This program executes the OpenFAST regression test suite through the use of
CTest and other custom scripts. The test data is contained in a git submodule,
r-test, which must be initialized prior to running. r-test can be initialized
with `git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python executeRegressionTestCase.py openfast_executable`  
Example: `python executeRegressionTestCase.py /path/to/openfast`

#### executeOpenfastCase.py
This program executes a single OpenFAST case.

Usage: `python executeOpenfastCase.py input_file openfast_executable`  
- `openfast_executable` is an optional argument pointing to the OpenFAST executable of choice.
- if `openfast_executable` is not given, an attempt will be made to find one in $PATH

Example: `python executeRegressionTestCase.py CaseDir/case01.fst`  
Example: `python executeRegressionTestCase.py CaseDir/case01.fst openfast`  
Example: `python executeRegressionTestCase.py CaseDir/case01.fst openfast/install/bin/openfast`  

#### executeRegressionTestCase.py
This program executes OpenFAST and a regression test for a single test case.
The test case must be one of the CertTest cases. The test data is contained in a git submodule,
r-test, which must be initialized prior to running. r-test can be initialized
with `git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python executeRegressionTestCase.py testname openfast_executable tolerance system_name compiler_id`  
Example: `python executeRegressionTestCase.py Test02 openfast 0.000001 [Darwin,RHEL,Windows] [Intel,GNU]`

#### pass_fail.py
This program determines whether a new solution has regressed from the "gold standard"
solution. It reads two OpenFAST binary output files (.outb), and computes the variance
of the two solution files for each output channel. If the max variance is smaller than
the given tolerance, the test case passes.

Usage: `python pass_fail.py solution1 solution2 tolerance`  
Example: `python pass_fail.py output-local/Test01.outb gold-standard/Test01.outb 0.00000001`

#### fast_io.py
This program reads OpenFAST output files in binary or ascii format and returns the data in a Numpy array.
