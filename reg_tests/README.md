# OpenFAST/reg_tests

This directory contains the regression testing suite for OpenFAST. Its contents are list here and further described below:
- r-test: a standalone repository containing the input files necessary to run the regression test cases and "gold standard" files to determine pass/fail in the regression test
- CMake/CMake configuration files
- Python scripts

The automated regression test runs CTest and can be executed in two ways:
- `make test`
  - Requires OpenFAST to have been built with `make`. Specifically, it is assumed that an `openfast` binary executable exists at `openfast/build/glue-codes/fast/openfast`. This method creates a subdirectory in the CMake build directory called `reg_tests` which contains the inputs to run the test cases and the local outputs.


- `executeFullRegressionTest.py`
  - Runs CTest independently of CMake using a steering script at `openfast/ctest/steer.cmake`. A local build directory is created at `openfast/ctest-build` which contains the inputs to run the test cases and the local outputs.

## r-test
This repository serves as a container for regression test data for OpenFAST. The test cases are taken from the FAST V8 CertTests. The repository contains:
- input files for test execution
- directories with outputs for various machine and compiler combinations
The output files serve as "gold standard" solutions for the regression test suite.

r-test is brought into OpenFAST as a git submodule. r-test can be initialized with `git submodule update --init --recursive` or updated with `git submodule update`.

## CTest/CMake
The configuration files consist of
- CMakeLists.txt
- CTestList.cmake

#### CMakeLists.txt
This is CMake file which configures the regression test in the CMake build directory. It should be left untouched unless advanced configuration of CMake or CTest is required.

#### CTestList.txt
This is the CTest configuration file which lists the test cases that run in the automated test. The test list can be modified as needed by commenting lines with a `#`, but the full regression test consists of all the tests listed in this file.

## Python Scripts
The included Pythons scripts are used to execute various parts of the automated regression test, so they should remain in their current location with their current name. Each script can be executed independently.

#### executeFullRegressionTest.py
This program executes the openfast regression test suite through the use of
CTest and other custom scripts. The test data is contained in a git submodule,
r-test, which must be initialized prior to running. r-test can be initialized
with `git submodule update --init --recursive` or updated with `git submodule update`.

Required dependencies are:
- Python 2.7
- CTest

Usage: `python executeRegressionTestCase.py testname openfast_executable tolerance system_name compiler_id`
Example: `python executeRegressionTestCase.py Test02 openfast 0.000001 [Darwin,RHEL,Windows] [Intel,GNU]`

#### executeOpenfastCase.py
This program executes a single OpenFAST case.

Usage: `python executeOpenfastCase.py input_file openfast_executable`  
- openfast_executable is an optional argument pointing to the openfast executable of choice.
- if openfast_executable is not given, an attempt will be made to find one in $PATH

Example: `python executeRegressionTestCase.py .../CaseDir/case01.fst`  
Example: `python executeRegressionTestCase.py .../CaseDir/case01.fst openfast`  
Example: `python executeRegressionTestCase.py .../CaseDir/case01.fst openfast/install/bin/openfast`  

#### executeRegressionTestCase.py
This program executes OpenFAST and a regression test for a single test case.
The test case must be one of the CertTest cases. The test data is contained in a git submodule,
r-test, which must be initialized prior to running. r-test can be initialized
with `git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python executeRegressionTestCase.py testname openfast_executable tolerance system_name compiler_id`  
Example: `python executeRegressionTestCase.py Test02 openfast 0.000001 Darwin Intel`  
