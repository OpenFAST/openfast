# OpenFAST/reg_tests

This directory contains the regression testing suite for OpenFAST. Its contents are list here and further described below.
- [r-test](https://github.com/openfast/r-test), a standalone repository containing the regression test data
- CMake/CTest configuration files
- Python scripts

The automated regression test runs CTest and can be executed by running the command `make test` from the build directory. If
the entire OpenFAST package is to be built, CMake will configure CTest to find the new binary at
`openfast/build/glue-codes/fast/openfast`. However, if the intention is to build only the test suite, the OpenFAST binary
should be specified in the CMake configuration under the `OPENFAST_EXECUTABLE` flag.

The regression test creates a subdirectory in the build directory calledA build directory is created at `ctest-build` which contains
the inputs to run the test cases and the locally generated outputs.

Dependencies to run the regression test suite are
- Python 3.0+
- Numpy
- CMake and CTest

## r-test
This repository serves as a container for regression test data for OpenFAST. The test cases are taken from the [FAST V8 CertTests](https://github.com/NWTC/FAST/tree/master/CertTest). The repository contains:
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

#### executeOpenfastCase.py
This program executes a single OpenFAST case.

Usage: `python3 executeOpenfastCase.py input_file openfast_executable`  
- `openfast_executable` is an optional argument pointing to the OpenFAST executable of choice.
- if `openfast_executable` is not given, an attempt will be made to find one in $PATH

Example: `python3 executeRegressionTestCase.py CaseDir/case01.fst`  
Example: `python3 executeRegressionTestCase.py CaseDir/case01.fst openfast`  
Example: `python3 executeRegressionTestCase.py CaseDir/case01.fst openfast/install/bin/openfast`  

#### executeRegressionTestCase.py
This program executes OpenFAST and a regression test for a single test case.
The test case must be one of the CertTest cases. The test data is contained in a git submodule,
r-test, which must be initialized prior to running. r-test can be initialized
with `git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python3 executeRegressionTestCase.py testname openfast_executable source_directory tolerance system_name compiler_id`  
Example: `python3 executeRegressionTestCase.py Test02 openfast path/to/openfast_repo 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`

#### pass_fail.py
This program determines whether a new solution has regressed from the "gold standard"
solution. It reads two OpenFAST binary output files (.outb), and computes the L2 norm
of the two solution files for each output channel. If the max norm is smaller than
the given tolerance, the test case passes.

Usage: `python3 pass_fail.py solution1 solution2 tolerance`  
Example: `python3 pass_fail.py output-local/Test01.outb gold-standard/Test01.outb 0.00000001`

#### fast_io.py
This program reads OpenFAST output files in binary or ascii format and returns the data in a Numpy array.
