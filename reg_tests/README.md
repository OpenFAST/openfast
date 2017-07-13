# openfast/reg_tests

This directory contains the regression test suite for OpenFAST and its modules. Its contents are listed here and further described below.
- [r-test](https://github.com/openfast/r-test), a standalone repository containing the regression test data
- CMake/CTest configuration files
- Module specific regression test execution scripts
- A `lib` subdirectory with lower level python scripts

Dependencies required to run the regression test suite are
- Python 2.7/3+
- Numpy
- CMake and CTest

## Execution
The automated regression test runs CTest and can be executed by running either of the commands `make test` or `ctest` from the build directory. If the entire OpenFAST package is to be built, CMake will configure CTest to find the new binary at `openfast/build/glue-codes/fast/openfast`. However, if the intention is to build only the test suite, the OpenFAST binary should be specified in the CMake configuration under the `CTEST_OPENFAST_EXECUTABLE` flag. There is also a corresponding `CTEST_[MODULE]_NAME` flag for each module that is included in the regression test.

The regression test can be executed manually with the included driver `manualRegressionTest.py`.

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
The included Python scripts are used to execute various parts of the automated regression test, so they should remain in their current location with their current name. Each script can be executed independently.

#### executeOpenfastRegressionCase.py
This program executes OpenFAST and a regression test for a single test case.
The test data is contained in a git submodule, r-test, which must be initialized
prior to running. r-test can be initialized with
`git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python executeOpenfastRegressionCase.py testname openfast_executable source_directory build_directory tolerance system_name compiler_id`  
Example: `python executeOpenfastRegressionCase.py Test02 openfast path/to/openfast_repo path/to/openfast_repo/build 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`

#### executeBeamdynRegressionCase.py
This program executes BeamDyn and a regression test for a single test case.
The test data is contained in a git submodule, r-test, which must be initialized
prior to running. r-test can be initialized with
`git submodule update --init --recursive` or updated with `git submodule update`.

Usage: `python executeBeamdynRegressionCase.py testname beamdyn_driver source_directory build_directory tolerance system_name compiler_id`  
Example: `python executeBeamdynRegressionCase.py Test02 beamdyn_driver path/to/openfast_repo path/to/openfast_repo/build 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`

#### manualRegressionTest.py
This program executes OpenFAST on all of the CertTest cases. It mimics the
regression test execution through CMake/CTest. All generated data goes into
`openfast/build/reg_tests`.

Usage: `python manualRegressionTest.py path/to/openfast_executable [Darwin,Linux,Windows] [Intel,GNU]`

#### lib/executeOpenfastCase.py
This program executes a single OpenFAST case.

Usage: `python executeOpenfastCase.py input_file openfast_executable`
- `openfast_executable` is an optional argument pointing to the OpenFAST executable of choice.
- if `openfast_executable` is not given, an attempt will be made to find one in $PATH

Example: `python executeOpenfastCase.py CaseDir/case01.fst`  
Example: `python executeOpenfastCase.py CaseDir/case01.fst openfast`  
Example: `python executeOpenfastCase.py CaseDir/case01.fst openfast/install/bin/openfast`

#### lib/executeBeamdynCase.py
This program executes a single BeamDyn case.

Usage: `python executeBeamdynCase.py input_file beamdyn_executable`
- `beamdyn_executable` is an optional argument pointing to the BeamDyn executable of choice.
- if `beamdyn_executable` is not given, an attempt will be made to find one in $PATH

Example: `python executeBeamdynCase.py CaseDir/case01.fst`
Example: `python executeBeamdynCase.py CaseDir/case01.fst beamdyn`
Example: `python executeBeamdynCase.py CaseDir/case01.fst openfast/install/bin/beamdyn`

#### lib/pass_fail.py
This program determines whether a new solution has regressed from the "gold standard"
solution. It reads two OpenFAST binary output files (.outb), and computes the L2 norm
of the two solution files for each output channel. If the max norm is smaller than
the given tolerance, the test case passes.

Usage: `python pass_fail.py solution1 solution2 tolerance`  
Example: `python pass_fail.py output-local/Test01.outb gold-standard/Test01.outb 0.00000001`

#### lib/fast_io.py
This program reads OpenFAST output files in binary or ascii format and returns the data in a Numpy array.
