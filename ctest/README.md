# OpenFAST/reg_tests

This directory contains configuration files for the regression testing suite for OpenFAST. Its contents are list here and further described below. These files should be left untouched unless advanced configuration of CMake or CTest is required.
- CTestConfig.cmake
- CTestTestfile.cmake
- steer.cmake

## CTestConfig.cmake
Currently unused, but reserved for future use. It is included here to silence a CMake warning.

## CTestTestfile.cmake
This configuration file contains the list of tests that will be executed by CTest. Each line in the file corresponds to an individual case of the regression test. More tests can be added by simply adding additional calls to `add_test` and ensuring the corresponding "gold standard" files exist.

## steer.cmake
The steering file is the driver for CTest when it is not initiated through CMake. This script configures the test tolerances and the directory structure for the test.
