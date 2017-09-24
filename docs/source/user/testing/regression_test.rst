Regression test
===============

The regression test executes a series of test cases which fully describe
OpenFAST and its module's capabilities. Each locally computed result is compared
to a static set of baseline results. To account for machine and compiler
differences, the regression test attempts to match the current machine and
compiler type to the appropriate solution set from these combinations

- macOS with GNU compiler (default)
- Red Hat Enterprise Linux with Intel compiler
- Windows with Intel compiler

The automated regression test runs
`CTest <https://cmake.org/Wiki/CMake/Testing_With_CTest>`__ and can be executed
by running either of the commands ``make test`` or ``ctest`` from the build
directory. If the entire OpenFAST package is to be built, CMake will configure
CTest to find the new binary at ``openfast/build/glue-codes/fast/openfast``.
However, if the intention is to build only the test suite, the OpenFAST binary
should be specified in the CMake configuration under the ``CTEST_OPENFAST_EXECUTABLE``
flag. There is also a corresponding ``CTEST_[MODULE]_NAME`` flag for each module
included in the regression test.

The regression test can be executed manually with the included driver
``reg_tests/manualRegressionTest.py``.

In both modes of execution a subdirectory is created in the build directory
called ``reg_tests`` where all of the input files for the test cases are copied
and all of the locally generated outputs are stored.

Ultimately, both CTest and the manual execution script call a series of Python
scripts and libraries in ``reg_tests`` and ``reg_tests/lib``. One such script is
``lib/pass_fail.py`` which reads the output files and computes a norm on each 
channel reported. If the maximum norm is greater than a preset tolerance, that particular
test is reported as failed. The failure criteria is outlined in pseudocode below.

::
  
  difference = abs(testData-baselineData)
  for i in nChannels
     if channelRange < 1 {
        norm[i] = MaxNorm( difference[:,i] )
     } else {
        norm[i] = MaxNorm( difference[:,i] ) / channelRange
     }
     
  if max(norm) < tolerance:
    success


Manual driver configuration
---------------------------

The regression test can be executed manually with the included driver
``openfast/reg_tests/manualRegressionTest.py``. This program reads a case list file at 
``openfast/reg_tests/r-test/glue-codes/fast/CaseList.md``. Cases can be removed 
or ignored with a ``#``. This driver program includes multiple optional flags 
which can be obtained by executing with the help option:
``openfast/reg_tests/manualRegressionTest.py -h``

For the 5MW test cases, an external ServoDyn controller must be compiled and 
included in the appropriate directory or all 5MW cases will fail without starting.
Compiling these controllers is automated with CMake and make through
``openfast/reg_tests/r-test/glue-codes/fast/compileDISCON.py``. This should be
completed before running the regression test so that the compiled libraries are
included in the build directory when the files are initially copied. If the build
directory is populated without the external controllers, they can be compiled with 
``compileDISCON.py`` and copied manually into 
``openfast/build/reg_tests/glue-codes/fast/5MW_Baseline/ServoData``.
Running the regression test with CTest
--------------------------------------
When driven by CTest, the regression test can be executed by running various
forms of the command ``ctest`` from the build directory. The basic commands are

- ``ctest`` - Run the entire regression test
- ``ctest -V`` - Run the entire regression test with verbose output
- ``ctest -R [TestName]`` - Run a test by name
- ``ctest -j [N]`` - Run all tests with N tests executing in parallel

Each regression test case contains a series of labels associating all of the
modules used. The labeling can be seen in the test instantiation in
``reg_tests/CTestList.cmake`` and called directly with

- ``ctest -L [Label]``

These flags can be compounded making useful variations of ``ctest`` such as

- ``ctest -V -L aerodyn14`` - Runs all cases that use AeroDyn14 with verbose output
- ``ctest -j 16 -L aerodyn14`` - Runs all cases that use AeroDyn14 in 16 concurrent processes
- ``ctest -V -R 5MW_DLL_Potential_WTurb`` - Runs the case with name "5MW_DLL_Potential_WTurb"

Regression test from scratch
----------------------------

- Build OpenFAST and the test suite

::

  git clone --recursive https://github.com/openfast/openfast.git
  cd openfast/reg_tests/r-tests/openfast
  python compileDISCON.py
  cd ../../
  mkdir build && cd build
  # Configure CMake with openfast/CMakeLists.txt - BUILD_TESTING, CTEST_OPENFAST_EXECUTABLE, CTEST_[MODULE]_EXECUTABLE
  cmake ..
  make
  ctest


- Build only the test suite

::

  git clone --recursive https://github.com/openfast/openfast.git
  cd openfast/reg_tests/r-tests/openfast
  python compileDISCON.py
  cd ../../
  mkdir build && cd build
  # Configure CMake with openfast/reg_tests/CMakeLists.txt - CTEST_OPENFAST_EXECUTABLE, CTEST_[MODULE]_EXECUTABLE
  cmake ../reg_tests
  ctest

- `Windows with Visual Studio regression test <regression_test_windows.html>`__

Follow the link above for a detailed procedure. It is summarized below though
excluding the procedure to build OpenFAST itself.

::

  git clone --recursive https://github.com/openfast/openfast.git
  cd openfast

  ## Build the ServoDyn external controller libraries
  # Open the Visual Studio Solution (DISCON.sln) located in ``openfast\vs-build\DISCON``
  # Choose Release and x64 for the Solutions Configuration and Solutions Platform
  # Build Solution

  ## Execute the OpenFAST regression Tests
  # Open a command prompt which is configured for Python (like Anaconda)
  cd openfast\reg_tests
  python manualRegressionTest.py ..\build\bin\FAST_x64.exe Windows Intel
