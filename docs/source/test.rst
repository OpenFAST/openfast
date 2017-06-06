Testing OpenFAST
================

OpenFAST automated testing is accomplished through the use of `CTest <https://cmake.org/Wiki/CMake/Testing_With_CTest>`__ and customized with a set of Python scripts.

All of the files corresponding to automated testing are contained in the ``reg_tests``
directory of the OpenFAST repository. The associated files include

- input files
- "gold standard" outputs
- various Python programs used in the tests

Dependencies
------------
- Python 3.0+
- Numpy
- CMake and CTest

Regression test
---------------
The automated regression test executes a series of test cases which fully describe the OpenFAST capability. Each
locally computed result is compared to a static set of "gold standard" results. To account for machine
and compiler differences, three combinations of "gold standards" are included

- macOS with GNU compiler
- Red Hat Enterprise Linux with Intel compiler
- Windows with Intel compiler

CTest can automatically determine the appropriate solution set, but in case none match the default is macOS with GNU compiler.

The comparison script reads the OpenFAST binary output files (.outb) and computes a norm on each channel reported. If the maximum norm
is greater than a predetermined tolerance, that particular test is reported as failed. The failure criteria is outlined in pseudocode below.

::

  for j in range(nChannels)
     norm_diff[j] = L2norm(localSolution[j]-goldSolution[j])
     rms_gold[j] = L2norm(goldSolution[j])

  norm = norm_diff / rms_gold

  if max(norm) < tolerance:
    success

Configuring the automated test
------------------------------
A critical step in configuring the automated test is getting the input files
and "gold standards". These are brought into OpenFAST through the git submodule ``r-test``
and can be initialized with ``git submodule update --init --recursive`` or updated with
``git submodule update``.

If the test will be executed without building OpenFAST, set the ``OPENFAST_EXECUTABLE`` CMake
variable to the executable to test. CMake variables can be configured in the CMake
GUI or in the command line interface with the command ``ccmake ../reg_tests``.

Running the automated test
--------------------------
The automated regression test runs CTest and can be executed by running the command ``make test`` from the build directory. If
the entire OpenFAST package is to be built, CMake will configure CTest to find the new binary at
``openfast/build/glue-codes/fast/openfast``. However, if the intention is to build only the test suite, the OpenFAST binary
should be specified in the CMake configuration under the ``OPENFAST_EXECUTABLE`` flag.

Test procedure from scratch
---------------------------
- Building all of OpenFAST

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  mkdir build && cd build
  cmake ..
  make
  make test


- Building only the test

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  mkdir build && cd build
  # Configure CMake - OPENFAST_EXECUTABLE
  cmake ../reg_tests
  make test
