Testing OpenFAST
================

OpenFAST testing is accomplished through the use of `CTest <https://cmake.org/Wiki/CMake/Testing_With_CTest>`__ and customized with a set of Python scripts.

All of the files corresponding to testing suite are contained in the ``reg_tests`` directory of the OpenFAST repository. The associated files include

- input files
- baseline solutions
- various Python programs used in the tests

Some required files are linked in the OpenFAST repository through a git submodule. Be sure to obtain all of the necessary files for testing by running ``git submodule update --init --recursive`` or update the required files with ``git submodule update``.

Dependencies
------------
- Python 2.7/3+
- Numpy
- CMake and CTest

Configuring the test suite
--------------------------
The test suite is built with `CMake <https://cmake.org/>`__ similar to OpenFAST. The default CMake configuration is useful, but may need customization for particular build environments. CMake variables can be configured in the `CMake GUI <https://cmake.org/download/>`__ or through the command line interface with the command ``ccmake``. If the entire OpenFAST package is to be built, the test related CMake variables can be configured during the OpenFAST CMake configuration. However, if only the test suite will be built, configure CMake using the ``reg_tests`` project with ``ccmake path/to/reg_tests`` or selecting ``reg_tests`` as the source directory in the CMake GUI.

The test specific CMake variables are

- BUILD_TESTING
- OPENFAST_EXECUTABLE
- [MODULE]_EXECUTABLE

Look at the `Installing OpenFAST <install.html>`__ page for more details on configuring the CMake targets.

Unit test
---------
Coming soon


Regression test
---------------
The regression test executes a series of test cases which fully describe OpenFAST and some associated submodule capabilities. Each locally computed result is compared to a static set of baseline results. To account for machine and compiler differences, CTest attempts to match the current machine and compiler type to the appropriate solution set from these combinations

- macOS with GNU compiler (default)
- Red Hat Enterprise Linux with Intel compiler
- Windows with Intel compiler

The comparison script reads the output files and computes a norm on each channel reported. If the maximum norm is greater than a predetermined tolerance, that particular test is reported as failed. The failure criteria is outlined in pseudocode below.

::

  for j in range(nChannels)
     norm_diff[j] = L2norm(localSolution[j]-baselineSolution[j])
     rms_baseline[j] = L2norm(baselineSolution[j])

  norm = norm_diff / rms_baseline

  if max(norm) < tolerance:
    success


Each regression test case contains a series of labels associating all of the modules used. The labeling can be seen in the test instantiation in ``reg_tests/CTestList.cmake``

Running the regression test
---------------------------
The test suite is driven by CTest and can be executed by running various forms of the command ``ctest`` from the build directory.

Run a test by name: ``ctest -R TestName``

Run all tests with a particular label: ``ctest -L Label``

Parellel test execution with N processes: ``-j N``

Verbose output: ``-V``

Extra verbose output: ``-VV``

Some common uses of ``ctest`` are:

- ``ctest -j 16``
- ``ctest -VV -L aerodyn14``


Regression test from scratch
--------------------------------------
- Build OpenFAST and the test suite

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  mkdir build && cd build
  # Configure CMake - BUILD_TESTING, OPENFAST_EXECUTABLE, [MODULE]_EXECUTABLE
  cmake ..
  make
  ctest


- Build only the test suite

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  mkdir build && cd build
  # Configure CMake - OPENFAST_EXECUTABLE, [MODULE]_EXECUTABLE
  cmake ../reg_tests
  ctest
