Testing OpenFAST
================

OpenFAST automated testing is accomplished through the use of `CTest <https://cmake.org/Wiki/CMake/Testing_With_CTest>` and customized with a set of Python scripts.

All of the files corresponding to automated testing are contained in the ``reg_tests``
and ``ctest`` directories of the OpenFAST repository. ``ctest`` contains configuration
files and should generally be left untouched. ``reg_tests`` contains the input files,
"gold standard" outputs, and various Python programs used in the tests.

Dependencies
------------
- Python 3.0+
- Numpy
- CTest distributed through CMake

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
     norm_diff[j] = L2norm(dict1[:,j]-dict2[:,j])
     rms_gold[j] = L2norm(dict2[:,j])

  norm = norm_diff / rms_gold

  if max(norm) < tolerance:
    success

Configuring the automated test
------------------------------
The most critical step in configuring the automated test is getting the input files
and "gold standards". These are brought into OpenFAST through the git submodule ``r-test``
and can be initialized with ``git submodule update --init --recursive`` or updated with
``git submodule update``.

Running the automated test
--------------------------
The automated regression test uses CTest and can be executed in three ways:

- ``make test``

  Requires OpenFAST to have been built with ``make``. Specifically, it is
  assumed that a binary executable exists at ``openfast/build/glue-codes/fast/openfast``.
  This method creates a subdirectory in the CMake build directory called ``reg_tests``
  which contains the inputs to run the test cases and the locally generated outputs.


- ``executeFullRegressionTest.py``

  Runs CTest independently of CMake using a steering script at ``openfast/ctest/steer.cmake``.
  This method requires the user to specify an OpenFAST executable in the Python program call.
  A build directory is created at ``openfast/ctest-build`` which contains the inputs to run the
  test cases and the locally generated outputs.


- Calling CTest directly

  The CTest steering script can be executed directly, but this is not recommended as
  CMake and the Python program perform additional configuration checks before executing.

Test procedure from scratch
---------------------------
- ``make test`` method

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  mkdir build && cd build
  cmake -DENABLE_TESTS:BOOL=ON ..
  make
  make test

- ``executeFullRegressionTest.py`` method

::

  git clone https://github.com/openfast/openfast.git
  cd openfast/reg_tests
  git submodule update --init --recursive
  python3 executeFullRegressionTest.py path/to/openfast

- Calling CTest directly

::

  git clone https://github.com/openfast/openfast.git
  cd openfast
  git submodule update --init --recursive
  ctest -S ctest/steer.cmake -V -DEXECUTABLE=path/to/openfast
