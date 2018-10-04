.. _regression_test:

Regression test
===============

The regression test executes a series of test cases which intend to fully describe
OpenFAST and its module's capabilities.

Jump to :ref:`regression_test_ctest`, :ref:`regression_test_example`, or :ref:`regression_test_windows`
for instructions on running the regression tests locally.

Each locally computed result is compared
to a static set of baseline results. To account for system, hardware, and compiler
differences, the regression test attempts to match the current machine and
compiler type to the appropriate solution set from these combinations:

- macOS with GNU compiler (default)
- Red Hat Enterprise Linux with Intel compiler
- Windows with Intel compiler

The compiler versions, specific math libraries, and hardware used to generate these baseline
solutions are documented in the
`r-test repository documentation <https://github.com/openFAST/r-test#baselines>`__. Currently,
the regression test supports only double precision solutions, so it is required
to build OpenFAST in double precision for testing. All baseline solutions are generated
with a double precision build.

The regression test system can be executed with CMake and CTest or manually with
an included Python driver. Both systems provide similar functionality with respect
to testing, but CTest integration provides access to multithreading, automation,
and test reporting via CDash. Both modes of execution require some configuration
as outlined below.

In both modes of execution a subdirectory is created in the build directory
called ``reg_tests`` where all of the input files for the test cases are copied
and all of the locally generated outputs are stored.

Ultimately, both CTest and the manual execution program call a series of Python
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

Dependencies
------------
    - Python 3+
    - Numpy
    - CMake and CTest (Optional)
    - matplotlib (Optional)

Manual driver configuration
---------------------------

The regression test can be executed manually with the included driver
``openfast/reg_tests/manualRegressionTest.py``. This program reads a case list file at
``openfast/reg_tests/r-test/glue-codes/openfast/CaseList.md``. Cases can be removed
or ignored with a ``#``. This driver program includes multiple optional flags
which can be obtained by executing with the help option:
``openfast/reg_tests/manualRegressionTest.py -h``

For the NREL 5MW turbine test cases, an external ServoDyn controller must be compiled and
included in the appropriate directory or all NREL 5MW cases will fail without starting.
More information is available in the documentation for the
`r-test repository <https://github.com/openfast/r-test#note---servodyn-external-controllers-for-5mw_baseline-cases>`__.

CTest configuration
-------------------

CTest is included with CMake and is mostly a set of preconfigured targets and
commands. To use the CTest driver for the regression test, CMake must be run with
one of two ``CMakeLists.txt``'s:

- openfast/CMakeList.txt
- openfast/reg_tests/CMakeLists.txt

CMake variables can be configured in the `CMake
GUI <https://cmake.org/download/>`__ or through the command line interface with
``ccmake``.

The regression test specific CMake variables are

- BUILD_TESTING
- CTEST_OPENFAST_EXECUTABLE
- CTEST_[MODULE]_EXECUTABLE
- CTEST_PLOT_ERRORS
- CTEST_REGRESSION_TOL

**IT IS IMPORTANT** to verify that NREL 5MW turbine external controllers are compiled
and placed in the correct location. More information is available in the documentation for the
`r-test repository <https://github.com/openfast/r-test#note---servodyn-external-controllers-for-5mw_baseline-cases>`__,
but be aware that these three DISCON controllers must exist

.. code-block:: bash

  openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoDyn/DISCON.dll
  openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoDyn/DISCON_ITIBarge.dll
  openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoDyn/DISCON_OC3Hywind.dll

This can be accomplished manually with the CMake projects included with the DISCON source codes
at ``openfast/reg_tests/r-test/glue-codes/openfast/5MW_Baseline/ServoDyn/``
or during CMake configuration by setting the ``CMAKE_INSTALL_PREFIX`` CMake variable.
If using this method, the install prefix variable should point to an existing and appropriate
location for CMake to place the compiled binaries. This is important because the NREL 5MW turbine external
controller CMake projects are preconfigured to install themselves in the appropriate
location in the build directory. Then, it is important to execute ``make install``
rather than simply ``make``. If ``CMAKE_INSTALL_PREFIX`` is not appropriately configured,
the install step may fail or openfast binaries may be placed in some inappropriate default location.

After CMake configuration, the automated regression test can be executed
by running either of the commands ``make test`` or ``ctest`` from the build
directory. If the entire OpenFAST package is to be built, CMake will configure
CTest to find the new binary at ``openfast/build/glue-codes/openfast/openfast``.
However, if the intention is to build only the test suite, the OpenFAST binary
should be specified in the CMake configuration under the ``CTEST_OPENFAST_EXECUTABLE``
flag. There is also a corresponding ``CTEST_[MODULE]_NAME`` flag for each module
included in the regression test.

.. _regression_test_ctest:

Running the regression test with CTest
--------------------------------------

When driven by CTest, the regression test can be executed by running various
forms of the command ``ctest`` from the build directory. The basic commands are

- ``ctest`` - Run the entire regression test
- ``ctest -N`` - Disable actual execution of tests; this is helpful in formulating a particular ctest command
- ``ctest -V`` - Run the entire regression test with verbose output
- ``ctest -R [TestName]`` - Run a test by name where TestName is a regex to search
- ``ctest -j [N]`` - Run all tests with N tests executing in parallel

Each regression test case contains a series of labels associating all of the
modules used. The labeling can be seen in the test instantiation in
``reg_tests/CTestList.cmake`` or with the command

- ``ctest --print-labels`` - Print all available test labels

Labels can be called directly with

- ``ctest -L [Label]``

These flags can be compounded making useful variations of ``ctest`` such as

- ``ctest -V -L aerodyn14`` - Runs all cases that use AeroDyn14 with verbose output
- ``ctest -j 16 -L aerodyn14`` - Runs all cases that use AeroDyn14 in 16 concurrent processes
- ``ctest -V -R 5MW_DLL_Potential_WTurb`` - Runs the case with name "5MW_DLL_Potential_WTurb"
- ``ctest -N -L beamdyn`` - Lists all tests with the "beamdyn" label
- ``ctest -N -R bd --print-labels`` - Lists the labels included in all tests matching the regex "bd"

The automated regression test writes new files only into the build directory. Specifically,
all locally generated solutions are located in the corresponding glue-code or module within
``openfast/build/reg_tests``. The baseline solutions contained in ``openfast/reg_tests/r-test``
are strictly read not modified by the automated process.

.. _regression_test_example:

Regression test example
-----------------------

- Build OpenFAST and the test suite

.. code-block:: bash

  git clone --recursive https://github.com/openfast/openfast.git
  # The default git branch is 'master'. If necessary, switch to your target branch:
  # git checkout dev
  mkdir build install && cd build
  # Configure CMake with openfast/CMakeLists.txt
  # - BUILD_TESTING
  # - CTEST_OPENFAST_EXECUTABLE
  # - CTEST_[MODULE]_EXECUTABLE
  cmake .. -DBUILD_TESTING=ON
  make install
  ctest

- Build only the test suite if an openfast binary already exists

.. code-block:: bash

  git clone --recursive https://github.com/openfast/openfast.git
  # The default git branch is 'master'. If necessary, switch to your target branch:
  # git checkout dev
  mkdir build install && cd build
  # Configure CMake with openfast/reg_tests/CMakeLists.txt
  # - BUILD_TESTING
  # - CTEST_OPENFAST_EXECUTABLE
  # - CTEST_[MODULE]_EXECUTABLE
  cmake ../reg_tests
  make install
  ctest

- :ref:`regression_test_windows`

Follow the link above for a detailed procedure. It is summarized below though
excluding the procedure to build OpenFAST itself.

.. code-block:: bash

  git clone --recursive https://github.com/openfast/openfast.git
  cd openfast

  ## Build the ServoDyn external controller libraries
  # Open the Visual Studio Solution (DISCON.sln) located in 'openfast\vs-build\DISCON'
  # Choose Release and x64 for the Solutions Configuration and Solutions Platform
  # Build Solution

  ## Execute the OpenFAST regression Tests
  # Open a command prompt which is configured for Python (like Anaconda)
  cd openfast\reg_tests
  python manualRegressionTest.py ..\build\bin\openfast_x64.exe Windows Intel 1e-5
