.. _regression_test:

Regression tests
================
The regression test executes a series of test cases which intend to fully
describe OpenFAST and its module's capabilities. Jump to one of the following
sections for instructions on running the regression
tests:

- :ref:`python_driver`
- :ref:`ctest_driver`
- :ref:`regression_test_example`
- :ref:`regression_test_windows`

Each locally computed result is compared to a static set of baseline
results. To account for system, hardware, and compiler
differences, the regression test attempts to match the current machine and
compiler type to the appropriate solution set from these combinations:

================== ============== ============================
 Operating System   Compiler       Hardware
================== ============== ============================
 macOS 10.15        GNU 10.2       2020 MacbookPro
 Ubuntu 20.04       Intel oneAPI   Docker
 Ubuntu 20.04       GNU 10.2       Docker
 Windows 10         Intel oneAPI   Dell Precision 3530
================== ============== ============================

The compiler versions, specific math libraries, and more info on hardware used
to generate the baseline solutions are documented in the
`r-test repository documentation <https://github.com/openFAST/r-test>`__. Currently,
the regression test supports only double precision builds.

The regression test system can be executed with CMake (through its included
test driver, CTest) or manually with
a custom Python driver. Both systems provide similar functionality with
respect to testing, but CTest integration provides access to multithreading,
automation, and test reporting via CDash. Both modes of execution require some
configuration as described in the following sections.

In both modes of execution a directory is created in the build directory
called ``reg_tests`` where all of the input files for the test cases are copied
and all of the locally generated outputs are stored. Ultimately, both CTest and
the manual execution program call a series of Python scripts and libraries in
``reg_tests`` and ``reg_tests/lib``. One such script is ``lib/pass_fail.py``
which reads the output files and computes a norm on each channel reported. If
the maximum norm is greater than the given tolerance, that particular test is
reported as failed. The failure criteria is outlined below.

.. code-block:: python

    difference = abs(testData - baselineData)
    for i in nChannels:
        if channelRange < 1:
            norm[i] = MaxNorm( difference[:,i] )
        else:
            norm[i] = MaxNorm( difference[:,i] ) / channelRange

    if max(norm) < tolerance:
        pass = True
    else:
        pass = False

Dependencies
------------
The following packages are required for regression testing:

- Python 3.7+
- Numpy
- CMake and CTest (Optional)
- Bokeh 2.4+ (Optional)

.. _python_driver:

Executing with Python driver
----------------------------
The regression test can be executed manually with the included driver at
``openfast/reg_tests/manualRegressionTest.py``. This program reads a case list
file at ``openfast/reg_tests/r-test/glue-codes/openfast/CaseList.md``. Cases
can be removed or ignored by starting that line with a ``#``. The driver
program includes multiple optional flags which can be obtained by
executing with the help option:

::

    >>>$ python manualRegressionTest.py -h
    usage: manualRegressionTest.py [-h] [-p [Plotting-Flag]] [-n [No-Execution]]
                                [-v [Verbose-Flag]] [-case [Case-Name]] [-module [Module-Name]]
                                Executable-Name Relative-Tolerance Absolute-Tolerance

    Executes OpenFAST or driver and a regression test for a single test case.

    positional arguments:
    Executable-Name       path to the executable
    Relative-Tolerance    Relative tolerance to allow the solution to deviate; expressed as order of magnitudes less than baseline.
    Absolute-Tolerance    Absolute tolerance to allow small values to pass; expressed as order of magnitudes less than baseline.

    optional arguments:
    -h, --help            show this help message and exit
    -p [Plotting-Flag], -plot [Plotting-Flag]
                            bool to include plots in failed cases
    -n [No-Execution], -no-exec [No-Execution]
                            bool to prevent execution of the test cases
    -v [Verbose-Flag], -verbose [Verbose-Flag]
                            bool to include verbose system output
    -case [Case-Name]     single case name to execute
    -module [Module-Name], -mod [Module-Name]
                            name of module to execute

.. note::

    For the NREL 5MW turbine test cases, an external ServoDyn controller must
    be compiled and included in the appropriate directory or all NREL 5MW
    cases will fail without starting. More information is available in the
    documentation for the `r-test repository <https://github.com/openfast/r-test#note---servodyn-external-controllers-for-5mw_baseline-cases>`__,
    but be aware that these three DISCON controllers must exist

    .. code-block:: bash

        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON.dll
        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON_ITIBarge.dll
        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON_OC3Hywind.dll

.. _ctest_driver:

Executing with CTest
--------------------
CTest is included with CMake and is primarily a set of preconfigured targets
and commands. To use the CTest driver for the regression test, execute CMake as
described in :ref:`installation`, but with this additional flag:
``-DBUILD_TESTING=ON``.

The regression test specific CMake variables are

::

    BUILD_TESTING
    CTEST_OPENFAST_EXECUTABLE
    CTEST_[MODULE]_EXECUTABLE where [MODULE] is the module name
    CTEST_PLOT_ERRORS
    CTEST_REGRESSION_TOL

Some additional resources that are required for the full regression test suite
are included in the CMake project. Specifically, external ServoDyn controllers
must be compiled for a given system and placed in a particular location. Thus,
be sure to execute the build command with the ``install`` target:

.. code-block:: bash

    # Configure CMake with testing enabled and accept the default
    # values for all other test-specific CMake variables
    cmake .. -DBUILD_TESTING=ON

    # Build and install
    make install

.. note::

    REMINDER: For the NREL 5MW turbine test cases, an external ServoDyn controller must
    be compiled and included in the appropriate directory or all NREL 5MW
    cases will fail without starting. More information is available in the
    documentation for the `r-test repository <https://github.com/openfast/r-test#note---servodyn-external-controllers-for-5mw_baseline-cases>`__,
    but be aware that these three DISCON controllers must exist

    .. code-block:: bash

        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON.dll
        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON_ITIBarge.dll
        openfast/build/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData/DISCON_OC3Hywind.dll

After CMake configuration and compiling, the automated regression test can be
executed by running either of the commands ``make test`` or ``ctest`` from the
``build`` directory. If the entire OpenFAST package is to be built, CMake will
configure CTest to find the new binary at
``openfast/build/glue-codes/openfast/openfast``. However, if the intention is
to build only the test suite, the OpenFAST binary should be specified in the
CMake configuration under the ``CTEST_OPENFAST_EXECUTABLE`` flag. There is
also a corresponding ``CTEST_[MODULE]_NAME`` flag for each module included in
the regression test.

When driven by CTest, the regression test can be executed by running various
forms of the command ``ctest`` from the build directory. The basic commands
are:

.. code-block:: bash

    # Run the entire regression test
    ctest

    # Disable actual execution of tests;
    # this is helpful in formulating a particular ctest command
    ctest -N

    # Run the entire regression test with verbose output
    ctest -V

    # Run tests by name where TestName is a regular expression (regex)
    ctest -R [TestName]

    # Run all tests with N tests executing in parallel
    ctest -j [N]

Each regression test case contains a series of labels associating all of the
modules used. The labeling can be seen in the test instantiation in
``reg_tests/CTestList.cmake`` or with the command:

.. code-block:: bash

    # Print all available test labels
    ctest --print-labels

The test cases corresponding to a particular label can be executed with this
command:

.. code-block:: bash

    # Filter the test cases corresponding to a particular label
    ctest -L [Label]

Flags can be compounded making useful variations such as

.. code-block:: bash

    # Run all cases that use AeroDyn14 with verbose output
    ctest -V -L aerodyn14

    # Run all cases that use AeroDyn14 in 16 concurrent processes
    ctest -j 16 -L aerodyn14

    # Run the case with name "5MW_DLL_Potential_WTurb" with verbose output
    ctest -V -R 5MW_DLL_Potential_WTurb

    # List all tests with the "beamdyn" label
    ctest -N -L beamdyn

    # List the labels included in all tests matching the regex "bd"
    ctest -N -R bd --print-labels


The automated regression test writes new files only into the build directory.
Specifically, all locally generated solutions are located in the corresponding
glue-code or module within ``openfast/build/reg_tests``. The baseline solutions
contained in ``openfast/reg_tests/r-test`` are strictly read and are not
modified by the automated process.

.. _regression_test_example:

Regression test examples
------------------------
The following examples illustrate methods of running the regression tests
on Unix-based systems. However, similar procedures can be used
on Windows with CMake and CTest. An alternate method of running the
regression tests on Windows is given in :ref:`reg_test_windows`.

Compile OpenFAST and execute with CTest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following example assumes the user is starting completely from scratch.
The commands below download the source code, configure the OpenFAST project
with CMake, compile all executables, and execute the full regression test
suite.

.. code-block:: bash

    # Download the source code from GitHub
    #    Note: The default branch is 'main'
    git clone --recursive https://github.com/openfast/openfast.git
    cd openfast

    # If necessary, switch to another target branch and update r-test
    git checkout dev
    git submodule update

    # Create the build and install directories and move into build
    mkdir build install && cd build

    # Configure CMake for testing
    # - BUILD_TESTING - turn ON
    # - CTEST_OPENFAST_EXECUTABLE - accept the default
    # - CTEST_[MODULE]_EXECUTABLE - accept the default
    cmake .. -DBUILD_TESTING=ON

    # Compile and install
    make install

    # Execute the full test suite with 4 concurrent processes
    ctest -j4

Configure with CMake and a given executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example assumes the user has a fully functional OpenFAST executable
available along with any necessary libraries, but does not have the source
code repository downloaded. This might be the case when executables are
distributed within an organization or downloaded from an
`OpenFAST Release <https://github.com/openfast/openfast/releases>`__.
Here, nothing will be compiled, but the test suite will be configured
with CMake for use with the CTest command.

.. code-block:: bash

    # Download the source code from GitHub
    #    Note: The default branch is 'main'
    git clone --recursive https://github.com/openfast/openfast.git
    cd openfast

    # If necessary, switch to another target branch and update r-test
    git checkout dev
    git submodule update

    # Create the build directory and move into it
    mkdir build && cd build

    # Configure CMake with openfast/reg_tests/CMakeLists.txt for testing
    # - BUILD_TESTING - turn ON
    # - CTEST_OPENFAST_EXECUTABLE - provide a path
    # - CTEST_[MODULE]_EXECUTABLE - provide a path
    cmake ../reg_tests \
        -DBUILD_TESTING=ON \
        -DCTEST_OPENFAST_EXECUTABLE=/home/user/Desktop/openfast_executable \
        -DCTEST_BEAMDYN_EXECUTABLE=/home/user/Desktop/beamdyn_driver

    # Install required files
    make install

    # Execute the full test suite with 4 concurrent processes
    ctest -j4

.. _example_python_driver:

Python driver with a given executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example assumes the user has a fully functional OpenFAST executable
available along with any necessary libraries, but does not have the source
code repository downloaded. This might be the case when executables are
distributed within an organization or downloaded from an
`OpenFAST Release <https://github.com/openfast/openfast/releases>`__.
Nothing will be compiled, but the test suite will be executed with the
included Python driver.

.. code-block:: bash

    # Download the source code from GitHub
    #    Note: The default branch is 'main'
    git clone --recursive https://github.com/openfast/openfast.git
    cd openfast

    # If necessary, switch to another target branch and update r-test
    git checkout dev
    git submodule update

    # Execute the Python driver
    cd reg_tests
    python manualRegressionTest.py -h
    # usage: manualRegressionTest.py [-h] [-p [Plotting-Flag]] [-n [No-Execution]]
    #                                [-v [Verbose-Flag]] [-case [Case-Name]] [-module [Module-Name]]
    #                                Executable-Name Relative-Tolerance Absolute-Tolerance
    # 
    # Executes OpenFAST or driver and a regression test for a single test case.
    # 
    # positional arguments:
    # Executable-Name       path to the executable
    # Relative-Tolerance    Relative tolerance to allow the solution to deviate; expressed as order of magnitudes less than baseline.
    # Absolute-Tolerance    Absolute tolerance to allow small values to pass; expressed as order of magnitudes less than baseline.
    # 
    # optional arguments:
    #   -h, --help            show this help message and exit
    #   -p [Plotting-Flag], -plot [Plotting-Flag]
    #                         bool to include plots in failed cases
    #   -n [No-Execution], -no-exec [No-Execution]
    #                         bool to prevent execution of the test cases
    #   -v [Verbose-Flag], -verbose [Verbose-Flag]
    #                         bool to include verbose system output
    #   -case [Case-Name]     single case name to execute
    #   -module [Module-Name], -mod [Module-Name]
    #                         name of module to execute

    python manualRegressionTest.py ..\build\bin\openfast_x64_Double.exe 2.0 1.9

.. _reg_test_windows:

Detailed example of running on Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :ref:`example_python_driver` example can be used for running the
regression tests on a Windows computer. However, a more detailed, step-by-step
description is given in :ref:`regression_test_windows`.

.. toctree::
   :maxdepth: 1
   :hidden:

   regression_test_windows.rst

.. _new_regression_test_case:

Adding test cases
-----------------
In all modes of execution, the regression tests are ultimately driven by a
series of Python scripts located in the ``openfast/reg_tests`` directory
with the naming scheme ``execute<Module>RegressionTest.py``.
The first step to adding a new regression test case is to verify that
a script exists for the target module. If it does not, an issue
should be opened in `OpenFAST Issues <https://github.com/openfast/openfast/issues>`_
to coordinate with the NREL team on creating this script.

The next step is to add the test case in the appropriate location in
the `r-test` submodule. The directory structure in r-test mirrors the
directory structure in OpenFAST, so module-level tests should be placed
in their respective module directories and glue-code tests go in
``r-test/glue-codes/openfast``. Note the naming scheme of files for
existing tests and adapt the new test case files accordingly. Specifically,
the main input file and output file names may be expected in a particular
convention by the Python scripts. Also, consider that any relative paths
within the input deck for the new test case must work within the r-test
directory structure.

Once the test directory exists, the test case must be registered with
the appropriate drivers. For OpenFAST glue-code tests, this happens both in
CMake and a standalone list of test cases. For CMake, edit the file
``openfast/reg_tests/CTestList.cmake``. The additional test should be
added in the section corresponding to the module or driver at the
bottom of that file. For the Python driver, the new test case must
be added to ``openfast/reg_tests/r-test/glue-codes/openfast/CaseList.md``.
At this point, the registration with CTest can be verified:

.. code-block:: bash

    # Move into the build directory
    cd openfast/build

    # Run CMake to take the new changes to the test list
    cmake .. -DBUILD_TESTING=ON  # If the BUILD_TESTING flag was previously enabled, this can be left off

    # List the registered tests, but don't run them
    ctest -N

For module regression tests, the only option for execution is with the
CMake driver, so follow the instructions above to edit ``CTestList.cmake``.

Finally, the new test cases in the r-test submodule must be added to the
r-test repository. To do this, open a new issue in `r-test Issues <https://github.com/openfast/r-test/issues>`_
requesting for support from the NREL team to commit your test.
