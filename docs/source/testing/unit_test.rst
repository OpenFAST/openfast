.. _unit_test:

Unit tests
==========
In a software package as dynamic and collaborative as OpenFAST, confidence in
multiple layers of code is best accomplished with a strong system of unit
tests. Through robust testing practices, the entire OpenFAST community can
understand the intention behind code blocks and debug or expand functionality
quicker and with more confidence and stability.

Unit testing in OpenFAST modules is accomplished through `test-drive <https://github.com/fortran-lang/test-drive>`__.
test-drive is compiled along with OpenFAST through CMake when the CMake variable ``BUILD_TESTING`` is
turned on (default off) and the CMake variable ``BUILD_UNIT_TESTING`` is on
(turned on by default when ``BUILD_TEST`` is on).

The BeamDyn and NWTC Library modules contain some sample unit tests and should
serve as a reference for future development and testing.

Dependencies
------------
The following packages are required for unit testing:

- CMake
- test-drive - Included in OpenFAST repo in unit_test/test-drive

Compiling
---------
Compiling the unit tests is handled with CMake similar to compiling OpenFAST
in general. After configuring CMake with ``BUILD_TESTING`` turned on, new
build targets are created for each module included in the unit test
framework named ``[module]_utest``. Then, ``make`` the target to test:

.. code-block:: bash

    cmake .. -DBUILD_TESTING=ON
    make beamdyn_utest

This creates a unit test executable at
``openfast/build/unit_tests/beamdyn_utest``.

Executing
---------
To execute a module's unit test, simply run the unit test binary. For example:

.. code-block:: bash

    >>>$ ./openfast/build/unit_tests/beamdyn_utest
    All tests PASSED

the pass or fail status is provided for each test as it's run. An error message is output when the test fails.
Failure cases display the following output:

.. code-block:: bash

    >>>$ ./unit_tests/beamdyn_utest
    # Testing: Crv
    Starting test_BD_CheckRotMat ... (1/6)
        ... test_BD_CheckRotMat [PASSED]
    Starting test_BD_ComputeIniNodalCrv ... (2/6)
        ... test_BD_ComputeIniNodalCrv [PASSED]
    Starting test_BD_CrvCompose ... (3/6)
        ... test_BD_CrvCompose [PASSED]
    Starting test_BD_CrvExtractCrv ... (4/6)
        ... test_BD_CrvExtractCrv [PASSED]
    Starting test_BD_CrvMatrixH ... (5/6)
    [Fatal] Uncaught error
    Code: 1 Message: A(1,1) simple rotation with known parameters: Pi on xaxis:
    Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO
    ERROR STOP 

    Error termination. Backtrace:
    #0  0xffff9f70d08b in ???
    #1  0xffff9f70ddb3 in ???
    #2  0xffff9f70f333 in ???

Adding unit tests
-----------------
Unit tests should be included for each new, *testable* code block (subroutine
or function). What is testable is the discretion of the developer, but an
element of the pull request review process will be evaluating test coverage.

New unit tests can be added to a ``tests`` directory alongside the ``src``
directory included in each module. For example, a module directory may be
structured as

::

  openfast/
    └── modules/
        └── sampledyn/
            ├── src/
            │   ├── SampleDyn.f90
            │   └── SampleDyn_Subs.f90
            └── tests/
                ├── sampledyn_utest.F90
                ├── test_SampleDyn_Feature1.F90
                ├── test_SampleDyn_Feature2.F90
                └── test_SampleDyn_Feature3.F90

Each unit test file must contain a module that exports a function which populates
a list of unit tests in accordance with the ``test-drive`` documentation. These modules
contain subroutines which take an ``error`` argument that is populated by the ``check`` 
subroutine provided by ``test-drive``. The ``sampledyn_utest.F90`` collects all of the
unit tests lists from the adjacent modules and runs them. These programs are compiled
via the ``unit_tests/CMakeLists.txt`` file so all relevant modules and programs are 
specified there. 

Refer to existing unit tests for the ``BeamDyn`` or ``NWTC Library`` unit tests for examples
of how to structure and build the unit test drivers. Also review the ``test-drive`` documentation at
`test-drive <https://github.com/fortran-lang/test-drive>`__.

Some useful topics to consider when developing and testing for OpenFAST are:

- `Test driven development <https://en.wikipedia.org/wiki/Test-driven_development#Test-driven_development_cycle>`__
- `Separation of concerns <https://en.wikipedia.org/wiki/Separation_of_concerns>`__
- `pFUnit usage <http://pfunit.sourceforge.net/page_Usage.html>`__
