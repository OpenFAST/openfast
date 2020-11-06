.. _unit_test:

Unit test
=========
In a software package as dynamic and collaborative as OpenFAST, confidence in
multiple layers of code is best accomplished with a strong system of unit
tests. Through robust testing practices, the entire OpenFAST community can
understand the intention behind code blocks and debug or expand functionality
quicker and with more confidence and stability.

Unit testing in OpenFAST modules is accomplished through `pFUnit <https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git>`__.
This framework provides a Fortran abstraction to the popular
`xUnit <https://en.wikipedia.org/wiki/XUnit>`__ structure. pFUnit is compiled
along with OpenFAST through CMake when the CMake variable ``BUILD_TESTING`` is
turned on.

The BeamDyn and NWTC Library modules contain some sample unit tests and should
serve as a reference for future development and testing.

Dependencies
------------
The following packages are required for unit testing:

- Python 3.7+
- CMake
- pFUnit - Included in OpenFAST repo through a git-submodule

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
``openfast/build/unit_tests/beamdyn/beamdyn_utest``.

Executing
---------
To execute a module's unit test, simply run the unit test binary. For example:

.. code-block:: bash

    >>>$ ./openfast/build/unit_tests/beamdyn/beamdyn_utest
    .............
    Time:         0.018 seconds

     OK
     (14 tests)

pFUnit will display a ``.`` for each unit test successfully completed
and a ``F`` for each failing test. If any tests do fail, the failure
criteria will be displayed listing which particular value caused
the failure. Failure cases display the following output:

.. code-block:: bash

    >>>$ ./unit_tests/beamdyn/beamdyn_utest
    .....F.......
    Time:         0.008 seconds

    Failure
    in:
    test_BD_CrvMatrixH_suite.test_BD_CrvMatrixH
        Location:
    [test_BD_CrvMatrixH.F90:48]
    simple rotation with known parameters: Pi on xaxis expected +0.5000000 but found: +0.4554637;  difference: |+0.4453627E-01| > tolerance:+0.1000000E-13;  first difference at element [1, 1].

    FAILURES!!!
    Tests run: 13, Failures: 1, Errors: 0
    Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO
    ERROR STOP *** Encountered 1 or more failures/errors during testing. ***

    Error termination. Backtrace:
    #0  0x1073b958c
    #1  0x1073ba295
    #2  0x1073bb1b6
    #3  0x106ecdd4f
    #4  0x1063fabee
    #5  0x10706691e

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
                ├── test_SampleDyn_Subroutine1.F90
                ├── test_SampleDyn_Subroutine2.F90
                └── test_SampleDyn_Subroutine3.F90

Each unit test must be contained in a unique file called
``test_[SUBROUTINE].F90`` where ``[SUBROUTINE]`` is the code block being
tested. The new files should contain a Fortran `module` which itself
contains a Fortran `subroutine` for each specific test case. Generally,
multiple tests will be required to fully test one subroutine.

Finally, update the CMake configuration for building a module's unit
test executable by copying an existing unit test CMake configuration
into a new module directory:

.. code-block:: bash

    cp -r openfast/unit_tests/beamdyn openfast/unit_tests/[module]

Then, modify the new ``CMakeLists.txt`` with the appropriate list of test
subroutines and module name variables.

For reference, a template unit test file is included at
``openfast/unit_tests/test_SUBROUTINE.F90``. Each unit test should fully test
the target code block. If full test coverage is not easily achievable, it may
be an indication that refactoring would be beneficial.

Some useful topics to consider when developing and testing for OpenFAST are:

- `Test driven development <https://en.wikipedia.org/wiki/Test-driven_development#Test-driven_development_cycle>`__
- `Separation of concerns <https://en.wikipedia.org/wiki/Separation_of_concerns>`__
- `pFUnit usage <http://pfunit.sourceforge.net/page_Usage.html>`__
