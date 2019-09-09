.. _unit_test:

Unit test
=========

In a software package as dynamic and collaborative as OpenFAST, confidence in multiple
layers of code is best accomplished with a strong system of unit tests.
Through robust testing practices, the entire OpenFAST community can
understand the intention behind code blocks and debug or expand functionality
quicker and with more confidence and stability.

Unit testing in OpenFAST modules is accomplished through `pFUnit <http://pfunit.sourceforge.net>`__. 
This framework provides a Fortran abstraction to the popular `xUnit <https://en.wikipedia.org/wiki/XUnit>`__ 
structure. pFUnit is compiled along with OpenFAST through CMake when 
the CMake variable ``BUILD_TESTING`` is turned on.

The BeamDyn module has been unit tested and should serve as a reference for future 
development and testing.

Dependencies
------------
    - Python 3+
    - CMake and CTest
    - Numpy and matplotlib (Optional)
    - pFUnit (Included in unit test suite)

Compiling the unit tests
------------------------

Compiling the unit tests is handled with CMake similar to compiling OpenFAST in general.
After configuring CMake with ``BUILD_TESTING`` on, new make targets are created for each
module included in the unit test framework named ``[module]_utest``. Then, simply make the target to test

::
  
  cmake .. -DBUILD_TESTING=ON
  make beamdyn_utest
  
This creates a binary unit test executable at 
``openfast/build/unit_tests/[module]_utest``.
  

Executing the unit tests
------------------------

To execute a module's unit test, simply run the unit test binary. For example,
::
  
  >>>$ ./openfast/build/unit_tests/beamdyn_utest
  .............
  Time:         0.018 seconds
    
   OK
   (13 tests)

pFUnit will display a ``.`` for each unit test successfully completed
and a ``F`` for each failing test. If any tests do fail, the failure 
criteria will be displayed listing which particular value caused 
the failure. Failure cases display the following output

::
  
  >>>$ ./unit_tests/beamdyn_utest 
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

Unit tests should be included for each new, *testable* code block (subroutine or function).
What is testable is the discretion of the developer, but a portion 
of the pull request review process will be evaluating test coverage.

New unit tests can be added to a ``tests`` directory alongside the ``src``
directory included in each module. For example, the BeamDyn module directory is
structured as

::
  
  openfast/
  |-- modules/
    |-- beamdyn/
      |-- src/
        |-- BeamDyn.f90
        `-- BeamDyn_Subs.f90
      `-- tests/
        |-- test_BD_Subroutine1.F90
        |-- test_BD_Subroutine2.F90
        `-- test_BD_Subroutine3.F90
    
Each unit test must be contained in a unique file called ``test_[SUBROUTINE].F90`` where
``[SUBROUTINE]`` is the code block being tested. Finally, update the CMake configuration
for building a module's unit test executable with the appropriate list of test subroutines
in ``openfast/unit_tests/CMakeLists.txt`` using the following format

::
  
  set(testlist
     test_SUBROUTINE1
     test_SUBROUTINE2
     test_SUBROUTINE3
  )
  # it is important to keep the quotes around "${testlist}" in the call below
  build_utest("module_name" "${testlist}")
 
For reference, a template unit test file is included at ``openfast/unit_tests/test_SUBROUTINE.F90``.

Each unit test should fully test the target code block. If full test coverage
is not easily achievable, it may be an indication that refactoring would be beneficial.

Some useful topics to consider when developing and testing for OpenFAST are:

- `Test driven development <https://en.wikipedia.org/wiki/Test-driven_development#Test-driven_development_cycle>`__
- `Separation of concerns <https://en.wikipedia.org/wiki/Separation_of_concerns>`__
- `pFUnit usage <http://pfunit.sourceforge.net/page_Usage.html>`__
