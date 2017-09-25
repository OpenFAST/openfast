Unit test
=========

In a software package as dynamic and collaborative as OpenFAST, confidence in multiple
layers of code can only be accomplished with a strong system of unit tests.
Through robust testing practices, the entire OpenFAST community should be able to
understand the intention behind code blocks and help debug or expand functionality
quicker and with more confidence and stability.

Unit testing in OpenFAST modules is accomplished through `pFUnit <http://pfunit.sourceforge.net>`__. This framework provides a Fortran
abstraction to the popular `xUnit <https://en.wikipedia.org/wiki/XUnit>`__ structure.
This package is compiled along with OpenFAST through CMake when 
the CMake variable ``BUILD_TESTING`` is turned on.

Adding unit tests
-----------------

Unit tests should be included for each new testable code block (subroutine or function).
What is testable is the discretion of the developer, but a portion 
of the pull request review process will be evaluating test coverage.

New unit tests can be added to a ``tests`` directory alongside the ``src``
directory included in each module. For example, BeamDyn unit tests are at 
``openfast/modules-local/beamdyn/tests``. Each unit test subroutine
must be contained in a unique file called ``test_[SUBROUTINE].F90`` where
``[SUBROUTINE]`` is the code block being tested. Finally, add ``test_[SUBROUTINE]``
to the corresponding module section in the unit test CMake file at 
``openfast/unit_tests/CMakeLists.txt``.

Ideally, each unit test will fully test the target code block. If full test coverage
is not easily achievable, it may be an indication that refactoring is necessary.
Some useful topics to consider when developing and testing for OpenFAST are:

- `Test driven development <>`__
- `Separation of concerns <>`__
- `Functional programming <>`__
- `Object oriented program <>`__
- `xUnit <>`__
- `pFUnit <>`__


Executing the unit tests
------------------------

Upon successfully compiling OpenFAST, pFUnit, and the unit tests
themselves, a unit test binary is created at ``openfast/build/unit_tests/[module]_utest``.
To execute a module's unit test, simple run the unit test binary:

``./openfast/build/unit_tests/beamdyn_utest``

pFUnit will display a ``.`` for each unit test successfully completed
and a ``F`` for each failing test. If any tests do fail, the failure 
criteria will be displayed listing which particular value caused 
the failure.

Passing tests:

::
  
  >>>$ ./unit_tests/beamdyn_utest 
  .............
  Time:         0.018 seconds
    
   OK
   (13 tests)


Failing tests:

::
  
   >>>$ ./unit_tests/beamdyn_utest 
   ....F........
   Time:         0.018 seconds
     
    OK
    (13 tests)


