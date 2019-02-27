.. _testing:

Testing OpenFAST
================

The OpenFAST test suite consists of system and module level regression tests and
unit tests. The regression test compares locally generated solutions to a set of
baseline solutions. The unit tests ensure that individual subroutines are
functioning as intended.

All of the necessary files corresponding to the regression test are contained in
the ``reg_tests`` directory. The unit test framework is housed in ``unit_tests``
while the actual tests are contained in the directory corresponding to the tested
module.

Configuring the test suite
--------------------------
Portions of the test suite are linked to the OpenFAST repository through
``git submodule``. Specifically,

- `r-test <https://github.com/openfast/r-test>`__
- `pFUnit <http://github.com/openfast/pfunit>`__

Be sure to clone the repo with the ``--recursive`` flag or execute
``git submodule update --init --recursive`` after cloning.

The test suite can be built with `CMake <https://cmake.org/>`__ similar to
OpenFAST. The default CMake configuration is useful, but may need customization
for particular build environments. See the installation documentation at :numref:`installation`
for more details on configuring the CMake targets.

While the unit tests must be built with CMake due to its external dependencies,
the regression test may be executed without building with CMake. :numref:`unit_test` and :numref:`regression_test`
have more information on unit testing and regression testing, respectively.

Test specific documentation
---------------------------
.. toctree::
   :maxdepth: 1

   unit_test.rst
   regression_test.rst
   regression_test_windows.rst
