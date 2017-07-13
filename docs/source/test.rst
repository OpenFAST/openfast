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

Dependencies
------------
- Python 2.7/3+
- Numpy
- CMake and CTest (Optional in regression test)

Configuring the test suite
--------------------------
A portion of the test suite is linked to the OpenFAST repository through
``git submodule``. Specifically,

- `r-test <https://github.com/openfast/r-test>`__
- `pFUnit <http://pfunit.sourceforge.net>`__

Be sure to clone the repo with the ``--recursive`` flag or execute
``git submodule update --init --recursive`` after cloning.

The test suite can be built with `CMake <https://cmake.org/>`__ similar to
OpenFAST. The default CMake configuration is useful, but may need customization
for particular build environments. CMake variables can be configured in the `CMake
GUI <https://cmake.org/download/>`__ or through the command line interface with
the command ``ccmake``. If the entire OpenFAST package is to be built, the test
related CMake variables can be configured during the OpenFAST CMake
configuration. However, if only the test suite will be built, configure CMake
using the ``reg_tests`` project with ``ccmake path/to/reg_tests`` or selecting
``reg_tests`` as the source directory in the CMake GUI.

The test specific CMake variables are

- BUILD_TESTING
- CTEST_OPENFAST_EXECUTABLE
- CTEST_[MODULE]_EXECUTABLE

Look at the `Installing OpenFAST <install.html>`__ page for more details on
configuring the CMake targets.

While the unit tests must be built with CMake due to its external dependencies,
the regression test may be executed without building with CMake. See the
`regression test <test/regression_test.html>`__ section for more info.

Test specific documentation
---------------------------
.. toctree::
   :maxdepth: 2

   Unit Test <test/unit_test.rst>
   Regression Test <test/regression_test.rst>
