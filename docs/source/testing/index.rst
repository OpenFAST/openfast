.. _testing:

Testing OpenFAST
================

OpenFAST is a complex software with many moving parts. In order to
maintain stability in new and existing code, a test suite is included
directly in the source code. Two primary levels of tests exist:
regression tests at the highest level and unit tests at the lowest
level. The regression tests compare locally generated results with
stored "baseline" results. These tests give an indication of whether the
full-system or sub-system response has changed. The unit tests focus on
a single subroutine or code block. These tests need not be physically
realistic and focus on the mathematics and exersizing of an algorithm.
The objective of the included tests is to quickly catch bugs or unexpected
changes in results. Additionally, the tests can help programmers
design their module and subroutine interfaces in a sustainable and
maintainable manner.

All of the necessary files corresponding to the regression tests are contained
in the ``reg_tests`` directory. The unit test framework is housed in
``unit_tests`` while the actual tests are contained in the directory
corresponding to the tested module.

The OpenFAST GitHub repository uses `GitHub Actions <https://github.com/openfast/openfast/actions>`_
to automatically execute the test suite for new commits and pull
requests. This cloud computing resource is available to all
GitHub users and is highly recommended as part of the development
workflow. After enabling GitHub Actions in an OpenFAST repository, simply
pushing new commits will trigger the tests.

Test specific documentation
---------------------------
.. toctree::
   :maxdepth: 1

   unit_test.rst
   regression_test.rst

Obtaining and configuring the test suite
----------------------------------------
Portions of the test suite are linked to the OpenFAST repository through a
`git submodule`. Specifically, the following two repositories are included:

- `r-test <https://github.com/openfast/r-test>`__
- `pFUnit <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`__

.. tip::

    Be sure to clone the repo with the ``--recursive`` flag or execute
    ``git submodule update --init --recursive`` after cloning.

The test suite is configured with CMake similar to the general OpenFAST
build process with an additional CMake flag:

.. code-block:: bash

    # BUILD_TESTING     - Build the testing tree (Default: OFF)
    cmake .. -DBUILD_TESTING:BOOL=ON

Aside from this flag, the default CMake configuration is suitable for most systems.
See the :ref:`understanding_cmake` section for more details on configuring
the CMake targets. While the unit tests must be built with CMake due to its external
dependencies, the regression test may be executed without CMake, as described in
:ref:`python_driver`.
