.. _testing:

Testing OpenFAST
================

The OpenFAST test suite consists of glue code and module level regression tests
and unit tests. The regression tests compare locally generated solutions to
a set of baseline solutions. The unit tests ensure that individual subroutines
are functioning as intended.

All of the necessary files corresponding to the regression tests are contained
in the ``reg_tests`` directory. The unit test framework is housed in
``unit_tests`` while the actual tests are contained in the directory
corresponding to the tested module.

Configuring the test suite
--------------------------
Portions of the test suite are linked to the OpenFAST repository through a
`git submodule`. Specifically,

- `r-test <https://github.com/openfast/r-test>`__
- `pFUnit <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`__

.. tip::

    Be sure to clone the repo with the ``--recursive`` flag or execute
    ``git submodule update --init --recursive`` after cloning.

The test suite can be configured with CMake similar to OpenFAST. The default
CMake configuration is suitable for most systems, but may need customization
for particular build environments. See the :ref:`understanding_cmake` section
for more details on configuring the CMake targets. While the unit tests must
be built with CMake due to its external dependencies, the regression test
may be executed without CMake.

Test specific documentation
---------------------------
.. toctree::
   :maxdepth: 1

   unit_test.rst
   regression_test.rst

Continuous integration
----------------------
A TravisCI configuration file is included with the OpenFAST source code at ``openfast/.travis.yml``.
The continuous integration infrastructure is still under development, but the
status for all branches and pull requests can be found on the
`TravisCI OpenFAST page <https://travis-ci.org/OpenFAST>`_.

For development and testing purposes, a version of the TravisCI test can be run
locally with Docker. The code snippet below outlines starting a TravisCI image
on Docker.

.. code-block:: bash

    # Running a travis ci image on docker locally

    # Run this on your local machine's command line
    BUILDID="build-1"
    INSTANCE="travisci/ci-garnet:packer-1512502276-986baf0"
    docker run --name $BUILDID -dit $INSTANCE /sbin/init
    docker exec -it $BUILDID bash -l

    # Now you're inside your docker image
    sudo apt-get update
    sudo apt-get install python3-pip
    sudo -E apt-get -yq --no-install-suggests --no-install-recommends install gfortran libblas-dev liblapack-dev
    git clone --depth=50 https://github.com/OpenFAST/openfast.git OpenFAST/openfast
    cd OpenFAST/openfast

    # Modify this line for the commit or pull request to build
    git fetch origin +refs/pull/203/merge:

    git checkout -qf FETCH_HEAD
    git submodule update --init --recursive

    export FC=/usr/bin/gfortran-7
    export DOUBLE_PRECISION=ON
    export TRAVIS_BUILD_INTEL=YES
    export TRAVIS_COMPILER=gcc
    export CC=gcc
    gcc --version
    pyenv shell 3.6.3

    source ~/.bashrc
    pip3 install numpy
    mkdir build && cd build
    cmake .. -DBUILD_TESTING=ON -DDOUBLE_PRECISION=$DOUBLE_PRECISION -DBUILD_SHARED_LIBS=ON
    make -j 8 install
