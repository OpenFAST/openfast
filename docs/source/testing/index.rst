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

Continuous Integration
----------------------
A TravisCI configuration file is included with the OpenFAST source code at ``openfast/.travis.yml``.
The continuous integration infrastructure is still under development, but the status for all branches
and pull requests can be found on the `TravisCI OpenFAST page <https://travis-ci.org/OpenFAST>`_.

Note that if you use the included TravisCI configuration, you will need to add your own Intel compiler
license serial number to your TravisCI project. Otherwise, simply remove the ``ifort`` line from the
environment list.

For development and testing purposes, a version of the TravisCI test can be run locally with Docker.
Below is a guide which should be modified depending on the particular build being replicated

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

    export INTEL_SERIAL_NUMBER=VFGH-65656FTH
    export FC=ifort
    export DOUBLE_PRECISION=ON
    export TRAVIS_BUILD_INTEL=YES
    export TRAVIS_COMPILER=gcc
    export CC=gcc
    gcc --version
    pyenv shell 3.6.3

    wget 'https://raw.githubusercontent.com/nemequ/icc-travis/master/install-icc.sh' # installs from http://registrationcenter-download.intel.com/akdlm/irc_nas/9061/parallel_studio_xe_2016_update3_online.sh
    chmod 755 install-icc.sh
    ./install-icc.sh --components ifort,icc,mkl
    source ~/.bashrc
    pip3 install numpy
    mkdir build && cd build
    cmake .. -DBUILD_TESTING=ON -DDOUBLE_PRECISION=$DOUBLE_PRECISION -DBUILD_SHARED_LIBS=ON
    make -j 8 install
