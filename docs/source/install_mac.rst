Installing on Mac OS X
======================

Installing Mostly Automatically Using Spack
-------------------------------------------

The following describes how to build Nalu and its dependencies
automatically on your Mac using Spack. This can also be
used as a template to build Nalu on any system with Spack.

Step 1
~~~~~~

This assumes you have a Homebrew installation of GCC installed already
(I am using GCC 6.2.0).

Step 2
~~~~~~

Checkout the official Spack repo from github:

``cd ${HOME} && git clone https://github.com/LLNL/spack.git``

Step 3
~~~~~~

Add Spack shell support to your *.profile*:

::

    export SPACK_ROOT=${HOME}/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

Step 4
~~~~~~

Copy the ``nalu`` and ``nalu-trilinos`` package files from the NaluSpack repo to
your installation of Spack:

::

    cd ${HOME} && git clone https://github.com/jrood-nrel/NaluSpack.git
    cp -R ${HOME}/NaluSpack/nalu ${SPACK_ROOT}/var/spack/repos/builtin/packages/
    cp -R ${HOME}/NaluSpack/nalu-trilinos ${SPACK_ROOT}/var/spack/repos/builtin/packages/

Step 5
~~~~~~

Try ``spack info nalu`` to see if Spack works. If it does, check the
compilers you have available by:

::

    machine:~ user$ spack compilers
    ==> Available compilers
    -- gcc ----------------------------------------------------------
    gcc@6.2.0  gcc@4.2.1

    -- clang --------------------------------------------------------
    clang@8.0.0-apple  clang@7.3.0-apple

Step 6
~~~~~~

Install Nalu with whatever version of GCC (6.2.0 for me) you have
installed from Homebrew and force the build to use CMake 3.6.1 because
3.7.1 currently doesn't build on OS X:

::

    spack install nalu %gcc@6.2.0 ^cmake@3.6.1

That should be it!


Installing Manually
-------------------

The Mac platform installation is managed by the extensive usage of
"HomeBrew". This package provides many Mac OS X builds and installations
required for Nalu.

Homebrew
~~~~~~~~

Download and Install
`Homebrew <https://github.com/Homebrew/homebrew/wiki/Installation>`__ on
your local home terminal:

Now you can install homebrew; "brew doctor" will be the first line
comamnd required.

::

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew doctor

Packages to be obtained from homebrew
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenMpi
^^^^^^^

::

    brew install openmpi

Cmake
^^^^^

::

    brew search cmake  //cmake is there and will pop up on the next line
    brew install cmake

libxml2
^^^^^^^

::

    brew search libxml2
    brew install libxml2

boost, 1.60.0\_2
^^^^^^^^^^^^^^^^

::

    brew install boost

SuperLu, 4.3
^^^^^^^^^^^^

::

    brew tap homebrew/science
    brew search superlu 
    brew install superlu43

The latest version of SuperLu provided by homebrew is not compatible
with head Trilinos due to Trilinos' usage of some deprecated methods

Once done with the HB install, make sure everybody can read the files:

::

    sudo chmod -R u+rwX,go+rX,go-w /usr/local/Cellar

Non-Homebrew
~~~~~~~~~~~~

Other Nalu required libraries must be managed outside of the homebrew
environment

Helpful notes
^^^^^^^^^^^^^

How to specify install location and compilers with ``configure`` and
``cmake``

When building packages that use configure, do a

::

    ./configure PREFIX=/usr/local/pachages/install CC=mpicc CXX=mpicxx

When building packages that use cmake, do a ``mkdir build; cd build``
then

::

    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX:PATH=/myPath/install ..`

For all non-HomeBrew packages, myPath above will be /usr/local/packages

::

    mkdir /usr/local/packages
    mkdir /usr/local/packages/install

yaml-cpp
^^^^^^^^

For versions of Nalu after the v1.1.0-release, Yaml is provided under
`github <https://github.com/jbeder/yaml-cpp>`__

::

    cd $nalu_build_dir/packages
    git clone https://github.com/jbeder/yaml-cpp

Build yaml-cpp

::

    cd $nalu_build_dir/packages/yaml-cpp
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$nalu_build_dir/install ..
    make
    make install

Pre-v1.1.0-release; yaml-cpp, Version 0.3.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yaml is provided by
`code.google.com <https://code.google.com/p/yaml-cpp/downloads/detail?name=yaml-cpp-0.3.0.tar.gz&can=2&q=>`__

Follow the yaml installation directions. This process will put the
created files in /user/local/include/yaml-cpp. Below are some high level
points:

::

    mkdir /usr/local/packages
    cd packages/
    curl -o yaml-cpp-0.3.0.tar.gz https://yaml-cpp.googlecode.com/files/yaml-cpp-0.3.0.tar.gz 
    tar -zxvf yaml-cpp-0.3.0.tar.gz 
    mv yaml-cpp yaml-cpp-0.3.0

This series of commands will create
``/usr/local/packages/yaml-cpp-0.3.0``

Next, build yaml-cpp

::

    cd /usr/local/packages/yaml-cpp-0.3.0
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX:PATH=/usr/local/packages/install ..
    make
    make install

zlib, 1.2.8
^^^^^^^^^^^

zlib is provided by the `zlib <http://www.zlib.net>`__ project.

::

    cd /usr/local/packages/
    curl -o zlib-1.2.8.tar.gz http://zlib.net/zlib-1.2.8.tar.gz
    tar -zxvf zlib-1.2.8.tar.gz 

Build zlib

::

    cd /usr/local/packages/zlib-1.2.8
    CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --archs="-arch x86_64" --prefix=/usr/local/packages/install/
    make
    make install

hdf5, 1.8.12
^^^^^^^^^^^^

hdf5 1.8.12 is provided by the
`HDF <http://www.hdfgroup.org/downloads/index.html>`__ group

::

    cd /usr/local/packages/
    curl -o hdf5-1.8.12.tar.gz http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
    tar -zxvf hdf5-1.8.12.tar.gz 

This series of commands will create ``/usr/local/packages/hdf5-1.8.12``

Build

::

    cd /usr/local/packages/hdf5-1.8.12
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" --enable-parallel --with-zlib=/usr/local/packages/install --prefix=/usr/local/packages/install
    make
    make install
    make check

Full Parallel-Enabled Nalu using netCDF (V. 4.3.3.1) and Parallel netCDF (V. 1.6.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to support all aspects of Nalu's parallel models, this
combination of products is required.

Parallel netCDF, Version 1.6.1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parallel netCDF is provided on the `Argon Trac
Page <https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download>`__.

::

    cd $nalu_build_dir/packages/
    tar -zxvf parallel-netcdf-1.6.1.tar.gz

Configure, build and install:

::

    cd parallel-netcdf-1.6.1
    ./configure --prefix=/usr/local/packages/install CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS=-L/usr/local/packages/install/lib --disable-fortran
    make
    make install

Note that I have created an install directory that might look like:
$nalu\_build\_dir/install

netCDF Version 4.3.3.1
^^^^^^^^^^^^^^^^^^^^^^

netCDF is provided on
`github <https://github.com/Unidata/netcdf-c/releases>`__

::

    cd $nalu_build_dir/packages/
    curl -o netcdf-c-4.3.3.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.3.1
    tar -zxvf netcdf-c-4.3.3.1.tar.gz

Configure, build and install

::

    cd netcdf-c-4.3.3.1
    ./configure --prefix=$nalu_install_dir CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS=-L/usr/local/packages/install/lib --enable-pnetcdf --enable-parallel-tests --enable-netcdf-4 --disable-shared --disable-fsync --disable-cdmremote --disable-dap --disable-doxygen --disable-v2
    make -j 4 
    make install
    make check

Note that when using Parallel netCDF, the proper install directories
must be added to the Trilinos configuration file.

Partial Parallel-Enabled Nalu using netCDF, Version 4.3.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Parallel netCDF is omitted, follow the instructions below. On some
platforms, however, autodecompostion may fail.

netcdf is provided on
`github <https://github.com/Unidata/netcdf-c/releases>`__

Scroll down until you see "NetCDF-C 4.3.1.1 (Bugfix Release)" or similar
Click on the "Source (tar.gz)" button to download and then move the tar
file to: ``/usr/local/packages``

::

    cd /usr/local/packages
    tar -xvf netcdf-c-4.3.1.1.tar

This process will create ``/usr/local/packages/netcdf-c-4.3.1.1``

Build

::

    cd /usr/local/packages/netcdf-c-4.3.1.1
    ./configure --prefix=/usr/local/packages/install CC=mpicc FC=mpif90 CXX=mpicxx CFLAGS="-I/usr/local/packages/install/include -O3" LDFLAGS="-L/usr/local/packages/install/lib" --disable-fsync --disable-cdmremote --disable-dap --disable-shared --disable-doxygen
    make
    make install
    make check

Trilinos
^^^^^^^^

Trilinos is managed by the `Trilinos <http://www.trilinos.org>`__
project and can be found on github.

Clone the latest version of Trilinos within ``/packages``:

::

    cd /usr/local/packages/
    git clone https://github.com/trilinos/Trilinos.git

In some cases, the master Trilinos code base may have build issues. This
is a rare occurance, however, some aspects to Trilinos that Nalu
require, e.g., Tpetra, kokkos, STK and Muelu are in ``active``
development. If problems arise, one can revert back to a possible
successful SHA-1 using bisect. Again, this is hopefully going to be
mitigated by the strong SQA efforts at SNL.

Nalu Releases
~~~~~~~~~~~~~

Unfortunately, github does not allow for a "live" wiki for each of the
existing branches of Nalu.wiki. As such, instructions for the particular
releases have been embedded within this head wiki file.

Release v1.0.0-release
^^^^^^^^^^^^^^^^^^^^^^

For the formal v1.0.0-release, check-out the following Trilinos Version:

::

        git checkout trilinos-release-12-0-branch       

This version is the expected Trilinos code base for the v1.0.0-release
Nalu code base. Now proceed to the build section.

Head Code Base
^^^^^^^^^^^^^^

Proceed to the build section without checking out the Trilinos
12-0-branch.

Build
^^^^^

Create new folder in Trilinos called build,

::

    cd /usr/local/packages/Trilinos
    mkdir build

Place into build the script one of the ``do-configTrilinos_*`` files.

``do-configTrilinos_*`` will be used to run cmake to build trilinos
correctly for Nalu. Note that there are two files: one for ``release``
and the other ``debug``. The files can be found on the Nalu GitHub site
​here or copied from ``$nalu_build_dir/packages/Nalu/build``, which is
created in the Nalu build step documented below. For example:

Pull latest version of do-configTrilinos\_\* from Nalu's GitHub site:

::

    curl -o $nalu_build_dir/packages/Trilinos/build/do-configTrilinos_release https://raw.githubusercontent.com/NaluCFD/Nalu/master/build/do-configTrilinos_release

or if you create the Nalu directory as directed below, simply copy one
of the ``do-configTrilinos_*`` files from local copy of Nalu's git
repository:

::

    cp $nalu_build_dir/packages/Nalu/build/do-configTrilinos_release $nalu_build_dir/packages/Trilinos/build

Now edit the do-configTrilinos\_release to modify the defined paths as
follows:

::

    mpi_base_dir=/usr/local/Cellar/open-mpi/1.8.3
    nalu_build_dir=/usr/local/packages

Next, note that some packages, i.e., boost and superLu were provided by
HomeBrew. As such, make sure that boost\_dir and super\_lu point to
``/usr/local/Cellar``

Build

::

    cd /usr/local/packages/Trilinos/build
    ./do-configTrilinos_release
    make 
    make install

Nalu, the guest of honor
~~~~~~~~~~~~~~~~~~~~~~~~

Nalu is provided by `github <https://github.com/NaluCFD/Nalu>`__

No doubt, you already have cloned Nalu. If not, execute the following
command in the location that you want Nalu:

::

    git clone https://github.com/NaluCFD/Nalu.git

Nalu Releases
^^^^^^^^^^^^^

One may either build the released Nalu version, v1.0.0-release, or the
head code base.

Release v1.0.0-release
^^^^^^^^^^^^^^^^^^^^^^

For the formal Nalu v1.0.0-release, you should have already cloned
Trilinos and built the 12.0 release version of Trilinos. To obtain the
consistent Nalu version, after the clone, checkout the Nalu release,

::

    git checkout v1.0.0-release

Now proceed to the build section below.

Head Code Base
^^^^^^^^^^^^^^

Proceed to the build section without checking out the Nalu
v1.0.0-release code repository.

Build
^^^^^

In ``Nalu/build``, you will find the
`CMakeLists.txt <https://github.com/NaluCFD/Nalu/blob/master/CMakeLists.txt>`__
and
`do-configNalu\_release <https://github.com/NaluCFD/Nalu/blob/master/build/do-configNalu_release>`__.

Again, note that there is a debug version as well. Copy the
do-configNalu\_release to a new, non-tracked file,

::

    cp do-configNalu_release do-configNaluNonTracked

Edit the paths at the top of the files by defining the
``nalu_build_dir variable`` as:

::

    nalu_build_dir=/usr/local/packages

within ``Nalu/build``, execute the following commands

::

    ./do-configNaluNonTracked
    make 

This process will create ``naluX`` within the ``Nalu/build`` location.
Setting the DEBUG CMake option will create a naluXd executable.

Other useful tools from, e.g., seacas, are under
``/usr/local/packages/install/trilinos/bin``

Testing
~~~~~~~

After the ``naluX`` executable is created, please proceed with
regression testing to ensure a proper build.

Instructions for the regression testing can be found under the
`NaluRtest <https://github.com/NaluCFD/NaluRtest>`__ directory.

Please take care to ensure that the NaluRtest branch is consistent with
the Trilinos and Nalu version desired. One may also need to checkout the
v1.0.0-release code base.

