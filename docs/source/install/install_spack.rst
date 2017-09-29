.. _install_spack:

Building OpenFAST with Spack
============================

OpenFAST uses the `CMake <https://cmake.org>`_ build system. 
We recommend building OpenFAST using `Spack <https://spack.readthedocs.io/en/latest>`__. 
However, we also provide some sample scripts in ``openfast/share`` if you choose to proceed without `Spack`.

Dependencies
------------

OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and ``LAPACK_LIBRARIES`` appropriately for cmake if the library isn't found in standard paths. Use `BLASLIB` as an example when using Intel MKL.
- For the optional C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`__ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`__ (provided by ``YAML_ROOT``)
- For the optional testing framework, Python 3+

CMake Build Instructions
------------------------
::

    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    mkdir build && cd build
    cmake ../ 
    make 
    
A sample installation shell script is also provided in the ``openfast/share``. Update the location of ``openfast_dir`` in the ``fast-build*`` scripts and run the script from ``openfast/`` as:
::

    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    bash share/fast-install.sh

Current CMake Options
~~~~~~~~~~~~~~~~~~~~~

-  ``BUILD_DOCUMENTATION`` -  Build documentation (Default: OFF)
-  ``BUILD_FAST_CPP_API`` - Enable building OpenFAST - C++ API (Default: OFF)
-  ``BUILD_SHARED_LIBS`` - Enable building shared libraries (Default: OFF)
-  ``CMAKE_BUILD_TYPE`` - Choose the build type: Debug Release (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Install path prefix, prepended onto install directories.
-  ``DOUBLE_PRECISION`` - Treat REAL as double precision (Default: ON)
-  ``FPE_TRAP_ENABLED`` -  Enable Floating Point Exception (FPE) trap in compiler options (Default: OFF)
-  ``ORCA_DLL_LOAD`` - Enable OrcaFlex library load (Default: OFF)
-  ``USE_DLL_INTERFACE`` - Enable runtime loading of dynamic libraries (Default: ON)

Building OpenFAST Semi-Automatically Using Spack on macOS or Linux
---------------------------------------------------------------------

The following describes how to build OpenFAST and its dependencies
mostly automatically on your Mac using `Spack <https://spack.readthedocs.io/en/latest>`_. 
This can also be used as a template to build OpenFAST on any 
Linux system with Spack.

These instructions were developed on macOS 10.11 with the following tools installed via Homebrew:

- GCC 6.3.0
- CMake 3.6.1
- pkg-config 0.29.2

Step 1
~~~~~~

Checkout the official Spack repo from github (we will checkout into ``${HOME}``):

``cd ${HOME} && git clone https://github.com/LLNL/spack.git``

Step 2
~~~~~~

Add Spack shell support to your ``.profile`` by adding the lines:

::

    export SPACK_ROOT=${HOME}/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

Step 3
~~~~~~

Copy the https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py file
to your installation of Spack:

::
   
    mkdir ${SPACK_ROOT}/var/spack/repos/builtin/packages/openfast
    cd ${SPACK_ROOT}/var/spack/repos/builtin/packages/openfast
    wget --no-check-certificate https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py

Step 4
~~~~~~

Try ``spack info openfast`` to see if Spack works. If it does, check the
compilers you have available by:

::

    machine:~ user$ spack compilers
    ==> Available compilers
    -- gcc ----------------------------------------------------------
    gcc@6.3.0  gcc@4.2.1

    -- clang --------------------------------------------------------
    clang@8.0.0-apple  clang@7.3.0-apple

Step 5
~~~~~~

Install OpenFAST with your chosen version of GCC:

::

    spack install openfast %gcc@6.3.0

To install OpenFAST with the C++ API, do:

::

    spack install openfast+cxx %gcc@6.3.0
    
That should be it! Spack will automatically use the most up-to-date dependencies 
unless otherwise specified. For example to constrain OpenFAST to use some specific 
versions of dependencies you could issue the Spack install command:

::

    spack install openfast %gcc@6.3.0 ^hdf5@1.8.16 

The executables and libraries will be located at

::
   
    spack location -i openfast

    
Add the appropriate paths to your ``PATH`` and ``LD_LIBRARY_PATH`` to run OpenFAST.
