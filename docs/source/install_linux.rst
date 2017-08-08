Building OpenFAST
=================

OpenFAST uses the `CMake <https://cmake.org>`__ build system. We recommend installing OpenFAST using `spack <https://spack.readthedocs.io/en/latest>`__. However, we also provide some sample scripts in ``share`` folder if you choose to install without `spack`.

CMake Build Instructions
------------------------

::

    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    mkdir build && cd build
    cmake ../ 
    make 


Current CMake Options
~~~~~~~~~~~~~~~~~~~~~

-  ``BUILD_DOCUMENTATION`` -  Build documentation (Default: OFF)
-  ``BUILD_FAST_CPP_API`` - Enable building FAST - C++ API (Default: OFF)
-  ``BUILD_SHARED_LIBS`` - Enable building shared libraries (Default: OFF)
-  ``CMAKE_BUILD_TYPE`` - Choose the build type: Debug Release (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Install path prefix, prepended onto install directories.
-  ``DOUBLE_PRECISION`` - Treat REAL as double precision (Default: ON)
-  ``FPE_TRAP_ENABLED`` -  Enable Floating Point Exception (FPE) trap in compiler options (Default: OFF)
-  ``ORCA_DLL_LOAD`` - Enable OrcaFlex library load (Default: OFF)
-  ``USE_DLL_INTERFACE`` - Enable runtime loading of dynamic libraries (Default: ON)

Dependencies
~~~~~~~~~~~~

OpenFAST has the following dependencies:

- ``LAPACK`` libraries provided through the variable ``BLASLIB``
- for the C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`__ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`__ (provided by ``YAML_ROOT``). 



Building OpenFAST Semi-Automatically Using Spack on Mac OS X or Linux
---------------------------------------------------------------------

The following describes how to build OpenFAST and its dependencies
mostly automatically on your Mac using 
`Spack <https://spack.readthedocs.io/en/latest>`__. 
This can also be used as a template to build OpenFAST on any 
Linux system with Spack.

Step 1
~~~~~~

This assumes you have a Homebrew installation of GCC installed already 
(we are using GCC 6.3.0). These instructions have been tested on OSX 10.11. MacOS 10.12 
will not build CMake with GCC anymore, so these instructions won't work 
in that case, but we have built OpenFAST using Spack on MacOS Sierra by
using Homebrew to install ``cmake`` and ``pkg-config`` and defining these 
as external packages in Spack (see 
`packages.yaml.mac_sierra <https://github.com/NaluCFD/NaluSpack/blob/master/spack_config/packages.yaml.mac_sierra>`__).

Step 2
~~~~~~

Checkout the official Spack repo from github (we will checkout into ``${HOME}``):

``cd ${HOME} && git clone https://github.com/LLNL/spack.git``

Step 3
~~~~~~

Add Spack shell support to your ``.profile`` by adding the lines:

::

    export SPACK_ROOT=${HOME}/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

Step 4
~~~~~~

Copy the `https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py`__ file
to your installation of Spack:

::
   
    mkdir ${SPACK_ROOT}/etc/spack/openfast ; cd ${SPACK_ROOT}/etc/spack/openfast
    wget --no-check-certificate https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py

Step 5
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

Step 6
~~~~~~

Install OpenFAST with whatever version of GCC (6.3.0 for us) you have
installed from Homebrew and force the build to use CMake 3.6.1 because
newer versions currently don't build on OS X:

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


Building OpenFAST manually on Mac OS X or Linux
-----------------------------------------------

A sample installation shell script is provided in the `share` folder. Run the script from `openfast_dir` as:

::
   
    git clone https://github.com/OpenFAST/OpenFAST.git
    cd OpenFAST
    bash share/install.sh
