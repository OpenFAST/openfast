.. _install_cmake_linux:

Building OpenFAST with CMake on Linux and Mac
=============================================

We describe here how to install OpenFAST (or any of its modules) using the `CMake <https://cmake.org>`_ build system on Linux or Mac OS systems.  Separate `CMake <https://cmake.org>`_ documentation is provided for Windows Cygwin users: see :ref:`install_cmake_cygwin`.

Required software for building OpenFAST 
---------------------------------------

In order to build OpenFAST using CMake, one needs the following minimum set of packages installed:

- Fortran compiler (e.g., gcc version 6.1.0 or Intel)

- C/C++ compiler

- CMake (version 2.8.12 or later)

OpenFAST third-party-library (TPL) dependencies
-----------------------------------------------

OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and ``LAPACK_LIBRARIES`` appropriately for cmake if the library is not found in standard paths. Use `BLASLIB` as an example when using Intel MKL.

- **Optional:** For the C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`_ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_ (provided by ``YAML_ROOT``)

- **Optional:** For the testing framework, Python 3+

CMake build instructions
------------------------

If one has the appropriate TPLs, cmake, and git installed, obtaining and building OpenFAST can be accomplished as follows:

.. code-block:: bash

    # obtain the source code; e.g., from the command line using git:
    git clone https://github.com/OpenFAST/OpenFAST.git

    # go to the OpenFAST directory
    cd OpenFAST

    # create a directory called `build`
    mkdir build 

    # go to the build directory
    cd build

    # execute cmake with the default options, which will create a Makefile
    cmake ../ 

    # execute a make command (with no target provided, equivalent to `make all`
    make 

This will build the OpenFAST suite in the ``build`` directory, which can be deleted for a clean build.

There are many  ``Makefile`` targets (besides ``all``), which can be listed via ``help``:

.. code-block:: bash

    # list available make targets
    make help

    # make a specific target, e.g.
    make beamdyn_driver



Current CMake options
~~~~~~~~~~~~~~~~~~~~~

Below is a list of current CMake options including their default settings (which will effect, e.g., the targets in a resulting ``Makefile``.  

-  ``BUILD_DOCUMENTATION`` -  Build documentation (Default: OFF)
-  ``BUILD_FAST_CPP_API`` - Enable building FAST - C++ API (Default: OFF)
-  ``BUILD_SHARED_LIBS`` - Enable building shared libraries (Default: OFF)
-  ``CMAKE_BUILD_TYPE`` - Choose the build type: Debug Release (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Install path prefix, prepended onto install directories.
-  ``DOUBLE_PRECISION`` - Treat REAL as double precision (Default: ON)
-  ``FPE_TRAP_ENABLED`` -  Enable Floating Point Exception (FPE) trap in compiler options (Default: OFF)
-  ``ORCA_DLL_LOAD`` - Enable OrcaFlex library load (Default: OFF)
-  ``USE_DLL_INTERFACE`` - Enable runtime loading of dynamic libraries (Default: ON)

CMake options can be executed from the command line as, e.g., the cmake command above could be exectuted as

.. code-block:: bash

    # e.g., to enable Makefile for local building of sphinx-based documentation
    cmake -DBUILD_DOCUMENTATION:BOOL=ON ..

    # e.g., to compile OpenFAST in single precision
    cmake -D:DOUBLE_PRECISIONBOOL=OFF ..
 

