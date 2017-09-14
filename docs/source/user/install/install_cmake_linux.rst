.. _install_cmake_linux:

Building OpenFAST with CMake on Linux and Mac
=============================================

We describe here how to install OpenFAST (or any of its modules) using the `CMake <https://cmake.org>`_ build system on Linux or Mac OS systems.  Separate `CMake <https://cmake.org>`_ documentation is provided for Windows Cygwin users: see :ref:`install_cmake_cygwin`.

Dependencies
------------

OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and ``LAPACK_LIBRARIES`` appropriately for cmake if the library is not found in standard paths. Use `BLASLIB` as an example when using Intel MKL.

- **Optional:** For the C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`_ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_ (provided by ``YAML_ROOT``)

- **Optional:** For the testing framework, Python 3+

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
-  ``BUILD_FAST_CPP_API`` - Enable building FAST - C++ API (Default: OFF)
-  ``BUILD_SHARED_LIBS`` - Enable building shared libraries (Default: OFF)
-  ``CMAKE_BUILD_TYPE`` - Choose the build type: Debug Release (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Install path prefix, prepended onto install directories.
-  ``DOUBLE_PRECISION`` - Treat REAL as double precision (Default: ON)
-  ``FPE_TRAP_ENABLED`` -  Enable Floating Point Exception (FPE) trap in compiler options (Default: OFF)
-  ``ORCA_DLL_LOAD`` - Enable OrcaFlex library load (Default: OFF)
-  ``USE_DLL_INTERFACE`` - Enable runtime loading of dynamic libraries (Default: ON)

