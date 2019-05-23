.. _install_cmake_windows:

Building OpenFAST on Windows with CMake and Visual Studio
=========================================================

We describe here how to install OpenFAST (or any of its modules) using the `CMake <https://cmake.org>`_ 
build system on Windows systems. Separate CMake documentation is 
provided for Cygwin users at :numref:`install_cmake_cygwin` and Linux/Mac users at :numref:`install_cmake_linux`.
A standalone Visual Studio solution also exists at `openfast/vs-build` and documentation is at :numref:`install_vs_windows`.

Required software for building OpenFAST
---------------------------------------

In order to build OpenFAST using CMake, one needs the following minimum set of packages installed:

- Fortran compiler (GNU compiler version above 4.6.0 or Intel compiler version above 11)

- C/C++ compiler

- Visual Studio

- CMake (version 2.8.12 or later)

OpenFAST third party library dependencies
-----------------------------------------

OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and ``LAPACK_LIBRARIES`` appropriately for CMake if the library is not found in standard paths. Use `BLASLIB` as an example when using Intel MKL.

- **Optional:** For the C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`_ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_ (provided by ``YAML_ROOT``)

- **Optional:** For the testing framework, Python 3+

CMake build instructions
------------------------

If one has the appropriate third party libraries, CMake, and git installed, obtaining and building OpenFAST can be accomplished as follows:

.. code-block:: bash

    # Obtain the source code; e.g., from the command line using git:
    git clone https://github.com/OpenFAST/OpenFAST.git

    # Go to the OpenFAST directory
    cd OpenFAST

    # Create a directory called `build`
    mkdir build 

    # Go to the build directory
    cd build

    # Execute CMake with the default options and a specific Visual Studio version
    # and build architecture. For a list of available CMake generators, run
    # `cmake .. -G`
    cmake .. -G "Visual Studio 14 2015 Win64"

    # Open the generated Visual Studio solution
    start OpenFAST.sln

Visual Studio will open a solution containing all of the OpenFAST projects, and you
can then build any module library, module driver, or glue code. Note that any time 
CMake is rerun, the Visual Studio solution will be regenerated causing the Visual Studio
GUI to lag momentarily while it reloads the data.

**The CMake-generated Visual Studio build is currently damaged.** Some modules are compiled
before their associated registry type files are seen by Visual Studio so an initial build
will fail. However, a simple work around is to run the build command in Visual Studio
multiple times until it succeeds.


CMake options
~~~~~~~~~~~~~

Below is a list of current CMake options including their default settings (which will effect, e.g., the targets in a resulting ``Makefile``.  

-  ``BUILD_DOCUMENTATION`` -  Build documentation (Default: OFF)
-  ``BUILD_OPENFAST_CPP_API`` - Enable building OpenFAST - C++ API (Default: OFF)
-  ``BUILD_SHARED_LIBS`` - Enable building shared libraries (Default: OFF)
-  ``BUILD_TESTING`` - Build the testing tree (Default: OFF)
-  ``CMAKE_BUILD_TYPE`` - Choose the build type: Debug Release (Default: Release)
-  ``CMAKE_INSTALL_PREFIX`` - Install path prefix, prepended onto install directories.
-  ``DOUBLE_PRECISION`` - Treat REAL as double precision (Default: ON)
-  ``FPE_TRAP_ENABLED`` -  Enable Floating Point Exception (FPE) trap in compiler options (Default: OFF)
-  ``ORCA_DLL_LOAD`` - Enable OrcaFlex library load (Default: OFF)
-  ``USE_DLL_INTERFACE`` - Enable runtime loading of dynamic libraries (Default: ON)

CMake options can be configured through command line, e.g.

.. code-block:: bash

    # to enable Makefile for local building of sphinx-based documentation
    cmake .. -DBUILD_DOCUMENTATION:BOOL=ON

    # to compile OpenFAST in single precision
    cmake .. -DDOUBLE_PRECISION:BOOL=OFF
 

Custom CMake builds
~~~~~~~~~~~~~~~~~~~

The CMake configuration and resulting build can be customized easily by explicitly setting CMake variables. In general,
this is done by passing a flag in the CMake configuration command

.. code-block:: bash

    cmake .. -D<CMAKE_FLAG>=ON
    cmake .. -D<CMAKE_FLAG>=\home\user\Desktop\this_thing

This syntax is the same as in setting a CMake option and the result is used very similarly in the CMake configuration files.
Common customizations revolve around choosing a compiler or math library; for example

.. code-block:: bash

    cmake .. -DCMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-8 -DLAPACK_LIBRARIES=/System/Library/Frameworks/Accelerate.framework -DLAPACK_LIBRARIES=/System/Library/Frameworks/Accelerate.framework

**NOTE** Many CMake configurations can also be set through an environment variable.
For example, when using Intel's MKL, the math libraries can be discovered automatically by setting the ``MKLROOT``
environment variable. The Fortran compiler can also be set explicitly with the ``FC`` environment variable.

Here is a good resource for useful CMake variables: `GitLab useful cmake variables <https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/Useful-Variables>`_.
The `CMake documentation <https://cmake.org/cmake/help/latest/>`_ is also helpful for searching
through variables and determining the resulting action.
