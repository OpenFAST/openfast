.. _install_cmake_linux:

Building OpenFAST with CMake on Linux and Mac
=============================================

We describe here how to install OpenFAST (or any of its modules) using the `CMake <https://cmake.org>`_ 
build system on Linux or Mac OS systems. Separate CMake documentation is 
provided for Windows users at :numref:`install_cmake_windows` and Cygwin users at :numref:`install_cmake_cygwin`.
Also, some template build scripts are available in ``openfast/share``.

Required software for building OpenFAST 
---------------------------------------

In order to build OpenFAST using CMake, one needs the following minimum set of packages installed:

- Fortran compiler (GNU compiler version above 4.6.0 or Intel compiler version above 11)

- C/C++ compiler

- GNU Make (version 3.81 or later)

- CMake (version 2.8.12 or later)

OpenFAST third party library dependencies
-----------------------------------------

OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and ``LAPACK_LIBRARIES`` appropriately for CMake if the library is not found in standard paths. Use `BLASLIB` as an example when using Intel MKL.

- **Optional:** For the C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`_ (provided by ``HDF5_ROOT``) and `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_ (provided by ``YAML_ROOT``)

- **Optional:** For the testing framework, Python 3+

.. _cmake-build-instructions:

CMake build instructions
------------------------

If one has the appropriate third party libraries, CMake, and git installed, obtaining and building OpenFAST can be accomplished as follows:

.. code-block:: bash

    # obtain the source code; e.g., from the command line using git:
    git clone https://github.com/OpenFAST/OpenFAST.git

    # go to the OpenFAST directory
    cd OpenFAST

    # create a directory called `build`
    mkdir build 

    # go to the build directory
    cd build

    # execute CMake with the default options, which will create a series of Makefiles
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

The CMake configuration and resulting build can be customized easily through explicitly setting CMake variables. In general,
this is done by passing a flag in the CMake configuration command

.. code-block:: bash

    cmake .. -D<CMAKE_FLAG>=ON
    cmake .. -D<CMAKE_FLAG>=/usr/local/bin/this_thing

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


Parallel build
~~~~~~~~~~~~~~

GNU Make has a parellel build option with the ``-jobs`` or ``-j`` flag, and the OpenFAST
CMake configuration handles setting up the dependencies for Make so the build can be 
parallelized. However, it is important to note that the only parallel portion
of the build process is in compiling the modules. Due to some interdependency between
modules, the max parallel level is around 12. The remaining portion of the build,
mainly compiling the OpenFAST library itself, takes a considerable amount of time
and cannot be parallelized.

An example parallel build command is ``make -j 8``.

