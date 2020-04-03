.. _Running-OLAF:

Running OLAF
============

This section discusses how to obtain and execute OLAF from a personal
computer.

Advanced command line users looking for a quick guide may follow the
commands listed here; other users are recommended to follow the detailed
instructions in the coming paragraphs.

::

   # sudo apt-get install gfortran-8 cmake liblapack-dev libblas
   #    OR
   # module load cmake comp-intel mkl
   git clone https://github.com/openfast/openfast
   cd openfast
   mkdir build
   cd build
   cmake ..
   make
   ./modules/glue-codes/openfast/openfast [INPUTFILE]

Downloading the OLAF Software
-----------------------------

Download the OpenFAST archive, which includes OLAF, from the git
repository at https://github.com/OpenFAST/openfast.git.

Compiling OpenFAST with OLAF
----------------------------

This is the same process as compiling OpenFAST. The reader may obtain
detailed instructions at the following address:
https://openfast.readthedocs.io/. For completeness, basic instructions
are included here.

Compiling on a Microsoft Windows Machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three main options are available to Windows users:

-  Commercial software: install Intel Fortran Compiler and Microsoft
   Visual Studio. Open the visual studio solution file,
   ``vs-build/FAST/FAST.sln, with Visual Studio. From the Visual Studio menu, select: Build/Build solution.``

-  Linux subsystem: Install a Linux subsystem for Windows
   (https://docs.microsoft.com/en-us/windows/wsl/install-win10. You may
   be choose Ubuntu. Then, in a terminal, follow the instructions for
   Linux machines provided in . Note that the compiled code will only be
   run within the Linux subsystem environment, but the input and output
   files can be accessed seamlessly from Windows.

-  Free software: Compilation using free software is an involved process
   reserved to advanced users. Users interested in such options are
   referred to https://openfast.readthedocs.io/.

.. _sec:CompileLinux:

Compiling on Linux/Mac Machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec:linuxdep:

Dependencies/requirements
^^^^^^^^^^^^^^^^^^^^^^^^^

Compiling OpenFAST requires CMake, a Fortran compiler (Intel Fortran or
GFortran), and Linear Algebra Package (LAPACK) libraries (LAPACK or Math
Kernal Library [MKL]).

-  If you are in a Linux system such as Ubuntu, all the dependencies can
   be installed as follows:

   ``$ sudo apt-get install gfortran-8 cmake liblapack-dev libblas``

-  If you are in a cluster environment, such tools are usually already
   installed but need to be loaded into the environment. For example,
   the following command can load CMake, the Intel Fortran compiler, and
   the MKL libraries:

   ``$ module load cmake comp-intel mkl``

Additional resources are provided in .

Compilation
^^^^^^^^^^^

The compilation steps are as follows:

#. Optional: if you are in a cluster environment, make sure the
   dependencies are loaded (see )

#. In the OpenFAST directory, create a build directory and move into it:

   ``$ mkdir build``

   ``$ cd build``

#. In the build directory, run the following commands:

   ``$ cmake ..``

   ``$ make``

Advanced users may choose to customize the compilation step by providing
extra arguments to the CMake command. For instance, compiling the code
with the ``DEBUG option and single precision is done as follows:``

``$ cmake .. -DCMAKE_BUILD_Type=Debug -DDOUBLE_PRECISION:BOOL=OFF``

Compiling the code with optimization flags is done as follows:

``$ cmake .. -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -xhost" -DDOUBLE_PRECISION:BOOL=OFF``

.. _sec:cmdOptions:

Command/option explanations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``git clone`` :math:`\--` Git is a free and open-source, distributed,
version-control system designed to handle everything from small to very
large projects with speed and efficiency. :math:`\$` git clone clones a
remote repository into a new directory:

``mkdir`` :math:`\--` make a new directory

``cd`` :math:`\--` move into the specified directory

``module load`` :math:`\--` load specified modules (required on most
high-performance computing systems)

``cmake`` :math:`\--` CMake is an open-source, cross-platform family of
tools designed to build, test, and package software

``DCMAKE_Fortran_FLAGS_RELEASE`` :math:`\--` specify optional
compilation flags

``threads`` and ``qopenmp`` :math:`\--` required flags to run FAST.Farm
with OpenMP parallelization

``O2`` :math:`\--` optimization-level specification (recommended when
not working to debug code)

``xhost`` :math:`\--` additional optimization-level specification

``DDOUBLE_PRECISION:BOOL=OFF`` :math:`\--` compile code in single
precision (recommended)

``make`` :math:`\--` compile code using Makefile located in directory
(here, Makefile was generated with CMake).

.. _running-olaf-1:

Running OLAF
------------

#. Executable file will be located in:

   -  PC: ``~/OpenFAST/build/bin/openfast``

   -  Mac/Linux: ``~/OpenFAST/build/glue-codes/openfast/openfast``

#. Create a directory to all your OpenFAST and OLAF input files. You
   will need:

   -  ``<file>.fst`` and associated ``<file>.dat`` files :math:`\-–`
      turbine files

   -  Wind inflow files

#. Run code with
   ``$~/OpenFAST/build/glue-codes/openfast/openfast <file>.fst``.
