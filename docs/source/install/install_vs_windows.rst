.. _install_vs_windows:

Building OpenFAST on Windows with Visual Studio
===============================================

These instructions are specifically for the standalone Visual Studio project at `openfast/vs-build`.
Separate CMake documentation is provided for Windows users at :numref:`cmake_windows`.

Prerequisites
------------------------

1. A version of Visual Studio (VS).  

    -  NOTE: not all VS Studio versions are supported by the Intel compilers. In general, the Fortran compiler must be newer than Visual Studio.  A list of Intel Fortran compatible VS versions and specific installation notes are found `here <https://software.intel.com/en-us/intel-parallel-studio-xe-compilers-required-microsoft-visual-studio>`_.    

    -  Currently VS 2019 Community Edition, VS 2022 Professional, and VS 2022 Community Edition have been tested with OpenFAST.  Download VS 2022 Community `here <https://aka.ms/vs/17/release/vs_community.exe>`__.

    -  When installing Visual Studio, select the ``Desktop development with C++`` under ``Workloads``.

    -  Note: The included C/C++ project files for MAP++ and the Registry are compatible with VS 2019, but will upgrade seemlessly to a newer version of VS.

2. Intel Fortran Compiler

    -  We recommend compiling with the IFX compiler from Intel.  This is included in the ``Intel Fortran Essentials`` installation package.  Currently tested with version 2025.3

    -  You can download ``Intel Fortran Essentials`` `here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html>`__.  Note: do not install the ``oneAPI HPC Toolkit``

    -  Only install Intel Fortran after you have completed your Visual Studio installation. Note that Intel Fortran must be compatible with your version of Visual Studio.  See `here <https://www.intel.com/content/www/us/en/developer/articles/reference-implementation/intel-compilers-compatibility-with-microsoft-visual-studio-and-xcode.html>`__ for compatibility tables.

3. Git for Windows

    -  Download and install `git <https://git-scm.com/download/win>`__ for Windows.
    
4. Python 3.x for Windows (for regression/unit testing)

    -  The testing framework of OpenFAST requires the use of Python.  

    -  Please see :numref:`testing`  on testing OpenFAST for further information on this topic.

    -  We have been working with Continuum's `Anaconda <https://www.anaconda.com/download/#windows>`__ installation of Python 3.6 for Windows.


Compiling OpenFAST
------------------

1. Open ``A command prompt``, or ``git bash`` shell from the ``Start`` menu

2. Create a directory where you will clone OpenFAST repository (change
   ``code`` to your preferred name)

   ::

    mkdir code
    cd code

3. Clone the OpenFAST repository

   ::

    git clone https://github.com/openfast/openfast.git

This will create a directory called ``openfast`` within the ``code``
directory.

4. Using Windows Explorer, navigate to the directory ``openfast\vs-build\FAST``
   and double-click on the ``FAST.sln`` Visual Studio solution file.  This will 
   open Visual Studio with the FAST solution and its associated projects.
   
NOTE: If you are using Visual Studio 2015 or newer, you will be asked to upgrade
both the ``Fast_Registry.vcxproj`` and the ``MAP_dll.vcxproj`` files to a newer
format.  Go ahead and accept the upgrade on those files.

5. Select the desired Solution Configuration, such as ``Release``, and the 
   desired Solution Platform, such as ``x64`` by using the drop down boxes 
   located below the menubar.
   
6. Build the solution using the ``Build->Build Solution`` menu option.

7. If the solution built without errors, the executable will be located under the ``openfast\build\bin`` folder.
