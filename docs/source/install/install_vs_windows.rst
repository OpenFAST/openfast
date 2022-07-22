.. _install_vs_windows:

Building OpenFAST on Windows with Visual Studio
===============================================

These instructions are specifically for the standalone Visual Studio project at `openfast/vs-build`.
Separate CMake documentation is provided for Windows users at :numref:`cmake_windows`.

Prerequisites
------------------------

1. A version of Visual Studio (VS).  

    -  Currently VS 2013 Professional and VS 2015 Community Edition have been tested with OpenFAST.

    -  A list of Intel Fortran compatible VS versions and specific installation notes are found `here <https://software.intel.com/en-us/intel-parallel-studio-xe-compilers-required-microsoft-visual-studio>`_.    

    -  The included C/C++ project files for MAP++ and the Registry are compatible with VS 2013, but will upgrade seemlessly to a newer version of VS.

    -  If you download and install `Visual Studio 2015 Community Edition <https://go.microsoft.com/fwlink/?LinkId=691978&clcid=0x409>`__, you will need to be sure and select the ``C/C++ component`` using the ``Customize`` option.

2. Intel Fortran Compiler

    -  Currently only version 2017.1 has been tested with OpenFAST, but any newer version should be compatible.

    -  You can download an Intel Fortran compiler `here <https://software.intel.com/en-us/fortran-compilers>`__.

    -  Only install Intel Fortran after you have completed your Visual Studio installation.

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
