.. _install_cmake_cygwin:

Building OpenFAST on Windows with CMake and Cygwin 64-bit
=========================================================

Installing prerequisites
------------------------

1. Download and install `Cygwin
   64-bit <https://cygwin.com/setup-x86_64.exe>`__. You will need to
   ``Run as Administrator`` to complete the installation process.

    -  Choose ``Install from internet``
    -  Choose the default install location
    -  Choose the default package download location
    -  Choose ``Direct connection``
    -  Choose ``http://www.gtlib.gatech.edu`` as the download site
    -  See next step for ``select packages``. Alternately, you can skip this
       step and run ``setup-x86_64.exe`` anytime later to select and install
       required software.

2. Select packages necessary for compiling ``OpenFAST``. Choose
   ``binary`` packages and not the source option.

    -  Choose ``Category`` view, we will be installing packages from
       ``Devel`` and ``Math``
    -  From ``Devel`` mark the following packages for installation

         -  ``cmake``
         -  ``cmake-doc``
         -  ``cmake-gui``
         -  ``cygwin-devel``
         -  ``gcc-core``
         -  ``gcc-fortran``
         -  ``gcc-g++``
         -  ``git``
         -  ``make``
         -  ``makedepend``

    -  From ``Math`` mark the following packages for installation

         -  ``liblapack-devel``
         -  ``libopenblas``

    -  Click ``Next`` and accept all additional packages that the setup
       process requests to install to satisfy dependencies

3. It is *recommended* that you reboot the machine after installing
   ``Cygwin`` and all the necessary packages.

Compiling OpenFAST
------------------

1. Open ``Cygwin64 Terminal`` from the ``Start`` menu

2. Create a directory where you will clone OpenFAST repository (change
   ``code`` to your preferred name)

::

   mkdir code
   cd code

3. Clone the OpenFAST repository

::

    git clone https://github.com/OpenFAST/OpenFAST.git

This will create a directory called ``OpenFAST`` within the ``code``
directory.

4. Create a build directory

::

    cd OpenFAST
    mkdir build
    cd build

5. Run ``cmake``. Note that this step is necessary only if you change
   compiler settings, or add new files to any of the ``CMakeLists.txt``.
   Modification of ``.f90`` files do not require you to run ``cmake``
   again, just re-run ``make`` command (see next item) to recompile with
   latest source code modifications.

::

    FC=gfortran cmake ../

Sample output is shown below:

::

    $ FC=gfortran cmake ../    
    -- The Fortran compiler identification is GNU 5.4.0    
    -- The C compiler identification is GNU 5.4.0    
    -- Check for working Fortran compiler: /usr/bin/gfortran.exe    
    -- Check for working Fortran compiler: /usr/bin/gfortran.exe  -- works    
    -- Detecting Fortran compiler ABI info    
    -- Detecting Fortran compiler ABI info - done    
    -- Checking whether /usr/bin/gfortran.exe supports Fortran 90    
    -- Checking whether /usr/bin/gfortran.exe supports Fortran 90 -- yes    
    -- Check for working C compiler: /usr/bin/cc    
    -- Check for working C compiler: /usr/bin/cc -- works    
    -- Detecting C compiler ABI info    
    -- Detecting C compiler ABI info - done    
    -- Detecting C compile features    
    -- Detecting C compile features - done    
    -- Looking for Fortran sgemm    
    -- Looking for Fortran sgemm - found    
    -- Looking for pthread.h    
    -- Looking for pthread.h - found    
    -- Looking for pthread_create    
    -- Looking for pthread_create - found    
    -- Found Threads: TRUE    
    -- A library with BLAS API found.    
    -- A library with BLAS API found.    
    -- Looking for Fortran cheev    
    -- Looking for Fortran cheev - found    
    -- A library with LAPACK API found.    
    -- Setting system file as: src/SysGnuLinux.f90    
    -- Configuring done    
    -- Generating done    
    -- Build files have been written to: /home/sanantha/code/OpenFAST/build

6. Compile ``OpenFAST``

::

    make

Grab a cup of coffee as this takes a while on Cygwin. Once the
compilation is completed, the ``OpenFAST`` executable is present in
``OpenFAST/build/glue-codes/fast/openfast.exe``

7. Test the executable

::

    $ glue-codes/fast/openfast.exe -h


    **************************************************************************************************
    FAST (v8.17.00a-bjj, 27-Aug-2016)

    Copyright (C) 2016 National Renewable Energy Laboratory

    This program comes with ABSOLUTELY NO WARRANTY. See the "license.txt" file distributed with this
    software for details.
    **************************************************************************************************

     Running FAST (v8.17.00a-bjj, 27-Aug-2016), compiled as a 64-bit application using double
     precision
     linked with NWTC Subroutine Library (v2.11.00, 12-Nov-2016)


     Syntax is:

        FAST_x64.exe [-h] <InputFile>

     where:

        -h generates this help message.
        <InputFile> is the name of the required primary input file.

     Note: values enclosed in square brackets [] are optional. Do not enter the brackets.


    FAST_InitializeAll:The required input file was not specified on the command line.

    FAST encountered an error during module initialization.
     Simulation error level: FATAL ERROR

     Aborting FAST.

\`\`\`

Other tips
----------

-  You can specify an installation location during your ``cmake``
   process so that the executable, libraries, and headers (e.g., ``MAP``
   and ``OpenFOAM`` headers) are installed in a common location that you
   can use to update your environment variables.

::

    # 1. Create an installation location mkdir -p ~/software

    # 2. Instruct CMake to use the custom install location FC=gfortran cmake
    -DCMAKE\_INSTALL\_PREFIX:PATH=$HOME/software ../

    # 3. Compile OpenFAST executable make

    # 4. Install OpenFAST to custom install location make install \`\`\`

With this step, you can execute ``make install`` after ``make`` (see
step 6 above). Now the ``openfast.exe`` and other executables (e.g.,
``aerodyn.exe``) are available in ``~/software/bin/`` directory.

-  If you desire to be able to run ``openfast.exe`` from the ``cmd``
   window, then you must add the ``C:\cygwin64\lib\lapack`` and
   ``C:\cygwin64\home\<USERNAME>\software\bin`` to your ``%PATH%``
   variable in environment setting. Replace ``<USERNAME>`` with your
   account name on windows system.

-  In addition to ``openfast.exe``, the current CMake setup also allows
   the user to compile other executables or libraries without compiling
   the entire codebase. Use ``make help`` to see what targets are
   available and then do ``make <TARGET>`` to choose your desired
   target. For example, ``make aerodyn`` will compile only the
   ``aerodyn.exe`` executable and its dependencies without compiling the
   remaining targets.
