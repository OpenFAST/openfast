.. _install_cygwin:

Building OpenFAST on Windows with CMake and Cygwin 64-bit
=========================================================
WARNING: This build process takes a significantly long amount of time.
If GNU tools are not required, it is recommended that Windows users see one
of the following sections:

- :ref:`download_binaries`
- :ref:`cmake_windows`
- :ref:`install_vs_windows`.

Installing prerequisites
------------------------

1. Download and install `Cygwin 64-bit <https://cygwin.com/setup-x86_64.exe>`__.
   You will need to ``Run as Administrator`` to complete the installation
   process.

    - Choose ``Install from internet``
    - Choose the default install location
    - Choose the default package download location
    - Choose ``Direct connection``
    - Choose a download site
    - See next step for ``select packages``. Alternately, you can skip this
      step and run ``setup-x86_64.exe`` anytime later to select and install
      required software.

2. Select packages necessary for compiling ``OpenFAST``. Choose ``binary``
   packages and not the source option.

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

    -  To run the test suite, install these optional packages from ``Python``:

         -  ``python3``
         -  ``Python3-numpy``

    -  Click ``Next`` and accept all additional packages that the setup
       process requests to install to satisfy dependencies

3. It is *recommended* that you reboot the machine after installing
   ``Cygwin`` and all the necessary packages.

Compiling OpenFAST
------------------
From here, pick up from the Linux with CMake instructions at
:ref:`cmake_unix`.

Other tips
----------
- If you would like to run ``openfast.exe`` from the ``cmd`` terminal, then you
  must add the ``C:\cygwin64\lib\lapack`` and
  ``C:\cygwin64\home\<USERNAME>\software\bin`` to your ``%PATH%`` variable in
  environment setting. Replace ``<USERNAME>`` with your account name on Windows
  system.

- It is suggested to compile with optimization level 2 for Cygwin. Do this by
  changing the build mode in the cmake command

    .. code-block:: bash

        cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
