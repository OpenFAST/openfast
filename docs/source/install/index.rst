.. _installation:

Installing OpenFAST
===================

The following pages provide instructions for building OpenFAST and/or its modules from source code.  
The developer team is moving towards a CMake-only approach that well supports Window Visual Studio users,
but at this time we provide a separate build path for those users.

Obtaining the OpenFAST source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenFAST can be cloned (i.e., downloaded) from its `Github Repository <https:// github.com/openfast/openfast>`_.
For example, from a command line:

::

    git clone https://github.com/OpenFAST/OpenFAST.git

It can also be downloaded directly from https://github.com/OpenFAST/OpenFAST.

Linux and Mac
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   install_cmake_linux.rst
   install_spack.rst

Windows
~~~~~~~

.. toctree::
   :maxdepth: 1

   install_cmake_windows.rst
   install_vs_windows.rst
   install_cmake_cygwin.rst
