.. _install_spack:

Building OpenFAST with Spack
============================

The process to build and install OpenFAST with
`Spack <https://spack.readthedocs.io/en/latest>`__  on Linux or macOS is
described here.

Dependencies
------------
OpenFAST has the following dependencies:

- LAPACK libraries. Users should set ``BLAS_LIBRARIES`` and
  ``LAPACK_LIBRARIES`` appropriately for CMake if the library isn't found
  in standard paths. Use `BLASLIB` as an example when using Intel MKL.
- For the optional C++ API, `HDF5 <https://support.hdfgroup.org/HDF5/>`__
  (provided by ``HDF5_ROOT``) and
  `yaml-cpp <https://github.com/jbeder/yaml-cpp>`__ (provided by ``YAML_ROOT``)
- For the optional testing framework, Python 3+ and Numpy

Building OpenFAST Semi-Automatically Using Spack on macOS or Linux
------------------------------------------------------------------

The following describes how to build OpenFAST and its dependencies
mostly automatically on macOS using
`Spack <https://spack.readthedocs.io/en/latest>`_. This can also be used as a
template to build OpenFAST on any Linux system with Spack.

These instructions were developed on macOS 10.11 with the following tools
installed via Homebrew:

- GCC 6.3.0
- CMake 3.6.1
- pkg-config 0.29.2

Step 1
~~~~~~
Checkout the official Spack repo from github (we will checkout into
``${HOME}``):

.. code-block:: bash

    cd ${HOME} && git clone https://github.com/LLNL/spack.git

Step 2
~~~~~~
Add Spack shell support to your ``.profile`` by adding the lines:

.. code-block:: bash

    export SPACK_ROOT=${HOME}/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

Step 3
~~~~~~
Copy the https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py file
to your installation of Spack:

.. code-block:: bash

    mkdir ${SPACK_ROOT}/var/spack/repos/builtin/packages/openfast
    cd ${SPACK_ROOT}/var/spack/repos/builtin/packages/openfast
    wget --no-check-certificate https://raw.githubusercontent.com/OpenFAST/openfast/dev/share/spack/package.py

Step 4
~~~~~~
Try ``spack info openfast`` to see if Spack works. If it does, check the
compilers you have available by:

.. code-block:: bash

    machine:~ user$ spack compilers
    ==> Available compilers
    -- gcc ----------------------------------------------------------
    gcc@6.3.0  gcc@4.2.1

    -- clang --------------------------------------------------------
    clang@8.0.0-apple  clang@7.3.0-apple

Step 5
~~~~~~
Install OpenFAST with your chosen version of GCC:

.. code-block:: bash

    spack install openfast %gcc@6.3.0

To install OpenFAST with the C++ API, do:

.. code-block:: bash

    spack install openfast+cxx %gcc@6.3.0

That should be it! Spack will automatically use the most up-to-date
dependencies unless otherwise specified. For example to constrain OpenFAST
to use some specific versions of dependencies you could issue the Spack
install command:

.. code-block:: bash

    spack install openfast %gcc@6.3.0 ^hdf5@1.8.16

The executables and libraries will be located at

.. code-block:: bash

    spack location -i openfast

Add the appropriate paths to your ``PATH`` and ``LD_LIBRARY_PATH`` to run
OpenFAST.
