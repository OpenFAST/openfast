.. _get_started:

Getting Started
===============

**Get the code:**

OpenFAST can be cloned (i.e., downloaded) from its `Github Repository <https:// github.com/OpenFAST/OpenFAST>`_, e.g., from the command line:
::

    git clone https://github.com/OpenFAST/OpenFAST.git

It can also be downloaded directly from https://github.com/OpenFAST/OpenFAST.git.

**Compile the code:**

See :ref:`installation`, for installation instructions (including dependency    requirements) for CMake, Spack, and Visual Studio and for multiple platforms    (Linux, Mac, Windows).
As an example, from the command line in a Mac or Linux environment:
::

    cd OpenFAST
    mkdir build && cd build
    cmake ../
    make

Note that one can see all of the `make` targets via
::

    make help


**Use the code:**

See :ref:`user_guide`, which is under construction.
In the interim, users may refer to the FAST v8 documentation at https://nwtc.nrel.gov/.

**Develop the code:**

See :ref:`dev_guide`, which is under construction.
In the interim, developers may consult the FAST v8 `Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`_.

