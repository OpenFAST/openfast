Building the Documentation
==========================

This document describes how to build Nalu's documentation.
The documentation is based on the use of Doxygen, Sphinx,
and Doxylink. Therefore we will need to install these tools
as well as many extensions of Sphinx that are utilized.

Install the Tools
-----------------

Install CMake, Doxygen, Sphinx, Doxylink, and the
extensions used. Doxygen uses the ``dot`` application
installed with GraphViz. Sphinx uses a combination
of extensions installed with ``pip install`` as well as some
that come with Nalu located in the ``_extensions``
directory. Using Homebrew on Mac OS X, 
this would look something like:

::

  brew install cmake
  brew install python
  brew install doxygen
  brew install graphviz
  pip install sphinx
  pip install sphinxcontrib-doxylink
  pip install sphinxcontrib-bibtex

Run CMake Configure
-------------------

In the Nalu repository checkout, create a ``mybuild`` directory.
Change to the build directory and run CMake with ``BUILD_DOCUMENTATION``
on. For example:

::

  cmake -DTrilinos_DIR:PATH=`spack location -i nalu-trilinos` \
  -DYAML_DIR:PATH=`spack location -i yaml-cpp` \
  -DENABLE_INSTALL:BOOL=ON -DCMAKE_BUILD_TYPE=RELEASE \
  -DBUILD_DOCUMENTATION:BOOL=ON \
  ..

If all of the main tools are found successfully, CMake should configure with the ability
to build the documentation. If Sphinx or Doxygen aren't found, the configure will skip
the documentation.


Make the Docs
-------------

Issue the command ``make docs`` which should first build the Doxygen documentation and
then the Sphinx documentation. If this completes successfully, the entry point to
the documentation should be in ``mybuild/docs/html/index.html``.
