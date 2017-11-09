.. _build_doc:

Building this documentation locally
===================================

This document describes how to build the OpenFAST documentation on your local   machine.  Documentation is automatically built and updated on readthedocs when  new material is pushed to the github repo. However, while developing            documentation, one should build locally to see changes quickly, and without the need to push your changes to see them on readthedocs.

The documentation is based on the use of Doxygen, Sphinx,
and Doxylink. Therefore users will need to install these tools
as well as several extensions of Sphinx that are utilized.


Install the Tools
-----------------

Install CMake, Doxygen, Sphinx, Doxylink, and the
extensions used. Doxygen uses the ``dot`` application
installed with GraphViz. Sphinx uses a combination
of extensions installed with ``pip install`` as well as some
that come with OpenFAST located in the ``_extensions``
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
  pip install sphinx_rtd_theme

Run CMake Configure and Make the Docs
-------------------------------------

In the OpenFAST repository checkout, if it has not been created yet,
create a ``build`` directory.  Change
to the build directory and run CMake with ``BUILD_DOCUMENTATION`` on.  If all
of the main tools are found successfully, CMake should configure with the
ability to build the documentation. If Sphinx or Doxygen aren't found, the
configure will skip the documentation.

Issue the command ``make docs`` which should first build the Doxygen
documentation and then the Sphinx documentation. If this completes
successfully, the entry point to the documentation should be in
``build/docs/html/index.html``.

For example, from the OpenFAST directory:

::

    mkdir build
    cd build
    cmake -DBUILD_DOCUMENTATION:BOOL=ON ..
    make docs

If you modify document source files in ``OpenFAST/docs/source``, you can simply update the html files through another ``make docs`` in ``OpenFAST/build``:

::

    make docs

Documentation Output
--------------------

After building the documentation, it can be access by opening the output in a browser.
Open the high level html file generated at ``openfast/build/docs/html/index.html``
and begin using the page as any other web page.


Additional Build Targets
------------------------

The html portion of the documentation can be built with ``make sphinx-html``, and
the output is available at ``openfast/build/docs/html/index.html``.

If LaTeX is installed, a pdf version of the documentation can be built with
``make sphinx-pdf``, and the output is available at ``openfast/build/docs/latex/Openfast.pdf``.
