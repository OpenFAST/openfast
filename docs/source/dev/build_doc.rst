.. _build_doc:

Building this documentation locally
===================================

This document describes how to build the OpenFAST documentation on your local machine. The official documentation is automatically built
and updated on `readthedocs <http://openfast.readthedocs.io/en/latest/>`__ when new material is pushed to the github repo.
However, while developing documentation, it is helpful to build locally in order to see changes quickly and without needing
to publish your changes to the public facing site.

Dependencies
------------
The documentation is built in `Sphinx <http://www.sphinx-doc.org/en/master/>`__ with optional support for 
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`__, `Doxylink <https://pythonhosted.org/sphinxcontrib-doxylink/>`__, and
`Graphviz <http://www.graphviz.org>`__. Therefore users will need to install these tools as well as several extensions of Sphinx that are utilized.

Doxygen and Graphviz can be installed directly from their website or with a package manager like ``brew``, ``yum``, or ``apt``.

The remaining tools are Python based and should be installed with `pip` using the requirements file at
``docs/requirements.txt``.

With CMake and Make
-------------------
In the OpenFAST repository checkout, if it has not been created yet,
create a ``build`` directory.  Change
to the build directory and run CMake with ``BUILD_DOCUMENTATION`` on.  If all
of the required tools are found successfully, CMake will configure with the
ability to build the documentation.

Issue the command ``make docs`` which will first build the Doxygen
documentation and then the Sphinx documentation. If this completes
successfully, the entry point to the documentation will be in
``build/docs/html/index.html``.

For example, from the OpenFAST directory:

::

    mkdir build
    cd build
    cmake .. -DBUILD_DOCUMENTATION=ON
    make docs

If you modify document source files in ``OpenFAST/docs/source``, you can simply update the html files through another ``make docs`` in ``OpenFAST/build``:

::

    make docs

Pure python
-----------
If CMake and Make are not available on your system, the documentation can be generated directly
with `sphinx`. **Note: This method does not generate the API documentation through Doxygen.**

First, align your build structure to the standard OpenFAST build by creating a directory 
at ``openfast/build``.

If all tools are available, move into ``openfast/build`` and run the `sphinx` command:

::

    # sphinx-build -b <builder-name> <source-directory> <output-directory>
    sphinx-build -b html ../docs ./docs/html



Documentation Output
--------------------

After building the documentation, it can be accessed by opening the output in a browser.
Open the high level html file generated at ``openfast/build/docs/html/index.html``
and begin using the page as any other web page.


Additional Build Targets
------------------------

The html portion of the documentation can be built with ``make sphinx-html``, and
the output is available at ``openfast/build/docs/html/index.html``.

If LaTeX is installed, a pdf version of the documentation can be built with
``make sphinx-pdf``, and the output is available at ``openfast/build/docs/latex/Openfast.pdf``.
