.. _build_doc:

Developing Documentation
========================
OpenFAST documentation is hosted on
`readthedocs <http://openfast.readthedocs.io/>`_. It is automatically generated
through the ``readthedocs`` build system from both the ``main`` and ``dev``
branches whenever new commits are added. This documentation uses the
`restructured text <http://www.sphinx-doc.org/en/main/usage/restructuredtext/basics.html>`_
markup language.

Building this documentation locally
-----------------------------------
The documentation is compiled with `Sphinx <http://sphinx-doc.org>`__, which is
a Python based tool. Install it and the other required Python packages listed
in ``openfast/docs/requirements.txt`` with ``pip`` or another Python package
manager.

These additional packages are optional and are not included in the requirements
file:

- `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`__
- `Doxylink <https://pythonhosted.org/sphinxcontrib-doxylink/>`__
- `Graphviz <http://www.graphviz.org>`__
- `LaTeX <https://www.latex-project.org>`__

Doxygen and Graphviz can be installed directly from their website or with a
package manager like ``brew``, ``yum``, or ``apt``.

The result of building the documentation locally will be a set of
HTML files and their accompanying required files. The main HTML file
will exist ``openfast/build/docs/html/index.html``. This file can
be opened with any browser to view and navigate the locally-generated
documentation as if it were any other web site.

Pure python build
~~~~~~~~~~~~~~~~~
If CMake and Make are not available on your system, the documentation can
be generated directly with `sphinx`.

.. note::

    This method does not generate the API documentation through Doxygen.

First, align your directory structure to the standard OpenFAST build by
creating a directory  at ``openfast/build``. Then, move into
``openfast/build`` and run this command:

.. code-block:: bash

    # sphinx-build -b <builder-name> <source-directory> <output-directory>
    sphinx-build -b html ../docs ./docs/html

If this completes successfully, an html file will be created at
``build/docs/html/index.html`` which can be opened with any web browser.

Building with CMake and Make
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the OpenFAST directory, create a ``build`` directory and move into it.
Then, run CMake with this flag: ``-DBUILD_DOCUMENTATION=ON``. CMake will
configure the build system with the necessary files for building
the documentation.

Next, run the command to compile the docs:

.. code-block:: bash

    make docs

This will first build the Doxygen API documentation and then the Sphinx
documentation. If this completes successfully, a html file will be
created at ``build/docs/html/index.html`` which can be opened with any web
browser.

The full procedure for configuring and building the documentation is:

.. code-block:: bash

    mkdir build
    cd build
    cmake .. -DBUILD_DOCUMENTATION=ON
    make docs

If any modifications are made to the source files in ``openfast/docs/source``,
you can simply update the html files by executing the ``make`` command again.

The table below lists make-targets related to the documentation.

======================= ================== ========================================
 Target                  Command            Output location
======================= ================== ========================================
 Full docs               make docs          openfast/build/docs/html/index.html
 Full docs               make sphinx        openfast/build/docs/html/index.html
 Doxygen API Reference   make doxygen
 HTML only               make sphinx-html   openfast/build/docs/html/index.html
 PDF only                make sphinx-pdf    openfast/build/docs/latex/Openfast.pdf
======================= ================== ========================================

Adding documentation
--------------------

Coming soon. Feel like contributing? Start here!
