.. _build_doc:

Developing Documentation
========================
OpenFAST documentation is hosted on
`readthedocs <http://openfast.readthedocs.io/>`_. It is automatically generated
through the readthedocs build system from both the ``master`` and ``dev``
branches whenever new commits are added. This documentation uses the
`restructured text <http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
markup language.

Building this documentation locally
-----------------------------------
The documentation is compiled with Sphinx, which is a Python based tool.
Install it and the other required Python packages listed in
``docs/requirements.txt`` with pip or another python package manager.

These packages are optional:
- `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`__
- `Doxylink <https://pythonhosted.org/sphinxcontrib-doxylink/>`__
- `Graphviz <http://www.graphviz.org>`__

Doxygen and Graphviz can be installed directly from their website or with a
package manager like ``brew``, ``yum``, or ``apt``.

Pure python build
-----------------
If CMake and Make are not available on your system, the documentation can
be generated directly with `sphinx`.
**Note: This method does not generate the API documentation through Doxygen.**

First, align your build structure to the standard OpenFAST build by creating
a directory  at ``openfast/build``.

If all tools are available, move into ``openfast/build`` and run the `sphinx`
command:

::

    # sphinx-build -b <builder-name> <source-directory> <output-directory>
    sphinx-build -b html ../docs ./docs/html

If this completes successfully, a html file will be created at
``build/docs/html/index.html`` which can be opened with any web browser.

Building with CMake and Make
----------------------------
In the OpenFAST directory, create a ``build`` directory and move into it.
Then, run CMake with this flag: ``-DBUILD_DOCUMENTATION=ON``.  If all
of the required tools are found successfully, CMake will configure with the
ability to build the documentation.

Next, run the command ``make docs`` which will first build the Doxygen
documentation and then the Sphinx documentation. If this completes
successfully, a html file will be created at ``build/docs/html/index.html``
which can be opened with any web browser.

For example, from the OpenFAST directory:

::

    mkdir build
    cd build
    cmake .. -DBUILD_DOCUMENTATION=ON
    make docs

If you modify document source files in ``openfast/docs/source``, you can simply
update the html files through another ``make docs`` in ``openfast/build``:

::

    make docs

Building only the html
~~~~~~~~~~~~~~~~~~~~~~
Generating the only html is much faster than compiling all components of the
documentation. This can be done with ``make sphinx-html``.

Building the PDF
~~~~~~~~~~~~~~~~
If LaTeX is installed, a pdf version of the documentation can be built with
``make sphinx-pdf``, and the output is available at
``OpenFAT/build/docs/latex/Openfast.pdf``.
