.. _dev_guide:

Developer Documentation
=======================

**Our goal as developers of OpenFAST is to ensure that it is well tested, well
documented, and self-sustaining software.** To that end, we
continually work to improve the documentation and test coverage along with
feature additions and improvements. This section of the documentation outlines
the processes and procedures we have established for external developers
to work with the NREL OpenFAST team on code development.

If you'd like to help with general OpenFAST development or work on a particular
feature, then first install OpenFAST following the
:doc:`installation instructions <../install/index>` for your machine. Next,
verify that your installation is valid by running the test suite following the
:doc:`testing instructions <../testing/index>`. While OpenFAST is compiling, we
encourage reading through the :ref:`development_philosophy` section to
understand the general workflow for individual and coordinated development.
Finally, be sure to review the :doc:`GitHub workflow <github_workflow>` to
avoid any merge or code conflicts.

With development happening in parallel between NREL, industry partners, and
universities, NREL relies on GitHub to coordinate efforts:

- `GitHub Issues <https://github.com/openfast/openfast/issues>`_ is the place
  to ask usage or development questions, report bugs, and
  suggest code enhancements
- `GitHub Pull Requests <https://github.com/openfast/openfast/pulls>`_
  is the place for engaging with the OpenFAST team to have your new code
  merged into the main repository.

For other questions regarding OpenFAST, please contact
`Mike Sprague <mailto:michael.a.sprague@nrel.gov>`_.

.. tip::

    The following sections provide valuable guidance on workflow and
    development tips which make the process more efficient and
    effective:

    - :ref:`github_workflow`
    - :ref:`code_style`
    - :ref:`debugging`

API Reference
~~~~~~~~~~~~~
Some subroutines and derived types throughout the source code have in-source
documentation which is compiled with Doxygen. Though this portion of the
documentation is always under development, the existing API reference can
be found in the following pages:

- `Main Page <../../html/index.html>`_
- `Index of Types <../../html/classes.html>`_
- `Source Files <../../html/files.html>`_

.. _development_philosophy:

Development Philosophy
~~~~~~~~~~~~~~~~~~~~~~

OpenFAST is intended to be a self sustaining community developed software.
A couple of tenets of this goal are that the code should be reasonably
straightforward to comprehend and manageable to improve. With that in mind, we
expect that new capabilities will include adequate testing and documentation.

We have the following guidance for developers:

- When fixing a bug, first introduce a unit test that exposes the bug, fix the
  bug, and submit a Pull Request. See :ref:`testing` and
  :ref:`github_workflow` for more information.

- When adding a new feature, create appropriate automated unit and regression
  tests as described in :ref:`testing`. The objective is to create a GitHub
  pull request that provides adequate verification and validation so that the
  NREL OpenFAST developer team can merge the pull request with confidence that
  the new feature is "correct" and supports our goal of self-sustaining
  software. See :ref:`pull_requests` for more information on submitting
  a pull request.

- If a code modification affects regression test results in an expected manner,
  work with the NREL OpenFAST developer team to upgrade the regression test
  suite via a GitHub issue or pull request at the `openfast/r-test <https://github.com/openfast/r-test>`_
  repository.

Development Guidelines
~~~~~~~~~~~~~~~~~~~~~~
The following sections provide extended guidance on how to develop source code,
interacting with the NREL OpenFAST team and other community contributors, and
generally debugging and building out features.

.. toctree::
    :maxdepth: 1

    github_workflow.rst
    code_style.rst
    build_doc.rst
    types_files.rst
    debugging.rst
    performance.rst
    versioning.rst

Other Documentation
~~~~~~~~~~~~~~~~~~~
Additional documentation exists that may be useful for developers seeking deeper
understanding of the solver and mathematics.

- `NWTC Programmer’s Handbook <https://drive.google.com/file/d/1bDV1fBkiZUWs6Tkzb6nhCMUQvHpN_OtM/view?usp=sharing>`_
   This is an overview of programming guidelines for FAST 8. While some syntax and minor details have
   changed in OpenFAST, most of this guide is still relevant.
- :download:`OutListParameters.xlsx <../../OtherSupporting/OutListParameters.xlsx>`
   This Excel file contains the full list of outputs for each module.  It is used to generate the
   Fortran code for the output channel list handling for each module (this code is generally in
   the _IO.f90 files).  The MATLAB script available in the
   `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`__ repository at *Utilities/GetOutListParameters.m*.
