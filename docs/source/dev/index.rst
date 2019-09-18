.. _dev_guide:

Developer Documentation
=======================

**Our goal as developers is to ensure that OpenFAST is a sustainable open
source software that is well tested and well documented.** To that end, we
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
universities, NREL relies on these GitHub tools to coordinate efforts:

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
If the Doxygen documentation is available, it can be found at the following
locations:

- `Main Page <../../../doxygen/html/index.html>`_
- `Index of Classes <../../../doxygen/html/classes.html>`_
- `Files <../../../doxygen/html/files.html>`_

.. _development_philosophy:

Development Philosophy
~~~~~~~~~~~~~~~~~~~~~~

OpenFAST is intended to be a self sustaining community developed software.
A couple of tenets of this goal are that the code should be reasonably
straightforward to comprehend and manageable to improve. With that in mind, we
expect that new capabilities will include adequate testing and documentation.

We have the following guidance for developers:

- When fixing a bug, first introduce a unit test that exposes the bug, fix the
  bug, and submit a Pull Request. See :numref:`testing` and
  :numref:`github_workflow` for information on testing and the GitHub workflow.

- When adding a new feature, create appropriate automated unit and regression
  tests as described in :numref:`testing`. The objective is to create a GitHub
  pull request that provides adequate verification and validation such that the
  NREL OpenFAST developer team can merge the pull request with confidence that
  the new feature is "correct" and supports our goal of self-sustaining
  software. See :numref:`pull_requests` for detailed information on submitting
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
    debugging.rst
    performance.rst
    versioning.rst
