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
:doc:`testing instructions <../testing/index>`.

After a successful and validated build, we encourage reading through the
:doc:`OpenFAST development philosophy <dev_phil>` to understand the general
workflow for individual and coordinated development. Finally, be sure to review
the :doc:`GitHub workflow <github_workflow>` to avoid any merge or code
conflicts.

Please use `GitHub Issues <https://github.com/openfast/openfast/issues>`_ to:

- ask usage or development questions
- report bugs
- suggest code enhancements

The `GitHub Pull Requests <https://github.com/openfast/openfast/pulls>`_ page
is the place for engaging with the OpenFAST team to have your new code
merged into the main repository.

With development happening in parallel between NREL, industry partners, and
universities, NREL relies on these GitHub tools to coordinate efforts.

For other questions regarding OpenFAST, please contact
`Mike Sprague <mailto:michael.a.sprague@nrel.gov>`_.

Documentation
-------------
OpenFAST documentation is hosted on
`readthedocs <http://openfast.readthedocs.io/>`_. It is automatically generated
from both the ``master`` and ``dev`` branches whenever new commits are added.
The documentation can also be generated locally as described in
:doc:`build docs <build_doc>` while in development.

A PDF of the documentation can be retrieved from
`readthedocs <http://openfast.readthedocs.io/>`_ by clicking the arrow on the
lower left corner of the page next to ``v:master`` or ``v:dev``.

While OpenFAST developer documentation is being enhanced here,
developers are encouraged to consult the legacy FAST v8
`Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`_.

.. toctree::
   :maxdepth: 1

   dev_phil.rst
   github_workflow.rst
   build_doc.rst
   doxy_doc.rst
