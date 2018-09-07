.. _dev_guide:

Developer Documentation
=======================

**Our goal as developers is to ensure that OpenFAST is a sustainable open source software that is well tested and well documented.**
To that end, we continually work to improve the documentation and test coverage along with feature additions and improvements.
This section of the documentation outlines the processes and procedures we have established for external developers
to work with the NREL OpenFAST team on code development.

Getting in touch
----------------
Please use `GitHub Issues <https://github.com/openfast/openfast/issues>`_ to:

- ask usage or development questions
- report bugs
- suggest code enhancements

For other questions regarding OpenFAST, please contact `Mike Sprague <mailto:michael.a.sprague@nrel.gov>`_.

Contributing to OpenFAST
------------------------
If you'd like to help with general OpenFAST development or work on a particular feature, then first install OpenFAST
following the :doc:`installation instructions <../install/index>` for your machine. Next, verify that your installation is 
valid by running the test suite following the :doc:`testing instructions <../testing/index>`.

After a successful and validated build, we encourage reading through the :doc:`OpenFAST development philosophy <dev_phil>` to
understand the general workflow for individual and coordinated development. Finally, be sure to review the :doc:`GitHub workflow <github_workflow>`
to avoid any merge or code conflicts.

Coordination
------------
With development happening in parallel between NREL, industry partners, and universities, duplicated effort is likely without
proper communication. In that regard, the NREL OpenFAST team maintains the GitHub
`Issues <https://github.com/openfast/openfast/issues>`_ and `Pull Request <https://github.com/openfast/openfast/pulls>`_ pages.
Any suggested bug fixes, improvements, or new features should be documented there.

Issues and work assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~
Issues should be opened with proper documentation and data to fully describe the problem or feature gap. It is here that
communication and coordination should happen regarding ongoing work for new development, and developers should make clear
any intention to complete a task.

Pull requests and reviews
~~~~~~~~~~~~~~~~~~~~~~~~~
When a code modification is ready for review, a pull request should be submitted along with all appropriate documentation and tests
as described in the :doc:`GitHub workflow <github_workflow>`. An NREL OpenFAST team member will assign a reviewer and work with the 
developer to have the code merged into the main repository.

Documentation
-------------
OpenFAST documentation is hosted on `readthedocs <http://openfast.readthedocs.io/>`_. It
is automatically generated from both the ``master`` and ``dev`` branches whenever new commits are added.
The documentation can also be generated locally as described in :doc:`build docs <build_doc>` while in development.

Note that a PDF of the documentation can be retrieved from `readthedocs <http://openfast.readthedocs.io/>`_ by clicking the arrow on the 
lower left corner of the page next to ``v:master`` or ``v:dev``.

While OpenFAST developer documentation is being enhanced here, developers are encouraged to
consult the legacy FAST v8 `Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`_.

Section Map
-----------
.. toctree::
   :maxdepth: 1

   dev_phil.rst
   github_workflow.rst
   build_doc.rst
   doxy_doc.rst

