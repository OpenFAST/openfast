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

.. _development_philosophy:

Development Philosophy and Guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenFAST is intended to be a self-sustaining, community developed software.
While the NREL OpenFAST team serves as the gatekeeper of the repository, we
actively encourage the community to share new ideas and contribute code.
Considerations for contributing code are outlined here.

Engagement with NREL
--------------------

The process for community code contribution starts with engaging directly
with the NREL OpenFAST team to define the scope of the work and coordinate
development efforts. This is particularly important since many groups
work on OpenFAST simultaneously. By engaging early, all developers can
stay up to date and minimize conflicts during the code merge.
The prefered method of communication is `GitHub Issues <https://github.com/openfast/openfast/issues>`_.
An initial post should contain all relevant information about the planned
development work, the areas of the software that will be impacted,
and any model validation materials. See :ref:`development_plan`
for more information on describing the planned work.

The NREL OpenFAST team is always working on internal projects
that require the majority of our attention, but we will make every effort
to engage with the community and support development efforts in
a reasonable time frame. After posting an Issue, the NREL OpenFAST
team may reach out to schedule a meeting to talk through the details.

.. _development_plan:

Development Plan / Implementation Plan
--------------------------------------
Significant code development efforts at NREL begin with the development
of a detailed implementation plan, and a few such plans are available to
download for reference:

- :download:`Development Plan for the Aerodynamic Linearization of OpenFAST <../../OtherSupporting/AeroDyn/AeroLin_2019-12.pdf>`
- :download:`FAST.Farm Development Plan <../../OtherSupporting/FAST.Farm_Plan_Rev25.doc>`.
- :download:`Implementation Plan - 2nd-order Forces Within HydroDyn <../../OtherSupporting/HydroDyn/HydroDyn_2ndOrderForces_Plan.pdf>`
- :download:`Implementation Plan - 2nd-order Wave Kinematics Within HydroDyn <../../OtherSupporting/HydroDyn/WAVE2_document.pdf>`

A good plan within the modularization framework of OpenFAST will
follow the definitions and nomenclature used by the
:download:`NWTC Programmer's Handbook <../../OtherSupporting/NWTC_Programmers_Handbook.pdf>`.
It should communication the following information:

- State whether the module is intended for loose coupling,
  tight coupling for time marching, and/or linearization.
- Define the module's inputs (including initialization), 
  outputs (including initialization), states (continous,
  discrete, and constraint), and parameters, including units.
- Lay out an example input file for the module.
- Explain the module's mathemetical formulation, including
  Jacobians (for tight coupling and linearization), in the form
  required of the framework.
- Prescribe how the module's inputs are derived from the 
  outputs of other specific modules
- Identify any potential numerical problems and how to avoid 
  them in the code.
- Lay out the module's subroutines using pseudocode (as opposed 
  to actual code), including identifying which mathematical formulas
  are used by which subroutines, and describing the algorithms used
  in the solution process.

This information is very helpful since it is easier to review, iterate,
and agree on a plan before making changes to source code. Additionally,
an implementation plan will greatly aid in the programming effort and is
a useful starting point for writing the user and develop documentation.

Qualities of a good submission
------------------------------

Development efforts should include adequate testing throughout
the development process. When possible, new subroutines should include unit-level tests,
and the existing regression tests should be run periodically to ensure that
the full system behavior has not changed unintentionally. For new features,
additional regression tests should be added to cover the new code.
If the regression test results change in an expected manner, the baseline
results should be updated locally and in the `openfast/r-test <https://github.com/openfast/r-test>`_
repository. The `r-test README <https://github.com/openfast/r-test#updating-the-baselines>`_
describes updating the baselines and the :ref:`testing`
section in this documentation contains additional details on testing.

New code should consider robustness from both the developer and user
perspectives. Here are some questions to consider during
code development:

- Is it clear to other developers how to use your subroutine?
- Does your new code exhibit clear and predictable behavior?
- How will your code perform under different qualities of data?
- How does your code impact the performance of the simulation?

Additionally, user and developer documentation should be included
with new code. User documentation includes theory, modeling guidance,
and a description of any inputs and outputs. User documentation should be
included as part of the online documentation described in :ref:`build_doc`.
Developer documentation is typically included in comments in the source
code. This should describe subroutine API's (inputs and outputs) as well
as any algorithms or lines of code that are unclear. Ask yourself
what you would need to know to fully understand your code if you don't
see it again for two years.

Submit for review and NREL feedback
-----------------------------------

New code can be submitted for review from the NREL OpenFAST team by
opening a `pull request <https://github.com/openfast/openfast/pulls>`_
as described in :ref:`github_workflow`. We will review the code for
accuracy, validity, quality, and robustness. Reviewing open source
code contributions can be difficult, so it is worthwhile to review
your own code and consider what information would help you to
determine whether it is ready to merge.

The review process begins with simply ensuring that the automated
tests pass in `GitHub Actions <https://github.com/openfast/openfast/actions>`_.
**Please ensure that all automated tests pass prior to requesting a review.**
After that, the process will involve some communication between the
reviewer and the submitter, possibly requests for more information
on the background or validation, and comments in the pull request
to gain additional insight into specific lines of code.

After a consensus is reached between the submitter and reviewer,
the pull request will be merged into the target branch (typically
`dev`) and the pull request will be closed. You're done!
This change will be included in the subsequent release of OpenFAST
when the `dev` branch is merged into `main`.

Bug fixes
---------

If you've found a bug in the code, it is important to fully describe
it both in a `GitHub Issue <https://github.com/openfast/openfast/issues>`_
and through a minimal test. Before making a commit with the bug fix,
commit the new test that exposes the bug. This test should fail.
Then, commit the bug fix and show that the test passes. The git-commit
history should look something like this (progresses bottom to top):

.. mermaid::

  gitGraph BT:
  options
  {
    "nodeSpacing": 60,
    "nodeFillColor": "white",
    "nodeStrokeWidth": 2,
    "nodeStrokeColor": "#747474",
    "lineStrokeWidth": 2,
    "branchOffset": 30,
    "lineColor": "grey",
    "leftMargin": 20,
    "branchColors": ["#007bff", "#ff2d54"],
    "nodeRadius": 5,
    "nodeLabel": {
      "width": 75,
      "height": 100,
      "x": -25,
      "y": 0
    }
  }
  end

  commit
  branch dev
  checkout dev
  commit "Merge pull request #123"
  commit "Merge pull request #124"
  branch bugfix
  checkout bugfix
  commit "Add unit test exposing out of bounds error"
  commit "Fix out of bounds error in array"
  checkout dev
  commit "Merge pull request #125"
  merge bugfix

See :ref:`testing` and :ref:`github_workflow` for more information.

Additional guidance
-------------------

The following sections provide extended guidance on developing source code,
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




API Reference
~~~~~~~~~~~~~
Some subroutines and derived types throughout the source code have in-source
documentation which is compiled with Doxygen. Though this portion of the
documentation is always under development, the existing API reference can
be found in the following pages:

- `Main Page <../../html/index.html>`_
- `Index of Types <../../html/classes.html>`_
- `Source Files <../../html/files.html>`_

Other Documentation
~~~~~~~~~~~~~~~~~~~
Additional documentation exists that may be useful for developers seeking deeper
understanding of the solver and mathematics.

- :download:`NWTC Programmer's Handbook <../../OtherSupporting/NWTC_Programmers_Handbook.pdf>`
   This is an overview of programming guidelines for FAST 8. While some syntax and minor details have
   changed in OpenFAST, most of this guide is still relevant.
- :download:`OutListParameters.xlsx <../../OtherSupporting/OutListParameters.xlsx>`
   This Excel file contains the full list of outputs for each module.  It is used to generate the
   Fortran code for the output channel list handling for each module (this code is generally in
   the _IO.f90 files).  The MATLAB script available in the
   `matlab-toolbox <https://github.com/OpenFAST/matlab-toolbox>`__ repository at *Utilities/GetOutListParameters.m*.
