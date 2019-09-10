.. _dev_philosophy:

OpenFAST Development Philosophy
===============================

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

Versioning
----------
OpenFAST follows `semantic versioning <https://semver.org>`_. In summary, this
means that with a version number as MAJOR.MINOR.PATCH, the components will be
incremented as follows:

- MAJOR version when introducing incompatible API changes,
- MINOR version when adding functionality in a backwards-compatible manner, and
- PATCH version when making backwards-compatible bug fixes.

For example, ``OpenFAST-v1.0.0-123-gabcd1234-dirty`` describes OpenFAST as:

- v1.0.0 is the MAJOR.MINOR.PATCH numbering system and corresponds to a tagged
  commit made by NREL on GitHub
- 123-g is the number of additional commits after the most recent tag for a
  build [the ``-g`` is for ``git``]
- abcd1234 is the first 8 characters of the current commit hash
- dirty denotes that local changes have been made but not committed

Code Style
----------
OpenFAST and its underlying modules are mostly written in Fortran adhering to
the 2003 standard, but modules can be written in C or C++. Indentation is
typically three spaces and no tabs.

Generally, code should be written such that it is straightforward to read.
Syntactic sugar or brevity should not detract from readability. The exception
to this is in situations where performance dictates a poorly readable code.
Here, comment blocks should be used to describe what is not readily apparent
in the code.
