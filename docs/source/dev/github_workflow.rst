.. _github_workflow:

Working with OpenFAST on GitHub
===============================
The majority of the collaboration and development for OpenFAST takes place
on the `github repository <http://github.com/openfast/openfast>`__. There,
`issues <http://github.com/openfast/openfast/issues>`__ and
`pull requests <http://github.com/openfast/openfast/pulls>`__
are discussed and new versions are released. It is the best mechanism for
engaging with the NREL OpenFAST team and other developers throughout
the OpenFAST community.

Issues and work assignment
--------------------------
Issues should be opened with proper documentation and data to fully describe
the problem or feature gap. It is here that communication and coordination
should happen regarding ongoing work for new development, and developers should
make clear any intention to complete a task.

.. _pull_requests:

Pull Requests
-------------
When a code modification is ready for review, a pull request should be
submitted along with all appropriate documentation and tests. An NREL OpenFAST
team member will assign a reviewer and work with the  developer to have the
code merged into the main repository.

New pull requests should contain

- A description of the need for modifications

  - If the pull request fixes a bug,
    the accompanying GitHub issue should be referenced

- A highlight of the work implemented
- Regression test results

  - If all tests pass, the summary print out should be provided
  - If any tests fail, an explanation of the failing
    cases and supporting data like plots should be included

- Updated unit tests, if applicable
- Updated documentation in applicable sections ready for compilation and
  deployment to `readthedocs <http://openfast.readthedocs.io>`__.

Git workflow and interacting with the main repository
-----------------------------------------------------
OpenFAST development should follow "Git Flow" when interacting with the github
repository. Git Flow is a git workflow outlining safe methods of pushing and
pulling commits to a shared repository. Maintaining Git Flow is critical to
prevent remote changes from blocking your local development.

Git Flow
--------
The Git Flow process is well defined and adopted throughout the software
development community. It is detailed nicely
`here <http://nvie.com/posts/a-successful-git-branching-model>`__
and the chart below provides a high level perspective.

.. image:: ../../_static/GitFlowFeatureBranches.png
    :align: center

Reference: http://nvie.com/posts/a-successful-git-branching-model

OpenFAST Specific Git Flow
--------------------------
It is important to consider how your current work will be affected by other
developer's commits and how your commits will affect other developers.
On public branches, avoid using
`git rebase <https://git-scm.com/book/en/v2/Git-Branching-Rebasing>`__
and never `force push <https://git-scm.com/docs/git-push#git-push---force>`__.

In OpenFAST development, the typical workflow follows this procedure

1. Fork the OpenFAST/OpenFAST repository on GitHub

2. Clone your new fork: ``git clone https://github.com/<youruser>/OpenFAST``

3. Add OpenFAST/OpenFAST as a remote: ``git remote add upstream https://github.com/OpenFAST/OpenFAST``

4. Create a feature branch for active development:
``git branch feature/a_great_feature`` or
``git checkout -b feature/a_great_feature``

5. Add new development on `feature/a_great_feature`:
``git add a_file.f90 && git commit -m "A message" && git push``

5. Update your feature branch with OpenFAST/dev:
``git pull upstream dev && git push``

6. Create a GitHub pull request to merge
``youruser/OpenFAST/feature/a_great_feature`` into ``OpenFAST/OpenFAST/dev``

