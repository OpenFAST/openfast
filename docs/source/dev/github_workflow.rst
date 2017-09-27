.. _github_workflow:

Workflow for interacting with the OpenFAST github.com repo
==========================================================

OpenFAST development should follow "Git Flow" when interacting with the github repository.
Git Flow is a git workflow outlining safe methods of pushing and pulling commits
to a shared repository. Maintaining Git Flow is critical to prevent remote changes
from blocking your local development.

Git Flow
--------

The Git Flow process is strictly defined and adopted throughout the software development
community. It is detailed nicely `here <https://datasift.github.io/gitflow/IntroducingGitFlow.html>`__
and the chart below provides a high level perspective.

.. image:: ../../_static/GitFlowFeatureBranches.png
    :align: center
    
Reference: http://nvie.com/posts/a-successful-git-branching-model


OpenFAST Specific Git Flow
--------------------------

In OpenFAST development, the typical workflow follows this procedure

1. Fork OpenFAST/OpenFAST on GitHub

2. Clone your new fork: ``git clone https://github.com/youruser/OpenFAST``
  
3. Create a feature branch for active development: ``git branch feature/a_great_feature``
  
4. Add new development on feature/a_great_feature

5. Create a GitHub pull request to merge ``youruser/OpenFAST/feature/a_great_feature`` into ``OpenFAST/OpenFAST/dev``
  

.. _pull_requests:

Pull Requests
-------------

New pull requests should contain

- A description of the need for modifications

  - If the pull request fixes a bug, the accompanying GitHub issue should be referenced
 
- A highlight of the work implemented
- Regression test results

  - If all tests pass, the summary print out should be provided
  - If any tests fail, an explanation of the failing cases and supporting data like plots should be included 
  
- Updated unit tests, if applicable
- Updated documentation in applicable sections ready for compilation and deployment to `readthedocs <http://openfast.readthedocs.io>`__.
- A review request from a specific member of the NREL OpenFAST team



    