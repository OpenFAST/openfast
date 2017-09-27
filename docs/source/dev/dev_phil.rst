.. _dev_philosophy:

OpenFAST-development philosophy
================================

**Under construction**

OpenFAST and its underlying modules are mostly written in Fortran (adhering to the 2003 standard), but modules can be written in C/C++. 
OpenFAST was created with the goal of being a community model, with developers and users from research laboratories, academia, and industry. 

**Our goal as developers is to ensure that OpenFAST is sustainable software that is well tested and well documented.**
To that end, we are continually improving the documentation and test coverage for existing code, and we expect that new capabilities will include adequate testing and documentation.

Moving forward, we have the following guidance for developers:

- When fixing a bug, first introduce a unit test that exposes the bug, fix the bug, and submit a Pull Request.  
  See :numref:`testing` and :numref:`github_workflow`.

- When adding a new feature, create appropriate automated unit and regression tests as described in :numref:`testing`.  
  The objective is to create a github.com Pull Request that provides adequate automated testing such that the  OpenFAST core developer team can merge the pull request with confidence that the new feature is "correct" and supports our goal of self-sustaining software.
  Importantly, any new-feature Pull Request should have adequate documentation ready for compilation and posting on http://openfast.readthedocs.io.

- If a code modification affects regression-test results in an expected manner, work with the OpenFAST developer team to upgrade the regression-test suite via a ``New Issue`` on the `github.com <https://github.com/openfast/openfast/issues>`_ repository.
   
- While OpenFAST developer documentation is being enhanced here, developers are encouraged to consult the legacy FAST v8 `Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`_.


