.. _dev_philosophy:

OpenFAST development philosophy
================================

**Under construction**

OpenFAST and its underlying modules are mostly written in Fortran (adhering to the 2003 standard), but modules can be written in C/C++. 
OpenFAST was created with the goal of being a community model, with developers and users from research laboratories, academia, and industry. 

**Our goal as developers is to ensure that OpenFAST is sustainable software that is well tested and well documented.**
To that end, we are continually improving the documentation and test coverage for existing code, and we expect that new capabilities will include adequate testing and documentation.

Moving forward, we have the following guidance for developers:

- When fixing a bug, first introduce a unit test that exposes the bug, and then fix the bug.  
  See :numref:`testing`.
  

- When adding a new feature, 

- If a code modification changes regression-test results in an expected manner, work with the OpenFAST developer team to upgrade the regression-test suite via a ``New Issue`` on the `github.com <https://github.com/openfast/openfast/issues>`_ repository.
   

 
