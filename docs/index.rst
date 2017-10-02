.. OpenFAST documentation master file, created by
   sphinx-quickstart on Wed Jan 25 13:52:07 2017.

######################
OpenFAST Documentation
######################

.. only:: html

   :Version: |version|
   :Date: |today|

Welcome to the documentation site for OpenFAST. 

Overview
-----------------

OpenFAST is an open-source wind turbine simulation tool that was established in 2017 with the FAST v8 code as its starting point (see :ref:`fast_to_openfast`).  
OpenFAST is a multi-physics, multi-fidelity tool for simulating the coupled dynamic response of wind turbines.  
Practically speaking, OpenFAST is the framework (or glue code) that couples computational modules for
aerodynamics, hydrodynamics for offshore structures, control and electrical system (servo) dynamics, and structural dynamics to enable coupled nonlinear aero-hydro-servo-elastic simulation in the time domain. 
OpenFAST enables the analysis of a range of wind turbine configurations, including two- or three-blade horizontal-axis rotor, pitch or stall regulation, rigid or teetering hub, upwind or downwind rotor, and lattice or tubular tower. 
The wind turbine can be modeled on land or offshore on fixed-bottom or floating substructures. 

OpenFAST and its underlying modules are mostly written in Fortran (adhering to the 2003 standard), but modules can be written in C/C++.
OpenFAST was created with the goal of being a community model, with developers and users from research laboratories, academia, and industry. 
Our goal is also to ensure that OpenFAST is sustainable software that is well tested and well documented.   
To that end, we are continually improving the documentation and test coverage for existing code, and we expect that new capabilities will include adequate testing and documentation.

OpenFAST is under development; our team at NREL is now enhancing this documentation and automated unit/regression testing. 
During this transition period, users can find FAST v8 documentation at https://nwtc.nrel.gov/.

This documentation
------------------

OpenFAST documentation is built using `Sphinx <http://www.sphinx-doc.org>`_, which uses `reStructuredText <http://docutils.sourceforge.net/rst.html>`_ as its markup language. 
Online documentation is hosted on `readthedocs <http://openfast.readthedocs.io>`_, where one can choose between documentation generated from the OpenFAST `master` or `dev` github branches  (if viewing this on http://openfast.readthedocs.io click on  ``Read the Docs`` "box" on the lower left corner of the browser screen for options).  

This documentation is divided into two parts:

:ref:`user_guide`

   Directed towards end-users, this part provides detailed documentation
   regarding installation and usage of the OpenFAST and its underlying modules,
   as well as theory and verification documentation. 
   Also included are instructions for using the automated test suite, which 
   serves as a suite of examples.

:ref:`dev_guide`

   The developer guide is targeted towards users wishing to extend the
   functionality provided within OpenFAST. Here you will find details
   regarding the code structure, API supported by various classes, and links to
   source code documentation extracted using Doxygen.

**Note:** If viewing this on http://openfast.readthedocs.io,  one can get this documentation in PDF form via ``Read the Docs`` "box" on the lower left corner of the browser screen.  


.. toctree::
   :numbered:
   :maxdepth: 1

   source/get_started.rst
   source/install/index.rst
   source/testing/index.rst
   source/user/index.rst
   source/dev/index.rst
   source/links.rst
   source/license.rst
   source/help.rst
   source/acknowledgements.rst

   Nightly Testing Results <http://my.cdash.org/index.php?project=OpenFAST>
   github.com Repository <https://github.com/openfast/openfast>

