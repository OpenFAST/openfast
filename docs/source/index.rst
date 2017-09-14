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
Online documentation is hosted on `readthedocs <http://openfast.readthedocs.io>`_, where one can choose between documentation generated from the OpenFAST `master` or `dev` github branches.  

This documentation is divided into two parts:

:ref:`user_guide`

   Directed towards end-users, this part provides detailed documentation
   regarding installation and usage of the various utilities available within
   this library. Here you will find a comprehensive listing of all available
   utilties, and information regarding their usage and current
   limitations that the users must be aware of.

:ref:`dev_guide`

   The developer guide is targeted towards users wishing to extend the
   functionality provided within this library. Here you will find details
   regarding the code structure, API supported by various classes, and links to
   source code documentation extracted using Doxygen.


Getting started
---------------

**Get the code:** 
OpenFAST can be cloned (i.e., downloaded) from its `Github Repository <https://github.com/OpenFAST/OpenFAST>`_, e.g., from the command line:
::

    git clone https://github.com/OpenFAST/OpenFAST.git

**Compile the code:** 
See :ref:`installation`, for installation instructions (including dependency requirements) for cmake, spack, and Visual Studio and for multiple platforms (Linux, Mac, Windows).
As an example, from the command line in a Mac or Linux environment:
::

    cd OpenFAST
    mkdir build && cd build
    cmake ../
    make

Note that one can see all of the `make` targets via
::

    make help


**Use the code:**
See :ref:`user_guide`, which is under construction.  
In the interim, users may refer to the FAST v8 documentation at https://nwtc.nrel.gov/. 

**Develop the code:** 
See :ref:`dev_guide`, which is under construction.
In the interim, developers may consult the FAST v8 `Programmer's Handbook <https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf>`_.

Licensing
---------

The OpenFAST software, including its underlying modules, are licensed under `Apache
License Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_ open-source
license.


Getting Help
------------

For possible bugs, enhancement requests, or code questions, please submit an issue at
the `OpenFAST Github repository <https://github.com/OpenFAST/OpenFAST>`__.

Users may find the established FAST v8 through the NWTC Information Portal:
https://nwtc.nrel.gov/

Please contact `Michael.A.Sprague@NREL.gov <mailto:Michael.A.Sprague@NREL.gov>`_. with questions regarding the OpenFAST development plan or how to contribute.


Acknowledgements
----------------

This software is developed and maintained by researchers at the `National Renewable Energy Laboratory <https://www.nrel.gov>`_ with funding from U.S. Department of Energy (DOE) Wind Energy Technology Office through the 
`Atmosphere to electrons (A2e) <https://a2e.energy.gov>`_ research initiative.

NREL gratefully acknowledges development contributions from the following organizations:

* Envision Energy USA, Ltd

* Brigham Young University

NREL gratefully acknowledges additional development support through an 
`Intel® Parallel Computing Center (IPCC) <https://software.intel.com/en-us/ipcc>`_

Important Links
---------------------------

*  OpenFAST Github Organization Page <https://github.com/OpenFAST>
*  OpenFAST Github Repository <https://github.com/OpenFAST/OpenFAST>

Documentation-Site Contents
---------------------------

.. toctree::
   :maxdepth: 1

   User documentation <user/index.rst>
   Developer documentation <dev/index.rst>

