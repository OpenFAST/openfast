.. _overview:

Overview
========

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

