FAST v8 and the transition to OpenFAST 
======================================

This page describes the transition from FAST v8, a computer-aided engineering tool for simulating the coupled dynamic response of wind turbines, to OpenFAST. OpenFAST was established by researchers at the National Renewable Energy Laboratory (NREL) in 2017, who were supported by the U.S. Department of Energy Wind Energy Technology Office (DOE-WETO). 

FAST v8
-------

FAST v8 is a computer-aided engineering tool for simulating the coupled dynamic response of wind turbines. FAST joins aerodynamics models, hydrodynamics models for offshore structures, control and electrical system (servo) dynamics models, and structural (elastic) dynamics models to enable coupled nonlinear aero-hydro-servo-elastic simulation in the time domain. The FAST tool enables the analysis of a range of wind turbine configurations, including two- or three-blade horizontal-axis rotor, pitch or stall regulation, rigid or teetering hub, upwind or downwind rotor, and lattice or tubular tower. The wind turbine can be modeled on land or offshore on fixed-bottom or floating substructures. FAST is based on advanced engineering models derived from fundamental laws, but with appropriate simplifications and assumptions, and supplemented where applicable with computational solutions and test data.

The aerodynamic models use wind-inflow data and solve for the rotor-wake effects and blade-element aerodynamic loads, including dynamic stall. The hydrodynamics models simulate the regular or irregular incident waves and currents and solve for the hydrostatic, radiation, diffraction, and viscous loads on the offshore substructure. The control and electrical system models simulate the controller logic, sensors, and actuators of the blade-pitch, generator-torque, nacelle-yaw, and other control devices, as well as the generator and power-converter components of the electrical drive. The structural-dynamics models apply the control and electrical system reactions, apply the aerodynamic and hydrodynamic loads, adds gravitational loads, and simulate the elasticity of the rotor, drivetrain, and support structure. Coupling between all models is achieved through a modular interface and coupler.

Transition to OpenFAST
----------------------

The release of OpenFAST v0.1.0 represents a transition to better support an open-source developer community across research laboratories, industry, and academia around FAST-based aero-hydro-servo-elastic engineering models of wind-turbines and wind-plants. OpenFAST aims to provide a solid software-engineering framework for FAST development including well documented source code, extensive automated regression and unit testing, and a robust multi-platform and compiler build system.

Algorithmically, OpenFAST v0.1.0 is the release most closely related to the last release of FAST,  version 8.16.  OpenFAST v0.1.0 includes the following organizational changes relative to FAST v8.16:

* A new github organization (https://github.com/openfast) has been established;

* The FAST program has been renamed OpenFAST;

* The OpenFAST glue codes, modules, module drivers, and compiling tools are contained within one repository: https://github.com/openfast/openfast;

* The OpenFAST regression test baseline solutions (formerly the Certification Tests or CertTest) reside in a standalone repository: https://github.com/openfast/r-test;

* Scripts and documentation are now provided for compiling OpenFAST using CMake on Mac, Linux, and Cygwin (Windows) systems;

* Visual Studio Projects (VS-Build) are provided for compiling OpenFAST on Windows, but the development team is working to automate the generation of Visual Studio build files via CMake in future releases;

* Github issues (https://github.com/openfast/openfast/issues) has been introduced for developers to report and track bugs, request feature enhancements, and to ask questions related to the source code, compiling, and regression/unit testing. General user-related questions on OpenFAST theory and usage should still be handled through the forum at https://wind.nrel.gov/forum/wind;

* During the transition, most user-related documentation is still provided through the NWTC Information Portal, https://nwtc.nrel.gov;

* Version numbering has been updated for OpenFAST, e.g., OpenFAST-v1.0.0-123-gabcd1234-dirty, where

  - v1.0.0 is the major-minor-bugfix numbering sytem; this corresponds to a tagged commit made by NREL on github;

  - 123-g is the number of additional commits after the most recent tag for a build [the '-g' is for 'git'];

  - abcd1234 is the first 8 characters of the current commit hash;

  - dirty denotes that local changes have been made but not committed;

* Because all modules are contained in the same repository, the version numbers of each module have been eliminated and now use the OpenFAST version number; old documentation may still refer to old version numbers;

The AeroDyn v15 aerodynamics module has been significantly updated.  The blade-element/momentum theory (BEMT) solution algorithm has been improved as follows:

* BEMT now functions for the case where the undisturbed velocity along the x-direction of the local blade coordinate system (Vx) is less than zero;

*   BEMT no longer aborts when a valid value of the inflow angle (ϕ) cannot be found; in this case, the inflow angle is computed geometrically (without induction);

*   The inflow angle (ϕ) is now initialized on the first call instead of using ϕ= 0, giving better results during simulation start up;

*   When hub- and/or tip-loss are enabled (HubLoss = True and/or TipLoss = True), tangential induction (a’) is set to 0 instead of -1 at the root and/or tip, respectively (axial induction (a) is still set to 1 at the root and/or tip);

*   Made the BEMT solution more efficient.
In addition, several bugs in AeroDyn v15 have been fixed, including:

*   Fixed a bug whereby when hub- and/or tip-loss are enabled (HubLoss = True and/or TipLoss = True) along with the Pitt/Peters skewed-wake correction (SkewMod = 2), BEMT no longer modifies the induction factors at the hub and/or tip, respectively;

*   Fixed a bug whereby the time series was affected after the linearization analysis with AeroDyn coupled to OpenFAST when frozen wake is enabled (FrozenWake = True);

*   Fixed a bug in the calculation of wind loads on the tower whereby the tower displacement was used in place of the tower velocity; also, tower strikes detected by AeroDyn to calculate the influence of the tower on the wind local to the blade are now treated as fatal errors instead of severe errors;

*   Fixed minor bugs in the unsteady airfoil aerodynamics model.

Updates to BeamDyn:

* Extensive cleanup of the source code

* Bug fixes

  - Trapezoidal points are now correctly defined by blade stations instead of key points

  - An off-diagonal term in the structural damping-induced stiffness (i.e., representing a change in the damping force with beam displacement) has been corrected

Other updates:

- A new module for user-specified platform loading (ExtPtfm) has been introduced. ExtPtfm allows the user to specify 6x6 added mass, damping, and stiffness matrices, as well as a 6x1 load vector to define loads to be applied to ElastoDyn’s tower base/platform, e.g., to support the modeling of substructures or foundations through a super-element representation (with super-element derived from external software). ExtPtfm also provides the user with a module to customize with more advanced platform applied loads. Module ExtPtfm can be enabled by setting CompSub to 2 in the FAST primary input file (a new option) and setting SubFile to the name of the file containing the platform matrices and load time history, but setting CompSub to 2 requires one to disable hydrodynamics (by setting CompHydro to 0). Please note that the introduction of option 2 for CompSub represents a minor input file change (the only input file change in OpenFAST v0.1.0), but the MATLAB conversion scripts have not yet been updated.

- A bug has been fixed in the Line2-to-Point mapping of loads. Previously, the augmented mesh was being formed using an incorrect projection, thus causing strange transfer of loads in certain cases. This could cause issues in the coupling between ElastoDyn and AeroDyn v15 and/or in the coupling between HydroDyn and SubDyn.

- In the ServoDyn control and electrical-system module, the units and sign of output parameter YawMom have been corrected.

- Minor fixes were made to the error checking in ElastoDyn and ServoDyn.

- The interface between FAST and CFD wrappers, e.g., SOWFA has been modified.


OpenFAST: Looking forward
-------------------------

Our goal is to continually improve OpenFAST documentation and to increase the coverage of automated unit and regression testing. 
In order to increase testing coverage and to maintain robust software, we will require  that

* new modules be equipped by the module developer(s) with sufficient module-specific unit and regression testing along with appropriate OpenFAST regression tests; 

* bug fixes include appropriate unit tests;

* new features/capabilities include appropriate unit and regression tests.
We are in the process of better instrumenting the BeamDyn module with extensive testing as a demonstration of requirements for new modules.   

For unit testing, we will employ the pFUnit framework (https://sourceforge.net/projects/pfunit).

For the time being OpenFAST provides project and solution files to support users developing and compiling using Visual Studio. However, the team is continually working to automate the generation of Visual Studio build files via CMake in future releases. 

Please contact `Michael.A.Sprague@NREL.gov <mailto:Michael.A.Sprague@NREL.gov>`_ with questions regarding the OpenFAST
development plan.
