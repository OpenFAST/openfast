# OpenFAST

Updated on 2017-03-17

## Overview

OpenFAST is an open-source wind turbine simulation tool that was established with the FAST v8 code as its starting point (see below).  OpenFAST was created with the goal of being a community model, with developers and users from research laboratories, academia, and industry.  Our objective is to ensure that OpenFAST is sustainable software that is well tested and well documented.  **OpenFAST is under development**; our team at NREL is now establishing online documentation and automated unit/regression testing.  During this transition period, users can find FAST v8 documentation at <https://nwtc.nrel.gov/>.

## FAST v8

(modified from <https://nwtc.nrel.gov/FAST>)

FAST v8 is a computer-aided engineering tool for simulating the coupled dynamic response of wind turbines. FAST joins aerodynamics models, hydrodynamics models for offshore structures, control and electrical system (servo) dynamics models, and structural (elastic) dynamics models to enable coupled nonlinear aero-hydro-servo-elastic simulation in the time domain. The FAST tool enables the analysis of a range of wind turbine configurations, including two- or three-blade horizontal-axis rotor, pitch or stall regulation, rigid or teetering hub, upwind or downwind rotor, and lattice or tubular tower. The wind turbine can be modeled on land or offshore on fixed-bottom or floating substructures. FAST is based on advanced engineering models derived from fundamental laws, but with appropriate simplifications and assumptions, and supplemented where applicable with computational solutions and test data.

The aerodynamic models use wind-inflow data and solve for the rotor-wake effects and blade-element aerodynamic loads, including dynamic stall. The hydrodynamics models simulate the regular or irregular incident waves and currents and solve for the hydrostatic, radiation, diffraction, and viscous loads on the offshore substructure. The control and electrical system models simulate the controller logic, sensors, and actuators of the blade-pitch, generator-torque, nacelle-yaw, and other control devices, as well as the generator and power-converter components of the electrical drive. The structural-dynamics models apply the control and electrical system reactions, apply the aerodynamic and hydrodynamic loads, adds gravitational loads, and simulate the elasticity of the rotor, drivetrain, and support structure. Coupling between all models is achieved through a modular interface and coupler.

## OpenFAST Documentation 
We are creating a Sphinx-based documentation site at <http://openfast.readthedocs.io>.  

Documentation for FAST v8 and its modules may be found at <https://nwtc.nrel.gov/>, while we are building the new site.

## Obtaining OpenFAST

You are in the [right place](https://github.com/OpenFAST/OpenFAST)! For those not familiar with git and github, there are many resources, e.g.,

* <https://guides.github.com>
* <https://try.github.io>
* <https://help.github.com/categories/bootcamp/>
* <https://desktop.github.com/>

## Compiling, Using & Developing OpenFAST

Details for compiling, using, and developing OpenFAST on Linux-based and Windows machines are being established at <http://openfast.readthedocs.io>.

## OpenFAST Help

Please use [github issues](https://github.com/OpenFAST/OpenFAST/issues) to:

* ask usage questions,
* report bugs,
* request code enhancements.

For other questions regarding OpenFAST, please contact [Mike Sprague](mailto:michael.a.sprague@nrel.gov).

Users and developers may also be interested in the NREL National Wind Technology Center (NWTC) [phpBB Forum](https://wind.nrel.gov/forum/wind/).

## OpenFAST Support

OpenFAST is being maintained and developed by researchers and software engineers at the [National Renewable Energy Laboratory](http://www.nrel.gov/) (NREL), with support from the US Department of Energy's Wind Energy Technology Office.  NREL gratefully acknowledges development contributions from the following organizations:

* Envision Energy USA, Ltd
* Brigham Young University
* [Intel&reg; Parallel Computing Center (IPCC)](https://software.intel.com/en-us/ipcc)
