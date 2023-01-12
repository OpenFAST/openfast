# MoorDyn Module
This is an externally developed module with further information
available on the developer's documentation site:
[Matt Hall](http://www.matt-hall.ca/moordyn.html).

## Overview
MoorDyn is a lumped-mass mooring line model for simulating the dynamics of
moorings connected to floating offshore structures. It accounts for internal
axial stiffness and damping forces, weight and buoyancy forces, hydrodynamic
forces from Morison's equation (assuming quiescent water so far), and vertical
spring-damper forces from contact with the seabed. MoorDyn's input file format
is based on that of MAP. The model supports arbitrary line interconnections,
clump weights and floats, and different line properties.

The Fortran implementation of MoorDyn, which has been developed
following the FAST Modularization Framework, is included as a module in
OpenFAST.

For the C++ implementation of MoorDyn, see http://www.matt-hall.ca/moordyn.
"MoorDyn C" can be compiled as a dynamically-linked library and features
simpler functions for easy coupling with models or scripts coded in C/C++,
Fortran, Matlab/Simulink, etc. It has recently been integrated into WEC-Sim.

Both forms of MoorDyn feature the same underlying mooring model, use similar
input and output conventions, and are being updated and improved in parallel.
They follow the same version numbering, with a "C" or "F" suffix for
differentiation.
