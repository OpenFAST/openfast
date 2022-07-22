# FEAMooring Module
The legacy version of this module and additional documentation are available
the [NWTC Software Portal](https://nwtc.nrel.gov/FEAMooring/).

## Overview
A TAMU-led team, in partnership with ABS and NREL, was selected for funding
under topic area 1.3 of DE-FOA-0000415 to develop an open-source
mooring-dynamics module for the analysis of floating offshore wind turbines.
This project led to the development of FEAMooring (or FEAM), a
finite-element-based mooring-dynamics module that has been coupled into
OpenFAST. FEAMooring can also be driven as a standalone code to compute mooring
dynamics uncoupled from OpenFAST.

FEAMooring calculates the mooring-line reaction forces at the fairlead
positions of the floating platform considering the mooring dynamics, including
inertia and drag forces at each line element. The module can analyze several
kinds of mooring systems, including catenary mooring, taut mooring, vertical
tendons, etc. provided that the proper line geometry and properties are
specified. While powerful, FEAMooring has the following limitations:

- Only single uniform moorings lines between the fairlead and anchor are
  considered; it is not possible to model multi-segmented lines, line
  interconnections, or weights/tanks
- Mooring-line bending stiffness is not considered
- Seabed friction is not considered
- Interface does not support coupling to HydroDyn; the added mass, drag, and
  buoyancy calculations assume still water
