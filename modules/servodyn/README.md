# ServoDyn Module
The legacy version of TMD and additional documentation are available
at the [NWTC Software Portal](https://nwtc.nrel.gov/TMD/).

## Overview
ServoDyn is the control and Electrical Drive Dynamics Module for the
OpenFAST framework.

Included in ServoDyn is the tuned mass damper (TMD) module which adds
functionality to OpenFAST that simulates the addition of TMDs in the
nacelle and/or tower for structural control. The TMDs are two independent,
one-DOF, linear mass-spring-damping elements that act in the fore-aft and
side-side directions or one single omni-directional TMD. They can be placed
relative to the nacelle reference position or base of the undeflected tower
using the options in the input file. The TMD module is added as a sub-module
of ServoDyn.
