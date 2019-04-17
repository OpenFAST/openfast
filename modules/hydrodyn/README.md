# HydroDyn Module
The legacy version of this module and additional documentation are available
the [NWTC Software Portal](https://nwtc.nrel.gov/HydroDyn/).

## Overview
HydroDyn is a time-domain hydrodynamics module that has been coupled into the
OpenFAST wind turbine computer-aided engineering (CAE) tool to enable
aero-hydro-servo-elastic simulation of offshore wind turbines. HydroDyn is
applicable to both fixed-bottom and floating offshore substructures. It can
also be driven as a standalone code to compute hydrodynamic loading uncoupled
from OpenFAST.

HydroDyn allows for multiple approaches for calculating the hydrodynamic loads
on a structure: a potential-flow theory solution, a strip-theory solution, or a
combination of the two. Waves in HydroDyn can be regular (periodic) or
irregular (stochastic) and long-crested (unidirectional) or short-crested (with
wave energy spread across a range of directions). HydroDyn treats waves using
first-order (linear Airy) or first- plus second-order wave theory with the
option to include directional spreading, but no wave stretching or higher order
wave theories are included. To minimize computational expense, Fast Fourier
Transforms are applied in the summation of all wave frequency components.
