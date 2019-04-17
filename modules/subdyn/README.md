# SubDyn Module
The legacy version of this module and additional documentation are available
at the [NWTC Software Portal](https://nwtc.nrel.gov/SubDyn/).

## Overview
SubDyn is a time-domain structural-dynamics module for multi-member
fixed-bottom substructures that has been coupled into the OpenFAST
aero-hydro-servo-elastic computer-aided engineering (CAE) tool. Substructure
types supported by SubDyn include monopiles, tripods, jackets, and other
lattice-type substructures common for offshore wind installations in shallow
and transitional water depths. SubDyn can also be used to model lattice
support structures for land-based wind turbines.

SubDyn follows the requirements of the FAST modularization framework and
couples to OpenFAST. It can also be driven as a standalone code to compute
the mode shapes, natural frequencies, and time-domain responses of
substructures, uncoupled from OpenFAST and in the absence of external loading
other than gravity and interface motion.

SubDyn relies on two main engineering schematizations
1. a linear frame finite-element beam model (LFEB)
2. a dynamics system reduction via Craig-Bampton’s (C-B) method

together with a Static-Improvement method, greatly reducing the number of modes
needed to obtain an accurate solution.
