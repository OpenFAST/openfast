.. _Future-Work:

Future Work
===========

This first implementation phase focused on single-turbine capabilities,
fulfilling the basic requirements for the design of large and novel
rotor concepts. Future development work will turn toward the
implementation of features enabling multiple-turbine simulations on
medium-to-large-scale computational clusters. The reduction of the
computational time will also be of focus. This may be achieved using
tree techniques such as the fast multipole method. Further algorithmic
options, such as vortex amalgamation in the far wake, will be considered
to speed up the simulation. The framework presented in this manual is
compatible with grid-free or grid-based vortex particle formulations.
Such particle-based implementations will also be envisaged in the
future. Further validation of the code against measurements and
higher-order tools will be pursued. Applications to cases known to be
challenging for the BEM algorithm will also be investigated, such as
highly flexible rotors, offshore floating turbines, small-scale wind
farms, multiple-rotor turbines, or kites.

The following list contains future work on OLAF software:

-  Lagrangian particles

-  Multiple turbines, integration into FAST.Farm

-  Code speed-up

-  Dedicated dynamic stall model
