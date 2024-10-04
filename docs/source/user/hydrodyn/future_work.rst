Future Work
===========

This list contains features that could be implemented in future
releases:

*  Enable times-series input motions for Morison members in the
   standalone HydroDyn (**MorisonInputsMod** = 2)

*  Enable tight-coupling to FAST, including linearization.

*  Enable full support for floating platform force flags.

*  Enable joint overlap calculations (**JointOvrlp** = 1).

*  Enable auto-generation of all possible output channels (**OutAll** =
   TRUE).

*  Add outputs pertaining to the total hydrodynamic applied loads at
   nodes along members and at joints.

*  Ensure that the output channels are written in the order they are
   entered.

*  Allow **RdtnDT** to be independent from the FAST simulation time
   step.

*  Add rotational inertia terms for fluid-filled members and marine
   growth.

*  Add convective fluid acceleration terms.

*  Allow for wave directional spreading to include energy spectra that
   varies with direction (requires changing from the equal-energy
   method).

*  Add higher order regular wave kinematics models for fixed-bottom
   substructures.

*  Add breaking wave-impact loads for fixed-bottom substructures.

*  Add pressure mapping for floating platforms.

*  Added automated computation and use of hydrostatic restoring matrix
   for strip-theory members.
