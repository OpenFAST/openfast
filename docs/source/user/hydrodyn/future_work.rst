Future Work
===========

This list contains features that could be implemented in future
releases:

*  Enable times-series input motions for Morison members in the
   standalone HydroDyn (**MorisonInputsMod** = 2)

*  Enable tight-coupling to FAST, including linearization.

*  Enable wave stretching (**WaveStMod** > 0).

*  Enable full support for floating platform force flags.

*  Enable joint overlap calculations (**JointOvrlp** = 1).

*  Enable auto-generation of all possible output channels (**OutAll** =
   TRUE).

*  Add outputs pertaining to the total hydrodynamic applied loads at
   nodes along members and at joints.

*  Ensure that the output channels are written in the order they are
   entered.

*  Allow for a WAMIT reference point location other than (0,0,0).

*  Allow **RdtnDT** to be independent from the FAST simulation time
   step.

*  Add distributed axial viscous-drag loads on tapered members.

*  Add rotational inertia terms for fluid-filled members and marine
   growth.

*  Calculate the effective 6x6 added-mass matrix from strip-theory
   members and place in the HydroDyn summary file.

*  Add graphics/animation capability to visualize the substructure
   geometry and motion, wave elevation, and hydrodynamic loads.

*  Add convective fluid acceleration terms.

*  Allow for wave directional spreading to include energy spectra that
   varies with direction (requires changing from the equal-energy
   method).

*  Add higher order regular wave kinematics models for fixed-bottom
   substructures.

*  Add breaking wave-impact loads for fixed-bottom substructures.

*  Add floating platform hydro-elastics.

*  Add pressure mapping for floating platforms.

*  Added automated computation and use of hydrostatic restoring matrix
   for strip-theory members.
