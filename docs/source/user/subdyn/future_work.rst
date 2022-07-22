.. _sd_future-work: 

Known Limitations and Future Work
=================================


The following list contains known current limitations in the code:

-  Tight coupling is not yet supported.

-  Only nontapered two-node Euler-Bernoulli (E-B) or Timoshenko (T)
   element formulations are available. (In the future, tapered E-B and
   tapered Timoshenko element formulations will be implemented.)

-  Only straight circular members are permitted. (In the future, a
   generic cross section will be allowed.)

-  The number of elements per member (**NDiv**) is constant throughout
   the structure.

-  Internal matrices are not stored in sparse form, limiting the total
   number of possible nodes/DOFs to about 300/1800.

-  The dynamics system reduction is performed in the absence of external
   loading (e.g., hydrodynamic added mass).

-  Gravitational loading does not impact the global substructure
   stiffness.

-  Loads (gravitational, inertial, hydrodynamic) can only be applied as
   concentrated loads at element nodes; distributed loads (per unit
   length) are not yet supported.

-  The overlap of multiple members connected to a single joint is not
   modeled with super-elements.

-  Member-level outputs are only available for up to nine nodes of up to
   nine members (although the **OutAll** flag can generate further
   outputs at the member end nodes).

-  No graphics/animation capability is yet available to visualize the
   substructure geometry, modes, motion, and loads.
