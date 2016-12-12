# BeamDyn
BeamDyn is a time-domain structural-dynamics module for slender structures. The module has been coupled into the FAST aero-hydro-servo-elastic wind turbine multi-physics engineering tool where it used to model blade structural dynamics. BeamDyn is designed to analyze beams that are made of composite materials, initially curved and twisted, and subject to large displacement and rotation deformations. BeamDyn can also be used for static analysis of beams.

BeamDyn is based on the geometrically exact beam theory (GEBT) and is implemented using Legendre spectral finite elements (LSFEs). GEBT supports full geometric nonlinearity and large deflection, with bending, torsion, shear, and extensional degree-of-freedom (DOFs); anisotropic composite material couplings (using full 6x6 mass and stiffness matrices, including bend-twist coupling); and a reference axis that permits blades that are not straight (supporting built-in curve, sweep, and sectional offsets).

BeamDyn was developed by NREL under U.S. Department of Energy support.  NREL gratefully acknowledges the development contributions from the following organizations:
* Envision Energy USA, Ltd