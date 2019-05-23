# BeamDyn Module
The legacy version of this module and additional documentation are available
the [NWTC Software Portal](https://nwtc.nrel.gov/BeamDyn/).

## Overview
BeamDyn is a time-domain structural-dynamics module for slender structures.
The module has been coupled into the FAST aero-hydro-servo-elastic wind turbine
multi-physics engineering tool where it used to model blade structural
dynamics. BeamDyn is designed to analyze beams that are made of composite
materials, initially curved and twisted, and subject to large displacement and
rotation deformations. BeamDyn can also be used for static analysis of beams.

BeamDyn is based on the geometrically exact beam theory (GEBT) and is
implemented using Legendre spectral finite elements (LSFEs). GEBT supports full
geometric nonlinearity and large deflection, with bending, torsion, shear, and
extensional degree-of-freedom (DOFs); anisotropic composite material couplings
(using full 6x6 mass and stiffness matrices, including bend-twist coupling);
and a reference axis that permits blades that are not straight (supporting
built-in curve, sweep, and sectional offsets).

BeamDyn was originally developed by Qi Wang, Mike Sprague, and Jason Jonkman
under an NREL Laboratory Directed Research and Development (LDRD)
Project: High-fidelity computational modeling of wind-turbine structural
dynamics (2011-2013; PI: M.A. Sprague).
Further development was funded by the DOE Wind Energy Technology Office. NREL
gratefully acknowledges the development contributions from the following organizations:
* [Envision Energy USA, Ltd](http://www.envision-energy.com)

## Manual
BeamDyn documentation is available on the OpenFAST
[ReadTheDocs](https://openfast.readthedocs.io/en/master/source/user/beamdyn/index.html) site.

## Relevant Publications

* Wang, Q., M.A. Sprague, J. Jonkman, N. Johnson, and B. Jonkman, 2017, BeamDyn: An efficient high-fidelity wind turbine blade solver in the FAST modular framework.
*Wind Energy*, [DOI:10.1002/we.2101](http://onlinelibrary.wiley.com/doi/10.1002/we.2101/full).

* Wang, Q., M. Sprague, J. Jonkman, and B. Jonkman, Partitioned nonlinear structural analysis of wind turbines using BeamDyn.
Proceedings of *34th Wind Energy Symposium, AIAA Science and Technology  Forum and Exposition 2016*, San Diego, California, 4--8 January 2016.

* Guntur, S., J. Jonkman, S. Schreck, B. Jonkman, Q. Wang, M. Sprague, M. Hind, and R. Sievers,
FAST v8 Verification and Validation Using Experiments from Aeroelastically Tailored Megawatt-Scale Wind Turbine Blades.
Proceedings of *34th Wind Energy Symposium, AIAA Science and Technology  Forum and Exposition 2016*, San Diego, California, 4--8 January 2016.  Also published as [NREL/CP-5000-65389](http://www.nrel.gov/docs/fy16osti/65389.pdf).

* C. Pavese, T. Kim, Q. Wang, J. Jonkman, M.A. Sprague,
HAWC2 and BeamDyn: Comparison Between Beam Structural Models for Aero-Servo-Elastic Frameworks, Proceedings of the *European Wind Energy Association Annual Conference and Exhibition 2015* (EWEA 2015), 17-20 November 2015, Paris, France pp. 1193-1201. Also published as [NREL/CP-5000-65115](http://www.nrel.gov/docs/fy16osti/65115.pdf).

* Wang, Q., N. Johnson, M.A. Sprague, J. Jonkman, BeamDyn:
A high-fidelity wind turbine blade solver in the FAST modular framework.
Proceedings of the *AIAA Science and Technology Forum and Exposition,
33rd ASME Wind Energy Symposium*, Kissimmee, FL, 5--9 January 2015. Also published as [NREL/CP-5000-63165](http://www.nrel.gov/docs/fy15osti/63165.pdf).

* Wang, Q., M.A. Sprague, and J.M. Jonkman,  Nonlinear
Legendre spectral finite elements for wind turbine blade dynamics.
Proceedings of the *AIAA Science and Technology Forum and
Exposition, 32nd ASME Wind Energy Symposium*,  National Harbor, MD, January 13--17, 2014. Also published as [NREL/CP-2C00-60759](http://www.nrel.gov/docs/fy14osti/60759.pdf).
