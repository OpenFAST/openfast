.. _sd_appendix_F:

Appendix F. Major Changes in SubDyn
===================================

When first released, SubDyn (v0.4) was included as an undocumented
feature of FAST v8 and packaged as a stand-alone archive. Since v0.4,
SubDyn has been well integrated into FAST v8 and OpenFast, and the stand-alone form is
also available. This appendix outlines significant modifications to
SubDyn made since v0.4. Following are the main changes that the user may
notice, but for more information, refer to the *changelog.txt* text file
within the official archive and the GitHub log.



V1.05.00 (October 2021)
-----------------------

- Version 1.05.00 integrates with OpenFAST version 3.1+

- SubDyn driver supports loads at given nodes
  
- Outputs of Craig-Bampton/Guyan and FEM Modes to JSON format for visualization

- Streamlined yaml file output, with rigid body mass matrix at different points

- Bug fix for rigid assemblies.

- Bug fix for mass reported in summary file


V1.04.00 (September 2020)
-------------------------

- Version 1.04.00 integrates with OpenFAST version 2.4

- Member types: beam, rigid link, pretension cable  

- Joint types: cantilever, universal, pin, ball

- Input of all terms for concentrated mass

- Guyan damping matrix

- Extra lever arm

- Coupling sith SoilDyn

- Inclusion of soil-structure interaction (SSI) via flexible degrees of fixity at the restrained nodes and a new input file that allows for 6x6 stiffness and mass matrices that simulate boundary conditions at those nodes.

- Controllable pretension cable elements


V1.03.00a-rrd (September 2017)
------------------------------

- Version 1.03.00a-rrd integrates with the  `OpenFast software <https://github.com/OpenFAST/OpenFAST>`__.



V1.01.01a-rrd (September 2014)
------------------------------

Version 1.01.01a-rrd integrates with the `FAST v8
software <http://wind.nrel.gov/designcodes/simulators/fast8>`__
v8.09.00a-bjj.

-  Finite-element eigenvalue bug fixes: the full system eigenvalues were
   incorrectly reported in the summary file, although with no further
   consequences on the results. This bug is now fixed.

-  Shear area correction factor improvement: the shear area correction
   factor in the Timoshenko treatment is now aligned with Steinboeck et
   al. (2013).

-  The formulation for the TP reaction has been rearranged to adhere to
   the theory manual, with no consequences on the output results.


V1.01.00a-rrd (June 2014)
------------------------------

Version 1.00.01a-rrd integrates with the `FAST v8 software <http://wind.nrel.gov/designcodes/simulators/fast8>`__
v8.08.00c-bjj.

The new implementation has well-defined data exchange interfaces
(`following the FAST modularization
framework <http://wind.nrel.gov/designcodes/simulators/developers/>`__)
that should make integration of SubDyn into other multiphysics software
packages much simpler.

Several improvements and bug fixes have been implemented since version
v0.4 and the module has undergone an extensive verification effort with
good results.

-  Eigensolver bug fixes: the LAPACK solver proved to be unstable in
   single precision, and now it is solely run in double precision,
   regardless of the precision used in other parts of the code.

-  The input file format has changed. Please refer to the sample input
   file in Appendix A and the following notes:

   -  First header line has been removed.

   -  Simulation Control Section:

      -  **SDeltaT**: The "DEFAULT" keyword (in place of 0.0) is now
         used to indicate that the glue-code time step will be used for
         time integration by SubDyn.

      -  **IntMethod**: Allowed values are now 1-4 (in place of 0-3).

      -  **SttcSolve**: New flag introduced. If set to TRUE, the
         static-improvement method (SIM) will be used.

   -  FEA and Craig-Bampton Parameters Section:

      -  In v0.4, the damping coefficients had to be specified for all
         retained Craig-Bampton modes, or just once for all the modes
         (if **CBMod** = FALSE). In this version, the user can input
         any number of damping coefficients. In case the number of
         retained C-B modes (**NModes**) is larger than the input
         number of damping coefficients (**JDampings**), the last
         damping value will be replicated for all the remaining modes.

   -  Base Reaction Joints, Interface Joints, Member, and Member Cosine
      Matrices Sections:

      -  One line with units, below the headers, is expected in all the
         tables of the input file.

   -  Output: Summary and Outfile Section:

      -  This section now also contains the parameters previously
         assigned under the Section titled "Output: Fast/Subdyn
         Output-File Variables"

-  Some of the quantities in the summary file have been fixed. Some of
   the output matrices were, in fact, being output with wrong values
   because of an index mismatch. The new summary file is shorter and
   does not contain some of the CB method matrices, unless the compiler
   directive, DEBUG, is set.

-  SIM. This new implementation helps minimize the number of needed
   modes to capture the contribution of certain loads (such as static
   gravity and buoyancy loads or high-frequency loads transferred from
   the turbine). In the previous version, a large number of internal
   modes were needed to engage substructural modes excited by static and
   high-frequency forces. These modes are no longer needed and fewer
   modes can be retained while still achieving accurate results (see
   also :numref:`subdyn-theory`). With SIM enabled, all modes that are not considered
   by the Craig-Bampton reduction are treated quasi-statically.

-  There is now the possibility of retaining no internal C-B modes, thus
   relying solely only on SIM, in those cases where the substructure`s
   first eigenfrequencies are much higher than the expected
   energy-containing modes of the entire system.

-  The coupling of SubDyn within FAST now includes full hydro-elastic
   coupling with the HydroDyn hydrodynamics module.

