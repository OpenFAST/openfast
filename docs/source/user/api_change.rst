.. _api_change:

API changes between versions
============================

This page lists the main changes in the OpenFAST API (input files) between different versions.
For completeness, some changes in the previous versions of the FAST software are also included.

The changes are tabulated according to the module input file, line number, and flag name.
The line number corresponds to the resulting line number after all changes are implemented.
Thus, be sure to implement each in order so that subsequent line numbers are correct.

OpenFAST v2.0.0 to OpenFAST v2.1.0
----------------------------------

No changes required.

OpenFAST v1.0.0 to OpenFAST v2.0.0
----------------------------------

========= ==== =============== =====================================================================================================================================================================
Removed in OpenFAST v2.0.0
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module    Line  Flag Name        Example Value
========= ==== =============== =====================================================================================================================================================================
BeamDyn    5   analysis_type   analysis_type  - 1: Static analysis; 2: Dynamic analysis
========= ==== =============== =====================================================================================================================================================================


========= ==== ================== =====================================================================================================================================================================
 Added in OpenFAST v2.0.0
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module   Line  Flag Name          Example Value
========= ==== ================== =====================================================================================================================================================================
AeroDyn   22   SkewModFactor      "default"     SkewModFactor      - Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0]
AeroDyn   30   Section header     ======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]
AeroDyn   31   DBEMT_Mod          2   DBEMT_Mod          - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
AeroDyn   32   tau1_const         4   tau1_const         - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]
BeamDyn    5   QuasiStaticInit    True          QuasiStaticInit - Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]
BeamDyn   11   load_retries       DEFAULT       load_retries     - Number of factored load retries before quitting the aimulation
BeamDyn   14   tngt_stf_fd        DEFAULT       tngt_stf_fd      - Flag to use finite differenced tangent stiffness matrix (-)
BeamDyn   15   tngt_stf_comp      DEFAULT       tngt_stf_comp    - Flag to compare analytical finite differenced tangent stiffness matrix  (-)
BeamDyn   16   tngt_stf_pert      DEFAULT       tngt_stf_pert    - perturbation size for finite differencing (-)
BeamDyn   17   tngt_stf_difftol   DEFAULT       tngt_stf_difftol - Maximum allowable relative difference between analytical and fd tangent stiffness (-)
BeamDyn   18   RotStates          True          RotStates       - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]
========= ==== ================== =====================================================================================================================================================================

FAST v8.16 to OpenFAST v1.0.0
-----------------------------

The transition from FAST v8 to OpenFAST is described in detail at :ref:`fast_to_openfast`.

========== ==== =============== ====================================================================================================
Removed in OpenFAST v1.0.
------------------------------------------------------------------------------------------------------------------------------------
Module     Line  Flag Name       Example Value
========== ==== =============== ====================================================================================================
OpenFAST   18   CompSub         0 CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}
========== ==== =============== ====================================================================================================


========= ==== =============== ====================================================================================================
 Added in OpenFAST v1.0.
-----------------------------------------------------------------------------------------------------------------------------------
 Module   Line  Flag Name       Example Value
========= ==== =============== ====================================================================================================
OpenFAST  18   CompSub         0 CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=External Platform MCKF}
AeroDyn   12   CavityCheck     False         CavitCheck         - Perform cavitation check? (flag)
AeroDyn   17   Patm            9999.9   Patm               - Atmospheric pressure (Pa) [used only when CavitCheck=True]
AeroDyn   18   Pvap            9999.9   Pvap               - Vapour pressure of fluid (Pa) [used only when CavitCheck=True]
AeroDyn   19   FluidDepth      9999.9   FluidDepth         - Water depth above mid-hub height (m) [used only when CavitCheck=True]
========= ==== =============== ====================================================================================================
