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
 Module   Line  Flag Name        Example Value
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

========= ==== =============== ====================================================================================================
 Added in OpenFAST v1.0.
-----------------------------------------------------------------------------------------------------------------------------------
 Module   Line  Flag Name       Example Value
========= ==== =============== ====================================================================================================
 AeroDyn   12   CavityCheck     False         CavitCheck         - Perform cavitation check? (flag)
 AeroDyn   17   Patm            9999.9   Patm               - Atmospheric pressure (Pa) [used only when CavitCheck=True]
 AeroDyn   18   Pvap            9999.9   Pvap               - Vapour pressure of fluid (Pa) [used only when CavitCheck=True]       
 AeroDyn   19   FluidDepth      9999.9   FluidDepth         - Water depth above mid-hub height (m) [used only when CavitCheck=True]
========= ==== =============== ====================================================================================================

FAST v7 to FAST v8
------------------

A major restructuration occured between FAST7 and FAST8. The Matlab scripts `ConvertFAST7to8.m` in the folder `ConvertFASTversions` of the `matlab-toolbox` repository (https://github.com/OpenFAST/matlab-toolbox) performs part of the conversion.
In FAST7 a typical wind turbine model would only consist of 4 files: a main file (`.fst`), an aerodyn14 file, and elastodyn files for the blade and tower.
In FAST8, the files are splitted with: a main file (`.fst`) and one input file per module (plus the files included by these modules such as the elastodyn files for the blade and tower).

**Changes in AeroDyn**
The format of the main aerodyn file has changed significantly. The script `ConvertFAST7to8.m` should perform some these conversions.

The aerodynamic polar files also need to be converted from Aerodyn14 to AeroDyn15 format. The matlab script `WritePolarAD15.m`  from the `matlab-toolbox` can be used to convert polar data to AeroDyn15 format.

One major difference is in the definition of the blade radius in the Aerodyn file. In FAST7, the variable `RNodes` indicates the location of a blade station point along the blade starting from the rotor apex. In AeroDyn15 the variable `BlSpn` starts at the root of the blade and not at the rotor apex.

**Changes in ElastoDyn Blade and Tower files**

* The parameters `CalcTMode` and `CalcBMode` on line 5 of the ElastoDyn tower and blade file have been removed. 

* The distributed blade properties were modified in the ElastoDyn Blade file.
In FAST7 the properties (lines 16-end) were:

::

    --------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------
    BlFract  AeroCent  StrcTwst  BMassDen  FlpStff       EdgStff       GJStff        EAStff        Alpha  FlpIner  EdgIner  PrecrvRef  PreswpRef  FlpcgOf  EdgcgOf  FlpEAOf  EdgEAOf
    (-)      (-)       (deg)     (kg/m)    (Nm^2)        (Nm^2)        (Nm^2)        (N)           (-)    (kg m)   (kg m)   (m)        (m)        (m)      (m)      (m)      (m)

In FAST8 the properties (lines 14-end) are:

::

    --------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------
    BlFract      PitchAxis      StrcTwst       BMassDen        FlpStff        EdgStff
      (-)           (-)          (deg)          (kg/m)         (Nm^2)         (Nm^2)

The `PitchAxis` column has no effect on the aerodynamic calculations done by AeroDyn15 so far.
