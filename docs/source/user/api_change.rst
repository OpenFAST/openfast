.. _api_change:

API changes between versions
============================

This page lists the main chages in the OpenFAST API between versions of OpenFAST.
These API changes only manifest in input file changes.
For completeness, some of the changes in the previous versions of the FAST program are also included.

Migration from OpenFAST v2.0.0 to OpenFAST v2.1.0
-------------------------------------------------

No changes required.


Migration from OpenFAST v1.0.0 to OpenFAST v2.0.0
-------------------------------------------------

========= ==== ===============  =====================================================================================================================================================================
Added in v2.0.0
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module   Line  Flag Name        Example Value
========= ==== ===============  =====================================================================================================================================================================
 AeroDyn   22   SkewMod          "default"     SkewModFactor      - Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0]
 AeroDyn   30   DBEMT section    ======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]
 AeroDyn   31   DBEMT_Mod        2   DBEMT_Mod          - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
 AeroDyn   32   tau1_const       4   tau1_const         - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1] 
========= ==== ===============  =====================================================================================================================================================================

Removed in v2.0.0
None

Migration from FAST v8.16 to OpenFAST v1.0.0
--------------------------------------------

The transition from FAST v8 to OpenFAST is described here :ref:`fast_to_openfast`. 

The only difference in input files between FAST v8 and OpenFAST1.0 lays in the AeroDyn input file, the version of AeroDyn being 15.03 and 15.04 respectively.

The AeroDyn file for OpenFast 1.0 requires the following additions:

* the parameter `CavityCheck` is inserted on line 12, after the line for the parameter `FrozenWake`:

::

    False         CavitCheck         - Perform cavitation check? (flag)

* three lines for the `Patm`, `Pvap` and `FluidDepth` are inserted on line 17-19, after the line for the parameter `SpdSound`:

::

       9999.9   Patm               - Atmospheric pressure (Pa) [used only when CavitCheck=True]
       9999.9   Pvap               - Vapour pressure of fluid (Pa) [used only when CavitCheck=True]            
       9999.9   FluidDepth         - Water depth above mid-hub height (m) [used only when CavitCheck=True]


Migration from FAST v7 to FAST v8
---------------------------------

A major restructuration occured between FAST7 and FAST8. The Matlab scripts `ConvertFAST7to8.m` in the folder `ConvertFASTversions` of the `matlab-toolbox` repository (https://github.com/OpenFAST/matlab-toolbox) performs part of the conversion.
In FAST7 a typical wind turbine model would only consist of 4 files: a main file (`.fst`), an aerodyn14 file, and elastodyn files for the blade and tower.
In FAST8, the files are splitted with: a main file (`.fst`) and one input file per module (plus the files included by these modules such as the elastodyn files for the blade and tower).

**Changes in AeroDyn**
The format of the main aerodyn file has changed significantly. The script `ConvertFAST7to8.m` should perform some these conversions.

The aerodynamic polar files also need to be converted from Aerodyn14 to AeroDyn15 format. The matlab script `WritePolarAD15.m`, which converts the polars will be added in the near future in the `matlab-toolbox` repository.

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
