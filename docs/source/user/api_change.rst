.. _api_change:

API changes between versions
============================

This page lists the main changes in the OpenFAST API (input files) between different versions.

The changes are tabulated according to the module input file, line number, and flag name.
The line number corresponds to the resulting line number after all changes are implemented.
Thus, be sure to implement each in order so that subsequent line numbers are correct.

OpenFAST v2.4.0 to OpenFAST `dev`
---------------------------------

============== ==== ================== =============================================================================================================================================================================
Added in OpenFAST `dev`
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module        Line  Flag Name          Example Value
============== ==== ================== =============================================================================================================================================================================
HydroDyn       53   ExctnMod                0   ExctnMod   - Wave Excitation model {0: None, 1: DFT, 2: state-space} (-) 
OpenFAST       44   CalcSteady         true     CalcSteady - Calculate a steady-state periodic operating point before linearization? [unused if Linearize=False] (flag)
OpenFAST       45   TrimCase                3   TrimCase   - Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only if CalcSteady=True] (-)
OpenFAST       46   TrimTol            0.0001   TrimTol    - Tolerance for the rotational speed convergence [used only if CalcSteady=True] (-)
OpenFAST       47   TrimGain            0.001   TrimGain   - Proportional gain for the rotational speed error (>0) [used only if CalcSteady=True] (rad/(rad/s) for yaw or pitch; Nm/(rad/s) for torque)
OpenFAST       48   Twr_Kdmp                0   Twr_Kdmp   - Damping factor for the tower [used only if CalcSteady=True] (N/(m/s))
OpenFAST       49   Bld_Kdmp                0   Bld_Kdmp   - Damping factor for the blades [used only if CalcSteady=True] (N/(m/s))
InflowWind     48   InitPosition(x)       0.0   InitPosition(x) - Initial offset in +x direction (shift of wind box) [Only used with WindType = 5] (m)
AeroDyn        13   CompAA             False                   CompAA             - Flag to compute AeroAcoustics calculation [only used when WakeMod=1 or 2]
AeroDyn        14   AA_InputFile       "unused"                AA_InputFile       - Aeroacoustics input file
AeroDyn        35   [separator line]   ======  OLAF cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when WakeMod=3]
AeroDyn        36   OLAFInputFileName  "Elliptic_OLAF.dat"     OLAFInputFileName - Input file for OLAF [used only when WakeMod=3]
AirFoilTables  4\*  BL_file            "unused"                BL_file           - The file name including the boundary layer characteristics of the profile. Ignored if the aeroacoustic module is not called.
IfW driver     6    [separator line]   ===================== File Conversion Options =================================
IfW driver     7    WrHAWC             false    WrHAWC    - Convert all data to HAWC2 format? (flag)
IfW driver     8    WrBladed           false    WrBladed  - Convert all data to Bladed format? (flag)
IfW driver     9    WrVTK              false    WrVTK     - Convert all data to VTK format? (flag)
InflowWind     7    VFlowAng                0   VFlowAng  - Upflow angle (degrees) (not used for native Bladed format WindType=7)
============== ==== ================== =============================================================================================================================================================================

Modified in OpenFAST `dev`
--------------------------

============== ==== ================== =============================================================================================================================================================================
 Module        Line  Flag Name          Example Value
============== ==== ================== =============================================================================================================================================================================
AirFoilTables  40\* filtCutOff         "DEFAULT"               filtCutOff        - Reduced frequency cut-off for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (-) [default = 0.5] 
============== ==== ================== =============================================================================================================================================================================

\*non-comment line count, excluding lines contained if NumCoords is not 0.

Additional nodal output channels added for :ref:`AeroDyn15<AD-Nodal-Outputs>`,
:ref:`BeamDyn<BD-Nodal-Outputs>`, and :ref:`ElastoDyn<ED-Nodal-Outputs>`.

OpenFAST v2.2.0 to OpenFAST v2.3.0
----------------------------------

============================================= ==== =============== ========================================================================================================================================================================================================
Removed in OpenFAST v2.3.0
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                                        Line  Flag Name        Example Value
============================================= ==== =============== ========================================================================================================================================================================================================
AeroDyn Airfoil Input File - Airfoil Tables   2    Ctrl            0   Ctrl              ! Control setting (must be 0 for current AirfoilInfo)
============================================= ==== =============== ========================================================================================================================================================================================================


============================================= ==== =============== ========================================================================================================================================================================================================
Added in OpenFAST v2.3.0
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                                        Line  Flag Name        Example Value
============================================= ==== =============== ========================================================================================================================================================================================================
AeroDyn Airfoil Input File - Airfoil Tables   2    UserProp         0   UserProp          ! User property (control) setting
AeroDyn                                       37   AFTabMod         1   AFTabMod          - Interpolation method for multiple airfoil tables {1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp} (-)
============================================= ==== =============== ========================================================================================================================================================================================================


OpenFAST v2.1.0 to OpenFAST v2.2.0
----------------------------------

No changes required.


OpenFAST v2.0.0 to OpenFAST v2.1.0
----------------------------------

============== ==== ================== =====================================================================================================================================================================
 Added in OpenFAST v2.1.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module        Line  Flag Name          Example Value
============== ==== ================== =====================================================================================================================================================================
BeamDyn driver 21   GlbRotBladeT0      True   GlbRotBladeT0 - Reference orientation for BeamDyn calculations is aligned with initial blade root?
============== ==== ================== =====================================================================================================================================================================

OpenFAST v1.0.0 to OpenFAST v2.0.0
----------------------------------

========= ==== ================== =====================================================================================================================================================================
Removed in OpenFAST v2.0.0
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module    Line Flag Name          Example Value
========= ==== ================== =====================================================================================================================================================================
BeamDyn    5   analysis_type      analysis_type  - 1: Static analysis; 2: Dynamic analysis
========= ==== ================== =====================================================================================================================================================================


========= ==== ================== =====================================================================================================================================================================
Added in OpenFAST v2.0.0
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module    Line Flag Name          Example Value
========= ==== ================== =====================================================================================================================================================================
AeroDyn   22   SkewModFactor      "default"     SkewModFactor    - Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0]
AeroDyn   30   Section header     ======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]
AeroDyn   31   DBEMT_Mod          2             DBEMT_Mod        - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
AeroDyn   32   tau1_const         4             tau1_const       - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]
BeamDyn    5   QuasiStaticInit    True          QuasiStaticInit  - Use quasi-static pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]
BeamDyn   11   load_retries       DEFAULT       load_retries     - Number of factored load retries before quitting the simulation
BeamDyn   14   tngt_stf_fd        DEFAULT       tngt_stf_fd      - Flag to use finite differenced tangent stiffness matrix (-)
BeamDyn   15   tngt_stf_comp      DEFAULT       tngt_stf_comp    - Flag to compare analytical finite differenced tangent stiffness matrix  (-)
BeamDyn   16   tngt_stf_pert      DEFAULT       tngt_stf_pert    - perturbation size for finite differencing (-)
BeamDyn   17   tngt_stf_difftol   DEFAULT       tngt_stf_difftol - Maximum allowable relative difference between analytical and fd tangent stiffness (-)
BeamDyn   18   RotStates          True          RotStates        - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]
========= ==== ================== =====================================================================================================================================================================



FAST v8.16 to OpenFAST v1.0.0
-----------------------------

The transition from FAST v8 to OpenFAST is described in detail at :ref:`fast_to_openfast`.

========== ==== =============== ====================================================================================================
Removed in OpenFAST v1.0.0
------------------------------------------------------------------------------------------------------------------------------------
Module     Line  Flag Name       Example Value
========== ==== =============== ====================================================================================================
OpenFAST   18   CompSub         0 CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}
========== ==== =============== ====================================================================================================


========== ==== =============== ====================================================================================================
Added in OpenFAST v1.0.0
------------------------------------------------------------------------------------------------------------------------------------
Module     Line  Flag Name       Example Value
========== ==== =============== ====================================================================================================
OpenFAST   18   CompSub         0 CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=External Platform MCKF}
AeroDyn    12   CavityCheck     False         CavitCheck         - Perform cavitation check? (flag)
AeroDyn    17   Patm            9999.9   Patm               - Atmospheric pressure (Pa) [used only when CavitCheck=True]
AeroDyn    18   Pvap            9999.9   Pvap               - Vapor pressure of fluid (Pa) [used only when CavitCheck=True]
AeroDyn    19   FluidDepth      9999.9   FluidDepth         - Water depth above mid-hub height (m) [used only when CavitCheck=True]
========== ==== =============== ====================================================================================================
