.. _api_change:

API changes between versions
============================

This page lists the main changes in the OpenFAST API (input files) between different versions.

The changes are tabulated according to the module input file, line number, and flag name.
The line number corresponds to the resulting line number after all changes are implemented.
Thus, be sure to implement each in order so that subsequent line numbers are correct.


OpenFAST v2.5.0 to OpenFAST `dev`
---------------------------------

Many changes were applied to SubDyn input file format. You may consult the following example:
:download:`(SubDyn's Input File) <./subdyn/examples/OC4_Jacket_SD_Input.dat>`: 
and the online SubDyn documentation.

-  ServoDyn

   -  The input file parser is updated to a keyword/value pair based input.
      Each entry must have a corresponding keyword with the same spelling as
      expected
   -  The TMD submodule of ServoDyn is replaced by an updated Structural Control
      module (StC) with updated capabilities and input file.

============================================= ==== =============== ========================================================================================================================================================================================================
OpenFAST v2.5.0 to OpenFAST dev 
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                                        Line  Flag Name        Example Value
============================================= ==== =============== ========================================================================================================================================================================================================
AeroDyn 15                                         TwrTi               0.0000000E+00  6.0000000E+00  1.0000000E+00  1.0000000E-01                 [additional column in *Tower Influence and Aerodynamics* table]
SubDyn                                         8   GuyanLoadCorr.      False   GuyanLoadCorection  - Include extra moment from lever arm at interface and rotate FEM for floating
SubDyn                                        15   GuyanDampMod        0       GuyanDampMod - Guyan damping {0=none, 1=Rayleigh Damping, 2=user specified 6x6 matrix}
SubDyn                                        16   RayleighDamp        0.001, 0.003   RayleighDamp - Mass and stiffness proportional damping  coefficients (Rayleigh Damping) [only if GuyanDampMod=1]
SubDyn                                        17   GuyanDampSize       6       GuyanDampSize - Guyan damping matrix size (square, 6x6) [only if GuyanDampMod=2]
SubDyn                                        18   GuyanDampMat        0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00 
SubDyn                                        -23  GuyanDampMat        0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00 
SubDyn                                        na   CablesSection       -------------------------- CABLE PROPERTIES  -------------------------------------
SubDyn                                        na   CablesSection       0   NCablePropSets   - Number of cable cable properties
SubDyn                                        na   CablesSection       PropSetID     EA          MatDens       T0 
SubDyn                                        na   CablesSection          (-)        (N)         (kg/m)        (N) 
SubDyn                                        na   RigidSection        ---------------------- RIGID LINK PROPERTIES ------------------------------------
SubDyn                                        na   RigidSection        0   NRigidPropSets - Number of rigid link properties
SubDyn                                        na   RigidSection        PropSetID   MatDens   
SubDyn                                        na   RigidSection          (-)       (kg/m)
HydroDyn                                      52   NBody              1   NBody          - Number of WAMIT bodies to be used (-) [>=1; only used when PotMod=1. If NBodyMod=1, the WAMIT data contains a vector of size 6*NBody x 1 and matrices of size 6*NBody x 6*NBody; if NBodyMod>1, there are NBody sets of WAMIT data each with a vector of size 6 x 1 and matrices of size 6 x 6]
HydroDyn                                      53   NBodyMod           1   NBodyMod       - Body coupling model {1: include coupling terms between each body and NBody in HydroDyn equals NBODY in WAMIT, 2: neglect coupling terms between each body and NBODY=1 with XBODY=0 in WAMIT, 3: Neglect coupling terms between each body and NBODY=1 with XBODY=/0 in WAMIT} (switch) [only used when PotMod=1]
ServoDyn                                      60   AeroControlSec     ---------------------- AERODYNAMIC FLOW CONTROL --------------------------------
ServoDyn                                      61   AfCmode                   0   AfCmode      - Airfoil control mode {0: none, 1: cosine wave cycle, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)
ServoDyn                                      61   AfC_Mean                  0   AfC_Mean     - Mean level for cosine cycling or steady value (-) [used only with AfCmode==1]
ServoDyn                                      61   AfC_Amp                   0   AfC_Amp      - Amplitude for for cosine cycling of flap signal (-) [used only with AfCmode==1]
ServoDyn                                      61   AfC_Phase                 0   AfC_Phase    - Phase relative to the blade azimuth (0 is vertical) for for cosine cycling of flap signal (deg) [used only with AfCmode==1]
ServoDyn                                      65   StCSection         ---------------------- STRUCTURAL CONTROL ---------------------------------------
ServoDyn                                      66   NumBStC                   0   NumBStC      - Number of blade structural controllers (integer)
ServoDyn                                      67   BStCfiles          "unused"   BStCfiles    - Name of the files for blade structural controllers (quoted strings) [unused when NumBStC==0]
ServoDyn                                      68   NumNStC                   0   NumNStC      - Number of nacelle structural controllers (integer)
ServoDyn                                      69   NStCfiles          "unused"   NStCfiles    - Name of the files for nacelle structural controllers (quoted strings) [unused when NumNStC==0]
ServoDyn                                      70   NumTStC                   0   NumTStC      - Number of tower structural controllers (integer)
ServoDyn                                      71   TStCfiles          "unused"   TStCfiles    - Name of the files for tower structural controllers (quoted strings) [unused when NumTStC==0]
ServoDyn                                      72   NumSStC                   0   NumSStC      - Number of substructure structural controllers (integer)
ServoDyn                                      73   SStCfiles          "unused"   SStCfiles    - Name of the files for substructure structural controllers (quoted strings) [unused when NumSStC==0]
ServoDyn                                      74   CablesSection      ---------------------- CABLE CONTROL -------------------------------------------
ServoDyn                                      75   CCmode                    0   CCmode       - Cable control mode {0: none, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)
============================================= ==== =============== ========================================================================================================================================================================================================

============================================= ====== =============== ======================================================================================================================================================================================================
Modified in OpenFAST dev
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                                        Line    Flag Name        Example Value
============================================= ====== =============== ======================================================================================================================================================================================================
AeroDyn 15                                    9      TwrShadow        0   TwrShadow          - Calculate tower influence on wind based on downstream tower shadow (switch) {0=none, 1=Powles model, 2=Eames model}
SubDyn                                        26     Joints           JointID JointXss JointYss JointZss JointType JointDirX  JointDirY JointDirZ JointStiff
SubDyn                                        27     Joints             (-)     (m)      (m)      (m)      (-)        (-)       (-)       (-)      (Nm/rad) 
SubDyn                                        na     Members          MemberID MJointID1 MJointID2 MPropSetID1 MPropSetID2 MType COSMID
SubDyn                                        na     Members            (-)       (-)       (-)        (-)         (-)      (-)   (-)
SubDyn                                        na     ConcentratedM    CMJointID  JMass    JMXX      JMYY      JMZZ       JMXY     JMXZ     JMYZ    MCGX  MCGY MCGZ
SubDyn                                        na     ConcentratedM      (-)      (kg)    (kg*m^2)  (kg*m^2)  (kg*m^2)  (kg*m^2)  (kg*m^2) (kg*m^2)  (m)  (m)   (m)
HydroDyn                                      48     ExtnMod              1   ExctnMod       - Wave-excitation model {0: no wave-excitation calculation, 1: DFT, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES \*.ssexctn INPUT FILE]
HydroDyn                                      49     RdtnMod              2   RdtnMod        - Radiation memory-effect model {0: no memory-effect calculation, 1: convolution, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES \*.ss INPUT FILE]
HydroDyn                                      50     RdtnTMax            60   RdtnTMax       - Analysis time for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod>0; determines RdtnDOmega=Pi/RdtnTMax in the cosine transform; MAKE SURE THIS IS LONG ENOUGH FOR THE RADIATION IMPULSE RESPONSE FUNCTIONS TO DECAY TO NEAR-ZERO FOR THE GIVEN PLATFORM!]
HydroDyn                                      51     RdtnDT          0.0125   RdtnDT         - Time step for wave radiation kernel calculations (sec) [only used when PotMod=1 and ExctnMod>0 or RdtnMod>0; DT<=RdtnDT<=0.1 recommended; determines RdtnOmegaMax=Pi/RdtnDT in the cosine transform]
HydroDyn                                      54     PotFile         "Barge"  PotFile        - Root name of potential-flow model data; WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3) (quoted string) [1 to NBody if NBodyMod>1] [MAKE SURE THE FREQUENCIES INHERENT IN THESE WAMIT FILES SPAN THE PHYSICALLY-SIGNIFICANT RANGE OF FREQUENCIES FOR THE GIVEN PLATFORM; THEY MUST CONTAIN THE ZERO- AND INFINITE-FREQUENCY LIMITS!]
HydroDyn                                      55     WAMITULEN            1   WAMITULEN      - Characteristic body length scale used to redimensionalize WAMIT output (meters) [1 to NBody if NBodyMod>1] [only used when PotMod=1]
HydroDyn                                      56     PtfmRefxt          0.0   PtfmRefxt      - The xt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]
HydroDyn                                      57     PtfmRefyt          0.0   PtfmRefyt      - The yt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]
HydroDyn                                      58     PtfmRefzt          0.0   PtfmRefzt      - The zt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1. If NBodyMod=2,PtfmRefzt=0.0]
HydroDyn                                      59     PtfmRefztRot       0.0   PtfmRefztRot   - The rotation about zt of the body reference frame(s) from xt/yt (degrees) [1 to NBody] [only used when PotMod=1]
HydroDyn                                      60     PtfmVol0          6000   PtfmVol0       - Displaced volume of water when the body is in its undisplaced position (m^3) [1 to NBody] [only used when PotMod=1; USE THE SAME VALUE COMPUTED BY WAMIT AS OUTPUT IN THE .OUT FILE!]
HydroDyn                                      61     PtfmCOBxt          0.0   PtfmCOBxt      - The xt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]
HydroDyn                                      62     PtfmCOByt          0.0   PtfmCOByt      - The yt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]
HydroDyn                                      69-74  AddF0                0   AddF0    - Additional preload (N, N-m) [If NBodyMod=1, one size 6*NBody x 1 vector; if NBodyMod>1, NBody size 6 x 1 vectors]
HydroDyn                                      75-80  AddCLin          0 0 0 0 0 0   AddCLin  - Additional linear stiffness (N/m, N/rad, N-m/m, N-m/rad)                     [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
HydroDyn                                      81-86  AddBLin          0 0 0 0 0 0   AddBLin  - Additional linear damping(N/(m/s), N/(rad/s), N-m/(m/s), N-m/(rad/s))        [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
HydroDyn                                      87-92  AddBQuad         0 0 0 0 0 0   AddBQuad - Additional quadratic drag(N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2) [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
HydroDyn                                      na     Simple Coef Tab  SimplCd    SimplCdMG    SimplCa    SimplCaMG    SimplCp    SimplCpMG   SimplAxCa  SimplAxCaMG  SimplAxCa  SimplAxCaMG  SimplAxCp   SimplAxCpMG
HydroDyn                                      na                        (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)
HydroDyn                                      na     Depth Coef Tab   Dpth      DpthCd   DpthCdMG   DpthCa   DpthCaMG       DpthCp   DpthCpMG   DpthAxCa   DpthAxCaMG    DpthAxCa   DpthAxCaMG       DpthAxCp   DpthAxCpMG
HydroDyn                                      na                       (m)       (-)      (-)        (-)      (-)            (-)      (-)          (-)        (-)           (-)        (-)              (-)         (-)
HydroDyn                                      na     Member Coef Tab  MemberID    MemberCd1     MemberCd2    MemberCdMG1   MemberCdMG2    MemberCa1     MemberCa2    MemberCaMG1   MemberCaMG2    MemberCp1     MemberCp2    MemberCpMG1   MemberCpMG2   MemberAxCd1   MemberAxCd2  MemberAxCdMG1 MemberAxCdMG2  MemberAxCa1   MemberAxCa2  MemberAxCaMG1 MemberAxCaMG2  MemberAxCp1  MemberAxCp2   MemberAxCpMG1   MemberAxCpMG2
HydroDyn                                      na                        (-)         (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)
HydroDyn                                      na     OutList names    *see OutlistParameters.xlsx for new and revised output channel names*
============================================= ====== =============== ======================================================================================================================================================================================================



============================================= ==== =============== ========================================================================================================================================================================================================
Removed in OpenFAST dev
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                                        Line  Flag Name        Example Value
============================================= ==== =============== ========================================================================================================================================================================================================
HydroDyn                                      68   na              ---------------------- FLOATING PLATFORM FORCE FLAGS  -------------------------- [unused with WaveMod=6]
HydroDyn                                      69   PtfmSgF           True             PtfmSgF        - Platform horizontal surge translation force (flag) or DEFAULT
HydroDyn                                      70   PtfmSwF           True             PtfmSwF        - Platform horizontal sway translation force (flag) or DEFAULT
HydroDyn                                      71   PtfmHvF           True             PtfmHvF        - Platform vertical heave translation force (flag) or DEFAULT
HydroDyn                                      72   PtfmRF            True             PtfmRF         - Platform roll tilt rotation force (flag) or DEFAULT
HydroDyn                                      73   PtfmPF            True             PtfmPF         - Platform pitch tilt rotation force (flag) or DEFAULT
HydroDyn                                      74   PtfmYF            True             PtfmYF         - Platform yaw rotation force (flag) or DEFAULT
============================================= ==== =============== ========================================================================================================================================================================================================




OpenFAST v2.4.0 to OpenFAST v2.5.0
----------------------------------

-  InflowWind

   -  The input file parser is updated to a keyword/value pair based input.
      Each entry must have a corresponding keyword with the same spelling as
      expected. See :numref:`input_file_overview` for an overview.
   -  Driver code includes ability to convert between wind types

============== ==== ================== =============================================================================================================================================================================
Added in OpenFAST v2.5.0
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module        Line  Flag Name          Example Value
============== ==== ================== =============================================================================================================================================================================
IfW driver     6    [separator line]   ===================== File Conversion Options =================================
IfW driver     7    WrHAWC               false    WrHAWC      - Convert all data to HAWC2 format? (flag)
IfW driver     8    WrBladed             false    WrBladed    - Convert all data to Bladed format? (flag)
IfW driver     9    WrVTK                false    WrVTK       - Convert all data to VTK format? (flag)
InflowWind     7    VFlowAng                  0   VFlowAng    - Upflow angle (degrees) (not used for native Bladed format WindType=7)
============== ==== ================== =============================================================================================================================================================================


============================ ====== ================================================ ====================================================================================
Modified in OpenFAST v2.5.0
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Module                       Line    Flag Name / section                              Example Value
============================ ====== ================================================ ====================================================================================
MoorDyn                        na    added CtrlChan column in LINE PROPERTIES table  .. code-block:: none

                                                                                        Line    LineType  UnstrLen  NumSegs   NodeAnch  NodeFair  Outputs  CtrlChan
                                                                                        (-)       (-)       (m)       (-)       (-)       (-)       (-)      (-)
                                                                                        1         main     835.35      20        1         4         -        0
============================ ====== ================================================ ====================================================================================


OpenFAST v2.3.0 to OpenFAST v2.4.0
----------------------------------

Additional nodal output channels added for :ref:`AeroDyn15<AD-Nodal-Outputs>`, :ref:`BeamDyn<BD-Nodal-Outputs>`, and :ref:`ElastoDyn<ED-Nodal-Outputs>`.

============== ==== ================== =============================================================================================================================================================================
Added in OpenFAST v2.4.0
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
============== ==== ================== =============================================================================================================================================================================

============== ==== ================== ======================================================================================================================================================= =========================
Modified in OpenFAST v2.4.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Module        Line  New Flag Name      Example Value                                                                                                                                           Previous Flag Name/Value
============== ==== ================== ======================================================================================================================================================= =========================
AirFoilTables  40\* filtCutOff         "DEFAULT"  filtCutOff   - Reduced frequency cut-off for low-pass filtering the AoA input to UA, as well as the 1st and 2nd deriv (-) [default = 0.5]     [default = 20]
InflowWind     17   Filename_Uni        "unused"  Filename_Uni - Filename of time series data for uniform wind field.      (-)                                                                  Filename
InflowWind     18   RefHt_Uni                 90  RefHt_Uni    - Reference height for horizontal wind speed                (m)                                                                  RefHt
InflowWind     35   RefHt_Hawc                90  RefHt_Hawc   - reference height; the height (in meters) of the vertical center of the grid (m)                                                RefHt
InflowWind     47   PLExp_Hawc               0.2  PLExp_Hawc   - Power law exponent (-) (used for PL wind profile type only)                                                                    PLExp
InflowWind     49   XOffset                    0  XOffset      - Initial offset in +x direction (shift of wind box)                                                                             InitPosition(x)
============== ==== ================== ======================================================================================================================================================= =========================

\*non-comment line count, excluding lines contained if NumCoords is not 0.




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
