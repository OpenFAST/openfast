removed variables:
------------------------------------------------------------------------------------------------------------------------------------

9999.9      TiDynBrk    - Time to initiate deployment of the dynamic generator brake [CURRENTLY IGNORED] (s)
   1        PSpnElN     - Number of the innermost blade element which is still part of the pitchable portion of the blade for partial-span pitch control [1 to BldNodes] [CURRENTLY IGNORED] (-)
            DynBrkFi    - File containing a mech-gen-torque vs HSS-speed curve for a dynamic brake [CURRENTLY IGNORED] (quoted string)
False       CalcTMode   - Calculate tower mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)
False       CalcBMode   - Calculate blade mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)
"NRELOffshrBsline5MW_Platform_Monopile_RF.dat"   PtfmFile    - Name of file containing platform properties (quoted string) [unused when PtfmModel=0]
False       Platform    - Read in additional data for platform (flag)
   3        ADAMSPrep   - ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create an ADAMS model, 3: do both} (switch)
   1        AnalMode    - Analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch)
"NRELOffshrBsline5MW_ADAMSSpecific.dat"          ADAMSFile   - Name of file containing ADAMS-specific input parameters (quoted string) [unused when ADAMSPrep=1]
"NRELOffshrBsline5MW_Linear.dat"                 LinFile     - Name of file containing FAST linearization parameters (quoted string) [unused when AnalMode=1]
            NoiseFile   - Name of file containing aerodynamic noise input parameters (quoted string) [used only when CompNoise=True]
False       CompNoise   - Compute aerodynamic noise (flag)


PtfmLdMod is now CompUserPtfmLd
TwrLdMod is now CompUserTwrLd
%BJJ ==> Check bugzilla bug 195


variables still to be removed:
------------------------------------------------------------------------------------------------------------------------------------


ADAMS columns in Tower and Blade files removed (keep code for now)
AeroCent column in blade file removed... look for repurcussions... should be in AeroDyn

---------------------- TIP-BRAKE -----------------------------------------------
   0.0      TBDrConN    - Tip-brake drag constant during normal operation, Cd*Area (m^2)
   0.0      TBDrConD    - Tip-brake drag constant during fully-deployed operation, Cd*Area (m^2)
   0.0      TpBrDT      - Time for tip-brake to reach full deployment once released (sec)

9999.9      TTpBrDp(1)  - Time to initiate deployment of tip brake 1 (s)
9999.9      TTpBrDp(2)  - Time to initiate deployment of tip brake 2 (s)
9999.9      TTpBrDp(3)  - Time to initiate deployment of tip brake 3 (s) [unused for 2 blades]
9999.9      TBDepISp(1) - Deployment-initiation speed for the tip brake on blade 1 (rpm)
9999.9      TBDepISp(2) - Deployment-initiation speed for the tip brake on blade 2 (rpm)
9999.9      TBDepISp(3) - Deployment-initiation speed for the tip brake on blade 3 (rpm) [unused for 2 blades]



tip-brake aerodynamics: let's not support it in this version
furling: let's not support it in this version
adams: let's not support it in this version


controls:
add simple yaw control: fixed yaw with this variable:
   0.0      YawNeut     - Neutral yaw position--yaw spring force is zero at this yaw (degrees)  <--- RENAME THIS VARIABLE



ed: Outputs that are invalid now:
      ! RotCq  -> uses airdensity
      ! RotCp  -> uses airdensity
      ! RotCt  -> uses airdensity
      ! HSShftCq -> uses airdensity
      ! HSShftCp -> uses airdensity
      ! GenCq -> uses airdensity
      ! GenCp -> uses airdensity

      ! NacYawErr -> only computed with CompAero.... is this a SrvD output now?

      HSSBrTq
      GenTq
      GenPwr

        WindVxi
        WindVyi
        WindVzi
       TotWindV
       HorWindV
      HorWndDir
      VerWndDir



      TFinAlpha
      TFinCLift
      TFinCDrag
      TFinDnPrs
      TFinCPFx
      TFinCPFy