
.. _ss-primary-input_example:

Appendix A: Example SeaState primary input file
===============================================

The following is a SeaState primary input file generating irregular (JONSWAP) waves internally
structure::

      ------- SeaState v1.00.* Input File --------------------------------------------
      Example SeaState primary input file
      False            Echo           - Echo the input file data (flag)
      ---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
           "default"   WtrDens        - Water density (kg/m^3)
           "default"   WtrDpth        - Water depth (meters) relative to MSL
           "default"   MSL2SWL        - Offset between still-water level and mean sea level (meters) [positive upward; unused when WaveMod = 6; must be zero if PotMod=1 or 2]
      ---------------------- SPATIAL DISCRETIZATION ---------------------------------------------------
                30.0   X_HalfWidth    – Half-width of the domain in the X direction (m) [>0, NOTE: X[nX] = nX*dX, where nX = {-NX+1,-NX+2,…,NX-1} and dX = X_HalfWidth/(NX-1)]
                30.0   Y_HalfWidth    – Half-width of the domain in the Y direction (m) [>0, NOTE: Y[nY] = nY*dY, where nY = {-NY+1,-NY+2,…,NY-1} and dY = Y_HalfWidth/(NY-1)]
                25.0   Z_Depth        – Depth of the domain the Z direction (m) relative to SWL [0 < Z_Depth <= WtrDpth+MSL2SWL; "default": Z_Depth = WtrDpth+MSL2SWL; Z[nZ] = ( COS( nZ*dthetaZ ) – 1 )*Z_Depth, where nZ = {0,1,…NZ-1} and dthetaZ = pi/( 2*(NZ-1) )]
                  10   NX             – Number of nodes in half of the X-direction domain (-) [>=2]
                  10   NY             – Number of nodes in half of the Y-direction domain (-) [>=2]
                  10   NZ             – Number of nodes in the Z direction (-) [>=2]
      ---------------------- WAVES ---------------------------------------------------
                   2   WaveMod        - Incident wave kinematics model {0: none=still water, 1: regular (periodic), 1P#: regular with user-specified phase, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: White noise spectrum (irregular), 4: user-defined spectrum from routine UserWaveSpctrm (irregular), 5: Externally generated wave-elevation time series, 6: Externally generated full wave-kinematics time series [option 6 is invalid for PotMod/=0], 7: User-defined wave frequency components} (switch)
                   1   WaveStMod      - Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching} (switch) [unused when WaveMod=0 or when PotMod/=0]
                 600   WaveTMax       - Analysis time for incident wave calculations (sec) [unused when WaveMod=0; determines WaveDOmega=2Pi/WaveTMax in the IFFT]
                 0.2   WaveDT         - Time step for incident wave calculations     (sec) [unused when WaveMod=0 or 7; 0.1<=WaveDT<=1.0 recommended; determines WaveOmegaMax=Pi/WaveDT in the IFFT]
                 2.0   WaveHs         - Significant wave height of incident waves (meters) [used only when WaveMod=1, 2, or 3]
                  10   WaveTp         - Peak-spectral period of incident waves       (sec) [used only when WaveMod=1 or 2]
           "DEFAULT"   WavePkShp      - Peak-shape parameter of incident wave spectrum (-) or DEFAULT (string) [used only when WaveMod=2; use 1.0 for Pierson-Moskowitz]
            0.314159   WvLowCOff      - Low  cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]
            1.570796   WvHiCOff       - High cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]
                   0   WaveDir        - Incident wave propagation heading direction                         (degrees) [unused when WaveMod=0 or 6]
                   0   WaveDirMod     - Directional spreading function {0: none, 1: COS2S}                  (-)       [only used when WaveMod=2,3, or 4]
                   1   WaveDirSpread  - Wave direction spreading coefficient ( > 0 )                        (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1]
                   1   WaveNDir       - Number of wave directions                                           (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1; odd number only]
                   0   WaveDirRange   - Range of wave directions (full range: WaveDir +/- 1/2*WaveDirRange) (degrees) [only used when WaveMod=2,3,or 4 and WaveDirMod=1]
           123456789   WaveSeed(1)    - First  random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]
              RANLUX   WaveSeed(2)    - Second random seed of incident waves [-2147483648 to 2147483647] for intrinsic pRNG, or an alternative pRNG: "RanLux"    (-)       [unused when WaveMod=0, 5, or 6]
               FALSE   WaveNDAmp      - Flag for normally distributed amplitudes                            (flag)    [only used when WaveMod=2, 3, or 4]
            "unused"   WvKinFile      - Root name of externally generated wave data file(s)        (quoted string)    [used only when WaveMod=5 or 6]
      ---------------------- 2ND-ORDER WAVES ----------------------------------------- [unused with WaveMod=0 or 6]
               FALSE   WvDiffQTF      - Full difference-frequency 2nd-order wave kinematics (flag)
               FALSE   WvSumQTF       - Full summation-frequency  2nd-order wave kinematics (flag)
                   0   WvLowCOffD     - Low  frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]
            1.256637   WvHiCOffD      - High frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]
            0.618319   WvLowCOffS     - Low  frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]
            3.141593   WvHiCOffS      - High frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]
      ---------------------- CONSTRAINED WAVES --------------------------------------- 
                   0   ConstWaveMod   - Constrained wave model: 0=none; 1=Constrained wave with specified crest elevation, alpha; 2=Constrained wave with guaranteed peak-to-trough crest height, HCrest (flag)
                   3   CrestHmax      - Crest height (2*alpha for ConstWaveMod=1 or HCrest for ConstWaveMod=2), must be larger than WaveHs (m) [unused when ConstWaveMod=0]
                  60   CrestTime      - Time at which the crest appears (s) [unused when ConstWaveMod=0]
                   0   CrestXi        - X-position of the crest (m) [unused when ConstWaveMod=0]
                   0   CrestYi        - Y-position of the crest (m) [unused when ConstWaveMod=0]
      ---------------------- CURRENT ------------------------------------------------- [unused with WaveMod=6]
                   0   CurrMod        - Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent} (switch)
                   0   CurrSSV0       - Sub-surface current velocity at still water level  (m/s) [used only when CurrMod=1]
           "DEFAULT"   CurrSSDir      - Sub-surface current heading direction (degrees) or DEFAULT (string) [used only when CurrMod=1]
                  20   CurrNSRef      - Near-surface current reference depth            (meters) [used only when CurrMod=1]
                   0   CurrNSV0       - Near-surface current velocity at still water level (m/s) [used only when CurrMod=1]
                   0   CurrNSDir      - Near-surface current heading direction         (degrees) [used only when CurrMod=1]
                   0   CurrDIV        - Depth-independent current velocity                 (m/s) [used only when CurrMod=1]
                   0   CurrDIDir      - Depth-independent current heading direction    (degrees) [used only when CurrMod=1]
      ---------------------- MacCamy-Fuchs diffraction model -------------------------
                   0   MCFD           - MacCamy-Fuchs member radius (ignored if radius <= 0) [must be 0 when WaveMod 0 or 6] 
      ---------------------- OUTPUT --------------------------------------------------
      False            SeaStSum       - Output a summary file [flag]
                   3   OutSwtch       - Output requested channels to: [1=SeaState.out, 2=GlueCode.out, 3=both files]
      "E15.7e2"        OutFmt         - Output format for numerical results (quoted string) [not checked for validity!]
      "A15"            OutSFmt        - Output format for header strings (quoted string) [not checked for validity!]
                   2   NWaveElev      - Number of points where the incident wave elevations can be computed (-)       [maximum of 9 output locations]
            0.0, 5.0   WaveElevxi     - List of xi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]
            0.0, 0.0   WaveElevyi     - List of yi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]
                   2   NWaveKin       - Number of points where the wave kinematics can be output (-)       [maximum of 9 output locations]
          0.0,   0.0   WaveKinxi - List of xi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]
          0.0,   5.0   WaveKinyi - List of yi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]
        -14.0, -17.0   WaveKinzi - List of zi-coordinates for points where the wave kinematics can be output (meters) [NWaveKin points, separated by commas or white space; usused if NWaveKin = 0]
      ---------------------- OUTPUT CHANNELS -----------------------------------------
      "Wave1Elev, Wave1Elv1, Wave1Elv2"               - Wave elevation 
      "Wave2Elev, Wave2Elv1, Wave2Elv2"
      "FVel1xi, FVel1yi, FVel1zi"  - fluid velocity      at location 1
      "FAcc1xi, FAcc1yi, FAcc1zi"  - fluid accelerations at location 1
      "FDynP1"                     - fluid dynamic pressure at location 1
      "FVel2xi, FVel2yi, FVel2zi"  - fluid velocity      at location 2
      "FAcc2xi, FAcc2yi, FAcc2zi"  - fluid accelerations at location 2
      "FDynP2" 
      END

Appendix B: Example SeaState driver input file
==============================================
The following is a SeaState driver input file
structure::

      Seastate driver file
      Compatible with SeaState v1.00
      FALSE            Echo               - Echo the input file data (flag)
      ---------------------- ENVIRONMENTAL CONDITIONS -------------------------------
      9.80665          Gravity            - Gravity (m/s^2)
      1025             WtrDens            - Water density (kg/m^3)
      200              WtrDpth            - Water depth (m)
      0                MSL2SWL            - Offset between still-water level and mean sea level (m) [positive upward]
      ---------------------- SEASTATE -----------------------------------------------
      "./seastate_input.dat" SeaStateInputFile  - Primary SeaState input file name (quoted string)
      "./seastate.SeaSt"     OutRootName        - The name which prefixes all SeaState generated files (quoted string)
          0                  WrWvKinMod         - Write wave kinematics? [0: Do not write any kinematics to file, 1: Write only the (0,0) wave elevations to file, 2: Write the complete wave kinematics to files, no files written if WaveMod=6]
       5001                  NSteps             - Number of time steps in the simulations (-)
        0.1                  TimeInterval       - Time step for the simulation (sec)
      ---------------------- Waves multipoint elevation output ----------------------                                                                                                                
      False                  WaveElevSeriesFlag - T/F flag to output the wave elevation field (for movies)
      END of driver input file

.. _sea-output-channels:

Appendix C. List of Output Channels
===================================

This is a list of all possible output channels for the SeaState module. 
The names are grouped by meaning, but can be ordered in the OUTPUT 
CHANNELS section of the primary SeaState input file as you see fit. 
α refers to the output position for either wave elevation or wave 
kinematics specified in the OUTPUT section of the primary SeaState input 
file, where α is a number in the range [1,NWaveElev] for wave elevation 
outputs and in the range [1,NWaveKin] for wave kinematics outputs. 
Setting α > NWaveElev or α > NWaveKin yields invalid output. All outputs 
are in the global inertial-frame coordinate system.

================================================================ ========================================================================================================== ==========================================================================================
Channel Name(s)                                                  Units                                                                                                      Description
================================================================ ========================================================================================================== ==========================================================================================
**Wave Elevation**
WaveαElev                                                        (m)                                                                                                        Total (first- plus second-order) wave elevations (up to 9 designated locations)
WaveαElv1                                                        (m)                                                                                                        First-order wave elevations (up to 9 designated locations)
WaveαElv2                                                        (m)                                                                                                        Second-order wave elevations (up to 9 designated locations)
**Wave and Current Kinematics**                                                                                                                                 
FVelαxi, FVelαyi, FVelαzi                                        (m/s), (m/s), (m/s)                                                                                        Total (first- plus second-order waves and current) fluid velocities at α
FAccαxi, FAccαyi, FAccαzi                                        (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Total (first- plus second-order waves) fluid accelerations at α
FDynPα                                                           (Pa)                                                                                                       Total (first- plus second-order waves) fluid dynamic pressure at α
FAccMCFαxi, FAccMCFαyi, FAccMCFαzi                               (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Scaled first-order-wave fluid accelerations for the MacCamy-Fuchs members in HydroDyn at α
================================================================ ========================================================================================================== ==========================================================================================
