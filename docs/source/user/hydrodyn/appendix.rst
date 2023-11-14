
.. _hd-primary-input_example:

Appendix A: OC4 Semi-submersible Input File
===========================================

The following is a HydroDyn primary input file for OC4 semi-submersible
structure::

      ------- HydroDyn Input File ----------------------------------------------------
      NREL 5.0 MW offshore baseline floating platform HydroDyn input properties for the OC4 Semi-submersible.
      False            Echo           - Echo the input file data (flag)
      ---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
      "DEFAULT"   WtrDens        - Water density (kg/m^3)
      "DEFAULT"   WtrDpth        - Water depth (meters)
      "DEFAULT"   MSL2SWL        - Offset between still-water level and mean sea level (meters) [positive upward; unused when WaveMod = 6; must be zero if PotMod=1 or 2]
      ---------------------- WAVES ---------------------------------------------------
                  3   WaveMod        - Incident wave kinematics model {0: none=still water, 1: regular (periodic), 1P#: regular with user-specified phase, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: White noise spectrum (irregular), 4: user-defined spectrum from routine UserWaveSpctrm (irregular), 5: Externally generated wave-elevation time series, 6: Externally generated full wave-kinematics time series [option 6 is invalid for PotMod/=0]} (switch)
                  0   WaveStMod      - Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching} (switch) [unused when WaveMod=0 or when PotMod/=0]
            4600   WaveTMax       - Analysis time for incident wave calculations (sec) [unused when WaveMod=0; determines WaveDOmega=2Pi/WaveTMax in the IFFT]
            0.2   WaveDT         - Time step for incident wave calculations     (sec) [unused when WaveMod=0; 0.1<=WaveDT<=1.0 recommended; determines WaveOmegaMax=Pi/WaveDT in the IFFT]
            1.2646   WaveHs         - Significant wave height of incident waves (meters) [used only when WaveMod=1, 2, or 3]
                  10   WaveTp         - Peak-spectral period of incident waves       (sec) [used only when WaveMod=1 or 2]
      "DEFAULT"        WavePkShp      - Peak-shape parameter of incident wave spectrum (-) or DEFAULT (string) [used only when WaveMod=2; use 1.0 for Pierson-Moskowitz]
            0.314159   WvLowCOff      - Low  cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]
            1.570796   WvHiCOff       - High cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]
                  0   WaveDir        - Incident wave propagation heading direction                         (degrees) [unused when WaveMod=0 or 6]
                  0   WaveDirMod     - Directional spreading function {0: none, 1: COS2S}                  (-)       [only used when WaveMod=2,3, or 4]
                  1   WaveDirSpread  - Wave direction spreading coefficient ( > 0 )                        (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1]
                  1   WaveNDir       - Number of wave directions                                           (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1; odd number only]
                  0   WaveDirRange   - Range of wave directions (full range: WaveDir +/- 1/2*WaveDirRange) (degrees) [only used when WaveMod=2,3,or 4 and WaveDirMod=1]
      123456789   WaveSeed(1)    - First  random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]
      1011121314   WaveSeed(2)    - Second random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]
      FALSE            WaveNDAmp      - Flag for normally distributed amplitudes                            (flag)    [only used when WaveMod=2, 3, or 4]
      ""               WvKinFile      - Root name of externally generated wave data file(s)        (quoted string)    [used only when WaveMod=5 or 6]
                  1   NWaveElev      - Number of points where the incident wave elevations can be computed (-)       [maximum of 9 output locations]
                  0   WaveElevxi     - List of xi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]
                  0   WaveElevyi     - List of yi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]
      ---------------------- 2ND-ORDER WAVES ----------------------------------------- [unused with WaveMod=0 or 6]
      FALSE            WvDiffQTF      - Full difference-frequency 2nd-order wave kinematics (flag)
      FALSE            WvSumQTF       - Full summation-frequency  2nd-order wave kinematics (flag)
                  0   WvLowCOffD     - Low  frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]
            1.256637   WvHiCOffD      - High frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]
            0.618319   WvLowCOffS     - Low  frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]
            3.141593   WvHiCOffS      - High frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]
      ---------------------- CURRENT ------------------------------------------------- [unused with WaveMod=6]
                  0   CurrMod        - Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent} (switch)
                  0   CurrSSV0       - Sub-surface current velocity at still water level  (m/s) [used only when CurrMod=1]
      "DEFAULT"        CurrSSDir      - Sub-surface current heading direction (degrees) or DEFAULT (string) [used only when CurrMod=1]
                  20   CurrNSRef      - Near-surface current reference depth            (meters) [used only when CurrMod=1]
                  0   CurrNSV0       - Near-surface current velocity at still water level (m/s) [used only when CurrMod=1]
                  0   CurrNSDir      - Near-surface current heading direction         (degrees) [used only when CurrMod=1]
                  0   CurrDIV        - Depth-independent current velocity                 (m/s) [used only when CurrMod=1]
                  0   CurrDIDir      - Depth-independent current heading direction    (degrees) [used only when CurrMod=1]
      ---------------------- FLOATING PLATFORM --------------------------------------- [unused with WaveMod=6]
                  1   PotMod         - Potential-flow model {0: none=no potential flow, 1: frequency-to-time-domain transforms based on WAMIT output, 2: fluid-impulse theory (FIT)} (switch)
      "HydroData/marin_semi"    PotFile        - Root name of potential-flow model data; WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3) (quoted string) [MAKE SURE THE FREQUENCIES INHERENT IN THESE WAMIT FILES SPAN THE PHYSICALLY-SIGNIFICANT RANGE OF FREQUENCIES FOR THE GIVEN PLATFORM; THEY MUST CONTAIN THE ZERO- AND INFINITE-FREQUENCY LIMITS!]
                  1   WAMITULEN      - Characteristic body length scale used to redimensionalize WAMIT output (meters) [only used when PotMod=1]
            13917   PtfmVol0       - Displaced volume of water when the platform is in its undisplaced position (m^3) [only used when PotMod=1; USE THE SAME VALUE COMPUTED BY WAMIT AS OUTPUT IN THE .OUT FILE!]
                  0   PtfmCOBxt      - The xt offset of the center of buoyancy (COB) from the platform reference point (meters)  [only used when PotMod=1]
                  0   PtfmCOByt      - The yt offset of the center of buoyancy (COB) from the platform reference point (meters)  [only used when PotMod=1]
                  1   ExctnMod       - Wave Excitation model {0: None, 1: DFT, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ssexctn INPUT FILE]   
                  1   RdtnMod        - Radiation memory-effect model {0: no memory-effect calculation, 1: convolution, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ss INPUT FILE]
                  60   RdtnTMax       - Analysis time for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod>0; determines RdtnDOmega=Pi/RdtnTMax in the cosine transform; MAKE SURE THIS IS LONG ENOUGH FOR THE RADIATION IMPULSE RESPONSE FUNCTIONS TO DECAY TO NEAR-ZERO FOR THE GIVEN PLATFORM!]
      "DEFAULT"        RdtnDT         - Time step for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod=1; DT<=RdtnDT<=0.1 recommended; determines RdtnOmegaMax=Pi/RdtnDT in the cosine transform]
      ---------------------- 2ND-ORDER FLOATING PLATFORM FORCES ---------------------- [unused with WaveMod=0 or 6, or PotMod=0 or 2]
                  0   MnDrift        - Mean-drift 2nd-order forces computed                                       {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero]
                  0   NewmanApp      - Mean- and slow-drift 2nd-order forces computed with Newman's approximation {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. Used only when WaveDirMod=0]
                  0   DiffQTF        - Full difference-frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}          [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero]
                  0   SumQTF         - Full summation -frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}
      ---------------------- FLOATING PLATFORM FORCE FLAGS  -------------------------- [unused with WaveMod=6]
      True             PtfmSgF        - Platform horizontal surge translation force (flag) or DEFAULT
      True             PtfmSwF        - Platform horizontal sway translation force (flag) or DEFAULT
      True             PtfmHvF        - Platform vertical heave translation force (flag) or DEFAULT
      True             PtfmRF         - Platform roll tilt rotation force (flag) or DEFAULT
      True             PtfmPF         - Platform pitch tilt rotation force (flag) or DEFAULT
      True             PtfmYF         - Platform yaw rotation force (flag) or DEFAULT
      ---------------------- PLATFORM ADDITIONAL STIFFNESS AND DAMPING  --------------
                  0             0             0             0             0             0   AddF0    - Additional preload (N, N-m)
                  0             0             0             0             0             0   AddCLin  - Additional linear stiffness (N/m, N/rad, N-m/m, N-m/rad)
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0    1451298897             0             0
                  0             0             0             0    1451298897             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0   AddBLin  - Additional linear damping(N/(m/s), N/(rad/s), N-m/(m/s), N-m/(rad/s))
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0   AddBQuad - Additional quadratic drag(N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2)
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
                  0             0             0             0             0             0
      ---------------------- AXIAL COEFFICIENTS --------------------------------------
                  2   NAxCoef        - Number of axial coefficients (-)
      AxCoefID  AxCd     AxCa     AxCp
      (-)    (-)      (-)      (-)
      1     0.00     0.00     1.00
      2     9.60     0.00     1.00
      ---------------------- MEMBER JOINTS -------------------------------------------
                  44   NJoints        - Number of joints (-)   [must be exactly 0 or at least 2]
      JointID   Jointxi     Jointyi     Jointzi  JointAxID   JointOvrlp   [JointOvrlp= 0: do nothing at joint, 1: eliminate overlaps by calculating super member]
      (-)     (m)         (m)         (m)        (-)       (switch)
      1     0.00000     0.00000   -20.00000      1            0
      2     0.00000     0.00000    10.00000      1            0
      3    14.43376    25.00000   -14.00000      1            0
      4    14.43376    25.00000    12.00000      1            0
      5   -28.86751     0.00000   -14.00000      1            0
      6   -28.86751     0.00000    12.00000      1            0
      7    14.43376   -25.00000   -14.00000      1            0
      8    14.43376   -25.00000    12.00000      1            0
      9    14.43375    25.00000   -20.00000      2            0
      10   -28.86750     0.00000   -20.00000      2            0
      11    14.43375   -25.00000   -20.00000      2            0
      12     9.23760    22.00000    10.00000      1            0
      13   -23.67130     3.00000    10.00000      1            0
      14   -23.67130    -3.00000    10.00000      1            0
      15     9.23760   -22.00000    10.00000      1            0
      16    14.43375   -19.00000    10.00000      1            0
      17    14.43375    19.00000    10.00000      1            0
      18     4.04145    19.00000   -17.00000      1            0
      19   -18.47520     6.00000   -17.00000      1            0
      20   -18.47520    -6.00000   -17.00000      1            0
      21     4.04145   -19.00000   -17.00000      1            0
      22    14.43375   -13.00000   -17.00000      1            0
      23    14.43375    13.00000   -17.00000      1            0
      24     1.62500     2.81500    10.00000      1            0
      25    11.43376    19.80385    10.00000      1            0
      26    -3.25000     0.00000    10.00000      1            0
      27   -22.87000     0.00000    10.00000      1            0
      28     1.62500    -2.81500    10.00000      1            0
      29    11.43376   -19.80385    10.00000      1            0
      30     1.62500     2.81500   -17.00000      1            0
      31     8.43376    14.60770   -17.00000      1            0
      32    -3.25000     0.00000   -17.00000      1            0
      33   -16.87000     0.00000   -17.00000      1            0
      34     1.62500    -2.81500   -17.00000      1            0
      35     8.43376   -14.60770   -17.00000      1            0
      36     1.62500     2.81500   -16.20000      1            0
      37    11.43376    19.80385     9.13000      1            0
      38    -3.25000     0.00000   -16.20000      1            0
      39   -22.87000     0.00000     9.13000      1            0
      40     1.62500    -2.81500   -16.20000      1            0
      41    11.43376   -19.80385     9.13000      1            0
      42    14.43376    25.00000   -19.94000      1            0
      43   -28.86751     0.00000   -19.94000      1            0
      44    14.43376   -25.00000   -19.94000      1            0
      ---------------------- MEMBER CROSS-SECTION PROPERTIES -------------------------
                  4   NPropSets      - Number of member property sets (-)
      PropSetID    PropD         PropThck
      (-)        (m)            (m)
      1        6.50000        0.03000          ! Main Column
      2       12.00000        0.06000          ! Upper Columns
      3       24.00000        0.06000          ! Base Columns
      4        1.60000        0.01750          ! Pontoons
      ---------------------- SIMPLE HYDRODYNAMIC COEFFICIENTS (model 1) --------------
      SimplCd    SimplCdMG    SimplCa    SimplCaMG    SimplCp    SimplCpMG   SimplAxCa  SimplAxCaMG  SimplAxCp   SimplAxCpMG
            (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)
            0.00        0.00        0.00        0.00        1.00        1.00        0.00        0.00        1.00        1.00
      ---------------------- DEPTH-BASED HYDRODYNAMIC COEFFICIENTS (model 2) ---------
                  0   NCoefDpth       - Number of depth-dependent coefficients (-)
      Dpth      DpthCd   DpthCdMG   DpthCa   DpthCaMG       DpthCp   DpthCpMG   DpthAxCa   DpthAxCaMG       DpthAxCp   DpthAxCpMG
      (m)       (-)      (-)        (-)      (-)            (-)      (-)          (-)        (-)              (-)         (-)
      ---------------------- MEMBER-BASED HYDRODYNAMIC COEFFICIENTS (model 3) --------
                  25   NCoefMembers       - Number of member-based coefficients (-)
      MemberID    MemberCd1     MemberCd2    MemberCdMG1   MemberCdMG2    MemberCa1     MemberCa2    MemberCaMG1   MemberCaMG2    MemberCp1     MemberCp2    MemberCpMG1   MemberCpMG2   MemberAxCa1   MemberAxCa2  MemberAxCaMG1 MemberAxCaMG2  MemberAxCp1  MemberAxCp2   MemberAxCpMG1   MemberAxCpMG2
      (-)         (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)          ! Main Column
      1          0.56          0.56          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 1
      2          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 2
      3          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 3
      4          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 1
      5          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 2
      6          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 3
      7          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 1
      23          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 2
      24          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 3
      25          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 1
      8          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 2
      9          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 3
      10          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 1
      11          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 2
      12          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 3
      13          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 1
      14          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 2
      15          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 3
      16          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 1
      17          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 2
      18          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 3
      19          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 1
      20          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 2
      21          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 3
      22          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00
      -------------------- MEMBERS -------------------------------------------------
                  25   NMembers       - Number of members (-)
      MemberID  MJointID1  MJointID2  MPropSetID1  MPropSetID2  MDivSize   MCoefMod  PropPot   [MCoefMod=1: use simple coeff table, 2: use depth-based coeff table, 3: use member-based coeff table] [ PropPot/=0 if member is modeled with potential-flow theory]
      (-)        (-)        (-)         (-)          (-)        (m)      (switch)   (flag)
      1         1          2           1            1         1.0000      3        TRUE           ! Main Column
      2         3          4           2            2         1.0000      3        TRUE           ! Upper Column 1
      3         5          6           2            2         1.0000      3        TRUE           ! Upper Column 2
      4         7          8           2            2         1.0000      3        TRUE           ! Upper Column 3
      5        42          3           3            3         1.0000      3        TRUE           ! Base Column 1
      6        43          5           3            3         1.0000      3        TRUE           ! Base Column 2
      7        44          7           3            3         1.0000      3        TRUE           ! Base Column 3
      23         9         42           3            3         1.0000      3        TRUE           ! Base column cap 1
      24        10         43           3            3         1.0000      3        TRUE           ! Base column cap 2
      25        11         44           3            3         1.0000      3        TRUE           ! Base column cap 3
      8        12         13           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 1
      9        14         15           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 2
      10        16         17           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 3
      11        18         19           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 1
      12        20         21           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 2
      13        22         23           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 3
      14        24         25           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 1
      15        26         27           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 2
      16        28         29           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 3
      17        30         31           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 1
      18        32         33           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 2
      19        34         35           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 3
      20        36         37           4            4         1.0000      3        TRUE           ! Cross Brace 1
      21        38         39           4            4         1.0000      3        TRUE           ! Cross Brace 2
      22        40         41           4            4         1.0000      3        TRUE           ! Cross Brace 3
      ---------------------- FILLED MEMBERS ------------------------------------------
                  2   NFillGroups     - Number of filled member groups (-) [If FillDens = DEFAULT, then FillDens = WtrDens; FillFSLoc is related to MSL2SWL]
      FillNumM FillMList             FillFSLoc     FillDens
      (-)      (-)                   (m)           (kg/m^3)
      3   2   3   4    -6.17           1025
      3   5   6   7   -14.89           1025
      ---------------------- MARINE GROWTH -------------------------------------------
                  0   NMGDepths      - Number of marine-growth depths specified (-)
      MGDpth     MGThck       MGDens
      (m)        (m)         (kg/m^3)
      ---------------------- MEMBER OUTPUT LIST --------------------------------------
                  0   NMOutputs      - Number of member outputs (-) [must be < 10]
      MemberID   NOutLoc    NodeLocs [NOutLoc < 10; node locations are normalized distance from the start of the member, and must be >=0 and <= 1] [unused if NMOutputs=0]
      (-)        (-)        (-)
      ---------------------- JOINT OUTPUT LIST ---------------------------------------
                  0   NJOutputs      - Number of joint outputs [Must be < 10]
      0           JOutLst        - List of JointIDs which are to be output (-)[unused if NJOutputs=0]
      ---------------------- OUTPUT --------------------------------------------------
      True             HDSum          - Output a summary file [flag]
      False            OutAll         - Output all user-specified member and joint loads (only at each member end, not interior locations) [flag]
                  1   OutSwtch        - Output requested channels to: [1=Hydrodyn.out, 2=GlueCode.out, 3=both files]
      "ES11.4e2"       OutFmt         - Output format for numerical results (quoted string) [not checked for validity!]
      "A11"            OutSFmt        - Output format for header strings (quoted string) [not checked for validity!]
      ---------------------- OUTPUT CHANNELS -----------------------------------------
      "Wave1Elev"               - Wave elevation at the platform reference point (0,  0)
      END of output channels and end of file. (the word "END" must appear in the first 3 columns of this line)

Appendix B: OC4 Semi-submersible Input File
===========================================
The following is a HydroDyn driver input file for OC4 semi-submersible
structure::

      HydroDyn Driver file for OC4 Semi-submersible.  
      Compatible with HydroDyn v2.03.*
      TRUE             Echo                - Echo the input file data (flag)
      ---------------------- ENVIRONMENTAL CONDITIONS -------------------------------
      9.80665          Gravity             - Gravity (m/s^2)
      1025             WtrDens             - Water density (kg/m^3)
      200              WtrDpth             - Water depth (meters)
      0                MSL2SWL             - Offset between still-water level and mean sea level (meters) [positive upward; unused when WaveMod = 6; must be zero if PotMod=1 or 2]
      ---------------------- HYDRODYN -----------------------------------------------
      "./OC4Semi.dat"  HDInputFile         - Primary HydroDyn input file name (quoted string)
      "./OC4Semi"      OutRootName         - The name which prefixes all HydroDyn generated files (quoted string)
      1                NSteps              - Number of time steps in the simulations (-)
      0.025            TimeInterval        - TimeInterval for the simulation (sec)
      ---------------------- WAMIT INPUTS -------------------------------------------
      1                WAMITInputsMod      - Inputs model {0: all inputs are zero for every timestep, 1: steadystate inputs, 2: read inputs from a file (InputsFile)} (switch)
      ""               WAMITInputsFile     - Name of the inputs file if InputsMod = 2 (quoted string)
      ---------------------- WAMIT STEADY STATE INPUTS  -----------------------------
      1.0   2.0   3.0   4.0   5.0   6.0    uWAMITInSteady         - input displacements and rotations at the platform reference point (m, rads)
      7.0   8.0   9.0  10.0  11.0  12.0    uDotWAMITInSteady      - input translational and rotational velocities at the platform reference point (m/s, rads/s)
      13.0 14.0  15.0  16.0  17.0  18.0    uDotDotWAMITInSteady   - input translational and rotational acccelerations at the platform reference point (m/s^2, rads/s^2)
      ---------------------- MORISON INPUTS -----------------------------------------
      0                MorisonInputsMod    - Inputs model {0: all inputs are zero for every timestep, 1: steadystate inputs, 2: read inputs from a file (InputsFile)} (switch)
      " "              MorisonInputsFile   - Name of the inputs file if InputsMod = 2 (quoted string)
      ---------------------- MORISON STEADY STATE INPUTS  ---------------------------
      1.0   2.0   3.0   4.0   5.0   6.0    uMorisonInSteady       - input displacements and rotations for the morison elements (m, rads)
      7.0   8.0   9.0  10.0  11.0  12.0    uDotMorisonInSteady    - input translational and rotational velocities for the morison elements (m/s, rads/s)
      13.0 14.0  15.0  16.0  17.0  18.0    uDotDotMorisonInSteady - input translational and rotational acccelerations for the morison elements (m/s^2, rads/s^2)
      END of driver input file

.. _hd-output-channels:

Appendix C. List of Output Channels
===================================

This is a list of all possible output parameters for the HydroDyn
module. The names are grouped by meaning, but can be ordered in the
OUTPUT CHANNELS section of the HydroDyn input file as you see fit. MαNβ,
refers to output node β of output member α, where α is a number in the
range [1,9] and corresponds to row α in the MEMBER OUTPUT LIST table and
β is a number in the range [1,9] and corresponds to location β in the
**NodeLocs** list of that table entry. Jα refers to output joint α,
where α is a number in the range [1,9] and corresponds to row α in the
JOINT OUTPUT LIST table. Bα refers to body α, where α is a number in
the range [1,9]. Setting α > NBody yields invalid output; if NBody > 9,
only the first 9 bodies can be output. Waveα refers to point α where
wave elevations can be output, where α is a number in the range [1,9].
Setting α > NWaveElev yields invalid output. All outputs are in the
global inertial-frame coordinate.

================================================================ ========================================================================================================== ========================================================================================
Channel Name(s)                                                  Units                                                                                                      Description
================================================================ ========================================================================================================== ========================================================================================
**Wave and Current Kinematics**                                                                                                                                 
WaveαElev                                                        (m)                                                                                                        Total (first- plus second-order) wave elevations (up to 9 designated locations)
WaveαElv1                                                        (m)                                                                                                        First-order wave elevations (up to 9 designated locations)
WaveαElv2                                                        (m)                                                                                                        Second-order wave elevations (up to 9 designated locations)
MαNβVxi, MαNβVyi, MαNβVzi                                        (m/s), (m/s), (m/s)                                                                                        Total (first- plus second-order) fluid particle velocities at MαNβ
MαNβAxi, MαNβAyi, MαNβAzi                                        (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                                Total (first- plus second-order) fluid particle accelerations at MαNβ
MαNβDynP                                                         (Pa)                                                                                                       Total (first- plus second-order) fluid particle dynamic pressure at MαNβ
JαVxi, JαVyi, JαVzi                                              (m/s), (m/s), (m/s)                                                                                        Total (first- plus second-order) fluid particle velocities at Jα
JαAxi, JαAyi, JαAzi                                              (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                                Total (first- plus second-order) fluid particle accelerations at Jα
JαDynP                                                           (Pa)                                                                                                       Total (first- plus second-order) fluid particle dynamic pressure at Jα
**Total and Additional Loads**                                                                                                                                              
BαAddFxi, BαAddFyi, BαAddFzi, BαAddMxi, BαAddMyi, BαAddMzi       (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Loads due to additional preload, stiffness, and damping at Bα
HydroFxi, HydroFyi, HydroFzi, HydroMxi, HydroMyi, HydroMzi       (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Total integrated hydrodynamic loads from both potential flow and strip theory at (0,0,0)
**Loads from Potential-Flow Solution**                                                                                                                                      
BαWvsFxi, BαWvsFyi, BαWvsFzi, BαWvsMxi, BαWvsMyi, BαWvsMzi       (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Total (first- plus second-order) wave-excitation loads from diffraction at Bα
BαWvsF1xi, BαWvsF1yi, BαWvsF1zi, BαWvsM1xi, BαWvsM1yi, BαWvsM1zi (N), (N), (N), (N·m), (N·m), (N·m)                                                                         First-order wave-excitation loads from diffraction at Bα
BαWvsF2xi, BαWvsF2yi, BαWvsF2zi, BαWvsM2xi, BαWvsM2yi, BαWvsM2zi (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Second-order wave-excitation loads from diffraction at Bα
BαHdSFxi, BαHdSFyi, BαHdSFzi, BαHdSMxi, BαHdSMyi, BαHdSMzi       (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Hydrostatic loads at Bα
BαRdtFxi, BαRdtFyi, BαRdtFzi, BαRdtMxi, BαRdtMyi, BαRdtMzi       (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Wave-radiation loads at Bα
**Structural Motions**                                                                                                                                                      
BαSurge, BαSway, BαHeave, BαRoll, BαPitch BαYaw                  (m), (m), (m), (rad), (rad), (rad)                                                                         Displacements and rotations at Bα
BαTVxi, BαTVyi, BαTVzi, BαRVxi, BαRVyi, BαRVzi                   (m/s), (m/s), (m/s), (rad/s), (rad/s), (rad/s)                                                             Translational and rotational velocities at Bα
BαTAxi, BαTAyi, BαTAzi, BαRAxi, BαRAyi, BαRAzi                   (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`) Translational and rotational accelerations at Bα
MαNβSTVxi, MαNβSTVyi, MαNβSTVzi                                  (m/s), (m/s), (m/s)                                                                                        Structural translational velocities at MαNβ
MαNβSTAxi, MαNβSTAyi, MαNβSTAzi                                  (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                                Structural translational accelerations at MαNβ
JαSTVxi, JαSTVyi, JαSTVzi                                        (m/s), (m/s), (m/s)                                                                                        Structural translational velocities at Jα
JαSTAxi, JαSTAyi, JαSTAzi                                        (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                                Structural translational accelerations at Jα
**Distributed Loads (Per Unit Length) on Members**                                                                                                                          
MαNβFDxi, MαNβFDyi, MαNβFDzi                                     (N/m), (N/m), (N/m)                                                                                        Viscous-drag forces at MαNβ
MαNβFIxi, MαNβFIyi, MαNβFIzi                                     (N/m), (N/m), (N/m)                                                                                        Fluid-inertia forces at MαNβ
MαNβFBxi, MαNβFByi, MαNβFBzi, MαNβMBxi, MαNβMByi, MαNβMBzi       (N/m), (N/m), (N/m), (N·m/m), (N·m/m), (N·m/m)                                                             Buoyancy loads at MαNβ
MαNβFBFxi, MαNβFBFyi, MαNβFBFzi, MαNβMBFxi, MαNβMBFyi, MαNβMBFzi (N/m), (N/m), (N/m), (N·m/m), (N·m/m), (N·m/m)                                                             Negative buoyancy loads due to flooding/ballasting at MαNβ
MαNβFMGxi, MαNβFMGyi, MαNβFMGzi, MαNβMMGxi, MαNβMMGyi, MαNβMMGzi (N/m), (N/m), (N/m), (N·m/m), (N·m/m), (N·m/m)                                                             Loads due to marine growth weight at MαNβ
MαNβFAMxi, MαNβFAMyi, MαNβFAMzi                                  (N/m), (N/m), (N/m)                                                                                        Hydrodynamic added-mass forces at MαNβ
MαNβFAGxi, MαNβFAGyi, MαNβFAGzi, MαNβMAGxi, MαNβMAGyi, MαNβMAGzi (N/m), (N/m), (N/m), (N·m/m), (N·m/m), (N·m/m)                                                             Marine growth mass inertia loads at MαNβ
MαNβFAFxi, MαNβFAFyi, MαNβFAFzi, MαNβMAFxi, MαNβMAFyi, MαNβMAFzi (N/m), (N/m), (N/m), (N·m/m), (N·m/m), (N·m/m)                                                             Flooding/ballasting mass inertia loads at MαNβ
**Lumped Loads at Joints**                                                                                                                                                  
JαFDxi, JαFDyi, JαFDzi                                           (N), (N), (N)                                                                                              Viscous-drag forces at Jα
JαFIxi, JαFIyi, JαFIzi                                           (N), (N), (N)                                                                                              Fluid-inertia forces at Jα
JαFBxi, JαFByi, JαFBzi, JαMBxi, JαMByi, JαMBzi                   (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Buoyancy loads at Jα
JαFBFxi, JαFBFyi, JαFBFzi, JαMBFxi, JαMBFyi, JαMBFzi             (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Negative buoyancy loads due to flooding/ballasting at Jα
JαFMGxi, JαFMGyi, JαFMGzi                                        (N), (N), (N)                                                                                              Forces due to marine growth weight at Jα
JαFAMxi, JαFAMyi, JαFAMzi                                        (N), (N), (N)                                                                                              Hydrodynamic added-mass forces at Jα
JαFAGxi, JαFAGyi, JαFAGzi, JαMAGxi, JαMAGyi, JαMAGzi             (N), (N), (N), (N·m), (N·m), (N·m)                                                                         Marine growth mass inertia loads at Jα
================================================================ ========================================================================================================== ========================================================================================
