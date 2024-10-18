
.. _hd-primary-input_example:

Appendix A: OC4 Semi-submersible Input File
===========================================

The following is a HydroDyn primary input file for OC4 semi-submersible
structure::

      ------- HydroDyn Input File ----------------------------------------------------
      NREL 5.0 MW offshore baseline floating platform HydroDyn input properties for the OC4 Semi-submersible.
      False            Echo           - Echo the input file data (flag)
      ---------------------- FLOATING PLATFORM --------------------------------------- [unused with WaveMod=6]
                   1   PotMod         - Potential-flow model {0: none=no potential flow, 1: frequency-to-time-domain transforms based on WAMIT output, 2: fluid-impulse theory (FIT)} (switch)
                   1   ExctnMod       - Wave-excitation model {0: no wave-excitation calculation, 1: DFT, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ssexctn INPUT FILE; if PtfmYMod=1, need ExctnMod=0 or 1]
                   0   ExctnDisp      - Method of computing Wave Excitation {0: use undisplaced position, 1: use displaced position, 2: use low-pass filtered displaced position) [only used when PotMod=1 and ExctnMod>0 and SeaState's WaveMod>0]} (switch)
                  10   ExctnCutOff    - Cutoff (corner) frequency of the low-pass time-filtered displaced position (Hz) [>0.0] [used only when PotMod=1, ExctnMod>0, and ExctnDisp=2]) [only used when PotMod=1 and ExctnMod>0 and SeaState's WaveMod>0]} (switch)
                   0   PtfmYMod       - Model for large platform yaw offset {0: Static reference yaw offset based on PtfmRefY, 1: dynamic reference yaw offset based on low-pass filtering the PRP yaw motion with cutoff frequency PtfmYCutOff} (switch)
                   0   PtfmRefY       - Constant (if PtfmYMod=0) or initial (if PtfmYMod=1) platform reference yaw offset (deg)
                0.01   PtfmYCutOff    - Cutoff frequency for the low-pass filtering of PRP yaw motion when PtfmYMod=1 [>0.0; unused when PtfmYMod=0] (Hz)
                  36   NExctnHdg      - Number of evenly distributed platform yaw/heading angles over the range of [-180, 180) deg for which the wave excitation shall be computed [>=2; unused when PtfmYMod=0] (-)
                   1   RdtnMod        - Radiation memory-effect model {0: no memory-effect calculation, 1: convolution, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ss INPUT FILE]
                  60   RdtnTMax       - Analysis time for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod>0; determines RdtnDOmega=Pi/RdtnTMax in the cosine transform; MAKE SURE THIS IS LONG ENOUGH FOR THE RADIATION IMPULSE RESPONSE FUNCTIONS TO DECAY TO NEAR-ZERO FOR THE GIVEN PLATFORM!]
              0.0125   RdtnDT         - Time step for wave radiation kernel calculations (sec) [only used when PotMod=1 and ExctnMod>0 or RdtnMod>0; DT<=RdtnDT<=0.1 recommended; determines RdtnOmegaMax=Pi/RdtnDT in the cosine transform]
                   1   NBody          - Number of WAMIT bodies to be used (-) [>=1; only used when PotMod=1. If NBodyMod=1, the WAMIT data contains a vector of size 6*NBody x 1 and matrices of size 6*NBody x 6*NBody; if NBodyMod>1, there are NBody sets of WAMIT data each with a vector of size 6 x 1 and matrices of size 6 x 6]
                   1   NBodyMod       - Body coupling model {1: include coupling terms between each body and NBody in HydroDyn equals NBODY in WAMIT, 2: neglect coupling terms between each body and NBODY=1 with XBODY=0 in WAMIT, 3: Neglect coupling terms between each body and NBODY=1 with XBODY=/0 in WAMIT} (switch) [only used when PotMod=1]
        "marin_semi"   PotFile        - Root name of potential-flow model data; WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3) (quoted string) [1 to NBody if NBodyMod>1] [MAKE SURE THE FREQUENCIES INHERENT IN THESE WAMIT FILES SPAN THE PHYSICALLY-SIGNIFICANT RANGE OF FREQUENCIES FOR THE GIVEN PLATFORM; THEY MUST CONTAIN THE ZERO- AND INFINITE-FREQUENCY LIMITS!]
                   1   WAMITULEN      - Characteristic body length scale used to redimensionalize WAMIT output (meters) [1 to NBody if NBodyMod>1] [only used when PotMod=1]
                   0   PtfmRefxt      - The xt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]
                   0   PtfmRefyt      - The yt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]
                   0   PtfmRefzt      - The zt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1. If NBodyMod=2,PtfmRefzt=0.0]
                   0   PtfmRefztRot   - The rotation about zt of the body reference frame(s) from xt/yt (degrees) [1 to NBody] [only used when PotMod=1]
               13917   PtfmVol0       - Displaced volume of water when the body is in its undisplaced position (m^3) [1 to NBody] [only used when PotMod=1; USE THE SAME VALUE COMPUTED BY WAMIT AS OUTPUT IN THE .OUT FILE!]
                   0   PtfmCOBxt      - The xt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]
                   0   PtfmCOByt      - The yt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]
      ---------------------- 2ND-ORDER FLOATING PLATFORM FORCES ---------------------- [unused with WaveMod=0 or 6, or PotMod=0 or 2]
                   0   MnDrift        - Mean-drift 2nd-order forces computed                                       {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. If NBody>1, MnDrift  /=8]
                   0   NewmanApp      - Mean- and slow-drift 2nd-order forces computed with Newman's approximation {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. If NBody>1, NewmanApp/=8. Used only when WaveDirMod=0]
                   0   DiffQTF        - Full difference-frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}          [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. If PtfmYMod=1, need DiffQTF=0]
                   0   SumQTF         - Full summation -frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}          [If PtfmYMod=1, need SumQTF=0]
      ---------------------- PLATFORM ADDITIONAL STIFFNESS AND DAMPING  -------------- [unused with PotMod=0 or 2]
                   0   AddF0    - Additional preload (N, N-m)  [If NBodyMod=1, one size 6*NBody x 1 vector; if NBodyMod>1, NBody size 6 x 1 vectors]
                   0
                   0
                   0
                   0
                   0
                   0             0             0             0             0             0   AddCLin  - Additional linear stiffness (N/m, N/rad, N-m/m, N-m/rad)  [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0   AddBLin  - Additional linear damping (N/(m/s), N/(rad/s), N-m/(m/s), N-m/(rad/s))  [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0   AddBQuad - Additional quadratic damping (N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2)  [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
                   0             0             0             0             0             0
      ---------------------- STRIP THEORY OPTIONS --------------------------------------
                   0   WaveDisp       - Method of computing Wave Kinematics {0: use undisplaced position, 1: use displaced position) } (switch) [If PtfmYMod=1, need WaveDisp=1]
                   0   AMMod          - Method of computing distributed added-mass force. (0: Only and always on nodes below SWL at the undisplaced position. 2: Up to the instantaneous free surface) [overwrite to 0 when WaveMod = 0 or 6 or when WaveStMod = 0 in SeaState]
      ---------------------- AXIAL COEFFICIENTS --------------------------------------
                   2   NAxCoef        - Number of axial coefficients (-)
      AxCoefID  AxCd     AxCa     AxCp    AxFDMod   AxVnCOff  AxFDLoFSc
         (-)    (-)      (-)      (-)      (-)       (-)        (-)
          1     0.00     0.00     1.00      0        0.00       1.00
          2     9.60     0.00     1.00      0        0.00       1.00
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
      SimplCd    SimplCdMG    SimplCa    SimplCaMG    SimplCp    SimplCpMG   SimplAxCd  SimplAxCdMG   SimplAxCa  SimplAxCaMG  SimplAxCp   SimplAxCpMG    SimplCb    SimplCbMG
         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)            (-)         (-)
         0.00        0.00        0.00        0.00        1.00        1.00        0.00        0.00        0.00        0.00        1.00        1.00           1.00        1.00
      ---------------------- DEPTH-BASED HYDRODYNAMIC COEFFICIENTS (model 2) ---------
                   0   NCoefDpth       - Number of depth-dependent coefficients (-)
      Dpth      DpthCd   DpthCdMG   DpthCa   DpthCaMG       DpthCp   DpthCpMG   DpthAxCd   DpthAxCdMG   DpthAxCa   DpthAxCaMG   DpthAxCp   DpthAxCpMG   DpthCb   DpthCbMG
      (m)       (-)      (-)        (-)      (-)            (-)      (-)        (-)        (-)          (-)        (-)          (-)        (-)           (-)      (-)
      ---------------------- MEMBER-BASED HYDRODYNAMIC COEFFICIENTS (model 3) --------
                  25   NCoefMembers       - Number of member-based coefficients (-)
      MemberID    MemberCd1     MemberCd2    MemberCdMG1   MemberCdMG2    MemberCa1     MemberCa2    MemberCaMG1   MemberCaMG2    MemberCp1     MemberCp2    MemberCpMG1   MemberCpMG2   MemberAxCd1   MemberAxCd2  MemberAxCdMG1 MemberAxCdMG2  MemberAxCa1   MemberAxCa2  MemberAxCaMG1 MemberAxCaMG2  MemberAxCp1  MemberAxCp2   MemberAxCpMG1   MemberAxCpMG2    MemberCb1     MemberCb2    MemberCbMG1   MemberCbMG2
         (-)         (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)              (-)           (-)           (-)           (-)
          1          0.56          0.56          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Main Column
          2          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Upper Column 1
          3          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Upper Column 2
          4          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Upper Column 3
          5          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base Column 1
          6          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base Column 2
          7          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base Column 3
         23          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base column cap 1
         24          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base column cap 2
         25          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Base column cap 3
          8          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Upper 1
          9          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Upper 2
         10          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Upper 3
         11          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Lower 1
         12          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Lower 2
         13          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Delta Pontoon, Lower 3
         14          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Upper 1
         15          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Upper 2
         16          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Upper 3
         17          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Lower 1
         18          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Lower 2
         19          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Y Pontoon, Lower 3
         20          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Cross Brace 1
         21          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Cross Brace 2
         22          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00             1.00          1.00          1.00          1.00          ! Cross Brace 3
      -------------------- MEMBERS -------------------------------------------------
                  25   NMembers       - Number of members (-)
      MemberID  MJointID1  MJointID2  MPropSetID1  MPropSetID2  MDivSize   MCoefMod   MHstLMod  PropPot   [MCoefMod=1: use simple coeff table, 2: use depth-based coeff table, 3: use member-based coeff table] [ PropPot/=0 if member is modeled with potential-flow theory]
        (-)        (-)        (-)         (-)          (-)        (m)      (switch)   (switch)  (flag)
          1         1          2           1            1         1.0000      3          1       TRUE           ! Main Column
          2         3          4           2            2         1.0000      3          1       TRUE           ! Upper Column 1
          3         5          6           2            2         1.0000      3          1       TRUE           ! Upper Column 2
          4         7          8           2            2         1.0000      3          1       TRUE           ! Upper Column 3
          5        42          3           3            3         1.0000      3          1       TRUE           ! Base Column 1
          6        43          5           3            3         1.0000      3          1       TRUE           ! Base Column 2
          7        44          7           3            3         1.0000      3          1       TRUE           ! Base Column 3
         23         9         42           3            3         1.0000      3          1       TRUE           ! Base column cap 1
         24        10         43           3            3         1.0000      3          1       TRUE           ! Base column cap 2
         25        11         44           3            3         1.0000      3          1       TRUE           ! Base column cap 3
          8        12         13           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Upper 1
          9        14         15           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Upper 2
         10        16         17           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Upper 3
         11        18         19           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Lower 1
         12        20         21           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Lower 2
         13        22         23           4            4         1.0000      3          1       TRUE           ! Delta Pontoon, Lower 3
         14        24         25           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Upper 1
         15        26         27           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Upper 2
         16        28         29           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Upper 3
         17        30         31           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Lower 1
         18        32         33           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Lower 2
         19        34         35           4            4         1.0000      3          1       TRUE           ! Y Pontoon, Lower 3
         20        36         37           4            4         1.0000      3          1       TRUE           ! Cross Brace 1
         21        38         39           4            4         1.0000      3          1       TRUE           ! Cross Brace 2
         22        40         41           4            4         1.0000      3          1       TRUE           ! Cross Brace 3
      ---------------------- FILLED MEMBERS ------------------------------------------
                   2   NFillGroups     - Number of filled member groups (-) [If FillDens = DEFAULT, then FillDens = WtrDens; FillFSLoc is related to MSL2SWL]
      FillNumM FillMList FillFSLoc     FillDens
      (-)      (-)       (m)           (kg/m^3)
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
                       JOutLst        - List of JointIDs which are to be output (-)[unused if NJOutputs=0]
      ---------------------- OUTPUT --------------------------------------------------
      True             HDSum          - Output a summary file [flag]
      False            OutAll         - Output all user-specified member and joint loads (only at each member end, not interior locations) [flag]
                   2   OutSwtch       - Output requested channels to: [1=Hydrodyn.out, 2=GlueCode.out, 3=both files]
      "E16.8e2"        OutFmt         - Output format for numerical results (quoted string) [not checked for validity!]
      "A11"            OutSFmt        - Output format for header strings (quoted string) [not checked for validity!]
      ---------------------- OUTPUT CHANNELS -----------------------------------------
      HydroFxi                   
      HydroFyi                   
      HydroFzi                   
      HydroMxi                   
      HydroMyi                   
      HydroMzi                   
      END of output channels and end of file. (the word "END" must appear in the first 3 columns of this line)

Appendix B: OC4 Semi-submersible Input File
===========================================
The following is a HydroDyn driver input file for OC4 semi-submersible
structure::

      ------- HydroDyn Driver Input File --------------------------------------------
      HydroDyn Driver file for OC4 Semi-submersible.
            FALSE   Echo                - Echo the input file data (flag)
      ---------------------- ENVIRONMENTAL CONDITIONS -------------------------------
          9.80665   Gravity             - Gravity (m/s^2)
             1025   WtrDens             - Water density (kg/m^3)
              200   WtrDpth             - Water depth (m)
                0   MSL2SWL             - Offset between still-water level and mean sea level (m) [positive upward]
      ---------------------- HYDRODYN -----------------------------------------------
      "./OC4Semi.dat"    HDInputFile       - Primary HydroDyn input file name (quoted string)
      "./SeaState.dat"   SeaStateInputFile - Primary SeaState input file name (quoted string)
      "./OC4Semi"        OutRootName       - The name which prefixes all HydroDyn generated files (quoted string)
            FALSE        Linearize         - Flag to enable linearization
             4801        NSteps            - Number of time steps in the simulation (-)   [60 seconds total]
           0.0125        TimeInterval      - Time step for the simulation (sec)
      ---------------------- PRP INPUTS (Platform Reference Point) ------------------
                0   PRPInputsMod      - Model for the PRP (platform reference point) inputs {0: all inputs are zero for every timestep, 1: steady-state inputs, 2: read inputs from a file (InputsFile)} (switch)
                0   PtfmRefzt         - Vertical distance from the ground level to the platform reference point (m)
      "not_used"    PRPInputsFile     - Filename for the PRP HydroDyn input InputsMod = 2 (quoted string)
      ---------------------- PRP STEADY STATE INPUTS  -------------------------------
                0,          0,          0,          0,          0,          0    uPRPInSteady         - PRP Steady-state (3) displacements and (3) rotations at the platform reference point (m, m, m, rad, rad, rad)
                0,          0,          0,          0,          0,          0    uDotPRPInSteady      - PRP Steady-state (3) translational and (3) rotational velocities at the platform reference point (m/s, rads/s)
                0,          0,          0,          0,          0,          0    uDotDotPRPInSteady   - PRP Steady-state (3) translational and (3) rotational accelerations at the platform reference point (m/s^2, rads/s^2)

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
MαNβVxi, MαNβVyi, MαNβVzi                                        (m/s), (m/s), (m/s)                                                                                        Total (first- plus second-order) fluid particle velocities at MαNβ
MαNβAxi, MαNβAyi, MαNβAzi                                        (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Total (first- plus second-order) fluid particle accelerations at MαNβ
MαNβDynP                                                         (Pa)                                                                                                       Total (first- plus second-order) fluid particle dynamic pressure at MαNβ
JαVxi, JαVyi, JαVzi                                              (m/s), (m/s), (m/s)                                                                                        Total (first- plus second-order) fluid particle velocities at Jα
JαAxi, JαAyi, JαAzi                                              (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Total (first- plus second-order) fluid particle accelerations at Jα
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
PRPSurge, PRPSway, PRPHeave, PRPRoll, PRPPitch, PRPYaw           (m), (m), (m), (rad), (rad), (rad)                                                                         Displacements and rotations at platform reference point (PRP)
PRPTVxi, PRPTVyi, PRPTVzi, PRPRVxi, PRPRVyi, PRPRVzi             (m/s), (m/s), (m/s), (rad/s), (rad/s), (rad/s)                                                             Translational and rotational velocities of the PRP
PRPTAxi, PRPTAyi, PRPTAzi, PRPRAxi, PRPRAyi, PRPRAzi             (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`) Translational and rotational accelerations of the PRP
BαSurge, BαSway, BαHeave, BαRoll, BαPitch BαYaw                  (m), (m), (m), (rad), (rad), (rad)                                                                         Displacements and rotations at Bα
BαTVxi, BαTVyi, BαTVzi, BαRVxi, BαRVyi, BαRVzi                   (m/s), (m/s), (m/s), (rad/s), (rad/s), (rad/s)                                                             Translational and rotational velocities at Bα
BαTAxi, BαTAyi, BαTAzi, BαRAxi, BαRAyi, BαRAzi                   (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`), (rad/s\ :sup:`2`) Translational and rotational accelerations at Bα
MαNβSTVxi, MαNβSTVyi, MαNβSTVzi                                  (m/s), (m/s), (m/s)                                                                                        Structural translational velocities at MαNβ
MαNβSTAxi, MαNβSTAyi, MαNβSTAzi                                  (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Structural translational accelerations at MαNβ
JαSTVxi, JαSTVyi, JαSTVzi                                        (m/s), (m/s), (m/s)                                                                                        Structural translational velocities at Jα
JαSTAxi, JαSTAyi, JαSTAzi                                        (m/s\ :sup:`2`), (m/s\ :sup:`2`), (m/s\ :sup:`2`)                                                          Structural translational accelerations at Jα
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
