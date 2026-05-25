import pytest
import os.path as osp
import platform
from pathlib import Path

# looking up OS for the correct executable extension
mactype = platform.system().lower()
if mactype in ["linux", "linux2", "darwin"]:
    exeExt = ""
elif mactype in ["win32", "windows", "cygwin"]: #NOTE: platform.system()='Windows', sys.platform='win32'
    libext = '.exe'
else:
    raise ValueError('Unknown platform type: '+mactype)

REPOSITORY_ROOT = osp.dirname(osp.dirname(osp.dirname(osp.dirname(__file__))))
BUILD_DIR = osp.join(REPOSITORY_ROOT, "build/reg_tests")

# Path to the OpenFAST executable
OF_PATH = osp.join(REPOSITORY_ROOT,"build/glue-codes/openfast",f"openfast{exeExt}")

def pytest_addoption(parser):
    parser.addoption("--executable", action="store", default=OF_PATH, help="Path to the OpenFAST executable")
    parser.addoption("--source_dir", action="store", default=REPOSITORY_ROOT, help="Path to the openfast repository")
    parser.addoption("--build_dir", action="store", default=BUILD_DIR, help="Path to the test data directory")


# ── Fixtures for new IO/Driver tests ──

@pytest.fixture
def r_test_dir():
    """Path to r-test root. Skip if not available."""
    candidate = Path(REPOSITORY_ROOT) / "reg_tests" / "r-test"
    if not candidate.exists():
        pytest.skip("r-test repo not available")
    return candidate


@pytest.fixture
def r_test_5mw_dir(r_test_dir):
    """Path to r-test 5MW_Land_DLL_WTurb case."""
    candidate = r_test_dir / "glue-codes" / "openfast" / "5MW_Land_DLL_WTurb"
    if not candidate.exists():
        pytest.skip("5MW_Land_DLL_WTurb r-test case not found")
    return candidate


@pytest.fixture
def r_test_fastfarm_dir(r_test_dir):
    """Path to r-test FAST.Farm cases."""
    candidate = r_test_dir / "glue-codes" / "fast-farm"
    if not candidate.exists():
        pytest.skip("r-test FAST.Farm cases not available")
    return candidate


@pytest.fixture
def sample_ed_file(tmp_path):
    """Minimal valid ElastoDyn .dat content for testing."""
    # Create minimal blade file
    blade_content = """\
------- ELASTODYN INDIVIDUAL BLADE INPUT FILE --------------------------
Test blade file
---------------------- BLADE PARAMETERS ----------------------------------------
6                      NBlInpSt    - Number of blade input stations (-)
1.0                    BldFlDmp(1) - Blade flap mode #1 structural damping (%)
1.0                    BldFlDmp(2) - Blade flap mode #2 structural damping (%)
1.0                    BldEdDmp(1) - Blade edge mode #1 structural damping (%)
---------------------- BLADE ADJUSTMENT FACTORS --------------------------------
1.0                    FlStTunr(1) - Blade flapwise modal stiffness tuner, 1st mode (-)
1.0                    FlStTunr(2) - Blade flapwise modal stiffness tuner, 2nd mode (-)
1.0                    AdjBlMs     - Factor to adjust blade mass density (-)
1.0                    AdjFlSt     - Factor to adjust blade flap stiffness (-)
1.0                    AdjEdSt     - Factor to adjust blade edge stiffness (-)
---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------
    BlFract      StrcTwst       BMassDen        FlpStff        EdgStff
      (-)         (deg)          (kg/m)         (Nm^2)         (Nm^2)
 0.000000000000000e+00  1.330000000000000e+01  6.780000000000000e+02  1.810000000000000e+10  1.810000000000000e+10
 2.000000000000000e-01  1.330000000000000e+01  6.780000000000000e+02  1.810000000000000e+10  1.810000000000000e+10
 4.000000000000000e-01  1.000000000000000e+01  4.000000000000000e+02  9.000000000000000e+09  9.000000000000000e+09
 6.000000000000000e-01  7.000000000000000e+00  2.500000000000000e+02  5.000000000000000e+09  5.000000000000000e+09
 8.000000000000000e-01  3.000000000000000e+00  1.200000000000000e+02  2.000000000000000e+09  2.000000000000000e+09
 1.000000000000000e+00  1.000000000000000e+00  1.000000000000000e+01  1.000000000000000e+08  1.000000000000000e+08
---------------------- BLADE MODE SHAPES ---------------------------------------
 0.0622                BldFl1Sh(2) - Flap mode 1, coeff of x^2
 1.7254                BldFl1Sh(3) -            , coeff of x^3
-3.2452                BldFl1Sh(4) -            , coeff of x^4
 4.7131                BldFl1Sh(5) -            , coeff of x^5
-2.2555                BldFl1Sh(6) -            , coeff of x^6
-0.5809                BldFl2Sh(2) - Flap mode 2, coeff of x^2
 1.2067                BldFl2Sh(3) -            , coeff of x^3
-15.5349               BldFl2Sh(4) -            , coeff of x^4
 29.7347               BldFl2Sh(5) -            , coeff of x^5
-13.8255               BldFl2Sh(6) -            , coeff of x^6
 0.3877                BldEdgSh(2) - Edge mode 1, coeff of x^2
 2.1179                BldEdgSh(3) -            , coeff of x^3
-4.4402                BldEdgSh(4) -            , coeff of x^4
 5.8456                BldEdgSh(5) -            , coeff of x^5
-2.9110                BldEdgSh(6) -            , coeff of x^6
"""
    blade_path = tmp_path / "test_blade.dat"
    blade_path.write_text(blade_content)

    # Create minimal tower file
    tower_content = """\
------- ELASTODYN TOWER INPUT FILE -------------------------------------
Test tower file
---------------------- TOWER PARAMETERS ----------------------------------------
3                      NTwInpSt    - Number of input stations to specify tower geometry
1.0                    TwrFADmp(1) - Tower 1st fore-aft mode structural damping ratio (%)
1.0                    TwrFADmp(2) - Tower 2nd fore-aft mode structural damping ratio (%)
1.0                    TwrSSDmp(1) - Tower 1st side-to-side mode structural damping ratio (%)
1.0                    TwrSSDmp(2) - Tower 2nd side-to-side mode structural damping ratio (%)
---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------
1.0                    FAStTunr(1) - Tower fore-aft modal stiffness tuner, 1st mode (-)
1.0                    FAStTunr(2) - Tower fore-aft modal stiffness tuner, 2nd mode (-)
1.0                    SSStTunr(1) - Tower side-to-side stiffness tuner, 1st mode (-)
1.0                    SSStTunr(2) - Tower side-to-side stiffness tuner, 2nd mode (-)
1.0                    AdjTwMa     - Factor to adjust tower mass density (-)
1.0                    AdjFASt     - Factor to adjust tower fore-aft stiffness (-)
1.0                    AdjSSSt     - Factor to adjust tower side-to-side stiffness (-)
---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------
  HtFract       TMassDen         TwFAStif       TwSSStif
   (-)           (kg/m)           (Nm^2)         (Nm^2)
 0.000000000000000e+00  8.857000000000000e+03  6.144000000000000e+11  6.144000000000000e+11
 5.000000000000000e-01  5.000000000000000e+03  3.000000000000000e+11  3.000000000000000e+11
 1.000000000000000e+00  3.475000000000000e+03  1.858000000000000e+11  1.858000000000000e+11
---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------
 0.7004                TwFAM1Sh(2) - Mode 1, coefficient of x^2 term
 2.1963                TwFAM1Sh(3) -       , coefficient of x^3 term
-5.6202                TwFAM1Sh(4) -       , coefficient of x^4 term
 6.2275                TwFAM1Sh(5) -       , coefficient of x^5 term
-2.5040                TwFAM1Sh(6) -       , coefficient of x^6 term
-26.0840               TwFAM2Sh(2) - Mode 2, coefficient of x^2 term
 70.5765               TwFAM2Sh(3) -       , coefficient of x^3 term
-79.6498               TwFAM2Sh(4) -       , coefficient of x^4 term
 51.0700               TwFAM2Sh(5) -       , coefficient of x^5 term
-14.9127               TwFAM2Sh(6) -       , coefficient of x^6 term
---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------
 0.7004                TwSSM1Sh(2) - Mode 1, coefficient of x^2 term
 2.1963                TwSSM1Sh(3) -       , coefficient of x^3 term
-5.6202                TwSSM1Sh(4) -       , coefficient of x^4 term
 6.2275                TwSSM1Sh(5) -       , coefficient of x^5 term
-2.5040                TwSSM1Sh(6) -       , coefficient of x^6 term
-26.0840               TwSSM2Sh(2) - Mode 2, coefficient of x^2 term
 70.5765               TwSSM2Sh(3) -       , coefficient of x^3 term
-79.6498               TwSSM2Sh(4) -       , coefficient of x^4 term
 51.0700               TwSSM2Sh(5) -       , coefficient of x^5 term
-14.9127               TwSSM2Sh(6) -       , coefficient of x^6 term
"""
    tower_path = tmp_path / "test_tower.dat"
    tower_path.write_text(tower_content)

    # Create minimal ElastoDyn main file
    ed_content = f"""\
------- ELASTODYN INPUT FILE -------------------------------------------
Test ElastoDyn file
---------------------- SIMULATION CONTROL --------------------------------------
False                  Echo        - Echo input data to "<RootName>.ech" (flag)
          3            Method      - Integration method {{1: RK4, 2: AB4, 3: ABM4}} (-)
    0.00625            DT          - Integration time step (s)
---------------------- DEGREES OF FREEDOM --------------------------------------
True                   FlapDOF1    - First flapwise blade mode DOF (flag)
True                   FlapDOF2    - Second flapwise blade mode DOF (flag)
True                   EdgeDOF     - First edgewise blade mode DOF (flag)
True                   PitchDOF    - Blade pitch DOF (flag)
False                  TeetDOF     - Rotor-teeter DOF (flag)
False                  DrTrDOF     - Drivetrain rotational-flexibility DOF (flag)
True                   GenDOF      - Generator DOF (flag)
False                  YawDOF      - Nacelle-yaw DOF (flag)
True                   TwFADOF1    - First fore-aft tower bending-mode DOF (flag)
True                   TwFADOF2    - Second fore-aft tower bending-mode DOF (flag)
True                   TwSSDOF1    - First side-to-side tower bending-mode DOF (flag)
True                   TwSSDOF2    - Second side-to-side tower bending-mode DOF (flag)
False                  PtfmSgDOF   - Platform horizontal surge translation DOF (flag)
False                  PtfmSwDOF   - Platform horizontal sway translation DOF (flag)
False                  PtfmHvDOF   - Platform vertical heave translation DOF (flag)
False                  PtfmRDOF    - Platform roll tilt rotation DOF (flag)
False                  PtfmPDOF    - Platform pitch tilt rotation DOF (flag)
False                  PtfmYDOF    - Platform yaw rotation DOF (flag)
---------------------- INITIAL CONDITIONS --------------------------------------
          0            OoPDefl     - Initial out-of-plane blade-tip displacement (m)
          0            IPDefl      - Initial in-plane blade-tip displacement (m)
          0            BlPitch(1)  - Blade 1 initial pitch (degrees)
          0            BlPitch(2)  - Blade 2 initial pitch (degrees)
          0            BlPitch(3)  - Blade 3 initial pitch (degrees)
          0            TeetDefl    - Initial or fixed teeter angle (degrees)
          0            Azimuth     - Initial azimuth angle for blade 1 (degrees)
       12.1            RotSpeed    - Initial or fixed rotor speed (rpm)
          0            NacYaw      - Initial or fixed nacelle-yaw angle (degrees)
          0            TTDspFA     - Initial fore-aft tower-top displacement (m)
          0            TTDspSS     - Initial side-to-side tower-top displacement (m)
          0            PtfmSurge   - Initial platform surge (m)
          0            PtfmSway    - Initial platform sway (m)
          0            PtfmHeave   - Initial platform heave (m)
          0            PtfmRoll    - Initial platform roll (deg)
          0            PtfmPitch   - Initial platform pitch (deg)
          0            PtfmYaw     - Initial platform yaw (deg)
---------------------- TURBINE CONFIGURATION -----------------------------------
          3            NumBl       - Number of blades (-)
       63.0            TipRad      - The distance from the rotor apex to the blade tip (meters)
        1.5            HubRad      - The distance from the rotor apex to the blade root (meters)
       -2.5            PreCone(1)  - Blade 1 cone angle (degrees)
       -2.5            PreCone(2)  - Blade 2 cone angle (degrees)
       -2.5            PreCone(3)  - Blade 3 cone angle (degrees)
          0            HubCM       - Distance from rotor apex to hub mass (meters)
          0            UndSling    - Undersling length (meters)
          0            Delta3      - Delta-3 angle for teetering rotors (degrees)
          0            AzimB1Up    - Azimuth value to use for I/O when blade 1 points up (degrees)
       -5.0            OverHang    - Distance from yaw axis to rotor apex (meters)
        1.9            ShftGagL    - Distance from rotor apex to shaft strain gages (meters)
       -5.0            ShftTilt    - Rotor shaft tilt angle (degrees)
       -3.09528        NacCMxn     - Downwind distance from tower-top to nacelle CM (meters)
          0            NacCMyn     - Lateral distance from tower-top to nacelle CM (meters)
        1.75           NacCMzn     - Vertical distance from tower-top to nacelle CM (meters)
       -3.09528        NcIMUxn     - Downwind distance from tower-top to nacelle IMU (meters)
          0            NcIMUyn     - Lateral distance from tower-top to nacelle IMU (meters)
        1.75           NcIMUzn     - Vertical distance from tower-top to nacelle IMU (meters)
        1.96256        Twr2Shft    - Vertical distance from tower-top to rotor shaft (meters)
       87.6            TowerHt     - Height of tower above ground level (meters)
       10.0            TowerBsHt   - Height of tower base above ground level (meters)
          0            PtfmCMxt    - Downwind distance from ground to platform CM (meters)
          0            PtfmCMyt    - Lateral distance from ground to platform CM (meters)
          0            PtfmCMzt    - Vertical distance from ground to platform CM (meters)
          0            PtfmRefxt   - Downwind distance from ground to platform ref point (meters)
          0            PtfmRefyt   - Lateral distance from ground to platform ref point (meters)
          0            PtfmRefzt   - Vertical distance from ground to platform ref point (meters)
---------------------- MASS AND INERTIA ----------------------------------------
          0            TipMass(1)  - Tip-brake mass, blade 1 (kg)
          0            TipMass(2)  - Tip-brake mass, blade 2 (kg)
          0            TipMass(3)  - Tip-brake mass, blade 3 (kg)
          0            PBrIner(1)  - Pitch bearing inertia, blade 1 (kg m^2)
          0            PBrIner(2)  - Pitch bearing inertia, blade 2 (kg m^2)
          0            PBrIner(3)  - Pitch bearing inertia, blade 3 (kg m^2)
          0            BlPIner(1)  - Blade pitch inertia, blade 1 (kg m^2)
          0            BlPIner(2)  - Blade pitch inertia, blade 2 (kg m^2)
          0            BlPIner(3)  - Blade pitch inertia, blade 3 (kg m^2)
      56780            HubMass     - Hub mass (kg)
     115926            HubIner     - Hub inertia about rotor axis (kg m^2)
          0            HubIner_Teeter - Hub inertia about teeter axis (kg m^2)
        534.116        GenIner     - Generator inertia about HSS (kg m^2)
     240000            NacMass     - Nacelle mass (kg)
    2607890            NacYIner    - Nacelle inertia about yaw axis (kg m^2)
          0            YawBrMass   - Yaw bearing mass (kg)
          0            PtfmMass    - Platform mass (kg)
          0            PtfmRIner   - Platform inertia for roll tilt rotation (kg m^2)
          0            PtfmPIner   - Platform inertia for pitch tilt rotation (kg m^2)
          0            PtfmYIner   - Platform inertia for yaw rotation (kg m^2)
          0            PtfmXYIner  - Platform xy moment of inertia (kg m^2)
          0            PtfmYZIner  - Platform yz moment of inertia (kg m^2)
          0            PtfmXZIner  - Platform xz moment of inertia (kg m^2)
---------------------- BLADE ---------------------------------------------------
         49            BldNodes    - Number of blade nodes (per blade) used for analysis (-)
"test_blade.dat"       BldFile(1)  - Name of file containing properties for blade 1
"test_blade.dat"       BldFile(2)  - Name of file containing properties for blade 2
"test_blade.dat"       BldFile(3)  - Name of file containing properties for blade 3
---------------------- ROTOR-TEETER --------------------------------------------
          0            TeetMod     - Rotor-teeter spring/damper model (switch)
          0            TeetDmpP    - Rotor-teeter damper position (degrees)
          0            TeetDmp     - Rotor-teeter damping constant (N-m/(rad/s))
          0            TeetCDmp    - Rotor-teeter rate-independent Coulomb-damping moment (N-m)
          0            TeetSStP    - Rotor-teeter soft-stop position (degrees)
          0            TeetHStP    - Rotor-teeter hard-stop position (degrees)
          0            TeetSSSp    - Rotor-teeter soft-stop linear-spring constant (N-m/rad)
          0            TeetHSSp    - Rotor-teeter hard-stop linear-spring constant (N-m/rad)
---------------------- YAW-FRICTION --------------------------------------------
          0            YawFrctMod  - Yaw-friction model (switch)
          0            M_CSmax     - Maximum static Coulomb friction torque
          0            M_FCSmax    - Maximum static Coulomb friction torque proportional to shear force
          0            M_MCSmax    - Maximum static Coulomb friction torque proportional to bending moment
          0            M_CD        - Dynamic Coulomb friction moment
          0            M_FCD       - Dynamic Coulomb friction moment proportional to shear force
          0            M_MCD       - Dynamic Coulomb friction moment proportional to bending moment
          0            sig_v       - Linear viscous friction coefficient
          0            sig_v2      - Quadratic viscous friction coefficient
          0            OmgCut      - Yaw angular velocity cutoff
---------------------- DRIVETRAIN ----------------------------------------------
         97.0          GBoxEff     - Gearbox efficiency (%)
       97.028          GBRatio     - Gearbox ratio (-)
    8.676e8            DTTorSpr    - Drivetrain torsional spring (N-m/rad)
    6.215e6            DTTorDmp    - Drivetrain torsional damper (N-m/(rad/s))
---------------------- FURLING -------------------------------------------------
False                  Furling     - Read in additional model properties for furling turbine (flag)
"unused"               FurlFile    - Name of file containing furling properties
---------------------- TOWER ---------------------------------------------------
         20            TwrNodes    - Number of tower nodes used for analysis (-)
"test_tower.dat"       TwrFile     - Name of file containing tower properties
---------------------- OUTPUT --------------------------------------------------
True                   SumPrint    - Print summary data to "<RootName>.sum" (flag)
          1            OutFile     - Switch to determine where output will be placed
False                  TabDelim    - Use tab delimiters in text tabular output file? (flag)
"ES10.3E2"             OutFmt      - Format used for text tabular output
          0            TStart      - Time to begin tabular output (s)
          1            DecFact     - Decimation factor for tabular output (-)
          0            NTwGages    - Number of tower nodes that have strain gages for output
          0            TwrGagNd    - List of tower nodes that have strain gages
          0            NBlGages    - Number of blade nodes that have strain gages for output
          0            BldGagNd    - List of blade nodes that have strain gages
                   OutList             - The next line(s) contains a list of output parameters.
END of OutList section (the word "END" must appear in the first 3 columns of the last OutList line)
---------------------------------------------------------------------------------------
"""
    p = tmp_path / "ElastoDyn.dat"
    p.write_text(ed_content)
    return p

