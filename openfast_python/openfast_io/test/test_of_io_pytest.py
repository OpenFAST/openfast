import pytest
import os.path as osp
import subprocess, sys
import platform

from openfast_io.FAST_reader import InputReader_OpenFAST
from openfast_io.FAST_writer import InputWriter_OpenFAST
from openfast_io.FAST_output_reader import FASTOutputFile

from openfast_io.FileTools import check_rtest_cloned
from pathlib import Path

REPOSITORY_ROOT = osp.dirname(osp.dirname(osp.dirname(osp.dirname(__file__))))
RTESTS_DIR = osp.join(REPOSITORY_ROOT, "reg_tests","r-test")
TEST_DATA_DIR = osp.join(RTESTS_DIR, "glue-codes", "openfast")

RUN_DIR = osp.join(REPOSITORY_ROOT, "build_ofio", "testSuite")
Path(RUN_DIR).mkdir(parents=True, exist_ok=True)


# Exercising the  various OpenFAST modules
FOLDERS_TO_RUN = [
    "5MW_Land_DLL_WTurb"                    ,       # "openfast;elastodyn;aerodyn15;servodyn")
    "5MW_OC3Mnpl_DLL_WTurb_WavesIrr"        ,       # "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
    "5MW_OC3Mnpl_DLL_WTurb_WavesIrr_Restart",       # "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore;restart")
    "5MW_ITIBarge_DLL_WTurb_WavesIrr"       ,       # "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
    "5MW_OC4Semi_WSt_WavesWN"               ,       # "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn;offshore")
    "5MW_Land_BD_DLL_WTurb"                 ,       # "openfast;beamdyn;aerodyn15;servodyn")
    "5MW_OC4Jckt_ExtPtfm"                   ,       # "openfast;elastodyn;extptfm")
    "HelicalWake_OLAF"                      ,       # "openfast;aerodyn15;olaf")
    "StC_test_OC4Semi"                      ,       # "openfast;servodyn;hydrodyn;moordyn;offshore;stc")
    "MHK_RM1_Fixed"                         ,       # "openfast;elastodyn;aerodyn15;mhk")
    "MHK_RM1_Floating"                      ,       # "openfast;elastodyn;aerodyn15;hydrodyn;moordyn;mhk")
    "Tailfin_FreeYaw1DOF_PolarBased"        ,       # "openfast;elastodyn;aerodyn15")
]

# looking up OS for the correct executable extension
mactype = platform.system().lower()
if mactype in ["linux", "linux2", "darwin"]:
    exeExt = ""
elif mactype in ["win32", "windows", "cygwin"]: #NOTE: platform.system()='Windows', sys.platform='win32'
    libext = '.exe'
else:
    raise ValueError('Unknown platform type: '+mactype)

# Path to the OpenFAST executable
of_path = osp.join(REPOSITORY_ROOT,"build/glue-codes/openfast",f"openfast{exeExt}")


def read_action(folder):
    print(f"Reading from {folder}")

    # Read input deck
    fast_reader = InputReader_OpenFAST()
    fast_reader.FAST_InputFile =  f'{folder}.fst'   # FAST input file (ext=.fst)
    fast_reader.FAST_directory = osp.join(TEST_DATA_DIR, folder)   # Path to fst directory files
    fast_reader.execute()

    return fast_reader.fst_vt



def write_action(folder, fst_vt):
    print(f"Writing to {folder}, with TMax = 2.0")

    fast_writer = InputWriter_OpenFAST()
    fast_writer.FAST_runDirectory = osp.join(RUN_DIR,folder)
    Path(fast_writer.FAST_runDirectory).mkdir(parents=True, exist_ok=True)
    fast_writer.FAST_namingOut    = folder

    fast_writer.fst_vt = dict(fst_vt)
    fst_vt = {}
    fst_vt['Fst', 'TMax'] = 2.
    fst_vt['Fst','OutFileFmt'] = 3
    fast_writer.update(fst_update=fst_vt)
    fast_writer.execute()

def run_action(folder):
    # Placeholder for the actual run action
    print(f"Running simulation for {folder}")
    subprocess.run([of_path, str(osp.join(RUN_DIR, folder, f"{folder}.fst"))], check=True)

def check_ascii_out(folder):
    # Placeholder for the actual check action
    print(f"Checking ASCII output for {folder}")
    asciiOutput = osp.join(RUN_DIR, folder, f"{folder}.out")
    fast_outout = FASTOutputFile(filename=asciiOutput)

def check_binary_out(folder):
    # Placeholder for the actual check action
    print(f"Checking binary output for {folder}")
    binaryOutput = osp.join(RUN_DIR, folder, f"{folder}.outb")
    fast_outout = FASTOutputFile(filename=binaryOutput)

# Begining of the test
def test_rtest_cloned():
    if check_rtest_cloned(TEST_DATA_DIR):
        assert True, "R-tests cloned properly"
    else:# stop the test if the r-tests are not cloned properly
        print("R-tests not cloned properly")
        sys.exit(1)

def test_openfast_executable_exists():
    if osp.exists(of_path):
        assert True, f"OpenFAST executable found at {of_path}"
    else: # stop the test if the OpenFAST executable is not found
        print(f"OpenFAST executable not found at {of_path}. Please build OpenFAST and try again.")
        sys.exit(1)

# # Define a list of action functions for parameterization
# actions = [
#     ("read", read_action),
#     ("write", write_action),
#     ("run", run_action),
#     ("check", check_action),
# ]

# Parameterize the test function to run for each folder and action
@pytest.mark.parametrize("folder", FOLDERS_TO_RUN)
# @pytest.mark.parametrize("action_name, action_func", actions)
def test_openfast_io_with_detailed_reporting(folder):
    try:
        # action_func(folder)
        action_name = "read"
        fst_vt = read_action(folder)

        action_name = "write"
        write_action(folder, fst_vt)

        action_name = "run"
        run_action(folder)

        action_name = "check ASCII"
        check_ascii_out(folder)

        action_name = "check binary"
        check_binary_out(folder)

    except Exception as e:
        pytest.fail(f"Action '{action_name}' for folder '{folder}' failed with exception: {e}")

def main():
    # Initialize any necessary setup here

    for folder in FOLDERS_TO_RUN:
        print(f"Processing folder: {folder}")

        # Assuming read_action, write_action, run_action, and check_action are defined elsewhere
        data = read_action(folder)
        write_action(folder, data)
        run_action(folder)
        check_ascii_out(folder)
        check_binary_out(folder)
        print(f"Successfully processed folder: {folder}")

if __name__ == "__main__":
    main()