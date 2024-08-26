import pytest
import os.path as osp
import platform

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

