import unittest
import os.path as osp
import subprocess, sys
import platform

from openfast.FAST_reader import InputReader_OpenFAST
from openfast.FAST_writer import InputWriter_OpenFAST
from openfast.FAST_output_reader import FASTOutputFile

REPOSITORY_ROOT = osp.dirname(osp.dirname(osp.dirname(osp.dirname(__file__))))
RTESTS_DIR = osp.join(REPOSITORY_ROOT, "reg_tests","r-test")
TEST_DATA_DIR = osp.join(RTESTS_DIR, "glue-codes", "openfast", "5MW_Land_DLL_WTurb")

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

class TestOFutils(unittest.TestCase):
    def testOF_Inputs(self):

        # Check if r-tests are available
        if not osp.isdir(RTESTS_DIR):
            Exception("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(RTESTS_DIR))


        # Read input deck
        fast_reader = InputReader_OpenFAST()
        fast_writer = InputWriter_OpenFAST()
        fast_reader.FAST_InputFile =  '5MW_Land_DLL_WTurb.fst'   # FAST input file (ext=.fst)
        fast_reader.FAST_directory = TEST_DATA_DIR   # Path to fst directory files
        fast_writer.FAST_runDirectory = osp.join('temp','OpenFAST')
        fast_writer.FAST_namingOut    = 'nrel5mw'

        with self.subTest('Reading', i=0):
            try:
                fast_reader.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Reading','Success')

        # Test the OF writer
        fast_writer.fst_vt = dict(fast_reader.fst_vt)
        fst_vt = {}
        fst_vt['Fst', 'TMax'] = 2.
        fst_vt['AeroDyn15', 'TwrAero'] = False
        fst_vt['Fst','CompMooring'] = 0
        fst_vt['Fst','CompServo'] = 0
        fst_vt['Fst','OutFileFmt'] = 3
        fast_writer.update(fst_update=fst_vt)

        with self.subTest('Writing', i=1):
            try:
                fast_writer.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Writing','Success')

        with self.subTest('Running', i=2):
            try:
                subprocess.run([of_path, str(osp.join(fast_writer.FAST_runDirectory, fast_writer.FAST_namingOut+'.fst'))])
                self.assertTrue(True)
            except:
                self.assertEqual('Running','Success')

    def testOF_Outputs(self):
        # Read output deck
        asciiOutput = osp.join('temp','OpenFAST',f'nrel5mw.out')
        binaryOutput = osp.join('temp','OpenFAST',f'nrel5mw.outb')
        
        with self.subTest('Reading ASCII output', i=0):
            try:
                fast_outout = FASTOutputFile(filename=asciiOutput)
                self.assertTrue(True)
            except:
                self.assertEqual('Writing','Success')

        with self.subTest('Reading Binary output', i=1):
            try:
                fast_outout = FASTOutputFile(filename=binaryOutput)
                self.assertTrue(True)
            except:
                self.assertEqual('Writing','Success')


if __name__ == "__main__":
    unittest.main()
