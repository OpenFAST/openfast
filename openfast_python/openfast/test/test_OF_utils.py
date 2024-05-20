import unittest
import os.path as osp
import subprocess
import platform

from openfast.FAST_reader import InputReader_OpenFAST
from openfast.FAST_writer import InputWriter_OpenFAST

REPOSITORY_ROOT = osp.dirname(osp.dirname(osp.dirname(osp.dirname(__file__))))
TEST_DATA_DIR = osp.join(REPOSITORY_ROOT, "openfast_python", "openfast", "test", "test_data")

mactype = platform.system().lower()
if mactype in ["linux", "linux2", "darwin"]:
    exeExt = ""
elif mactype in ["win32", "windows", "cygwin"]: #NOTE: platform.system()='Windows', sys.platform='win32'
    libext = '.exe'
else:
    raise ValueError('Unknown platform type: '+mactype)

of_path = osp.join("../../../build/glue-codes/openfast",f"openfast{exeExt}")
bin_dir  = osp.dirname(of_path)

class TestOFutils(unittest.TestCase):

    def testOF_utils(self):
        # Read input deck
        fast_reader = InputReader_OpenFAST()
        fast_writer = InputWriter_OpenFAST()
        fast_reader.FAST_InputFile =  'IEA-15-240-RWT-UMaineSemi.fst'   # FAST input file (ext=.fst)
        fast_reader.FAST_directory = osp.join(TEST_DATA_DIR, 'OpenFAST_models', 'IEA-15-240-RWT', 'IEA-15-240-RWT-UMaineSemi')   # Path to fst directory files
        fast_writer.FAST_runDirectory = osp.join('temp','OpenFAST')
        fast_writer.FAST_namingOut    = 'iea15'

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


if __name__ == "__main__":
    unittest.main()
