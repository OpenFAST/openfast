import unittest
from openfast.FAST_reader import InputReader_OpenFAST
from openfast.FAST_writer import InputWriter_OpenFAST
from openfast.FAST_wrapper import FAST_wrapper
from openfast.runFAST_pywrapper import runFAST_pywrapper
from openfast.LinearFAST import LinearFAST
from rosco import discon_lib_path
import os.path as osp
import shutil
import platform


REPOSITORY_ROOT = osp.dirname(osp.dirname(osp.dirname(osp.dirname(__file__))))
TEST_DATA_DIR = osp.join(REPOSITORY_ROOT, "openfast_python", "openfast", "test", "test_data")

of_path = shutil.which('openfast')
bin_dir  = osp.dirname(of_path)
lib_dir  = osp.abspath( osp.join(osp.dirname(bin_dir), 'lib') )


mactype = platform.system().lower()
if mactype in ["linux", "linux2"]:
    libext = ".so"
    staticext = ".a"
elif mactype in ["win32", "windows", "cygwin"]: #NOTE: platform.system()='Windows', sys.platform='win32'
    libext = '.dll'
    staticext = ".lib"
elif mactype == "darwin":
    libext = '.dylib'
    staticext = ".a"
else:
    raise ValueError('Unknown platform type: '+mactype)

found = False
for libname in ['libopenfastlib', 'openfastlib']:
    for d in [lib_dir, bin_dir]:
        lib_path = osp.join(d, libname+libext)
        if osp.exists(lib_path):
            found = True
            break
    if found:
        break

class TestOFutils(unittest.TestCase):

    def testOF_utils(self):
        # Read input deck
        fast_reader = InputReader_OpenFAST()
        fast_writer = InputWriter_OpenFAST()
        fast_wrap   = FAST_wrapper() 
        fast_obj    = runFAST_pywrapper()
        fast_reader.FAST_InputFile = fast_obj.FAST_InputFile = 'IEA-15-240-RWT-UMaineSemi.fst'   # FAST input file (ext=.fst)
        fast_reader.FAST_directory = fast_obj.FAST_directory = osp.join(TEST_DATA_DIR, 'OpenFAST_models', 'IEA-15-240-RWT', 'IEA-15-240-RWT-UMaineSemi')   # Path to fst directory files
        fast_writer.FAST_runDirectory = fast_wrap.FAST_directory = osp.join('temp','OpenFAST')
        fast_obj.FAST_runDirectory    = osp.join('temp2','OpenFAST')
        fast_obj.FAST_namingOut       = fast_writer.FAST_namingOut    = 'iea15'

        with self.subTest('Reading', i=0):
            try:
                fast_reader.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Reading','Success')

        # Test the OF writer
        fast_writer.fst_vt = dict(fast_reader.fst_vt)
        fst_vt = {}
        fst_vt['Fst', 'TMax'] = 20.
        fst_vt['AeroDyn15', 'TwrAero'] = False
        fst_vt['ServoDyn', 'DLL_FileName'] = discon_lib_path
        fst_vt['Fst','CompMooring'] = 0
        fast_writer.update(fst_update=fst_vt)
        with self.subTest('Writing', i=1):
            try:
                fast_writer.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Writing','Success')

        # Execute the written file
        fast_wrap.FAST_exe = of_path   # Path to executable
        fast_wrap.FAST_InputFile = osp.join(fast_writer.FAST_namingOut+'.fst')
        with self.subTest('Running', i=2):
            try:
                fast_wrap.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Running','Success')

        # Test the whole sequence in one go
        fast_obj.FAST_exe = None
        fast_obj.FAST_lib = lib_path
        fast_obj.case     = fst_vt
        with self.subTest('Batching', i=3):
            try:
                fast_obj.execute()
                self.assertTrue(True)
            except:
                self.assertEqual('Batching','Success')

    def testLinearFAST(self):

        lin_fast = LinearFAST(debug_level=2)
        lin_fast.FAST_exe = of_path
        lin_fast.FAST_InputFile           = 'IEA-15-240-RWT-Monopile.fst'   # FAST input file (ext=.fst)
        lin_fast.FAST_directory           = osp.join(TEST_DATA_DIR, "OpenFAST_models", "IEA-15-240-RWT", "IEA-15-240-RWT-Monopile")   # Path to fst directory files
        lin_fast.FAST_runDirectory        = osp.join(REPOSITORY_ROOT,'outputs','iea_mono_lin')
        
        # Read monopile input for fst_vt
        fast_read = InputReader_OpenFAST()
        fast_read.FAST_InputFile = lin_fast.FAST_InputFile
        fast_read.FAST_directory = lin_fast.FAST_directory
        try:
            fast_read.execute()
            self.assertTrue(True)
        except:
            self.assertEqual('Reading','Success')
            
        lin_fast.fst_vt = fast_read.fst_vt
        lin_fast.v_rated                    = 10.74  # needed as input from RotorSE or something, to determine TrimCase for linearization
        lin_fast.wind_speeds                 = [16]
        lin_fast.DOFs                       = ['GenDOF','TwFADOF1'] #,'PtfmPDOF']  # enable with 
        lin_fast.TMax                       = 50   # should be 1000-2000 sec or more with hydrodynamic states
        lin_fast.NLinTimes                  = 2
        
        try:
            lin_fast.gen_linear_cases(inputs={'U_init':lin_fast.wind_speeds,'pitch_init':[12.86]})
            #lin_fast.gen_linear_model()     #"Unable to find steady-state solution."
            self.assertTrue(True)
        except:
            self.assertEqual('Linear','Success')


if __name__ == "__main__":
    unittest.main()
