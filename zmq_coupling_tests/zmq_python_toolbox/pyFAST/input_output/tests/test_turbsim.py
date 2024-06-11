import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test
from pyFAST.input_output import TurbSimFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('TurbSim_*.*', TurbSimFile)

    def test_TurbSim(self):
        # --- Test without tower
        F = TurbSimFile(os.path.join(MyDir,'TurbSim_NoTwr.bts'))
        F.write(      os.path.join(MyDir,'TurbSim_NoTwr_TMP.bts'))
        F2= TurbSimFile(os.path.join(MyDir,'TurbSim_NoTwr_TMP.bts'))
        os.remove(    os.path.join(MyDir,'TurbSim_NoTwr_TMP.bts'))
        np.testing.assert_almost_equal(F['u'][0,:,:,:],F2['u'][0,:,:,:],4)
        np.testing.assert_almost_equal(F['u'][1,:,:,:],F2['u'][1,:,:,:],4)
        np.testing.assert_almost_equal(F['u'][2,:,:,:],F2['u'][2,:,:,:],4)
        # --- Test with tower
        F = TurbSimFile(os.path.join(MyDir,'TurbSim_WithTwr.bts'))
        np.testing.assert_almost_equal(F['u'][2,-1,1,3], 0.508036, 5)
        np.testing.assert_almost_equal(F['u'][0, 4,2,0], 7.4867466, 5)
        np.testing.assert_almost_equal(F['uTwr'][0, 4, :], [6.1509, 6.4063, 8.9555, 7.6943], 4)
        F.write(      os.path.join(MyDir,'TurbSim_WithTwr_TMP.bts'))
        F2= TurbSimFile(os.path.join(MyDir,'TurbSim_WithTwr_TMP.bts'))
        os.remove(    os.path.join(MyDir,'TurbSim_WithTwr_TMP.bts'))
        np.testing.assert_almost_equal(F['u'][0,:,:,:],F2['u'][0,:,:,:],3)
        np.testing.assert_almost_equal(F['u'][1,:,:,:],F2['u'][1,:,:,:],3)
        np.testing.assert_almost_equal(F['u'][2,:,:,:],F2['u'][2,:,:,:],3)

if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
