import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.mannbox_file import MannBoxFile


class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('MannBox_*.*', MannBoxFile)

    def test_MannBox(self):
        # --- Test read/write
        F = MannBoxFile(os.path.join(MyDir,'MannBox_2x4x8.bin'))
        F.write(        os.path.join(MyDir,'MannBox_2x4x8_TMP.bin'))
        F2= MannBoxFile(os.path.join(MyDir,'MannBox_2x4x8_TMP.bin'))
        os.remove(      os.path.join(MyDir,'MannBox_2x4x8_TMP.bin'))
        np.testing.assert_almost_equal(F['field'].shape ,[2,4,8])
        np.testing.assert_almost_equal(F['field'][:,:,:],F2['field'][:,:,:],8)
        np.testing.assert_almost_equal(F['field'][1,3,5], -3.6654968, 6) 

if __name__ == '__main__':
    unittest.main()
