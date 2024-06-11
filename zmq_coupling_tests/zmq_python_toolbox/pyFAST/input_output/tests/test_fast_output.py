import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output import FASTOutputFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTOut*.*', FASTOutputFile)

    def DF(self,FN):
        """ Reads a file and return a dataframe """ 
        return FASTOutputFile(os.path.join(MyDir,FN)).toDataFrame()
 
    def test_FASTOut(self):
        self.assertEqual(self.DF('FASTOut.out').values[-1,1],1036)
 
    def test_FASTOutBin(self):
        # --- Test reading
        F = FASTOutputFile(os.path.join(MyDir,'FASTOutBin.outb'))
        M = F.toDataFrame()
        self.assertAlmostEqual(M['GenPwr_[kW]'].values[-1],40.57663190807828)
        # --- Test writing
        tempFilename = '_FASTOutBin_out.outb'
        # Write to tempfile
        F.write(tempFilename)
        # Read written file
        F2= FASTOutputFile(tempFilename)
        # Test that read data match
        np.testing.assert_almost_equal(F.data,F2.data, 4)
        np.testing.assert_almost_equal(F.data[-1,-1] ,40.57663190807828, 10)
        np.testing.assert_almost_equal(F2.data[-1,-1],40.57663190807828, 10)
        self.assertEqual(F2.info['attribute_names'][-1],'GenPwr')
        self.assertEqual(F2.info['attribute_units'][-1],'kW')
        # cleanup
        try:
            os.remove(tempFilename)
        except:
            pass

if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
