import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test
from pyFAST.input_output import FASTLinearizationFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTLin*.*', FASTLinearizationFile)

    def test_FASTLin(self):

        # --- Test basic read
        F=FASTLinearizationFile(os.path.join(MyDir,'FASTLin.lin'))
        self.assertAlmostEqual(F['A'][3,1], 3.91159454E-04 )
        self.assertAlmostEqual(F['u'][7]   ,4.00176055E+04)

        # Test properties
        np.testing.assert_almost_equal(F.nx  , 4)
        np.testing.assert_almost_equal(F.nu  , 9)
        np.testing.assert_almost_equal(F.ny  , 16)
        np.testing.assert_almost_equal(F.nz  , 0 )

        # Test keys
        np.testing.assert_almost_equal(F['Azimuth']  ,  5.8684 , 4)
        np.testing.assert_almost_equal(F['RotSpeed'] ,  1.2367 , 4)
        self.assertEqual(F['WindSpeed'],  None  ) # NOTE: might become NaN in the future

        # --- Test methods
        dfs = F.toDataFrame() # Make sure this runs

        # Test EVA
        fd, zeta, Q, f0 = F.eva()
        np.testing.assert_almost_equal(f0  , [0.394858], 4)
        np.testing.assert_almost_equal(zeta, [0.06078], 4)

        # Test state removal
        F.removeStates(pattern='generator')
        fd, zeta, Q, f0 = F.eva()
        dfs = F.toDataFrame() # Make sure this runs
        np.testing.assert_almost_equal(F.nx  , 2)
        np.testing.assert_almost_equal(f0  , [0.394258], 4)
        np.testing.assert_almost_equal(zeta, [0.0603 ], 4)


        # --- Test lin file with M (only for EB's special branch...)
        F=FASTLinearizationFile(os.path.join(MyDir,'FASTLin_EDM.lin'))
        dfs=F.toDataFrame()
        M=dfs['M']
        self.assertAlmostEqual(M['7_TwFADOF1']['7_TwFADOF1'],0.436753E+06)
        self.assertAlmostEqual(M['13_GeAz']['13_GeAz']     , 0.437026E+08)



if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
