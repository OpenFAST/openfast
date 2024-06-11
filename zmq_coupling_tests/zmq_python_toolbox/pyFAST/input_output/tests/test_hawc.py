import unittest
import os
import numpy as np
import pyFAST.input_output as weio
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.hawc2_dat_file import HAWC2DatFile
from pyFAST.input_output.hawc2_ae_file import HAWC2AEFile
from pyFAST.input_output.hawc2_pc_file import HAWC2PCFile
from pyFAST.input_output.hawc2_st_file import HAWC2StFile
from pyFAST.input_output.hawcstab2_ind_file import HAWCStab2IndFile
from pyFAST.input_output.hawcstab2_pwr_file import HAWCStab2PwrFile


class Test(unittest.TestCase):
 
    def test_001_read_all(self):
        reading_test('HAWC*.*', weio.read)

    def DF(self,FN):
        """ Reads a file with weio and return a dataframe """ 
        return weio.read(os.path.join(MyDir,FN)).toDataFrame()

    def test_HAWC2(self):
        F=HAWC2DatFile(os.path.join(MyDir,'HAWC2_out_ascii.dat'))
        DF=F.toDataFrame()
        self.assertEqual(DF.values[-1,1],-1.72572E+03)
        self.assertEqual(DF.values[-1,-1], 3.63349E+03)
        self.assertEqual(DF.columns[0], 'Time_[s]')
        self.assertEqual(DF.columns[1], 'WSPgl.coo.,Vy_[m/s]')

        # Test that "exported dat files" are the same
        # NOTE: cannot do comparison of sel files since names are different
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        os.remove(os.path.join(MyDir,'HAWC2_out_ascii_TMP.sel'))
        os.remove(os.path.join(MyDir,'HAWC2_out_ascii_TMP2.sel'))

    def test_HAWC2_st(self):
        # --- not FPM
        F=HAWC2StFile(os.path.join(MyDir,'HAWC2_st.st'))
        dfs=F.toDataFrame()
        set11=dfs['1_1']
        set22=dfs['2_2']
        self.assertEqual(set11['m_[kg/m]'].values[-1], 2536.27)
        self.assertEqual(set11['A_[m^2]'].values[-1], 0.298)
        self.assertEqual(set22['r_[m]'].values[-1], 1.96256)
        self.assertEqual(set22['ri_x_[m]'].values[-1], 1.36)
        # --- FPM
        F=HAWC2StFile(os.path.join(MyDir,'HAWC2_st_fpm.st'))
        dfs=F.toDataFrame()
        set11=dfs['1_1']
        np.testing.assert_almost_equal(set11['m_[kg/m]'].values[-1], 5.6348074, 3)
        np.testing.assert_almost_equal(set11['K66'].values[-1], 8.41526513e04, 3)

    def test_HAWC2_pc(self):
        F=HAWC2PCFile(os.path.join(MyDir,'HAWC2_pc.dat'))
        self.assertEqual(len(F.data.pc_sets),1)
        thicknesses = F.data.pc_sets[1][0]
        firstPolar  = F.data.pc_sets[1][1][0]
        np.testing.assert_almost_equal(thicknesses, [24.1, 30.1, 36, 48, 60, 100])
        self.assertEqual(firstPolar.shape, (105,4))
        np.testing.assert_almost_equal(firstPolar[0,0], -180)
        np.testing.assert_almost_equal(firstPolar[-1,0], 180)

    def test_BHAWC(self):
        F=HAWC2DatFile(os.path.join(MyDir,'BHAWC_out_ascii.sel'))
        DF=F.toDataFrame()
        self.assertEqual(DF.values[-1,1], 147.85)
        self.assertEqual(DF.columns[0], 't_[s]')
        self.assertEqual(DF.columns[1], 'ang_azi_[deg]')

        # Testing that "exported" sel files are the same
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        os.remove(os.path.join(MyDir,'BHAWC_out_ascii_TMP.dat'))
        os.remove(os.path.join(MyDir,'BHAWC_out_ascii_TMP2.dat'))

        # Testing that "exported" dat files are the same
        F=HAWC2DatFile(os.path.join(MyDir,'BHAWC_out_ascii.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        os.remove(os.path.join(MyDir,'BHAWC_out_ascii_TMP.sel'))
        os.remove(os.path.join(MyDir,'BHAWC_out_ascii_TMP2.sel'))

    def test_HAWCStab2(self):
        # power file
        F=HAWCStab2PwrFile(os.path.join(MyDir,'HAWCStab2.pwr'))
        DF=F.toDataFrame()
        self.assertAlmostEqual(DF.values[-1,1],0.1553480512E+05)
        self.assertAlmostEqual(DF.values[-1,-1], 0.3181950053E+09)
        self.assertEqual(DF.columns[0], 'V_[m/s]')
        self.assertEqual(DF.columns[1], 'P_[kW]')
        # induction files
        F=HAWCStab2IndFile(os.path.join(MyDir,'HAWCStab2_u3000.ind'))  # normal .ind
        DF=F.toDataFrame()
        self.assertAlmostEqual(DF.values[-1,1],0.517961E+00)
        self.assertAlmostEqual(DF.values[-1,-1], 0.354614E-02)
        self.assertEqual(DF.columns[0], 's_[m]')
        self.assertEqual(DF.columns[1], 'A_[-]')
        F=HAWCStab2IndFile(os.path.join(MyDir,'HAWCStab2_defl_u3000.ind'))  # defl .ind
        DF=F.toDataFrame()
        self.assertAlmostEqual(DF.values[-1,1],19)
        self.assertAlmostEqual(DF.values[-1,-1], 0.242932E-05)
        self.assertEqual(DF.columns[0], 's_[m]')
        self.assertEqual(DF.columns[1], 'Element_no_[-]')
        F=HAWCStab2IndFile(os.path.join(MyDir,'HAWCStab2_fext_u3000.ind'))  # fext .ind
        DF=F.toDataFrame()
        self.assertAlmostEqual(DF.values[-1,1],20)
        self.assertAlmostEqual(DF.values[-1,-1], -0.170519E+03)
        self.assertEqual(DF.columns[0], 's_[m]')
        self.assertEqual(DF.columns[1], 'Node_[-]')

if __name__ == '__main__':
    #Test().test_HAWC2_st()
    unittest.main()
