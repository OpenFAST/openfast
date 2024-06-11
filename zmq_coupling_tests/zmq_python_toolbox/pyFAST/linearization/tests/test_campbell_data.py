import unittest
import os
import glob
import numpy as np
import pyFAST.linearization.mbc as mbc
import pyFAST.linearization.campbell as camp


MyDir=os.path.join(os.path.dirname(__file__))

class Test(unittest.TestCase):

    def mbc3_standstill(self, lin_file):
        # Script Parameters
        BladeLen     = 40.04                # Blade length, used to tune relative modal energy [m]
        TowerLen     = 55.59                # Tower length, used to tune relative modal energy [m]

        # Derived parameters
        lin_files = np.array([lin_file])

        # Performing MBC (NOTE: not stricly necessary without rotation)
        mbc_data, matData = mbc.fx_mbc3(lin_files, verbose=False)
        CD = camp.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)

        nModesMax = np.min([len(CD['Modes']),10])
        Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
        Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])
        LogDec = Damp*100*2*np.pi
        return Freq, Damp, LogDec

    def test_mbc3_standstill(self):
        lin_file     = os.path.join(MyDir,'../../../data/example_files/Standstill.1.lin') 
        Freq, Damp, LogDec = self.mbc3_standstill(lin_file)
        np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.450, 0.669], 3)
        np.testing.assert_almost_equal(LogDec[:3],[1.9505,2.1309,5.0649], 4)

    def test_mbc3_standstill_old(self):
        lin_file = os.path.join(MyDir, '../../../data/example_files/Standstill_old.1.lin')
        Freq, Damp, LogDec = self.mbc3_standstill(lin_file)
        np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.449, 0.667], 3)
        np.testing.assert_almost_equal(LogDec[:3],[1.9497,2.1162,5.0113], 4)

    def test_mbc3_ED_Azi3(self):
        # --- ElastoDyn, second order states only
        # Script Parameters
        BladeLen     = 61.5  # Blade length, used to tune relative modal energy [m] 
        TowerLen     = 87.6  # Tower length, used to tune relative modal energy [m] 
        lin_files = glob.glob(os.path.join(MyDir,'../../../data/linearization_outputs/ws03.0*.lin'))

        # Performing MBC
        mbc_data, matData = mbc.fx_mbc3(lin_files, verbose=False)
        CD = camp.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)

        nModesMax = np.min([len(CD['Modes']),10])
        Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
        Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])*100
        np.testing.assert_almost_equal(Freq[:3],[0.31402747,0.33140717,0.6263423], 3)
        np.testing.assert_almost_equal(Damp[:3],[0.43860177,6.03443072,2.48116647], 4)

    def test_mbc3_EDHD(self):
        # --- ElastoDyn and HydroDyn, mix of second order states and first order
        # Script Parameters
        BladeLen     = 61.5  # Blade length, used to tune relative modal energy [m]
        TowerLen     = 87.6  # Tower length, used to tune relative modal energy [m]
        lin_files = [os.path.join(MyDir,'../../../data/example_files/StandstillSemi_ForID_EDHD.1.lin')] 

        # Performing MBC
        mbc_data, matData = mbc.fx_mbc3(lin_files, verbose=False)
        CD = camp.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)

        nModesMax = np.min([len(CD['Modes']),10])
        Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
        Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])*100
        Freq_ref = [0.00852735,0.00865675,0.01333642 ,0.03147209 ,0.03544405,0.03550137,0.04915103,0.07296769 ,0.07801498 ,0.08493981]
        Damp_ref = [0.16886295,3.43729876,25.10525282,56.73819059,2.56872571,2.45772964,0.18473601,23.22413272,11.67420729,24.10895804]

        np.testing.assert_almost_equal(Freq, Freq_ref, 4)
        np.testing.assert_almost_equal(Damp, Damp_ref, 4)

    def test_mbc3_EDBD(self):
        # --- ElastoDyn and BeamDyn, second order states
        # Script Parameters
        BladeLen     = 100.0  # Blade length, used to tune relative modal energy [m]
        TowerLen     = 100.0  # Tower length, used to tune relative modal energy [m]
        lin_files = [os.path.join(MyDir,'../../../data/example_files/BAR_URC_EDBD.1.lin')] 

        # Performing MBC
        mbc_data, matData = mbc.fx_mbc3(lin_files, verbose=False)
        CD = camp.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)

        nModesMax = np.min([len(CD['Modes']),10])
        Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
        Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])*100
        Freq_ref = [0.18549539,0.18727956,1.12675789,1.25491926,9.03091505,9.05584944,9.31045707,9.36779001,11.30081115,13.78922866]
        Damp_ref = [0.40656628,0.40815288,0.64834856,0.75063003,6.02032061,8.60631293,6.29886588,8.37207306,9.96489352,9.70093435]

        np.testing.assert_almost_equal(Freq, Freq_ref, 4)
        np.testing.assert_almost_equal(Damp, Damp_ref, 4)


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
