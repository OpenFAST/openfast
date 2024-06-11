import unittest
import os
import numpy as np
import pyFAST
import pyFAST.linearization as lin

MyDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_ED(self):
        # Check that identify modes works as intended for a standstill case with ElastoDyn
        # NOTE: there are still some issues for the 2nd tower FA
        # 
        BladeLen=40
        TowerLen=55
        lin_files = [os.path.join(MyDir,'../../../data/example_files/Standstill_ForID_ED.1.lin')]

        mbc_data, matData = lin.fx_mbc3(lin_files, verbose=False)
        CD = lin.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)

        modeID_table,modesDesc = lin.IdentifyModes([CD])
        modeIDs= lin.IdentifiedModesDict([CD], modeID_table, modesDesc)[0]
        #for k,v in modeIDs.items():
        #    print(k,v)

        modeID_table_Manual = [0 , 2 , 1 , 4 , 6 , 5 , 7 , 8 , 3 ,     13 , 14 , 9 , 11 , 10 , 0 ] 
        modeID_table_Ref    = [0 , 2 , 1 , 4 , 6 , 5 , 7 , 3 , 8 , 12, 13 , 14 , 9 , 10 , 11 , 0 ] 

        np.testing.assert_equal(modeID_table[:,0],modeID_table_Ref)

        np.testing.assert_equal(modeIDs['1st Tower FA']['ID'],2)
        np.testing.assert_equal(modeIDs['1st Tower SS']['ID'],1)
        np.testing.assert_equal(modeIDs['1st Blade Flap (Regressive)']['ID'],4)
        np.testing.assert_equal(modeIDs['1st Blade Edge (Regressive)']['ID'],7)
        np.testing.assert_equal(modeIDs['2nd Tower FA']['ID'],13) # Should be 13
        np.testing.assert_equal(modeIDs['2nd Tower SS']['ID'],14)
        np.testing.assert_equal(modeIDs['2nd Blade Flap (Regressive)']['ID'],9)

        np.testing.assert_almost_equal(modeIDs['1st Tower FA']['f0'],0.44331, 4)
        np.testing.assert_almost_equal(modeIDs['1st Tower SS']['f0'],0.42219, 4)
        np.testing.assert_almost_equal(modeIDs['1st Blade Flap (Regressive)']['f0'],0.90737, 4)
        np.testing.assert_almost_equal(modeIDs['1st Blade Edge (Regressive)']['f0'],1.82162, 4)
        np.testing.assert_almost_equal(modeIDs['2nd Tower SS']['f0'],4.28604,4)
        np.testing.assert_almost_equal(modeIDs['2nd Blade Flap (Regressive)']['f0'],2.53595,4)

    def test_HD(self):
        # Check that identify modes works as intended for a standstill case with ElastoDyn and Hydrodyn
        # NOTE: there are still some issues for the 2nd tower FA
        # 
        BladeLen=40
        TowerLen=55
        lin_files = [os.path.join(MyDir,'../../../data/example_files/StandstillSemi_ForID_EDHD.1.lin')]

        mbc_data, matData = lin.fx_mbc3(lin_files, verbose=False)
        CD = lin.campbell_diagram_data_oneOP(mbc_data,BladeLen,TowerLen)
        sSummary = lin.campbellData2TXT(CD, nFreqOut=65)

        modeID_table,modesDesc = lin.IdentifyModes([CD])
        modeIDs= lin.IdentifiedModesDict([CD], modeID_table, modesDesc)[0]

#         print(sSummary)
# 
#         for k,v in modeIDs.items():
#             print(k,v)
# 
#         print(modeID_table)

        np.testing.assert_equal(modeIDs['1st Blade Flap (Regressive)' ]['ID'], 47)
        np.testing.assert_equal(modeIDs['1st Blade Flap (Progressive)']['ID'], 50)
        #np.testing.assert_equal(modeIDs['1st Blade Edge (Regressive)' ]['ID'], 55)
#         np.testing.assert_equal(modeIDs['1st Blade Edge (Progressive)']['ID'], 56)
        np.testing.assert_equal(modeIDs['2nd Blade Flap (Regressive)' ]['ID'], 58)
        np.testing.assert_equal(modeIDs['2nd Blade Flap (Progressive)']['ID'], 59)
        np.testing.assert_equal(modeIDs['2nd Blade Flap (Collective)' ]['ID'], 60)
        np.testing.assert_equal(modeIDs['1st Blade Edge (Collective)' ]['ID'], 51)

        np.testing.assert_equal(modeIDs['Platform sway' ]['ID'], 1)
        np.testing.assert_equal(modeIDs['Platform surge' ]['ID'], 2)
        np.testing.assert_equal(modeIDs['Platform heave' ]['ID'], 7)
        np.testing.assert_equal(modeIDs['Platform yaw' ]['ID'], 3)


        # Test that might need to updated after the indentification is improved
        # NOTE: the results below are probaly not good
#         np.testing.assert_equal(modeIDs['1st Tower FA']['ID'],52)
#         np.testing.assert_equal(modeIDs['1st Tower SS']['ID'],53)
#         np.testing.assert_equal(modeIDs['1st Blade Flap (Collective)' ]['ID'], 54) 
#         np.testing.assert_equal(modeIDs['2nd Tower SS' ]['ID'], 62)

if __name__ == '__main__':
    #Test().test_ED()
    #Test().test_HD()
    unittest.main()
