import unittest
import os
import numpy as np
from pyFAST.input_output.tests.helpers_for_test import MyDir, reading_test 
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_input_file import ExtPtfmFile
from pyFAST.input_output.fast_input_file import ADPolarFile
from pyFAST.input_output.fast_input_file import EDBladeFile
from pyFAST.input_output.fast_wind_file  import FASTWndFile


class Test(unittest.TestCase):
 
    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTIn*.*', FASTInputFile)

    def test_FASTIn(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_BD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['PitchK'],2.0e+07)
        self.assertAlmostEqual(F['MemberGeom'][-1,2],61.5)
        self.assertAlmostEqual(F['MemberGeom'][-2,3],0.023000)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_BD_bld.dat'))
        F.test_ascii(bCompareWritesOnly=False,bDelete=True)
        self.assertEqual(F['DampingCoeffs'][0][0],0.01)
        # TODO BeamDyn Blade properties are not really "user friendly"
        self.assertAlmostEqual(F['BeamProperties']['span'][1],1.0)
        self.assertAlmostEqual(F['BeamProperties']['K'][1][0,0],1.8e+08) # K11 @ section 2
        self.assertAlmostEqual(F['BeamProperties']['M'][1][0,0],1.2) # M11 @ section 2

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['RotSpeed'],0.2)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED_bld.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['BldEdgSh(6)'],-0.6952)
        F.comment = 'ElastoDyn file'

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED_twr.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['AdjFASt'],1)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertTrue(F['TipLoss'])

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ExtPtfm_SubSef.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['StiffnessMatrix'][2,2],1.96653266e+09)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_HD.dat'))
        #F.test_ascii(bCompareWritesOnly=True,bDelete=True) # TODO
        self.assertAlmostEqual(F['RdtnDT'],0.0125)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_IF_NoHead.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertAlmostEqual(F['Z0'],0.03)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_SbD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['Joints'][0,3],-100)
        self.assertEqual(int(F['Members'][0,1]),1)
        self.assertEqual(int(F['Members'][0,2]),2)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_SD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['PitManRat(1)'],2)
        
    def test_FASTADBld(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15_bld.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertTrue('NumBlNds' in F.keys())
        df = F.toDataFrame()
        self.assertEqual(df['BlChord_[m]'].values[-1], 1.419)
        self.assertTrue('c2_Swp_Approx_[m]' in df.keys())
#         import matplotlib.pyplot as plt
#         fig,axes = plt.subplots(1, 3, sharey=False, figsize=(10.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax= axes[0]
#         ax.plot( df['BlSpn_[m]'], df['BlCrvAC_[m]'], 'k-'    , label='Prebend AC')
#         #ax.plot( df['BlSpn_[m]'], df['x'], '--'    , label='x')
#         ax.plot( df['BlSpn_[m]'], df['c2_Crv_Approx_[m]'], ':'    , label='c2')
#         ax.legend()
#         ax= axes[1]
#         ax.plot( df['BlSpn_[m]'], df['BlSwpAC_[m]'], 'k-'    , label='Sweep AC')
#         #ax.plot( df['BlSpn_[m]'], df['y'], '--'    , label='y')
#         ax.plot( df['BlSpn_[m]'], df['c2_Swp_Approx_[m]'], ':'    , label='c2')
#         ax.legend()
#         ax= axes[2]
#         ax.plot( df['BlSpn_[m]'], df['AC_Approx_[-]'], '-'    , label='AC')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()

    def test_FASTADPol(self):
        #F=FASTInputFile(os.path.join(MyDir,'FASTIn_arf_coords.txt'))
        #print(F.keys())
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15_arfl.dat'))
        df = F.toDataFrame()
        self.assertTrue('Cn_pot_[-]' in df.keys())
        # --- Test Dedicated code
        F = ADPolarFile()
        F.read(os.path.join(MyDir,'FASTIn_AD15_arfl.dat'))

    def test_FASTADPolMulti(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15_arf_multitabs.dat'))
        F.test_ascii(bCompareWritesOnly=False,bDelete=True)

        dfs = F.toDataFrame()
        self.assertTrue('AFCoeff_2' in dfs.keys())

        df1 = dfs['AFCoeff_1']
        df2 = dfs['AFCoeff_2']
        self.assertTrue('Cn_pot_[-]' in df2.keys())

        self.assertEqual(df1.shape[0],23)
        self.assertEqual(df2.shape[0],24)

        F = ADPolarFile(numTabs=2)
        F.write('_DUMMY')

    def test_FASTEDBld(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED_bld.dat'))
        F.test_ascii(bCompareWritesOnly=True, bDelete=True)
        self.assertEqual(F['BldEdgSh(6)'],-0.6952)
        df = F.toDataFrame()
        self.assertAlmostEqual(df['ShapeFlap1_[-]'].values[40],0.8530735996)
        # --- Test Dedicated code
        F = EDBladeFile()
        F.read(os.path.join(MyDir,'FASTIn_ED_bld.dat'))

    def test_FASTExt(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ExtPtfm_SubSef.dat'))
        F.test_ascii(bCompareWritesOnly=False, bDelete=True)
        self.assertEqual(F['StiffnessMatrix'][2,2],1.96653266e+09)
        # --- Test Dedicated code
        F = ExtPtfmFile()
        F.read(os.path.join(MyDir,'FASTIn_ExtPtfm_SubSef.dat'))
        F.test_ascii(bCompareWritesOnly=False, bDelete=True)
        df=F.toDataFrame()
        self.assertAlmostEqual(df['InpF_Fx_[N]'].values[-1], 1660.749680)

    def test_FASTWnd(self):
        F=FASTWndFile(os.path.join(MyDir,'FASTWnd.wnd'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)

    def test_FASTInGraph(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_HD.dat'))
        #graph = F.toGraph()
        #print(graph)
        #self.assertEqual(len(graph.Nodes), 4)
        #self.assertEqual(len(graph.Elements), 3)
# 
        #F=FASTInputFile(os.path.join(MyDir,'FASTIn_SbD.dat'))
        #print(F)
        #graph = F.toGraph()
#         self.assertEqual(len(graph.Nodes), 2)
#         self.assertEqual(len(graph.Elements), 1)
    def test_FASTInMoorDyn(self):
        # MoorDyn version 1
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_MD-v1.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(float(F['LineTypes'][0,1]),0.02)

        # MoorDyn version 2
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_MD-v2.dat'))
        #F.write(os.path.join(MyDir,'FASTIn_MD-v2.dat---OUT'))
        self.assertTrue('Points'    in F.keys())
        self.assertTrue('LineTypes' in F.keys())
        self.assertTrue('LineProp'  in F.keys())
        self.assertEqual(F['LineProp'].shape   , (3,7))
        self.assertEqual(F['LineTypes'].shape  , (1,10))
        self.assertEqual(F['Points'].shape  , (6,9))
        self.assertEqual(len(F['Outlist'])  , 6)
        self.assertEqual(F['Outlist'][0]  , 'FairTen1')
        self.assertEqual(F['LineProp'][0,0] , '1')
        self.assertEqual(F['LineProp'][0,1] , 'main')
        self.assertEqual(F['LineProp'][0,6] , '-')

    def test_FASTInAirfoil(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15_arfl.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertTrue('InterpOrd'  in F.keys())
        self.assertTrue('AFCoeff'    in F.keys())
        self.assertEqual(F['AFCoeff'].shape, (30,4))

if __name__ == '__main__':
    #Test().test_FASTEDBld()
    #Test().test_FASTADBld()
    #Test().test_FASTADPol()
    #Test().test_FASTADPolMulti()
    #Test().test_FASTExt()
    #Test().test_FASTIn()
    unittest.main()
